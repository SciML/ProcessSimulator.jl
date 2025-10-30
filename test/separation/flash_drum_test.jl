using ProcessSimulator


# Load Clapeyron for thermodynamic properties
using Clapeyron



    
components = ["methane", "propane"]
eos = PR(components, idealmodel = ReidIdeal)  # Peng-Robinson equation of state

# =============================================================================
# Simulation Scenario
# Flash drum originally at 30bar with equimolar methane/propane holdup (global) and a feed mix
# at different composition and vaporized fraction comes in at 30bar 



# =============================================================================
# 2. FLASH DRUM SIZING
# =============================================================================
# Design criteria for vapor-liquid separator:
# - Residence time: 5-10 minutes for liquid
# - L/D ratio: typically 3-5 for horizontal drums, 2-4 for vertical
# - Vapor velocity: limited to avoid entrainment (typically < 1 m/s)

# Vessel sizing based on residence time and vapor velocity
initial_state = Clapeyron.qp_flash(eos, 0.5, 3e6, 100.0*[0.5, 0.5])  # 25°C, 30 bar, equimolar
volume(initial_state) #m^3


# Estimate flash conditions at lower pressure
flash_pressure = 20e5  # 20 bar
flash_temp = initial_state.data.T    # 263.15 K (-10°C) - colder to get more liquid

# Rough estimate: at these conditions, expect ~50% vaporization
# For 5 min residence time with 50 mol/s liquid:
# Liquid volume needed ≈ 5 min × 60 s/min × 50 mol/s / ρ_liquid
# Assuming ρ_liquid ≈ 15000 mol/m³ (typical for light hydrocarbons)
estimated_liquid_holdup = 5 * 60 * 50 / 15000  # ≈ 1 m³

# Design vertical drum with L/D = 3
drum_diameter = 1.5    # m
drum_height = 4.5      # m (L/D = 3)
drum_volume = π * (drum_diameter/2)^2 * drum_height  # ≈ 7.95 m³

geometry = CylindricalTank(D = drum_diameter, L = drum_height)

println("=== Flash Drum Design ===")
println("Diameter: $(drum_diameter) m")
println("Height: $(drum_height) m")
println("Total Volume: $(round(drum_volume, digits=2)) m³")
println("L/D ratio: $(round(drum_height/drum_diameter, digits=2))")
println()


    # =============================================================================
    # 4. VALVE SIZING
    # =============================================================================
    # Inlet valve: pressure drop from 30 bar to 20 bar
    # Cv = Q / √(ΔP) where Q is flow and ΔP is pressure drop
    # For gas: Cv = (Q_scfh / 963) × √(ρ/ΔP) × √(T/520)
    # Simplified: use moderate Cv values

    inlet_valve_cv = 10.0   # Flow coefficient (imperial units equivalent)
    liquid_valve_cv = 5.0   # Smaller valve for liquid control
    vapor_valve_cv = 15.0   # Larger valve for vapor (lower density)

    # Valve pressure drop functions
    # For simplicity: ΔP ∝ (flow)² / Cv²
    function valve_pressure_drop(cv)
        return (state, port) -> begin
            # Simple pressure drop: ΔP = k × ṅ² / Cv²
            # k is a constant, here we use 1e4 for reasonable pressure drops
            k = 1e4
            return k * (port.ṅ[1])^2 / cv^2
        end
    end

    # =============================================================================
    # 5. BUILD SYSTEM COMPONENTS
    # =============================================================================

    # Feed stream (high pressure)
    @named feed = FixedBoundary_pTzn_(
        medium = medium,
        p = feed_pressure,
        T = feed_temp,
        z = feed_composition,
        flowrate = feed_flowrate,
        flowbasis = :molar
    )

    # Inlet valve (throttling valve)
    inlet_valve_state = pTNVState(
        flash_pressure,  # Outlet pressure
        flash_temp,
        feed_composition .* feed_flowrate,
        nothing
    )

    @named inlet_valve = Valve(
        medium = medium,
        Cv = inlet_valve_cv,
        valve_equation = valve_pressure_drop(inlet_valve_cv),
        name = "inlet_valve",
        state = inlet_valve_state
    )

    # Flash drum
    @named flash_drum = SteadyStateFlashDrum(
        medium = medium,
        state = flash_state,
        Q = 0.0  # Adiabatic operation
    )

    # Liquid outlet valve
    liquid_valve_pressure = flash_pressure - 2e5  # Drain to 18 bar
    liquid_valve_state = pTNVState(
        liquid_valve_pressure,
        flash_temp,
        [0.3, 0.7] .* 40.0,  # Assuming liquid-rich in propane
        nothing
    )

    @named liquid_valve = Valve(
        medium = medium,
        Cv = liquid_valve_cv,
        valve_equation = valve_pressure_drop(liquid_valve_cv),
        name = "liquid_valve",
        state = liquid_valve_state
    )

    # Vapor outlet valve
    vapor_valve_pressure = flash_pressure - 1e5  # Vapor to 19 bar
    vapor_valve_state = pTNVState(
        vapor_valve_pressure,
        flash_temp,
        [0.8, 0.2] .* 60.0,  # Assuming vapor-rich in methane
        nothing
    )

    @named vapor_valve = Valve(
        medium = medium,
        Cv = vapor_valve_cv,
        valve_equation = valve_pressure_drop(vapor_valve_cv),
        name = "vapor_valve",
        state = vapor_valve_state
    )

    # Outlet pressure boundaries
    @named liquid_sink = ConstantPressure(medium = medium, p = liquid_valve_pressure)
    @named vapor_sink = ConstantPressure(medium = medium, p = vapor_valve_pressure)

    # =============================================================================
    # 6. CONNECT SYSTEM
    # =============================================================================
    @named flash_system = ODESystem([
        # Feed -> Inlet Valve -> Flash Drum
        connect(feed.OutPort, inlet_valve.InPort)
        connect(inlet_valve.OutPort, flash_drum.InPort)

        # Flash Drum -> Outlet Valves -> Sinks
        connect(flash_drum.LiquidOutPort, liquid_valve.InPort)
        connect(flash_drum.VaporOutPort, vapor_valve.InPort)
        connect(liquid_valve.OutPort, liquid_sink.Port)
        connect(vapor_valve.OutPort, vapor_sink.Port)
    ], t)

    # =============================================================================
    # 7. SIMPLIFY AND SOLVE
    # =============================================================================
    println("=== Simplifying system ===")
    sys_simplified = structural_simplify(flash_system)

    println("Number of equations: ", length(equations(sys_simplified)))
    println("Number of unknowns: ", length(unknowns(sys_simplified)))
    println()

    # Initial conditions
    u0 = [
        sys_simplified.flash_drum.Nᵢ[1] => 100.0,  # Initial methane holdup
        sys_simplified.flash_drum.Nᵢ[2] => 100.0,  # Initial propane holdup
        sys_simplified.flash_drum.ControlVolumeState.T => flash_temp
    ]

    # Guesses for algebraic variables
    guesses = [
        sys_simplified.flash_drum.V[2] => 2.0,  # Initial liquid volume (m³)
        sys_simplified.flash_drum.V[3] => drum_volume - 2.0,  # Initial vapor volume
        sys_simplified.flash_drum.nᴸⱽ[1] => 100.0,  # Initial liquid moles
        sys_simplified.flash_drum.nᴸⱽ[2] => 100.0,  # Initial vapor moles
        sys_simplified.flash_drum.h_liquid => 1.5,  # Initial liquid height (m)
    ]

    # =============================================================================
    # 8. CREATE AND SOLVE PROBLEM
    # =============================================================================
    println("=== Creating ODE problem ===")
    tspan = (0.0, 500.0)  # 500 seconds simulation
    prob = ODEProblem(sys_simplified, u0, tspan, guesses = guesses)

    println("=== Solving system ===")
    sol = solve(prob, FBDF(autodiff = false), abstol = 1e-8, reltol = 1e-6)

    # =============================================================================
    # 9. EXTRACT AND DISPLAY RESULTS
    # =============================================================================
    println("\n=== STEADY-STATE RESULTS ===")

    # Flash drum properties
    liquid_level = sol[sys_simplified.flash_drum.h_liquid][end]
    liquid_volume = sol[sys_simplified.flash_drum.V[2]][end]
    vapor_volume = sol[sys_simplified.flash_drum.V[3]][end]

    println("\nFlash Drum:")
    println("  Liquid level: $(round(liquid_level, digits=3)) m")
    println("  Liquid volume: $(round(liquid_volume, digits=3)) m³")
    println("  Vapor volume: $(round(vapor_volume, digits=3)) m³")
    println("  Liquid fraction: $(round(liquid_volume/drum_volume*100, digits=1))%")

    # Pressures
    drum_pressure = sol[sys_simplified.flash_drum.ControlVolumeState.p][end]
    liquid_outlet_p = sol[sys_simplified.flash_drum.LiquidOutPort.p][end]
    vapor_outlet_p = sol[sys_simplified.flash_drum.VaporOutPort.p][end]

    println("\nPressures:")
    println("  Drum vapor space: $(round(drum_pressure/1e5, digits=2)) bar")
    println("  Liquid outlet: $(round(liquid_outlet_p/1e5, digits=2)) bar")
    println("  Vapor outlet: $(round(vapor_outlet_p/1e5, digits=2)) bar")
    println("  Hydrostatic head: $(round((liquid_outlet_p - vapor_outlet_p)/1e5, digits=3)) bar")

    # Flow rates
    liquid_flow = sol[sys_simplified.flash_drum.LiquidOutPort.ṅ[1]][end]
    vapor_flow = sol[sys_simplified.flash_drum.VaporOutPort.ṅ[1]][end]

    println("\nFlow Rates:")
    println("  Feed: $(round(feed_flowrate, digits=2)) mol/s")
    println("  Liquid outlet: $(round(abs(liquid_flow), digits=2)) mol/s")
    println("  Vapor outlet: $(round(abs(vapor_flow), digits=2)) mol/s")
    println("  Material balance: $(round(feed_flowrate + liquid_flow + vapor_flow, digits=4)) mol/s (should be ≈ 0)")

    # Compositions
    println("\nCompositions (methane/propane):")
    println("  Feed: $(feed_composition)")
    liquid_comp_1 = sol[sys_simplified.flash_drum.LiquidOutPort.z[1,1]][end]
    liquid_comp_2 = sol[sys_simplified.flash_drum.LiquidOutPort.z[2,1]][end]
    println("  Liquid: [$(round(liquid_comp_1, digits=3)), $(round(liquid_comp_2, digits=3))]")

    vapor_comp_1 = sol[sys_simplified.flash_drum.VaporOutPort.z[1,1]][end]
    vapor_comp_2 = sol[sys_simplified.flash_drum.VaporOutPort.z[2,1]][end]
    println("  Vapor: [$(round(vapor_comp_1, digits=3)), $(round(vapor_comp_2, digits=3))]")

    # =============================================================================
    # 10. VERIFICATION TESTS
    # =============================================================================
    @test sol.retcode == ReturnCode.Success
    @test liquid_level > 0.0 && liquid_level < drum_height
    @test liquid_outlet_p > vapor_outlet_p  # Hydrostatic head check
    @test abs(feed_flowrate + liquid_flow + vapor_flow) < 1.0  # Material balance
    @test liquid_comp_2 > feed_composition[2]  # Liquid enriched in heavy component (propane)
    @test vapor_comp_1 > feed_composition[1]   # Vapor enriched in light component (methane)

    println("\n=== ALL TESTS PASSED ===")