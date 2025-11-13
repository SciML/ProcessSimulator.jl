@component function ρTz_ThermodynamicState_(;medium, name)

    pars = []

    vars = @variables begin
        ϕ(t)[1:medium.Constants.nphases - 1],                                                  [description = "phase fraction"]
        ρ(t)[1:medium.Constants.nphases],                                                      [description = "molar density"] 
        z(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                               [description = "mole fraction", irreducible = true] 
        p(t),                                                                                  [description = "pressure"]
        T(t),                                                                                  [description = "Temperature"]
    end

    eq = [
        [ρ[j] ~ PT_molar_density(medium.EoSModel, p, T, collect(z[:, j]), phase = medium.Constants.phase_names[j]) for j in 2:medium.Constants.nphases]...
        ρ[1] .~ 1.0/(sum([ϕ[j - 1]/ρ[j] for j in 2:medium.Constants.nphases]))...
    ]

    if medium.Constants.nphases > 2
        eq_extra = [sum(collect(ϕ)) ~ 1.0 
        ϕ[1] ~ flash_vaporized_fraction(medium.EoSModel, p, T, collect(z[:, 1]))[1]]
        guesses_ϕ = [ϕ[j - 1] => medium.Guesses.ϕ[j - 1] for j in 2:medium.Constants.nphases]

    else
        eq_extra = []
        guesses_ϕ = []
    end

    guesses_z = [z[i, j] => medium.Guesses.x[i, j] for i ∈ 1:medium.Constants.Nc, j ∈ 1:medium.Constants.nphases]

    guesses_p = [p => medium.Guesses.p]

    guesses_T = [T => medium.Guesses.T]


    ODESystem([eq...; eq_extra...], t, collect(Iterators.flatten(vars)), pars; name, guesses = [guesses_z...; guesses_p...; guesses_T...; guesses_ϕ...])

end


@connector function PhZConnector_(;medium, name)

    vars = @variables begin
        p(t),                                                                                 [description = "pressure"]                      
        h(t)[1:medium.Constants.nphases],                                                     [description = "molar enthalpy"]  
        z(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                              [description = "overall mole fraction"]
        ṅ(t)[1:medium.Constants.nphases],                                                     [description = "molar flow", connect = Flow]              
    end

    pars = []

    eqs = []

	eq = eqs==[] ? Equation[] : eqs

    ODESystem(eq, t, collect(Iterators.flatten(vars)), pars; name)

end

@component function ConstantFlowRate(; medium, name, flowrate, flowbasis = :molar)
    
    pars = @parameters begin
        flowrate = flowrate                                                           
    end

    vars = []

    systems = @named begin
        Port = PhZConnector_(medium = medium)
    end

    eqs = [
          Port.ṅ[1] ~ XtoMolar(flowrate, medium, nothing, flowbasis)
        ]

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; systems, name)

end

@component function ConstantPressure(; medium, name, p)
    
    pars = @parameters begin
        pressure = p                                                           
    end

    vars = []

    systems = @named begin
        Port = PhZConnector_(medium = medium)
    end

    eqs = [
          Port.p ~ pressure
        ]

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; systems, name)
end

@component function ConnHouse(; medium, name)
    
    pars = []

    vars = []

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
    end

    System(Equation[], t, vars, pars; systems, name)

end


@connector function SurfaceConnector_(;medium, name)

    vars = @variables begin
        T(t),                                                                                      [description = "Temperature", output = true]                      
        μ(t)[1:medium.Constants.Nc],                                                               [description = "Chemical potential or an estimate of it", output = true]  
        ϕₘ(t)[1:medium.Constants.Nc],                                                              [description = "Molar Flow Rate", output = true]
        ϕₕ(t),                                                                                     [description = "Heat Flow Rate", output = true]       
    end

    pars = []

    ODESystem(Equation[], t, collect(Iterators.flatten(vars)), pars; name)

end

@component function Surface(;medium, name)

    systems = @named begin
        InPort = SurfaceConnector_(medium = medium)
        OutPort = SurfaceConnector_(medium = medium)
    end

    vars = @variables begin
        kₘ(t)[1:medium.Constants.Nc],                                                                                  [description = "mass transfer coefficient"]
        kₕ(t),                                                                                                         [description = "heat transfer coefficient"]
    end

    eqs = [
            scalarize(kₘ .~ mass_transfer_coefficient(medium, InPort.T, InPort.μ))
            kₕ ~ heat_transfer_coefficient(medium, InPort.T, InPort.μ)
            InPort.ϕₕ ~ kₕ*(InPort.T - OutPort.T)
            scalarize(InPort.ϕₘ .~ kₘ.*(InPort.μ .- OutPort.μ))...
            scalarize(InPort.ϕₘ .+ OutPort.ϕₘ .~ 0)
            InPort.ϕₕ + OutPort.ϕₕ ~ 0
    ]

    pars = []

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; systems, name)

end


@component function TwoPortControlVolume_(;medium, name)

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    vars = @variables begin
        rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "mass source or sink - volumetric  basis"]
        rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "molar source or sink - through surface in contact with other phases"]      
        Nᵢ(t)[1:medium.Constants.Nc],                                                                  [description = "molar holdup"]
        nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                        [description = "molar holdup in each phase, excluding the overall phase"]              
        U(t),                                                                                          [description = "internal energy holdup"]                                               
        V(t)[1:medium.Constants.nphases],                                                              [description = "volume"]    
        Q(t),                                                                                          [description = "heat flux"]
        Wₛ(t),                                                                                         [description = "shaft work"]
    end

    pars = []


 
    eqs = [

        # Energy balance

        D(U) ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) + Q + Wₛ

        # Mole balances

        [D(Nᵢ[i]) ~ InPort.ṅ[1]*InPort.z[i, 1] + sum(dot(collect(OutPort.ṅ[2:end]), collect(OutPort.z[i:i, 2:end]))) + sum(collect(rᵥ[i, 2:end].*V[2:end])) + rₐ[i, 1] for i in 1:medium.Constants.Nc]...

        scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))...

        scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))...

        [Nᵢ[i] ~ sum(dot(nᴸⱽ, collect(ControlVolumeState.z[i, 2:end]))) for i ∈ 1:medium.Constants.Nc]...

        scalarize(sum(Nᵢ) ~ sum(nᴸⱽ))

        scalarize(ControlVolumeState.z[:, 1] .~ Nᵢ ./ sum(collect(Nᵢ)))...

        ControlVolumeState.ϕ[1] ~ nᴸⱽ[1]/sum(collect(nᴸⱽ))

        
        # Control Volume properties
        V[1] ~ sum(collect(V[2:end]))
        V[2]*ControlVolumeState.ρ[2] ~ nᴸⱽ[1]
        V[3]*ControlVolumeState.ρ[3] ~ nᴸⱽ[2]
        
  
        # Outlet port properties
        [OutPort.h[j] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.Constants.nphases]...
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        OutPort.p ~ ControlVolumeState.p

        scalarize(OutPort.z .~ ControlVolumeState.z)...

        OutPort.ṅ[2] ~ ControlVolumeState.ϕ[1]*OutPort.ṅ[1]
        OutPort.ṅ[3] ~ ControlVolumeState.ϕ[2]*OutPort.ṅ[1]


        ]

        return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])

end


@component function TwoPortControlVolume_SteadyState(;medium, name)

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    vars = @variables begin
        rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "mass source or sink - volumetric  basis"]
        rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "molar source or sink - through surface in contact with other phases"]      
        Nᵢ(t)[1:medium.Constants.Nc],                                                                  [description = "molar holdup"]
        nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                        [description = "molar holdup in each phase, excluding the overall phase"]              
        U(t),                                                                                          [description = "internal energy holdup"]                                               
        V(t)[1:medium.Constants.nphases],                                                              [description = "volume"]    
        Q(t),                                                                                          [description = "heat flux"]
        Wₛ(t),                                                                                         [description = "shaft work"]
    end

    pars = []
 
    eqs = [

        # Steady-state energy balance (no accumulation)
        0 ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) + Q + Wₛ

        # Steady-state mole balances (no accumulation)
        [0 ~ InPort.ṅ[1]*InPort.z[i, 1] + sum(dot(collect(OutPort.ṅ[2:end]), collect(OutPort.z[i:i, 2:end]))) + sum(collect(rᵥ[i, 2:end].*V[2:end])) + rₐ[i, 1] for i in 1:medium.Constants.Nc]...

        scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))...

        scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))...

        [Nᵢ[i] ~ sum(dot(nᴸⱽ, collect(ControlVolumeState.z[i, 2:end]))) for i ∈ 1:medium.Constants.Nc]...

        scalarize(sum(Nᵢ) ~ sum(nᴸⱽ))

        scalarize(ControlVolumeState.z[:, 1] .~ Nᵢ ./ sum(collect(Nᵢ)))...

        ControlVolumeState.ϕ[1] ~ nᴸⱽ[1]/sum(collect(nᴸⱽ))

        
        # Control Volume properties
        V[1] ~ sum(collect(V[2:end]))
        V[2]*ControlVolumeState.ρ[2] ~ nᴸⱽ[1]
        V[3]*ControlVolumeState.ρ[3] ~ nᴸⱽ[2]
        
  
        # Outlet port properties
        #[OutPort.h[j] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.Constants.nphases]...
        [OutPort.h[j] ~ pT_enthalpy(medium.EoSModel,
                                    ControlVolumeState.p, ControlVolumeState.T,
                                    collect(ControlVolumeState.z[:, j]),
                                    medium.Constants.phase_names[j]) for j in 2:medium.Constants.nphases]...
                                    
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        OutPort.p ~ ControlVolumeState.p

        scalarize(OutPort.z .~ ControlVolumeState.z)...

        OutPort.ṅ[2] ~ ControlVolumeState.ϕ[1]*OutPort.ṅ[1]
        OutPort.ṅ[3] ~ ControlVolumeState.ϕ[2]*OutPort.ṅ[1]


        ]

        return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])

end


@component function ThreePortControlVolume_(;medium, name)
    """
    Three-port control volume: 1 inlet, 2 separate outlets (liquid and vapor)
    Designed for flash drums with phase-separated exits
    """

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        LiquidOutPort = PhZConnector_(medium = medium)
        VaporOutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    vars = @variables begin
        rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "mass source or sink - volumetric basis"]
        rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "molar source or sink - through surface"]
        Nᵢ(t)[1:medium.Constants.Nc],                                                                  [description = "molar holdup"]
        nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                        [description = "molar holdup in each phase"]
        U(t),                                                                                          [description = "internal energy holdup"]
        V(t)[1:medium.Constants.nphases],                                                              [description = "volume"]
        Q(t),                                                                                          [description = "heat flux"]
        Wₛ(t),                                                                                         [description = "shaft work"]
    end

    pars = []

    eqs = [
        # Energy balance
        D(U) ~ InPort.h[1]*InPort.ṅ[1] +
               LiquidOutPort.h[1]*LiquidOutPort.ṅ[1] +
               VaporOutPort.h[1]*VaporOutPort.ṅ[1] + Q + Wₛ

        # Component mole balances
        [D(Nᵢ[i]) ~ InPort.ṅ[1]*InPort.z[i, 1] +
                     LiquidOutPort.ṅ[1]*LiquidOutPort.z[i, 1] +
                     VaporOutPort.ṅ[1]*VaporOutPort.z[i, 1] +
                     sum(collect(rᵥ[i, 2:end].*V[2:end])) + rₐ[i, 1] for i in 1:medium.Constants.Nc]...

        scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))...
        scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))...

        # Holdup relationships
        [Nᵢ[i] ~ sum(dot(nᴸⱽ, collect(ControlVolumeState.z[i, 2:end]))) for i ∈ 1:medium.Constants.Nc]...
        scalarize(sum(Nᵢ) ~ sum(nᴸⱽ))

        scalarize(ControlVolumeState.z[:, 1] .~ Nᵢ ./ (sum(collect(Nᵢ)) + 1e-10))...
        ControlVolumeState.ϕ[1] ~ nᴸⱽ[1]/(sum(collect(nᴸⱽ)) + 1e-10)

        # Control Volume properties
        V[1] ~ sum(collect(V[2:end]))
        V[2]*ControlVolumeState.ρ[2] ~ nᴸⱽ[1]
        V[3]*ControlVolumeState.ρ[3] ~ nᴸⱽ[2]

        # Liquid outlet - carries liquid phase (index 2)
        LiquidOutPort.h[2] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[2], ControlVolumeState.T, collect(ControlVolumeState.z[:, 2]))
        LiquidOutPort.h[3] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[3], ControlVolumeState.T, collect(ControlVolumeState.z[:, 3]))
        LiquidOutPort.h[1] ~ LiquidOutPort.h[2]  # Overall enthalpy = liquid enthalpy (pure liquid exit)

        # Liquid outlet composition - only liquid phase exits
        scalarize(LiquidOutPort.z[:, 2] .~ ControlVolumeState.z[:, 2])...
        scalarize(LiquidOutPort.z[:, 3] .~ ControlVolumeState.z[:, 3])...  # Keep for consistency
        scalarize(LiquidOutPort.z[:, 1] .~ ControlVolumeState.z[:, 2])...  # Overall = liquid composition

        # Liquid outlet flows - only liquid phase
        LiquidOutPort.ṅ[2] ~ LiquidOutPort.ṅ[1]  # Total flow = liquid flow
        LiquidOutPort.ṅ[3] ~ 1e-10  # Minimal vapor (numerically stable)

        # Vapor outlet - carries vapor phase (index 3)
        VaporOutPort.h[2] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[2], ControlVolumeState.T, collect(ControlVolumeState.z[:, 2]))
        VaporOutPort.h[3] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[3], ControlVolumeState.T, collect(ControlVolumeState.z[:, 3]))
        VaporOutPort.h[1] ~ VaporOutPort.h[3]  # Overall enthalpy = vapor enthalpy (pure vapor exit)

        VaporOutPort.p ~ ControlVolumeState.p  # Vapor pressure = tank pressure

        # Vapor outlet composition - only vapor phase exits
        scalarize(VaporOutPort.z[:, 2] .~ ControlVolumeState.z[:, 2])...  # Keep for consistency
        scalarize(VaporOutPort.z[:, 3] .~ ControlVolumeState.z[:, 3])...
        scalarize(VaporOutPort.z[:, 1] .~ ControlVolumeState.z[:, 3])...  # Overall = vapor composition

        # Vapor outlet flows - only vapor phase
        VaporOutPort.ṅ[3] ~ VaporOutPort.ṅ[1]  # Total flow = vapor flow
        VaporOutPort.ṅ[2] ~ 1e-10  # Minimal liquid (numerically stable)
    ]

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])
end


@component function ThreePortControlVolume_SteadyState(;medium, name)
    """
    Three-port control volume (steady-state): 1 inlet, 2 separate outlets (liquid and vapor)
    """

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        LiquidOutPort = PhZConnector_(medium = medium)
        VaporOutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    vars = @variables begin
        rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "mass source or sink - volumetric basis"]
        rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "molar source or sink - through surface"]
        Nᵢ(t)[1:medium.Constants.Nc],                                                                  [description = "molar holdup"]
        nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                        [description = "molar holdup in each phase"]
        U(t),                                                                                          [description = "internal energy holdup"]
        V(t)[1:medium.Constants.nphases],                                                              [description = "volume"]
        Q(t),                                                                                          [description = "heat flux"]
        Wₛ(t),                                                                                         [description = "shaft work"]
    end

    pars = []

    eqs = [
        # Steady-state energy balance
        0 ~ InPort.h[1]*InPort.ṅ[1] +
            LiquidOutPort.h[1]*LiquidOutPort.ṅ[1] +
            VaporOutPort.h[1]*VaporOutPort.ṅ[1] + Q + Wₛ

        # Steady-state component mole balances
        [0 ~ InPort.ṅ[1]*InPort.z[i, 1] +
             LiquidOutPort.ṅ[1]*LiquidOutPort.z[i, 1] +
             VaporOutPort.ṅ[1]*VaporOutPort.z[i, 1] +
             sum(collect(rᵥ[i, 2:end].*V[2:end])) + rₐ[i, 1] for i in 1:medium.Constants.Nc]...

        scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))...
        scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))...

        # Holdup relationships
        [Nᵢ[i] ~ sum(dot(nᴸⱽ, collect(ControlVolumeState.z[i, 2:end]))) for i ∈ 1:medium.Constants.Nc]...
        scalarize(sum(Nᵢ) ~ sum(nᴸⱽ))

        scalarize(ControlVolumeState.z[:, 1] .~ Nᵢ ./ (sum(collect(Nᵢ))))...
        ControlVolumeState.ϕ[1] ~ nᴸⱽ[1]/(sum(collect(nᴸⱽ)))

        # Control Volume properties
        V[1] ~ sum(collect(V[2:end]))
        V[2]*ControlVolumeState.ρ[2] ~ nᴸⱽ[1]
        V[3]*ControlVolumeState.ρ[3] ~ nᴸⱽ[2]

        # Liquid outlet properties
        LiquidOutPort.h[2] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[2], ControlVolumeState.T, collect(ControlVolumeState.z[:, 2]))
        LiquidOutPort.h[3] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[3], ControlVolumeState.T, collect(ControlVolumeState.z[:, 3]))
        LiquidOutPort.h[1] ~ LiquidOutPort.h[2]

        LiquidOutPort.p ~ ControlVolumeState.p
        scalarize(LiquidOutPort.z[:, 2] .~ ControlVolumeState.z[:, 2])...
        scalarize(LiquidOutPort.z[:, 3] .~ ControlVolumeState.z[:, 3])...
        scalarize(LiquidOutPort.z[:, 1] .~ ControlVolumeState.z[:, 2])...

        LiquidOutPort.ṅ[2] ~ LiquidOutPort.ṅ[1]
        LiquidOutPort.ṅ[3] ~ 1e-10

        # Vapor outlet properties
        VaporOutPort.h[2] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[2], ControlVolumeState.T, collect(ControlVolumeState.z[:, 2]))
        VaporOutPort.h[3] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[3], ControlVolumeState.T, collect(ControlVolumeState.z[:, 3]))
        VaporOutPort.h[1] ~ VaporOutPort.h[3]

        VaporOutPort.p ~ ControlVolumeState.p  # Vapor pressure = tank pressure
        scalarize(VaporOutPort.z[:, 2] .~ ControlVolumeState.z[:, 2])...
        scalarize(VaporOutPort.z[:, 3] .~ ControlVolumeState.z[:, 3])...
        scalarize(VaporOutPort.z[:, 1] .~ ControlVolumeState.z[:, 3])...

        VaporOutPort.ṅ[3] ~ VaporOutPort.ṅ[1]
        VaporOutPort.ṅ[2] ~ 1e-10
    ]

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])
end


@component function ClosedControlVolume_(;medium, name)

    systems = @named begin
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    vars = @variables begin
        rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "mass source or sink - volumetric  basis"]
        rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "molar source or sink - through surface in contact with other phases"]
        Nᵢ(t)[1:medium.Constants.Nc],                                                                  [description = "molar holdup"]
        nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                        [description = "molar holdup in each phase, excluding the overall phase"]
        U(t),                                                                                          [description = "internal energy holdup"]
        V(t)[1:medium.Constants.nphases],                                                              [description = "volume"]
        Q(t),                                                                                          [description = "heat flux"]
        Wₛ(t),                                                                                         [description = "shaft work"]
    end

    pars = []

    eqs = [

        # Energy balance

        D(U) ~  Q + Wₛ

        # Mole balance

        [D(Nᵢ[i]) ~ sum(collect(rᵥ[i, 2:end].*V[2:end])) + rₐ[i, 1] for i in 1:medium.Constants.Nc]...

        scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))...

        scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))...

        [Nᵢ[i] ~ sum(dot(nᴸⱽ, collect(ControlVolumeState.z[i, 2:end]))) for i ∈ 1:medium.Constants.Nc]...

        scalarize(sum(Nᵢ) ~ sum(nᴸⱽ))


        # Thermodynamic state equations

        scalarize(ControlVolumeState.z[:, 1] .~ Nᵢ ./ sum(collect(Nᵢ)))...

        ControlVolumeState.ϕ[1] ~ nᴸⱽ[1]/sum(collect(nᴸⱽ))


        # Control Volume properties

        V[1] ~ sum(collect(V[2:end]))


        ]

        return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])

end

"""
Struct wrapper for Boundary_pTzn component
"""
mutable struct Boundary_pTzn{M <: AbstractFluidMedium, S <: AbstractThermodynamicState}
    medium::M
    state::S
    flowrate::Real
    flowbasis::Symbol
    odesystem
    model_type::Symbol
end

"""
    Boundary_pTzn(; medium, state = nothing, p = nothing, T = nothing, z = nothing, flowrate, flowbasis = :molar, name)

Create a fixed boundary condition (inlet/outlet stream) with specified state and flowrate.

# Arguments
- `medium`: Medium specification
- `state`: Thermodynamic state (e.g., pTzState, pTNVState) [optional if p, T, z provided]
- `p`: Pressure [Pa] [optional if state provided]
- `T`: Temperature [K] [optional if state provided]
- `z`: Mole fraction vector [optional if state provided]
- `flowrate`: Flow rate (basis specified by `flowbasis`)
- `flowbasis`: Flow basis (`:molar`, `:mass`, or `:volumetric`) [default: `:molar`]
- `name`: Component name

# Returns
- `Boundary_pTzn` struct with embedded ODESystem and model_type = :Feed

# Examples
```julia
# Using state object
state = pTzState(10*101325.0, 350.15, [0.4, 0.6, 0.0])
@named S1 = Boundary_pTzn(medium = medium, state = state, flowrate = 100.0)

# Using p, T, z directly
@named S1 = Boundary_pTzn(medium = medium, p = 10*101325.0, T = 350.15,
                          z = [0.4, 0.6, 0.0], flowrate = 100.0)
```
"""
function Boundary_pTzn(; medium, state = nothing, p = nothing, T = nothing, z = nothing, flowrate, flowbasis = :molar, name)
    # Determine which form of input was provided
    if !isnothing(state)
        # Use the provided state
        if !isnothing(p) || !isnothing(T) || !isnothing(z)
            error("Cannot specify both 'state' and individual 'p', 'T', 'z' parameters")
        end
    elseif !isnothing(p) && !isnothing(T) && !isnothing(z)
        # Create state from p, T, z
        state = pTNVState(p, T, z, :Pressure)
    else
        error("Must provide either 'state' or all of 'p', 'T', and 'z'")
    end

    # Use resolve_guess! to update medium and state (just like CSTR)
    medium, state, _ = resolve_guess!(medium, state)

    # Extract p, T, z from state (works with both pTzState and pTNVState)
    p_val = state.p
    T_val = state.T
    z_val = if hasproperty(state, :z)
        state.z
    elseif hasproperty(state, :N)
        # Convert molar amounts to mole fractions
        state.N ./ sum(state.N)
    else
        error("State must have either 'z' (mole fractions) or 'N' (molar amounts)")
    end

    odesystem = Boundary_pTzn_Model(medium = medium, p = p_val, T = T_val, z = z_val, flowrate = flowrate, flowbasis = flowbasis, name = name)
    return Boundary_pTzn(medium, state, flowrate, flowbasis, odesystem, :Feed)
end

@component function Boundary_pTzn_Model(;medium, p, T, z, flowrate, name, flowbasis = :molar)

    systems = @named begin
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end
    
    vars = []

    pars = []

    eqs = [
        # Port equations
        OutPort.p ~ p
        scalarize(OutPort.z[:, 1] .~ z)...
        scalarize(OutPort.z[:, 2:end] .~ ControlVolumeState.z[:, 2:end])...

        OutPort.ṅ[1] ~ -XtoMolar(flowrate, medium, ControlVolumeState, flowbasis)
        OutPort.ṅ[2] ~ ControlVolumeState.ϕ[1]*OutPort.ṅ[1]
        OutPort.ṅ[3] ~ ControlVolumeState.ϕ[2]*OutPort.ṅ[1]
        #OutPort.ṅ[2] + OutPort.ṅ[3] ~ OutPort.ṅ[1] #Same as ∑ϕ = 1

        #[OutPort.h[i] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[i], ControlVolumeState.T, ControlVolumeState.z[:, i]) for i in 2:medium.Constants.nphases]...
        [OutPort.h[j] ~ pT_enthalpy(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, j]), medium.Constants.phase_names[j]) for j in 2:medium.Constants.nphases]...
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        # State equations
        scalarize(ControlVolumeState.z[:, 1] .~ z)...
        scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, p, T, z))...
        scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, p, T, z))...
        ControlVolumeState.T ~ T
        ControlVolumeState.p ~ p
        #scalarize(ControlVolumeState.ϕ[1] ~ flash_vaporized_fraction(medium.EoSModel, p, T, z)[1])
    ]

    return ODESystem(eqs, t, vars, collect(Iterators.flatten(pars)); name, systems = [systems...])


end

@component function FixedBoundary_pTz_(;medium, p, T, z, name)

    systems = @named begin
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end
    
    vars = []

    pars = []

    eqs = [
        # Port equations
        OutPort.p ~ p
        scalarize(OutPort.z[:, 1] .~ z)...
        scalarize(OutPort.z[:, 2:end] .~ ControlVolumeState.z[:, 2:end])...

        OutPort.ṅ[2] ~ ControlVolumeState.ϕ[1]*OutPort.ṅ[1]
        OutPort.ṅ[3] ~ ControlVolumeState.ϕ[2]*OutPort.ṅ[1]
        

        [OutPort.h[i] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[i], ControlVolumeState.T, collect(ControlVolumeState.z[:, i])) for i in 2:medium.Constants.nphases]...
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        # State equations
        scalarize(ControlVolumeState.z[:, 1] .~ z)...
        scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, p, T, z))...
        scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, p, T, z))...
        ControlVolumeState.T ~ T
        ControlVolumeState.p ~ p
        
    ]

    return ODESystem(eqs, t, vars, collect(Iterators.flatten(pars)); name, systems = [systems...])


end

export ρTz_ThermodynamicState_, PhZConnector_, ConstantFlowRate, ConstantPressure, ConnHouse, SurfaceConnector
export Surface, TwoPortControlVolume_, TwoPortControlVolume_SteadyState
export ThreePortControlVolume_, ThreePortControlVolume_SteadyState
export ClosedControlVolume_, Boundary_pTzn, Boundary_pTzn_Model, FixedBoundary_pTzn_, FixedBoundary_pTz_



