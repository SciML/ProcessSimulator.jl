
abstract type AbstractValve end

mutable struct Valve{M <: AbstractFluidMedium, C <: Real, F <: Function, S <: AbstractString, ST <: AbstractThermodynamicState} <: AbstractValve
    medium::M
    state::ST
    molar_flowrate_guess::C
    Cv::C
    opening_setpoint::C
    f::F
    phase::S
    odesystem
    model_type::Symbol
end


function Valve(; medium, state::S, Cv, f::Function, flowrate_guess = 1.0, flowbasis = :volume, name) where S <: AbstractThermodynamicState
    # Use resolve_guess! to update medium and state (consistent with CSTR and Boundary_pTzn)
    medium, state, phase = resolve_guess!(medium, state)

    # Extract p, T, z from state
    p = state.p
    T = state.T
    z = state.N/sum(state.N)

    flowstate = rhoTzState(medium.Guesses.ρ[1], T, z)
    opening_setpoint = 0.5
    molar_flowrate_guess = XtoMolar(flowrate_guess, medium, flowstate, flowbasis)
    odesystem = Valve__(medium = medium, Cv = Cv, ΔP_f = f, setpoint = opening_setpoint, phase = phase, name = name)
    return Valve(medium, state, molar_flowrate_guess, Cv, opening_setpoint, f, phase, odesystem, :Valve)
end


function Valve__(; medium, Cv, ΔP_f, setpoint, phase, name) 

    odesystem = TwoPortControlVolume0D(medium = medium, phase = phase, name = name)
    @unpack Q, OutPort, InPort, ControlVolumeState = odesystem
    
    extra_vars = @variables begin
        opening(t), [description = "Valve actual opening"]   
    end

    ΔP = InPort.p - OutPort.p

    extra_eqs = [Q ~ 0,
                OutPort.ṅ[1]/ControlVolumeState.ρ[1] ~ -opening*Cv*ΔP_f(ΔP), 
                D(opening) ~ -1.0*(opening - setpoint)]
    
    newsys = extend(System(extra_eqs, t, collect(Iterators.flatten(extra_vars)), []; name), odesystem)
    return newsys
end


# ──────────────────────────────────────────────────────────────────────────────
# ErgunDrop  – packed-bed pressure-drop component (extend-based, like Valve__)
# ──────────────────────────────────────────────────────────────────────────────

abstract type AbstractPorousDrop end

mutable struct ErgunDrop{M <: AbstractFluidMedium, SM, TK, S <: AbstractThermodynamicState} <: AbstractPorousDrop
    medium::M
    solidmedium::SM
    tank::TK
    state::S
    molar_flowrate_guess
    phase::String
    odesystem
    model_type::Symbol
end


function ErgunDrop(; medium, solidmedium, tank, state::S, flowrate_guess = 1.0, flowbasis = :volume, name) where S <: AbstractThermodynamicState
    
    medium, state, phase = resolve_guess!(medium, state)
    T = state.T
    z = state.N / sum(state.N)

    flowstate = rhoTzState(medium.Guesses.ρ[1], T, z)
    molar_flowrate_guess = XtoMolar(flowrate_guess, medium, flowstate, flowbasis)
    odesystem = ErgunDrop__(medium = medium, solidmedium = solidmedium, tank = tank, phase = phase, name = name)
    return ErgunDrop(medium, solidmedium, tank, state, molar_flowrate_guess, phase, odesystem, :ErgunDrop)
end

function ErgunDrop__(; medium, solidmedium, tank, phase = "vapor", name) 

    odesystem = TwoPortControlVolume0D(medium = medium, phase = phase, name = name)
    @unpack Q, OutPort, InPort, ControlVolumeState = odesystem

    extra_vars = @variables begin
        ΔP(t), [description = "Pressure drop across the bed segment"]   
    end
    
    #Imperative assignments for aliasing intermediate symbolic expressions in the Ergun equation (= sign)
    M̄ = molecular_weight(medium.EoSModel, collect(ControlVolumeState.z[:, 1]))
    ṁ = InPort.ṅ[1]/ControlVolumeState.ρ[1]/M̄
    vs = ṁ/(π*(tank.D/2)^2)  
    ε = tank.porosity
    ΔL = tank.H
    dp = solidmedium.Constants.particle_size*2.0
    μ = viscosity(medium, ControlVolumeState.p, abs(ControlVolumeState.T), collect(ControlVolumeState.z[:, 1]))
    _ΔP = InPort.p - OutPort.p
    ρ = ControlVolumeState.ρ[1]


    extra_eqs = [Q ~ 0,
                ΔP ~ _ΔP,
                _ΔP ~ 150.0*μ*ΔL/dp^2*(1-ε)^2/ε^3*vs + 1.75*ΔL*ρ/dp*(1-ε)/ε^3*abs(vs)*vs, 
                ]
    
    newsys = extend(System(extra_eqs, t, collect(Iterators.flatten(extra_vars)), []; name, guesses = [InPort.ṅ[1] => molar_flowrate_guess]), odesystem)
    return newsys
end



export Valve, Valve_, ErgunDrop