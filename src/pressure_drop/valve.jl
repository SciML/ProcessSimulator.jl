
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
end


function Valve(; medium, state_guess::S, Cv, f, flowrate_guess = 1.0, flowbasis = :volume, name) where S <: pTzState
    p, T, z = state_guess.p, state_guess.T, state_guess.z
    medium.Guesses = EosBasedGuesses(medium.EoSModel, p, T, z, Val(:Pressure))
    phase = ifelse(medium.Guesses.ϕ[2] ≈ 1.0, "vapor", "liquid")
    flowstate = rhoTzState(medium.Guesses.ρ[1], T, z)
    opening_setpoint = 0.5
    molar_flowrate_guess = XtoMolar(flowrate_guess, medium, flowstate, flowbasis)
    odesystem = Valve_(medium = medium, Cv = Cv, ΔP_f = f, setpoint = opening_setpoint, phase = phase, name = name)
    return Valve(medium, state_guess, molar_flowrate_guess, Cv, opening_setpoint, f, phase, odesystem)
end

@component function Valve_(;medium, Cv, ΔP_f, setpoint, phase, name)

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    pars = []

    vars = @variables begin
        opening_setpoint(t), [description = "Valve opening set point"]
        opening(t), [description = "Valve actual opening"]   
    end

    if phase == "vapor"

        phase_eq = [scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...]
    
    else

        phase_eq = [scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...]

    end

    eqs = [

        # Energy balance

        0.0 ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) 

        # Mole balance per component
        [0.0 ~ InPort.ṅ[1]*InPort.z[i, 1] + sum(dot(collect(OutPort.ṅ[2:end]), collect(OutPort.z[i:i, 2:end]))) for i in 1:medium.Constants.Nc]...
        [InPort.z[i, 1] ~ OutPort.z[i, 1] for i in 1:medium.Constants.Nc]...

        0.0 ~ InPort.ṅ[1] + OutPort.ṅ[1] #Mole balance
                 

        # Flow equation
        OutPort.ṅ[1]/ControlVolumeState.ρ[1] ~ -opening*Cv*ΔP_f(InPort.p - OutPort.p)
        OutPort.ṅ[2] ~ OutPort.ṅ[1]*ControlVolumeState.ϕ[1]
        OutPort.ṅ[3] ~ OutPort.ṅ[1]*ControlVolumeState.ϕ[2]


        # Outlet port properties
        [OutPort.h[j] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.Constants.nphases]...
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        # Port properties
        OutPort.p ~ ControlVolumeState.p
        scalarize(OutPort.z .~ ControlVolumeState.z)...


        D(opening) ~ -1.0*(opening - opening_setpoint)
        opening_setpoint ~ setpoint

        ]

        return System([phase_eq...; eqs...], t, collect(Iterators.flatten(vars)), pars; systems, name)

end


@component function ErgunDrop(;medium, solidmedium, tank, phase = "vapor", name)

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    pars = []

    vars = []

    if phase == "vapor"

        phase_eq = [scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...]
    
    else

        phase_eq = [scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...]

    end


    #Imperative assignments
    M̄ = molecular_weight(collect(ControlVolumeState.z[:, 1]), collect(ControlVolumeState.z[:, 1]))
    vs = (InPort.ṅ[1]/ControlVolumeState.ρ[1]/M̄)/(π*(tank.D/2)^2) #imperative assignment 
    ε = tank.porosity
    ΔL = tank.H
    dp = solidmedium.Constants.particle_size*2.0
    μ = viscosity(medium, ControlVolumeState.p, abs(ControlVolumeState.T), collect(ControlVolumeState.z[:, 1]))
    ΔP = (InPort.p - OutPort.p)
    ρ = ControlVolumeState.ρ[1]
    

    eqs = [

        # Energy balance

        0.0 ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) 

        # Mole balance per component
        [0.0 ~ InPort.ṅ[1]*InPort.z[i, 1] + sum(dot(collect(OutPort.ṅ[2:end]), collect(OutPort.z[i:i, 2:end]))) for i in 1:medium.Constants.Nc]...
        [InPort.z[i, 1] ~ OutPort.z[i, 1] for i in 1:medium.Constants.Nc]...

        0.0 ~ InPort.ṅ[1] + OutPort.ṅ[1] #Mole balance
                 

        # Flow equation
        ΔP ~ 150.0*μ*ΔL/dp^2*(1-ε)^2/ε^3*vs + 1.75*ΔL*ρ/dp*(1-ε)/ε^3*abs(vs)*vs
        OutPort.ṅ[2] ~ OutPort.ṅ[1]*ControlVolumeState.ϕ[1]
        OutPort.ṅ[3] ~ OutPort.ṅ[1]*ControlVolumeState.ϕ[2]


        # Outlet port properties
        [OutPort.h[j] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.Constants.nphases]...
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))
        OutPort.p ~ ControlVolumeState.p
        scalarize(OutPort.z .~ ControlVolumeState.z)...


        ]

        return System([phase_eq...; eqs...], t, collect(Iterators.flatten(vars)), pars; systems, name)

end

export Valve, Valve_