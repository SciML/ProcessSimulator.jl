mutable struct Splitter{M <: AbstractFluidMedium, S <: AbstractThermodynamicState} <: AbstractSeparator
    medium::M
    state::S
    odesystem
end

function Splitter(; medium, state, split_ratio, name)
    medium, state, phase = resolve_guess!(medium, state)
    odesystem = SplitterModel(medium = medium, phase = phase, state = state, split_ratio = split_ratio, name = name)
    return Splitter(medium, state, odesystem)
end

@component function SplitterModel(;medium, phase, state, split_ratio, name)

    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        OutPort2 = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    pars = @parameters begin
        α = split_ratio, [description = "split fraction"]
    end

    vars = []

    if phase == "vapor"

        phase_eq = [scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...]
    
    else

        phase_eq = [scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...]

    end

    hₗᵥ = [ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.Constants.nphases]
    h_all = dot(collect(hₗᵥ), collect(ControlVolumeState.ϕ))

    eqs = [

        # Energy balance
        0.0 ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) + OutPort2.h[1]*(OutPort2.ṅ[1]) 

        # Mole balance per component
        [0.0 ~ InPort.ṅ[1]*InPort.z[i, 1] + OutPort.ṅ[1]*OutPort.z[i, 1] + OutPort2.ṅ[1]*OutPort2.z[i, 1] for i in 1:medium.Constants.Nc]...

        [0.0 ~ OutPort.ṅ[1]*OutPort.z[i, 1] + OutPort.ṅ[2]*OutPort.z[i, 2] + OutPort.ṅ[3]*OutPort.z[i, 3] for i in 1:medium.Constants.Nc]...

        # Overall Mole balance
        0.0 ~ InPort.ṅ[1] + OutPort.ṅ[1] + OutPort2.ṅ[1]

        #Branch 1 mole balance
        0.0 ~ OutPort.ṅ[1] + InPort.ṅ[1]*α 

                 
        # Flow equation
        OutPort.ṅ[2] ~ OutPort.ṅ[1]*ControlVolumeState.ϕ[1]
        OutPort.ṅ[3] ~ OutPort.ṅ[1]*ControlVolumeState.ϕ[2]

        OutPort2.ṅ[2] ~ OutPort2.ṅ[1]*ControlVolumeState.ϕ[1]
        OutPort2.ṅ[3] ~ OutPort2.ṅ[1]*ControlVolumeState.ϕ[2]

        # Outlet port1 properties
        scalarize(OutPort.h[2:end] .~ hₗᵥ)
        OutPort.h[1] ~ h_all

        # Outlet port2 properties
        scalarize(OutPort2.h[2:end] .~ hₗᵥ)
        OutPort2.h[1] ~ h_all

        # ControlVolume properties
        OutPort.p ~ ControlVolumeState.p #Sets control volume pressure based on outlet port1 pressure
        scalarize(ControlVolumeState.z .~ OutPort.z)... #Sets control volume composition based on outlet port1 composition

        OutPort2.p ~ ControlVolumeState.p
        scalarize(ControlVolumeState.z .~ OutPort2.z)...

        state.p ~ ControlVolumeState.p

        ]

        return System([phase_eq...; eqs...], t, collect(Iterators.flatten(vars)), pars; systems, name)

end

export Splitter, SplitterModel