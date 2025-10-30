
# ==================== Dynamic Flash Drum ====================

mutable struct DynamicFlashDrum{M <: AbstractFluidMedium, S <: AbstractThermodynamicState, G <: AbstractTank} <: AbstractSeparator
    medium::M
    state::S
    gemotry::G
    odesystem
end

function DynamicFlashDrum(; medium, state, geometry, Q, W, name)
    medium, state, phase = resolve_guess!(medium, state)
    odesystem = DynamicFlashDrumModel(medium = medium, state = state, Q = Q, W = W, name = name)
    return DynamicFlashDrum(medium, state, geometry, odesystem)
end

@component function DynamicFlashDrumModel(; medium, state, geometry, Q = 0.0, W = 0.0, name)

    @named CV = ThreePortControlVolume_(medium = medium)
    @unpack Nᵢ, V, InPort, LiquidOutPort, VaporOutPort, ControlVolumeState, rₐ, rᵥ = CV

    vars = @variables begin
        h_liquid(t),    [description = "liquid level height (m)"]
        g(t)
    end

    pars = @parameters begin
        g = 9.81,   [description = "gravitational acceleration (m/s²)"]
    end

    # Basic equations
    eqs = [
        # No reactions or surface mass transfer
        scalarize(rₐ[:, 2:end] .~ 0.0)...
        scalarize(rᵥ[:, 2:end] .~ 0.0)...

        # Heat input
        CV.Q ~ Q

        # Shaft work
        CV.Wₛ ~ W

        scalarize(ControlVolumeState.z[:, end] .~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...

        # Internal energy definition
        CV.U ~ (ControlVolumeState.h[1] - ControlVolumeState.p/ControlVolumeState.ρ[1])*sum(collect(Nᵢ))

        # Total volume constraint (design specification)
        CV.V[1] ~ state.V
    ]

    # Flash drum specific equations
    eq_flash = [

        # Liquid level calculation from liquid volume
        h_liquid ~ V[2] / cross_section_area_(geometry)

        # Liquid outlet pressure includes hydrostatic head
        # P_liquid = P_vapor + ρ_liquid * g * h_liquid
        LiquidOutPort.p ~ ControlVolumeState.p + (ControlVolumeState.ρ[2] * molecular_weight(medium.EoSModel, collect(ControlVolumeState.z[:, 2]))) * g * h_liquid
    ]

    pars = []

    return extend(ODESystem([eqs...; eq_flash...], t, collect(Iterators.flatten(vars)), pars; name), CV)
end



export DynamicFlashDrum
