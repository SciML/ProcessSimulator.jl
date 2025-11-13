abstract type AbstractAdsorber end

@component function AdsorptionInterface(;fluidmedium, solidmedium, adsorbentmass, name)

    systems = @named begin
        SolidSurface = Surface(medium = solidmedium)
        FluidSurface = Surface(medium = fluidmedium)
    end

    eqs = [ domain_connect(FluidSurface.OutPort, SolidSurface.InPort)
            scalarize(SolidSurface.InPort.μ .~ adsorbentmass*loading(solidmedium.isotherm, collect(FluidSurface.OutPort.μ), FluidSurface.OutPort.T))...
            FluidSurface.OutPort.T ~ SolidSurface.InPort.T
            scalarize(FluidSurface.OutPort.ϕₘ + SolidSurface.InPort.ϕₘ .~ 0.0)...
            FluidSurface.OutPort.ϕₕ + SolidSurface.InPort.ϕₕ ~ 0.0
        ]

    vars = []

    pars = []

    guess_T = [FluidSurface.OutPort.T => solidmedium.Guesses.T]
    guess_μ = [FluidSurface.OutPort.μ[i] => fluidmedium.Guesses.x[i, 1]*fluidmedium.Guesses.p for i in 1:fluidmedium.Constants.Nc]

    System(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...], guesses = [guess_T..., guess_μ...])

end



mutable struct WellMixedAdsorber{FM <: AbstractFluidMedium, SM <: AbstractFluidMedium, S <: AbstractThermodynamicState} <: AbstractAdsorber
    fluidmedium::FM
    solidmedium::SM
    state::S
    phase::String
    odesystem
    model_type::Symbol
end

"""
    WellMixedAdsorber(; fluidmedium, solidmedium, state, porosity, V, name)

Creates a well-mixed adsorber with fluid and solid phases in equilibrium.

# Arguments
- `fluidmedium`: Fluid medium specification
- `solidmedium`: Solid/adsorbent medium specification
- `state`: Thermodynamic state for the fluid phase
- `porosity`: Bed porosity (void fraction)
- `V`: Total volume [m³]
- `name`: Component name

# Returns
- `WellMixedAdsorber` struct with embedded ODESystem and model_type = :Adsorber
"""
function WellMixedAdsorber(; fluidmedium, solidmedium, state::S, porosity, V, name) where S <: AbstractThermodynamicState
    # Use resolve_guess! to update fluidmedium and state
    fluidmedium, state, phase = resolve_guess!(fluidmedium, state)

    # Extract pressure from state
    p = state.p

    odesystem = WellMixedAdsorber_Model(fluidmedium = fluidmedium, solidmedium = solidmedium,
                                        porosity = porosity, V = V, p = p, phase = phase, name = name)

    return WellMixedAdsorber(fluidmedium, solidmedium, state, phase, odesystem, :Adsorber)
end

@component function WellMixedAdsorber_Model(;fluidmedium, solidmedium, porosity, V, p, phase, name)
 
    mass_of_adsorbent = V * solidmedium.EoSModel.ρ_T0 * porosity
    A  = V*(1.0 - porosity)*area_per_volume(solidmedium) #Interfacial area

    systems = @named begin

        mobilephase = TwoPortControlVolume_(medium = fluidmedium)
        stationaryphase = ClosedControlVolume_(medium = solidmedium)
        interface = AdsorptionInterface(fluidmedium = fluidmedium, solidmedium = solidmedium, adsorbentmass = mass_of_adsorbent)
        
    end

    eqs = [
        
        # Custom energy balance
        mobilephase.U ~ (mobilephase.OutPort.h[1] - mobilephase.ControlVolumeState.p/mobilephase.ControlVolumeState.ρ[1])*sum(collect(mobilephase.Nᵢ))
        stationaryphase.U ~ (ρT_enthalpy(solidmedium.EoSModel, stationaryphase.ControlVolumeState.ρ[1], stationaryphase.ControlVolumeState.T, collect(stationaryphase.ControlVolumeState.z[:, 1])))*mass_of_adsorbent #Pressure is not relevant for internal energy of the particle.

        # Fluid Surface ports equalities
        interface.FluidSurface.InPort.T ~ mobilephase.ControlVolumeState.T
        mobilephase.Q ~ interface.FluidSurface.OutPort.ϕₕ*A 

        # Solid Surface ports equalities
        scalarize(interface.SolidSurface.OutPort.μ .~ stationaryphase.Nᵢ)...
        scalarize(interface.SolidSurface.OutPort.T ~ stationaryphase.ControlVolumeState.T)

        # Heat and mass transfer
        stationaryphase.Q ~ interface.SolidSurface.InPort.ϕₕ*A + sum(-collect(isosteric_heat(solidmedium.isotherm, interface.FluidSurface.OutPort.μ, interface.FluidSurface.OutPort.T).*stationaryphase.rₐ[:, end]))
        scalarize(stationaryphase.rₐ[:, end] .~ interface.SolidSurface.InPort.ϕₘ)...

        # No Volumetric sink/source constraints
        scalarize(mobilephase.rᵥ[:, 2:end] .~ 0.0)...
        scalarize(stationaryphase.rᵥ[:, 2:end] .~ 0.0)...

        #Volume constraint
        stationaryphase.V[1] ~ V*(1.0 - porosity) #Volume of stationary phase

        #Control Volume Pressure Equality
        mobilephase.ControlVolumeState.p ~ stationaryphase.ControlVolumeState.p

        #No shaft work
        stationaryphase.Wₛ ~ 0.0
        mobilephase.Wₛ ~ 0.0

    ]

    if phase == "liquid" #MTK can't handle equation change in mid simulation, so pick one phase.
        
        phase_eqs = [

            scalarize(mobilephase.rₐ[:, 2] .~ interface.FluidSurface.OutPort.ϕₘ*A)...
            scalarize(mobilephase.rₐ[:, end] .~ 0.0)...

            scalarize(mobilephase.ControlVolumeState.z[:, end] .~ flash_mol_fractions_vapor(fluidmedium.EoSModel, mobilephase.ControlVolumeState.p, mobilephase.ControlVolumeState.T, collect(mobilephase.ControlVolumeState.z[:, 1])))...

            scalarize(interface.FluidSurface.InPort.μ .~ mobilephase.nᴸⱽ[1]/mobilephase.V[2])... #This is phase specific

            #Perfect pressure control (If you fix volume for liquid phase, it leads to problems as liquid phase is incompressible)
            mobilephase.ControlVolumeState.p ~ p # Only relevant for liquid phase

        ]

    elseif phase == "vapor" 

        phase_eqs = [

            scalarize(mobilephase.rₐ[:, 2] .~ 0.0)...
            scalarize(mobilephase.rₐ[:, end] .~ interface.FluidSurface.OutPort.ϕₘ*A)...
            
            scalarize(mobilephase.ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(fluidmedium.EoSModel, mobilephase.ControlVolumeState.p, mobilephase.ControlVolumeState.T, collect(mobilephase.ControlVolumeState.z[:, 1])))...

            scalarize(interface.FluidSurface.InPort.μ .~ mobilephase.ControlVolumeState.p*mobilephase.ControlVolumeState.z[:, end])... #Partial pressure in vapor phase/rigorously is fugacity

            #Volume Constraint          
            mobilephase.V[1] ~ V*(porosity) #Volume of mobile phase
        ]

    end

    vars = []

    pars = []

    return System([eqs...; phase_eqs...], t, collect(Iterators.flatten(vars)), pars; name, systems = systems)

end

export WellMixedAdsorber, WellMixedAdsorber_Model, AdsorptionInterface