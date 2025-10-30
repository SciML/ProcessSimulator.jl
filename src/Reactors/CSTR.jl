mutable struct DynamicCSTR{M <: AbstractFluidMedium, R <: AbstractReaction, S <: AbstractThermodynamicState} <: AbstractReactor
    medium::M
    state::S
    reactionset::R
    phase::S
    odesystem
    model_type::Symbol
end

function DynamicCSTR(; medium, reactionset, state::S, W, Q, name) where S <: pTNVState
    medium, state, phase = resolve_guess!(medium, state)
    odesystem = DynamicCSTRModel(medium = medium, reactions = reactionset, state = state, W = W, Q = Q, phase = phase, name = name)
    return DynamicCSTR(medium, state, reactionset, phase, odesystem, :CSTR)
end


@component function DynamicCSTRModel(;medium, reactions::PowerLawReaction, state, W = 0.0, Q = nothing, phase = "liquid", name)

    @named CV = TwoPortControlVolume_(medium = medium)
    @unpack Nᵢ, V, InPort, OutPort, ControlVolumeState, rₐ, rᵥ, Wₛ = CV

    vars = @variables begin
        cᵢ(t)[1:medium.FluidConstants.Nc],    [description = "bulk concentrations"]   
        X(t)[1:medium.FluidConstants.Nc],     [description = "conversion"]            
    end
    
    if !isnothing(Q)
        q_eq = [CV.Q ~ Q]
    else
        q_eq = []
    end

    eqs = [
        scalarize(cᵢ .~ Nᵢ/V[1])...
        Wₛ ~ W
        scalarize(rₐ[:, 2:end] .~ 0.0)...
        CV.U ~ (ControlVolumeState.h[1] - ControlVolumeState.p/ControlVolumeState.ρ[1])*sum(collect(Nᵢ));
        q_eq...
    ]

    

    if phase == "liquid"

        eq_reaction = [
            
            #Only liquid phase constraints
            scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...

            # Reaction kinetics
            scalarize(rᵥ[:, 2] .~ sum(Rate.(reactions, cᵢ, ControlVolumeState.T)))...
            scalarize(rᵥ[:, 3] .~ 0.0)...

            #Conversion definition
            scalarize(X .~ (InPort.ṅ[2].*InPort.z₂ .+ OutPort.ṅ[2].*OutPort.z₂)./(InPort.ṅ[2].*InPort.z₂ .+ 1e-8))...

            ControlVolumeState.p ~ state.p #Liquid phase do not fix volume since liquid is incompressible

        ]

    elseif phase == "vapor"

        eq_reaction = [
            #Only liquid phase constraints
            scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...

            # Reaction kinetics
            scalarize(rᵥ[:, 2] .~ 0.0)...
            scalarize(rᵥ[:, 3] .~ sum(Rate.(reactions, cᵢ, ControlVolumeState.T)))...
            
            scalarize(X .~ (InPort.ṅ[3].*InPort.z₃ .+ OutPort.ṅ[3].*OutPort.z₃)./(InPort.ṅ[3].*InPort.z₃ .+ 1e-8))...

            CV.V[1] ~ state.V
        ]
    end

    pars = []

    return extend(ODESystem([eqs...;eq_reaction...], t, collect(Iterators.flatten(vars)), pars; name), CV)

end


#Base homogeneous CSTR
mutable struct SteadyStateCSTR{M <: AbstractFluidMedium, R <: AbstractReaction, S <: AbstractString} <: AbstractReactor
    medium::M
    state
    reactionset::R
    phase::S
    odesystem
    model_type::Symbol
end

function SteadyStateCSTR(;medium, reactionset, limiting_reactant, state, W, Q, name)

    medium, state, phase = resolve_guess!(medium, state);

    if isnothing(limiting_reactant)
        limiting_reactant = medium.Constants.iupacName[1]
    end

    odesystem = SteadyStateCSTRModel(medium = medium, reactions = reactionset, 
    limiting_reactant = limiting_reactant, state = state, W = W, Q = Q, phase = phase, name = name)

    if !isnothing(Q) #If heat is given use, else fix temperature and calculate heat
        q_eq = [odesystem.Q ~ Q]
    else    
        @unpack ControlVolumeState = odesystem
        q_eq = [ControlVolumeState.T ~ state.T]
    end

    newsys = extend(System(q_eq, t, [], []; name), odesystem)
    return SteadyStateCSTR(medium, state, reactionset, phase, newsys, :CSTR)
end

function ConversionSteadyStateCSTR(; medium, reactionset, limiting_reactant, state, W, Q, conversion, name)
    cstr = SteadyStateCSTR(medium = medium, reactionset = reactionset, limiting_reactant = limiting_reactant, state = state, W = W, Q = Q, name = name)
    odesys = cstr.odesystem
    newsys = extend(System([odesys.X ~ conversion], t, [], []; name), odesys)
    cstr.odesystem = newsys
    return SteadyStateCSTR(cstr.medium, cstr.state, cstr.reactionset, cstr.phase, cstr.odesystem, :CSTR)
end

function FixedVolumeSteadyStateCSTR(; medium, reactionset, limiting_reactant, state, W, Q, volume, name) #Overwrites state volume and recalculate number of moles
    cstr = SteadyStateCSTR(medium = medium, reactionset = reactionset, limiting_reactant = limiting_reactant, state = state, W = W, Q = Q, name = name)
    odesys = cstr.odesystem
    @unpack V = odesys
    newsys = extend(System([V[1] ~ volume], t, [], []; name), odesys) #Fix overall volume to be V
    cstr.odesystem = newsys
    return SteadyStateCSTR(cstr.medium, cstr.state, cstr.reactionset, cstr.phase, cstr.odesystem, :CSTR)
end


@component function SteadyStateCSTRModel(;medium, reactions, limiting_reactant = nothing, state, W = 0.0, Q = nothing, phase = "liquid", name)

    @named CV = TwoPortControlVolume_SteadyState(medium = medium)
    @unpack U, Nᵢ, V, InPort, OutPort, ControlVolumeState, rₐ, rᵥ, Wₛ = CV

    vars = @variables begin
        cᵢ(t)[1:medium.Constants.Nc],         [description = "bulk concentrations"]   
        X(t),                                 [description = "Limiting reactant conversion"]            
    end


    eqs = [
        scalarize(cᵢ .~ Nᵢ/V[1])...
        Wₛ ~ W
        scalarize(rₐ[:, 2:end] .~ 0.0)...
        U ~ (OutPort.h[1] - ControlVolumeState.p/ControlVolumeState.ρ[1])*sum(collect(Nᵢ))
        ControlVolumeState.p ~ state.p
    ]

    limiting_index = findfirst(x -> x == limiting_reactant, medium.Constants.iupacName)

    if phase == "liquid"

        eq_reaction = [
            
            #Only liquid phase constraints
            scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...

            # Reaction kinetics
            scalarize(rᵥ[:, 2] .~ Rate(reactions, cᵢ, ControlVolumeState.T))...
            scalarize(rᵥ[:, 3] .~ 0.0)...

            #Conversion definition
            scalarize(X ~ (InPort.ṅ[2].*InPort.z[limiting_index, 2] .+ OutPort.ṅ[2].*OutPort.z[limiting_index, 2])./(InPort.ṅ[2].*InPort.z[limiting_index, 2] .+ 1e-8))
        
        ]

    elseif phase == "vapor"

        eq_reaction = [
            #Only vapor phase constraints
            scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...

            # Reaction kinetics
            scalarize(rᵥ[:, 2] .~ 0.0)...
            scalarize(rᵥ[:, 3] .~ sum(Rate.(reactions, cᵢ, ControlVolumeState.T)))...
            
            #Conversion definition
            scalarize(X ~ (InPort.ṅ[3].*InPort.z[limiting_index, 3] .+ OutPort.ṅ[3].*OutPort.z[limiting_index, 3])./(InPort.ṅ[3].*InPort.z[limiting_index, 3] .+ 1e-8))

        ]
        
    end

    pars = []

    return extend(System([eqs...;eq_reaction...], t, collect(Iterators.flatten(vars)), pars; name), CV)

end


export SteadyStateCSTR, DynamicCSTR, SteadyStateCSTRModel, DynamicCSTRModel, ConversionSteadyStateCSTR, FixedVolumeSteadyStateCSTR
