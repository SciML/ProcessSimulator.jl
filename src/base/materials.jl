abstract type AbstractFluidMedium end

abstract type AbstractEoSBased <: AbstractFluidMedium end

abstract type AbstractSinkSource end

abstract type AbstractReaction <: AbstractSinkSource end

const R = 8.31446261815324 # J/(mol K)


function XtoMolar(flowrate, medium, state, flowbasis)
    if flowbasis == :molar
        return flowrate
    elseif flowbasis == :mass
        return MasstoMolar(flowrate, medium, state)
    elseif flowbasis == :volume
        return VolumetoMolar(flowrate, medium, state)
    else
        error("Invalid flow basis: $flowbasis")
    end
end

function MasstoMolar(flowrate, medium::M, state::S) where {M <: AbstractFluidMedium}
    M̄ = sum(state.z[:, 1] .* medium.FluidConstants.molarMass) #kg/mol
    return flowrate / M̄ 
end

function VolumetoMolar(flowrate, medium::M, state::S) where {M <: AbstractFluidMedium}
    return flowrate*state.ρ[1] 
end

struct PowerLawReaction{T <: Real, Arr <: AbstractArray{T}} <: AbstractReaction
    species::Array{String}        # Reactants and Products names
    ν::Arr                        # Stoichiometry
    n::Arr                        # Reaction order
    A::T                          # Arrhenius constant
    Eₐ::T                         # Activation energy
end


function _Rate(SinkSource::PowerLawReaction, cᵢ, T)
    A, Eₐ, n, ν = SinkSource.A, SinkSource.Eₐ, SinkSource.n, SinkSource.ν
    r = A * exp(-Eₐ / (R * T)) * prod(cᵢ[i]^n[i] for i in eachindex(cᵢ))
    return r.*ν
end

Rate(SinkSource, cᵢ, T) = _Rate(SinkSource, cᵢ, T)

Broadcast.broadcasted(::typeof(Rate), reactions, cᵢ, T) = broadcast(_Rate, reactions, Ref(cᵢ), T)

struct LDFAdsorption{K <: AbstractArray} <: AbstractSinkSource
    k::K                            # Mass transfer coefficient in 1/s.
end

struct BasicFluidConstants{S <: Union{AbstractString, Nothing}, M <: Union{Real, AbstractVector{<: Real}}, N <: Int, V <: Union{Nothing, AbstractVector{<: AbstractString}}}
    iupacName::S           # "Complete IUPAC name (or common name, if non-existent)";
    casRegistryNumber::S   # "Chemical abstracts sequencing number (if it exists)";
    chemicalFormula::S     # "Chemical formula";
    structureFormula::S    # "Chemical structure formula";
    molarMass::M           # "Molar mass";
    Nc::N                  # "Number of components";
    nphases::N             # "Maximum number of phases";
    phaseNames::V          # "Label of the phases"
end

##That should be part of the PropertyModels package


function BasicFluidConstants(molarMass::M) where M <: AbstractVector{<: Real}
    Nc = length(molarMass)
    return BasicFluidConstants(nothing, nothing, nothing, nothing, molarMass, Nc, 3, ["overall", "liquid", "vapor"])
end

struct EosBasedGuesses{M <: Any, V <: Real, D <: AbstractArray{V}, F <: AbstractArray{V}}
    EoSModel::M
    p::V #Pressure
    T::V #Temperature
    ρ::D #Molar density per phase
    x::F #Mole fraction per phase
    h::D #Molar enthalpy per phase
    pᵇᵈ::D #Bubble and dew pressure
end

function EosBasedGuesses(EoSModel::M, p::V, T::V, z::D) where {M <: Any, V <: Real, D <: AbstractArray{ <: Real}}

    ## Bubble and dew pressure
    pᵇ = bubble_pressure(EoSModel, T, z)[1]
    pᵈ = dew_pressure(EoSModel, T, z)[1]
    pᵇᵈ = [pᵇ, pᵈ]

    Nc = length(z)
    nᵢⱼ = zeros(V, Nc, 2)
    x = zeros(V, Nc, 3)
    ρ = zeros(V, 3)
    h = zeros(V, 3)

    if p ≤ pᵇ && p ≥ pᵈ
        sol = TP_flash(EoSModel, p, T, z)
        ϕ = sol[1]
        x .= sol[2]
        ρₗ = PT_molar_density(EoSModel, p, T, x[:, 1], phase = "liquid") #Assumes only two phases
        ρᵥ = PT_molar_density(EoSModel, p, T, x[:, 2], phase = "vapor")  
        ρₒᵥ = 1.0/(ϕ[1]/ρₗ + ϕ[2]/ρᵥ)
        ρ .= [ρₒᵥ, ρₗ, ρᵥ]

        ## Enthalpy
        hₗ = ρT_enthalpy(EoSModel, ρ[2], T, x[:, 2])
        hᵥ = ρT_enthalpy(EoSModel, ρ[3], T, x[:, 3])
        hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
        h .= [hₒᵥ, hₗ, hᵥ]

    elseif p > pᵇ
        nᵢⱼ[:, 1] .= z
        x[:, 1:2] .= z
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
        println(ϕ)
        ρ[2] = PT_molar_density(EoSModel, p, T, z, phase = "liquid")
        ρ[3] = PT_molar_density(EoSModel, p, T, z, phase = "vapor")
        ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
        hₗ = ρT_enthalpy(EoSModel, ρ[2], T, z)
        hᵥ = ρT_enthalpy(EoSModel, ρ[3], T, z)
        hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
        h .= [hₒᵥ, hₗ, hᵥ]

    elseif p < pᵈ
        nᵢⱼ[:, 2] .= z
        x[:, [1, 3]] .= z
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
        ρ[2] = PT_molar_density(EoSModel, p, T, z, phase = "liquid")
        ρ[3] = PT_molar_density(EoSModel, p, T, z, phase = "vapor")
        ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
        hₗ = ρT_enthalpy(EoSModel, ρ[2], T, z)
        hᵥ = ρT_enthalpy(EoSModel, ρ[3], T, z)
        hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
        h .= [hₒᵥ, hₗ, hᵥ]
    end


    return EosBasedGuesses(EoSModel, p, T, ρ, x, h, pᵇᵈ)
end

struct EoSBased{F <: BasicFluidConstants, E <: Any, G <: EosBasedGuesses} <: AbstractEoSBased
    Constants::F
    EoSModel::E
    Guesses::G
end




