abstract type AbstractFluidMedium end

abstract type AbstractEoSBased <: AbstractFluidMedium end

abstract type AbstractThermodynamicState end

mutable struct pTNVState{T <: Union{Real, Nothing}, K <: AbstractArray{<:T}} <: AbstractThermodynamicState
    p::T
    T::T
    N::K
    V
end

function pTNVState(p, T, N; base = :Pressure)
    if base == :Pressure
    return pTNVState(p, T, N, nothing)
    else
    return pTNVState(nothing, T, N, V)
    end
end

mutable struct pTzState{T <: Real, V <: AbstractArray{<:T}} <: AbstractThermodynamicState
    p::T
    T::T
    z::V
end

mutable struct rhoTzState{T <: Real, V <: AbstractArray{<:T}} <: AbstractThermodynamicState
    ρ::T
    T::T
    z::V
end

export pTzState, pTNVState

function XtoMolar(flowrate, medium, state, flowbasis)
    if flowbasis == :molar
        return flowrate
    elseif flowbasis == :mass
        return MasstoMolar(flowrate, medium, state.z[:, 1])
    elseif flowbasis == :volume
        return VolumetoMolar(flowrate, state.ρ[1])
    else
        error("Invalid flow basis: $flowbasis")
    end
end

function MasstoMolar(flowrate, medium, x)
    M̄ = sum(x .* medium.FluidConstants.molarMass) #kg/mol
    return flowrate / M̄ 
end

function VolumetoMolar(flowrate, density)
    return flowrate./density 
end

struct BasicFluidConstants{S <: Union{AbstractString, Nothing, AbstractVector{<: AbstractString}}, M <: Union{Nothing, Real, AbstractVector{<: Real}}, N <: Int, V <: Union{Nothing, AbstractVector{<: AbstractString}}}
    iupacName::S           # "Complete IUPAC name (or common name, if non-existent)";
    molarMass::M           # "Molar mass";
    Nc::N                  # "Number of components";
    nphases::N             # "Allowed phases";
    phase_names::V          # "Label of the phases"
end


struct EosBasedGuesses{M <: Any, V <: Real, D <: AbstractArray{V}, F <: AbstractArray{V}}
    EoSModel::M
    p::V #Pressure
    T::V #Temperature
    ρ::D #Molar density per phase
    x::F #Mole fraction per phase
    h::D #Molar enthalpy per phase
    ϕ #vaporized fraction
end


"""
    EoSBased{F<:BasicFluidConstants, E<:Any, G<:EosBasedGuesses, T<:TransportModel} <: AbstractEoSBased

Represents a fluid medium characterized by an equation of state (EoS) model.

# Fields
- `Constants::F`: Physical and chemical constants of the fluid components (molecular weights, names, etc.)
- `EoSModel::E`: Equation of state model used for thermodynamic calculations (e.g., Peng-Robinson, GERG-2008)
- `TransportModel::T`: Model for transport properties (viscosity, thermal conductivity, etc.)
- `Guesses::G`: Current state variable guesses (pressure, temperature, densities, compositions)

# Examples
```julia
# Create a methane-ethane mixture with Peng-Robinson EoS
components = ["methane", "ethane"]
eos = PR(components)  # Peng-Robinson EoS from Clapeyron.jl
medium = EoSBased(components, eos)

# Create with a specific initial state (room temperature, slightly elevated pressure)
state = pTzState(5e5, 298.15, [0.7, 0.3])  # 5 bar, 25°C, 70% methane
medium_with_state = EoSBased(components, eos, state)

# Create with custom transport models and state
mt_model = ConstantMassTransferCoeff([0.5, 0.6])  # Component-specific coefficients
ht_model = ConstantHeatTransferCoeff(15.0)  # W/(m²·K)
visc_model = ChapmanEnskogModel(components)
transport = TransportModel(mt_model, ht_model, visc_model)
medium_full = EoSBased(components, eos, transport, state)

# Access thermodynamic properties
ρ_liquid = medium.Guesses.ρ[2]  # Liquid phase density (mol/m³)
h_vapor = medium.Guesses.h[3]  # Vapor phase enthalpy (J/mol)
vapor_frac = medium.Guesses.ϕ[2]  # Vapor fraction

"""
mutable struct EoSBased{F <: BasicFluidConstants, E <: Any, G <: EosBasedGuesses, T <: TransportModel} <: AbstractEoSBased
    Constants::F
    EoSModel::E
    TransportModel::T
    Guesses::G
end

# --- Main Outer Constructor ---
# This is the only constructor that should create the EoSBased instance directly.
function EoSBased(components::S, eosmodel, transportmodel::T, state::ST) where {S <: AbstractVector{<: AbstractString}, T <: TransportModel, ST <: pTzState}
    constants = BasicFluidConstants(components)
    _p, _T, _z = state.p, state.T, state.z
    guesses = EosBasedGuesses(eosmodel, _p, _T, _z, Val(:Pressure))
    return EoSBased(constants, eosmodel, transportmodel, guesses)
end


# --- Convenience Constructors ---
# These constructors provide default values and then call the main constructor.

# Constructor without a specified state (uses a default state)
function EoSBased(components::S, eosmodel, transportmodel::T) where {S <: AbstractVector{<: AbstractString}, T <: TransportModel}
    default_state = pTzState(101325.0, 298.15, ones(length(components))/length(components))
    return EoSBased(components, eosmodel, transportmodel, default_state)
end

# Constructor without a specified transport model (creates a default one)
function EoSBased(components::S, eosmodel, state::ST) where {S <: AbstractVector{<: AbstractString}, ST <: pTzState}
    masstransfermodel = ConstantMassTransferCoeff(0.5*ones(length(components)))
    heat_transfer_model = ConstantHeatTransferCoeff(10.0)
    # Assuming you have a default viscosity model or can pass nothing
    transportmodel = TransportModel(masstransfermodel, heat_transfer_model, nothing)
    return EoSBased(components, eosmodel, transportmodel, state)
end

# Constructor accepting pTNVState - converts N to z (mole fractions)
function EoSBased(components::S, eosmodel, state::ST) where {S <: AbstractVector{<: AbstractString}, ST <: pTNVState}
    # Convert molar amounts to mole fractions
    N_total = sum(state.N)
    z = state.N ./ N_total

    # Create pTzState from pTNVState
    ptz_state = pTzState(state.p, state.T, z)

    # Call the pTzState constructor
    return EoSBased(components, eosmodel, ptz_state)
end

# Constructor with pTNVState and transport model
function EoSBased(components::S, eosmodel, transportmodel::T, state::ST) where {S <: AbstractVector{<: AbstractString}, T <: TransportModel, ST <: pTNVState}
    # Convert molar amounts to mole fractions
    N_total = sum(state.N)
    z = state.N ./ N_total

    # Create pTzState from pTNVState
    ptz_state = pTzState(state.p, state.T, z)

    # Call the main constructor
    return EoSBased(components, eosmodel, transportmodel, ptz_state)
end

# Constructor with only components and EoS model (uses all defaults)
function EoSBased(components::S, eosmodel) where {S <: AbstractVector{<: AbstractString}}
    default_state = pTzState(101325.0, 298.15, ones(length(components))/length(components))
    # This now calls the constructor that creates the default transport model
    return EoSBased(components, eosmodel, default_state)
end


# --- Keyword Argument Constructors ---
# These simply forward the keyword arguments to the appropriate positional constructor.

function EoSBased(;
    components,
    eosmodel,
    transportmodel = nothing,
    state = nothing
)
    # This single function handles all keyword-based calls.
    # It decides which positional constructor to call based on the provided arguments.
    
    if transportmodel !== nothing && state !== nothing
        # All arguments provided
        return EoSBased(components, eosmodel, transportmodel, state)
    elseif transportmodel !== nothing && state === nothing
        # transportmodel provided, but no state
        return EoSBased(components, eosmodel, transportmodel)
    elseif transportmodel === nothing && state !== nothing
        # state provided, but no transportmodel
        return EoSBased(components, eosmodel, state)
    else # Both transportmodel and state are nothing
        # Only components and eosmodel provided
        return EoSBased(components, eosmodel)
    end
end

##That should be part of the PropertyModels package

function BasicFluidConstants(iupacName)
    Nc = length(iupacName)
    return BasicFluidConstants(iupacName, nothing, Nc, 3, ["overall", "liquid", "vapor"])
end


function EosBasedGuesses(EoSModel::M, p::V, T::V, z::D, ::Val{:Pressure}) where {M <: Any, V <: Real, D <: AbstractArray{ <: Real}}

    sol = TP_flash(EoSModel, p, T, z)
    ϕ = sol[1]
    x = sol[2]
    ρ = zeros(V, 3)
    h = zeros(V, 3)

    ## Enthalpy
    hₗ = enthalpy(EoSModel, p, T, x[:, 2], phase = "liquid")
    hᵥ = enthalpy(EoSModel, p, T, x[:, 3], phase = "vapor")
    ρ[2] = PT_molar_density(EoSModel, p, T, x[:, 2], phase = "liquid")
    ρ[3] = PT_molar_density(EoSModel, p, T, x[:, 3], phase = "vapor")
    ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
    hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
    h .= [hₒᵥ, hₗ, hᵥ]

    return EosBasedGuesses(EoSModel, p, T, ρ, x, h, ϕ)
end

function EosBasedGuesses(EoSModel::M, V::K, T::K, N::D, ::Val{:Volume}) where {M <: Any, K <: Real, D <: AbstractArray{ <: Real}}
    return nothing
end

function resolve_guess!(medium, state)
    phase = "unknown"
    p, T, N, V = state.p, state.T, state.N, state.V
    z = N./sum(N)
    if isnothing(V)
        medium.Guesses = EosBasedGuesses(medium.EoSModel, p, T, z, Val(:Pressure))
        phase = ifelse(medium.Guesses.ϕ[2] ≈ 1.0, "vapor", "liquid")
        state.V = sum(N)./medium.Guesses.ρ[1]
    elseif isnothing(p)
        medium.Guesses = EosBasedGuesses(medium.EoSModel, V, T, N, Val(:Volume))
        phase = ifelse(medium.Guesses.ϕ[2] ≈ 1.0, "vapor", "liquid")
        state.p = medium.Guesses.p
    else
        phase = ifelse(medium.Guesses.ϕ[2] ≈ 1.0, "vapor", "liquid")
    end
    return medium, state, phase
end


#Physicochemical properties

function PT_molar_density(EoSModel, p, T, x; phase = "unknown")
    return NaN   
end

function TP_flash(EoSModel, p, T, x)
    return NaN
end

function is_stable(EoSModel, p, T, x)
    return true
end

function is_VT_stable(EoSModel, v, T, x) 
    return true
end

function flash_mol_fractions(EoSModel, p, T, x)
    return NaN
end

function flash_mol_fractions_liquid(EoSModel, p, T, x) 
         z = NaN
         return z
end

function flash_mol_fractions_vapor(EoSModel, p, T, x)
    z = flash_mol_fractions(EoSModel, p, T, x)[:, 3] 
    return z
end

function flash_vaporized_fraction(EoSModel, p, T, x)
         ϕ = TP_flash(EoSModel, p, T, x)[1] 
         return ϕ
end


function ρT_enthalpy(EoSModel, ρ, T, x)
    return NaN
end

function pT_enthalpy(EoSModel, p, T, x; phase = :unknown)
    return NaN
end

function ρT_internal_energy(EoSModel, ρ, T, x)
    return NaN
end

function molecular_weight(model, z::AbstractVector)
    return NaN
end


#Transport functions
function mass_transfer_coefficient(fluid::A, T = 273.15, x = ones(fluid.Constants.Nc)/Nc) where A <: AbstractFluidMedium
    return mass_transfer_coefficient(fluid.TransportModel.MassTransferModel, fluid, T, x)
end

function mass_transfer_coefficient(model::M, fluid::A, T, x) where {M <: ConstantMassTransferCoeff, A <: AbstractFluidMedium}
    return model.k
end

#write specific dispatches here
function heat_transfer_coefficient(model::M, fluid::A, T, x) where {M <: ConstantHeatTransferCoeff, A <: AbstractFluidMedium}
    return model.k
end

function heat_transfer_coefficient(fluid::A, T = 273.15, x = ones(fluid.Constants.Nc)/Nc) where A <: AbstractFluidMedium
    return heat_transfer_coefficient(fluid.TransportModel.HeatTransferModel, fluid, T, x)
end


function viscosity(model::M, p, T, z) where M <: AbstractFluidMedium
    return viscosity(model.TransportModel.ViscosityModel, p, T, z)
end

function viscosity(model, p, T, z)
    return NaN
end


export EoSBased
export is_stable, is_VT_stable, TP_flash, flash_mol_fractions, flash_mol_fractions_liquid, flash_mol_fractions_vapor, flash_vaporized_fraction, PT_molar_density, ρT_enthalpy, ρT_internal_energy, molecular_weight, pT_enthalpy
export mass_transfer_coefficient, heat_transfer_coefficient, viscosity




