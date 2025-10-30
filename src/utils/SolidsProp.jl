abstract type AbstractSolidMedium end
abstract type AbstractSolidEoSModel end

struct BasicMaterialConstants{S <: Union{AbstractString, Nothing}, N <: Int, V <: Union{Nothing, AbstractVector{<: AbstractString}},
     T <: Union{Real, Nothing}, T2 <: Union{Real, Nothing}}
    adsorbent_name::S       # "Complete IUPAC name (or common name, if non-existent)";
    solute_names::V         # "Names of the solutes";
    Nc::N                   # "Number of components";
    nphases::N              # "Maximum number of phases";
    phase_names::V          # "Label of the phases"
    particle_size::T        # "Particle size of the adsorbent";
    pore_size::T2            # "Pore size of the adsorbent"; 
end

# This is the single, correct constructor. It replaces the two duplicate ones.
function BasicMaterialConstants(adsorbent_name, solute_names, phase_names, particle_size)
    Nc = length(solute_names)
    nphases = length(phase_names)
    # This calls the inner constructor, providing `nothing` for the optional `pore_size`.
    return BasicMaterialConstants(adsorbent_name, solute_names, Nc, nphases, phase_names, particle_size, nothing)
end

struct SolidEoSModel{C <: Union{AbstractVector{<: Real}, AbstractVector{<: Nothing}}, T <: Real} <: AbstractSolidEoSModel
    ρ0::T
    ρ_T0::T
    coeffs_ρ::C
    h0::T
    h_T0::T
    coeffs_h::C
end

function SolidEoSModel(;ρ0, ρ_T0, coeffs_ρ, h0, h_T0, coeffs_h)
    return SolidEoSModel(ρ0, ρ_T0, coeffs_ρ, h0, h_T0, coeffs_h)
end


struct Adsorbent{I <: Any,
     E <: Any, T <: Any, C <: Union{BasicMaterialConstants, Nothing}} <: AbstractSolidMedium
    isotherm::I
    EoSModel::E
    TransportModel::T
    Constants::C
    Guesses
end

struct SolidGuesses{S <: Real, X<:AbstractArray{<: S}}
    p::S
    T::S
    x::X
end

function Adsorbent(;adsorbent_name, components, particle_size, isotherm, EoSModel, transport_model)
    if length(isotherm.isotherms) != length(components)
         return "Isotherm and solute names must have the same length."

    else
        constants = BasicMaterialConstants(adsorbent_name, components, ["overall", "particle"], particle_size)
        guess = SolidGuesses(101325.0, 273.15, ones(constants.Nc, constants.nphases)./constants.Nc)
        return Adsorbent(isotherm, EoSModel, transport_model, constants, guess)
    end
end


function PT_molar_density(EoSModel::M, p, T, x; phase = "particle") where M <: SolidEoSModel
    # Assumes fourth order polynomial dependency with T and that the adsorbed fluid does not affect the density.
    coeffs = EoSModel.coeffs_ρ
    Tref = EoSModel.ρ_T0
    ρ_ref = EoSModel.ρ0
    ΔT = T - Tref
    return ρ_ref + coeffs[1]*ΔT+ coeffs[2]*ΔT^2 + coeffs[3]*ΔT^3 + coeffs[4]*ΔT^4 
end


function ρT_enthalpy(EoSModel::M, ρ, T, x; phase = "particle") where M <: SolidEoSModel
    # Assumes fourth order polynomial dependency with T and that the adsorbed fluid does not affect the enthalpy (mostly ok).
    coeffs = EoSModel.coeffs_h
    Tref = EoSModel.h_T0
    h0 = EoSModel.h0
    ΔT = T - Tref
    return h0 + coeffs[1]*ΔT + coeffs[2]*ΔT^2 + coeffs[3]*ΔT^3 + coeffs[4]*ΔT^4 
end

function area_per_volume(solid::A, sphericity = 1.0) where A <: AbstractSolidMedium
    # Returns the area per volume of the solid medium.
    # Assumes spherical particles.
    rₚ = solid.Constants.particle_size
    return 3.0 ./ rₚ*sphericity
end

function mass_transfer_coefficient(model::M, solid::A, T, x) where {M <: HomogeneousDiffusivityCoeff, A <: AbstractSolidMedium}
    rₚ² = solid.Constants.particle_size^2
    Dₕ = model.Dh
    return 15.0 .* Dₕ ./ rₚ²
end

function mass_transfer_coefficient(solid::A, T = 273.15, x = ones(adsorbent.Constants.Nc)/Nc) where A <: AbstractSolidMedium
    return mass_transfer_coefficient(solid.TransportModel.MassTransferModel, solid, T, x)
end

#write specific dispatches here
function heat_transfer_coefficient(model::M, solid::A, T, x) where {M <: ConstantHeatTransferCoeff, A <: AbstractSolidMedium}
    return model.k
end

function heat_transfer_coefficient(solid::A, T = 273.15, x = ones(adsorbent.Constants.Nc)/Nc) where A <: AbstractSolidMedium
    return heat_transfer_coefficient(solid.TransportModel.HeatTransferModel, solid, T, x)
end

export SolidEoSModel, Adsorbent


