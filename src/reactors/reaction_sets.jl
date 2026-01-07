"""
    Reaction sets and kinetics utilities for ProcessSimulator.jl

This module provides data structures and utilities for managing sets of
chemical reactions and computing reaction kinetics.
"""

"""
    ReactionSet

A collection of related reactions that occur together in a reactor.
Provides convenient access to stoichiometry, kinetics, and thermodynamic
properties for multiple reactions.

# Fields
- `reactions::Vector{Reaction}`: Vector of individual reactions
- `name::String`: Name/identifier for the reaction set
- `components::Vector{String}`: Names of all components involved
"""
struct ReactionSet
    reactions::Vector{Reaction}
    name::String
    components::Vector{String}

    function ReactionSet(reactions::Vector{Reaction}, components::Vector{String};
                         name::String = "ReactionSet")
        # Validate that all reactions have consistent stoichiometry length
        N_c = length(components)
        for (i, r) in enumerate(reactions)
            if length(r.ν) != N_c
                throw(ArgumentError(
                    "Reaction $i has stoichiometry of length $(length(r.ν)), " *
                    "expected $N_c (number of components)"))
            end
        end
        new(reactions, name, components)
    end
end

"""
    ReactionSet(; reactions, components, name="ReactionSet")

Keyword constructor for ReactionSet.
"""
function ReactionSet(; reactions::Vector{Reaction}, components::Vector{String},
                     name::String = "ReactionSet")
    ReactionSet(reactions, components; name = name)
end

"""
    nreactions(rs::ReactionSet)

Return the number of reactions in the set.
"""
nreactions(rs::ReactionSet) = length(rs.reactions)

"""
    ncomponents(rs::ReactionSet)

Return the number of components in the reaction set.
"""
ncomponents(rs::ReactionSet) = length(rs.components)

"""
    stoichiometry_matrix(rs::ReactionSet)

Return the stoichiometric matrix (N_reactions × N_components).
Each row corresponds to a reaction, each column to a component.
Negative values indicate reactants, positive values indicate products.
"""
function stoichiometry_matrix(rs::ReactionSet)
    N_r = nreactions(rs)
    N_c = ncomponents(rs)
    ν = zeros(N_r, N_c)
    for (i, r) in enumerate(rs.reactions)
        ν[i, :] .= r.ν
    end
    return ν
end

"""
    reaction_rates(rs::ReactionSet, p, T, x)

Compute the reaction rates for all reactions in the set.
Returns a vector of rates in mol/s.

# Arguments
- `rs::ReactionSet`: The reaction set
- `p`: Pressure (Pa)
- `T`: Temperature (K)
- `x`: Mole fractions (vector of length N_components)
"""
function reaction_rates(rs::ReactionSet, p, T, x)
    return [r.r(p, T, x) for r in rs.reactions]
end

"""
    species_production_rates(rs::ReactionSet, p, T, x)

Compute the net production rate of each species due to all reactions.
Returns a vector of rates in mol/s for each component.

# Arguments
- `rs::ReactionSet`: The reaction set
- `p`: Pressure (Pa)
- `T`: Temperature (K)
- `x`: Mole fractions (vector of length N_components)
"""
function species_production_rates(rs::ReactionSet, p, T, x)
    rates = reaction_rates(rs, p, T, x)
    ν = stoichiometry_matrix(rs)
    return ν' * rates  # N_c × 1 vector
end

"""
    reaction_enthalpies(rs::ReactionSet, T)

Compute the enthalpy of reaction for each reaction at temperature T.
Returns a vector of enthalpies in J/mol.
"""
function reaction_enthalpies(rs::ReactionSet, T)
    return [r.Δhᵣ(T) for r in rs.reactions]
end

"""
    total_heat_of_reaction(rs::ReactionSet, p, T, x)

Compute the total heat released/absorbed by all reactions.
Returns heat rate in J/s.

# Arguments
- `rs::ReactionSet`: The reaction set
- `p`: Pressure (Pa)
- `T`: Temperature (K)
- `x`: Mole fractions (vector of length N_components)
"""
function total_heat_of_reaction(rs::ReactionSet, p, T, x)
    rates = reaction_rates(rs, p, T, x)
    enthalpies = reaction_enthalpies(rs, T)
    return sum(rates .* enthalpies)
end

"""
    ArrheniusKinetics

Helper struct for Arrhenius-type reaction kinetics:
    k = A * exp(-Ea / (R * T))

# Fields
- `A`: Pre-exponential factor (units depend on reaction order)
- `Ea`: Activation energy (J/mol)
"""
struct ArrheniusKinetics
    A::Float64   # Pre-exponential factor
    Ea::Float64  # Activation energy (J/mol)
end

"""
    rate_constant(k::ArrheniusKinetics, T)

Compute the rate constant at temperature T (K) using Arrhenius equation.
"""
function rate_constant(k::ArrheniusKinetics, T)
    R = 8.314  # J/(mol·K)
    return k.A * exp(-k.Ea / (R * T))
end

"""
    PowerLawKinetics

Helper struct for power-law rate expressions:
    r = k(T) * prod(C_i^order_i)

# Fields
- `arrhenius::ArrheniusKinetics`: Temperature-dependent rate constant
- `orders::Vector{Float64}`: Reaction orders for each component
"""
struct PowerLawKinetics
    arrhenius::ArrheniusKinetics
    orders::Vector{Float64}
end

"""
    reaction_rate(k::PowerLawKinetics, T, C)

Compute the reaction rate given temperature T and concentrations C.
"""
function reaction_rate(k::PowerLawKinetics, T, C)
    k_val = rate_constant(k.arrhenius, T)
    return k_val * prod(C .^ k.orders)
end

export ReactionSet, nreactions, ncomponents, stoichiometry_matrix
export reaction_rates, species_production_rates, reaction_enthalpies, total_heat_of_reaction
export ArrheniusKinetics, PowerLawKinetics, rate_constant, reaction_rate
