abstract type AbstractSinkSource end

abstract type AbstractReaction <: AbstractSinkSource end

struct PowerLawReaction{T <: Real, Arr <: AbstractArray{T}} <: AbstractReaction
    species::Array{String}        # Reactants and Products names
    ν::Arr                        # Stoichiometry
    n::Arr                        # Reaction order
    A::T                          # Arrhenius constant
    Eₐ::T                         # Activation energy
end

"""
    PowerLawReaction(; stoichiometry::Dict, order::Dict, A, Eₐ, components=nothing)

Construct a PowerLawReaction using dictionaries for stoichiometry and reaction order.

# Arguments
- `stoichiometry::Dict{String, <:Real}`: Dictionary mapping component names to stoichiometric coefficients
  - Negative for reactants, positive for products
  - Example: `Dict("A" => -1, "B" => -1, "C" => 2)` for A + B → 2C
- `order::Dict{String, <:Real}`: Dictionary mapping component names to reaction orders
  - Example: `Dict("A" => 1.0, "B" => 1.0)` for first order in A and B
- `A::Real`: Arrhenius pre-exponential factor (units depend on overall order)
- `Eₐ::Real`: Activation energy (J/mol)
- `components::Union{Nothing, Vector{String}}`: Optional ordered list of component names
  - If provided, the reaction arrays will follow this exact order
  - If not provided, species are sorted alphabetically

# Examples
```julia
# With automatic ordering (alphabetical)
reaction = PowerLawReaction(
    stoichiometry = Dict("A" => -1.0, "B" => -1.0, "C" => 2.0),
    order = Dict("A" => 1.0, "B" => 1.0),
    A = 1e10,
    Eₐ = 50000.0
)

# With specified component order (matches medium definition)
components = ["ethylene oxide", "water", "ethylene glycol"]
reaction = PowerLawReaction(
    stoichiometry = Dict("ethylene oxide" => -1.0, "water" => -1.0, "ethylene glycol" => 1.0),
    order = Dict("ethylene oxide" => 1.0, "water" => 1.0),
    A = 1e10,
    Eₐ = 50000.0,
    components = components
)
```
"""
function PowerLawReaction(; stoichiometry::Dict{String, <:Real}, order::Dict{String, <:Real}, A::Real, Eₐ::Real, components::Union{Nothing, Vector{String}}=nothing)
    # Determine species order
    if components !== nothing
        # Use the provided component order
        species = components

        # Verify all species in stoichiometry/order are in components
        all_rxn_species = unique([keys(stoichiometry)..., keys(order)...])
        for sp in all_rxn_species
            if sp ∉ species
                error("Species '$sp' in reaction is not in the provided components list: $components")
            end
        end
    else
        # Get all unique species from both dictionaries and sort alphabetically
        all_species = unique([keys(stoichiometry)..., keys(order)...])
        species = sort(collect(all_species))
    end

    # Build stoichiometry array (default to 0 if not specified)
    ν = [get(stoichiometry, sp, 0.0) for sp in species]

    # Build reaction order array (default to 0 if not specified)
    n = [get(order, sp, 0.0) for sp in species]

    # Promote to common type
    T = promote_type(eltype(ν), eltype(n), typeof(A), typeof(Eₐ))
    ν_arr = convert(Vector{T}, ν)
    n_arr = convert(Vector{T}, n)
    A_val = convert(T, A)
    Eₐ_val = convert(T, Eₐ)

    return PowerLawReaction(species, ν_arr, n_arr, A_val, Eₐ_val)
end


function _Rate(SinkSource::PowerLawReaction, cᵢ, T)
    A, Eₐ, n, ν = SinkSource.A, SinkSource.Eₐ, SinkSource.n, SinkSource.ν
    r = A * exp(-Eₐ / (R * T)) * prod(cᵢ[i]^n[i] for i in eachindex(cᵢ))
    return r.*ν
end

Rate(SinkSource, cᵢ, T) = _Rate(SinkSource, cᵢ, T)

#= Broadcast.broadcasted(::typeof(Rate), reactions, cᵢ, T) = broadcast(_Rate, reactions, Ref(cᵢ), T) =#


# ============================================================================
# Reaction Network Composition
# ============================================================================

"""
    AbstractReactionSet

Abstract type for collections of reactions forming a reaction network.
"""
abstract type AbstractReactionSet end


"""
    PowerLawReactionSet{T <: Real}

A collection of PowerLawReactions forming a reaction network.

# Fields
- `species::Vector{String}`: All unique species in the network
- `reactions::Vector{PowerLawReaction}`: Individual reactions
- `ν::Matrix{T}`: Stoichiometric coefficient matrix (Nc × Nr)
  - Rows: species, Columns: reactions
  - ν[i,j] = stoichiometric coefficient of species i in reaction j
- `n::Matrix{T}`: Reaction order matrix (Nc × Nr)
  - Rows: species, Columns: reactions
  - n[i,j] = reaction order of species i in reaction j
- `A::Vector{T}`: Pre-exponential factors for each reaction
- `Eₐ::Vector{T}`: Activation energies for each reaction (J/mol)
"""
struct PowerLawReactionSet{T <: Real} <: AbstractReactionSet
    species::Vector{String}
    reactions::Vector{<:PowerLawReaction}
    ν::Matrix{T}   # Stoichiometric matrix (Nc × Nr)
    n::Matrix{T}   # Reaction order matrix (Nc × Nr)
    A::Vector{T}   # Pre-exponential factors
    Eₐ::Vector{T}  # Activation energies

    # Inner constructor to allow flexible reaction vector types
    function PowerLawReactionSet(species::Vector{String},
                                  reactions::Vector{<:PowerLawReaction},
                                  ν::Matrix{T},
                                  n::Matrix{T},
                                  A::Vector{T},
                                  Eₐ::Vector{T}) where T <: Real
        new{T}(species, reactions, ν, n, A, Eₐ)
    end
end

"""
    _Rate(network::PowerLawReactionSet, cᵢ, T)

Calculate net production rates for all species in a reaction network.

# Arguments
- `network::PowerLawReactionSet`: Reaction network with Nc species and Nr reactions
- `cᵢ::AbstractVector`: Concentration vector (length Nc)
- `T::Real`: Temperature (K)

# Returns
- `Vector`: Net production rate for each species (length Nc)
  - rᵢ = Σⱼ νᵢⱼ * rⱼ where rⱼ is the rate of reaction j

# Math
For each reaction j:
  rⱼ = Aⱼ * exp(-Eₐⱼ/(R*T)) * ∏ᵢ cᵢ^nᵢⱼ

Net rate for species i:
  rᵢ = Σⱼ νᵢⱼ * rⱼ

In matrix form:
  r = ν * r_rxn
where r_rxn is the vector of reaction rates
"""
function _Rate(network::RE, cᵢ, T) where RE <: PowerLawReactionSet
    Nc = length(network.species)  # Number of components
    Nr = length(network.reactions)  # Number of reactions

    # Calculate individual reaction rates
    r_rxn = zeros(eltype(cᵢ), Nr)

    for j in 1:Nr
        # Rate constant: k_j = A_j * exp(-Eₐ_j / (R*T))
        k_j = network.A[j] * exp(-network.Eₐ[j] / (R * T))

        # Concentration term: ∏ᵢ cᵢ^nᵢⱼ
        conc_term = one(eltype(cᵢ))
        for i in 1:Nc
            if network.n[i, j] != 0
                conc_term *= cᵢ[i]^network.n[i, j]
            end
        end

        # Reaction rate: r_j = k_j * ∏ᵢ cᵢ^nᵢⱼ
        r_rxn[j] = k_j * conc_term
    end

    # Net production rate for each species: rᵢ = Σⱼ νᵢⱼ * r_j
    # This is a matrix-vector multiplication: r = ν * r_rxn
    r_net = network.ν * r_rxn

    return r_net
end

Rate(network::PowerLawReactionSet, cᵢ, T) = _Rate(network, cᵢ, T)




"""
    PowerLawReactionSet(reactions::Vector{PowerLawReaction})

Compose multiple PowerLawReactions into a reaction network.

# Arguments
- `reactions::Vector{PowerLawReaction}`: Vector of individual reactions

# Returns
- `PowerLawReactionSet`: Unified reaction network with matrices

# Examples
```julia
# Define individual reactions
rxn1 = PowerLawReaction(
    stoichiometry = Dict("A" => -1.0, "B" => 1.0),
    order = Dict("A" => 1.0),
    A = 1e10, Eₐ = 50000.0
)

rxn2 = PowerLawReaction(
    stoichiometry = Dict("B" => -1.0, "C" => 1.0),
    order = Dict("B" => 1.0),
    A = 5e8, Eₐ = 60000.0
)

# Compose into network
network = PowerLawReactionSet([rxn1, rxn2])

# Access matrices
network.ν  # [Nc × Nr] stoichiometric matrix
network.n  # [Nc × Nr] reaction order matrix
```
"""
function PowerLawReactionSet(reactions::Vector{<:PowerLawReaction})
    if isempty(reactions)
        error("Cannot create empty reaction set")
    end

    # Collect all unique species across all reactions
    all_species = String[]
    for rxn in reactions
        append!(all_species, rxn.species)
    end
    species = sort(unique(all_species))

    Nc = length(species)  # Number of components
    Nr = length(reactions)  # Number of reactions

    # Build species index map
    species_idx = Dict(sp => i for (i, sp) in enumerate(species))

    # Determine output type from reactions
    T = promote_type([typeof(rxn.A) for rxn in reactions]...,
                     [typeof(rxn.Eₐ) for rxn in reactions]...,
                     [eltype(rxn.ν) for rxn in reactions]...,
                     [eltype(rxn.n) for rxn in reactions]...)

    # Initialize matrices
    ν_matrix = zeros(T, Nc, Nr)
    n_matrix = zeros(T, Nc, Nr)
    A_vec = zeros(T, Nr)
    Eₐ_vec = zeros(T, Nr)

    # Fill matrices from each reaction
    for (j, rxn) in enumerate(reactions)
        # Store kinetic parameters
        A_vec[j] = rxn.A
        Eₐ_vec[j] = rxn.Eₐ

        # Map species to global indices and fill matrices
        for (k, sp) in enumerate(rxn.species)
            i = species_idx[sp]
            ν_matrix[i, j] = rxn.ν[k]
            n_matrix[i, j] = rxn.n[k]
        end
    end

    return PowerLawReactionSet(species, reactions, ν_matrix, n_matrix, A_vec, Eₐ_vec)
end

"""
    PowerLawReactionSet(rxn1::PowerLawReaction, rxn2::PowerLawReaction, rxns::PowerLawReaction...)

Compose multiple PowerLawReactions using varargs.

# Examples
```julia
network = PowerLawReactionSet(rxn1, rxn2, rxn3)
```
"""
function PowerLawReactionSet(rxn1::PowerLawReaction, rxn2::PowerLawReaction, rxns::PowerLawReaction...)
    all_rxns = [rxn1, rxn2, rxns...]
    return PowerLawReactionSet(all_rxns)
end

"""
    Base.:+(rxn1::PowerLawReaction, rxn2::PowerLawReaction)

Compose two reactions using the + operator to form a reaction network.

# Examples
```julia
network = rxn1 + rxn2 + rxn3
```
"""
function Base.:+(rxn1::PowerLawReaction, rxn2::PowerLawReaction)
    return PowerLawReactionSet([rxn1, rxn2])
end

function Base.:+(set::PowerLawReactionSet, rxn::PowerLawReaction)
    return PowerLawReactionSet([set.reactions..., rxn])
end

export AbstractReaction, AbstractReactionSet, PowerLawReaction, PowerLawReactionSet, Rate