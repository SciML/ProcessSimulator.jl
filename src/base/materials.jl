abstract type AbstractMaterialSource end

struct Reaction{T<:Real}
    ν::Vector{T}                    # Stoichiometry
    r::Function                     # Reaction rate (defined as f(p,T,xᵢ)) in mol/s
    Δhᵣ::Function                   # Reaction enthalpy (defined as f(T)) in J/mol
    Reaction{T}(; ν, r, Δhᵣ) where {T} = new{T}(ν, r, Δhᵣ)
    Reaction(; ν::Vector{T}, r, Δhᵣ) where {T} = new{T}(ν, r, Δhᵣ)
end

struct MaterialSource{T<:Real, R<:Reaction} <: AbstractMaterialSource
    name::String                    # Name of the material source
    components::Vector{String}      # Component names
    N_c::Int                        # Number of components
    Mw::Vector{T}                   # Molar weight in kg/mol
    pressure::Function              # Pressure function (defined as f(ϱ,T,xᵢ;kwargs...)) in Pa
    molar_density::Function         # Molar density function (defined as f(p,T,xᵢ;kwargs...)) in mol/m³
    VT_internal_energy::Function    # Internal energy function (defined as f(ϱ,T,xᵢ;kwargs...)) in J/mol
    VT_enthalpy::Function           # Enthalpy function (defined as f(ϱ,T,xᵢ;kwargs...)) in J/mol
    VT_entropy::Function            # Entropy function (defined as f(ϱ,T,xᵢ;kwargs...)) in J/(mol K)
    tp_flash::Function              # Flash function (defined as f(p,T,xᵢ;kwargs...))
    reaction::Vector{R}             # Reaction struct
end

function MaterialSource(components::Union{String, Vector{String}}; kwargs...)
    components = components isa String ? [components] : components

    # Check for mandatory keyword arguments
    mandatory = [:Mw, :molar_density, :VT_enthalpy]
    [haskey(kwargs, k) || throw(ArgumentError("Missing keyword argument $k"))
     for k in mandatory]

    N_c = length(components)
    length(kwargs[:Mw]) == N_c ||
        throw(ArgumentError("Length of Mw must be equal to the number of components"))
    name = haskey(kwargs, :name) ? kwargs[:name] : join(components, "_")

    f_NA(field) = error("Function $field not defined in MaterialSource")

    Mw_vec = kwargs[:Mw] isa Number ? [kwargs[:Mw]] : collect(kwargs[:Mw])
    T = eltype(Mw_vec)
    reactions = get(kwargs, :reactions, Reaction{T}[])

    # Determine the reaction type - use Reaction{T} as default for empty vectors
    R = isempty(reactions) ? Reaction{T} : eltype(reactions)

    MaterialSource{T, R}(
        name,
        components,
        N_c,
        Mw_vec,
        get(kwargs, :pressure, (a, T, n; kws...) -> f_NA(:pressure)),
        kwargs[:molar_density],
        get(kwargs, :VT_internal_energy, (a, T, n; kws...) -> f_NA(:VT_internal_energy)),
        kwargs[:VT_enthalpy],
        get(kwargs, :VT_entropy, (a, T, n; kws...) -> f_NA(:VT_entropy)),
        get(kwargs, :tp_flash, (a, T, n; kws...) -> f_NA(:tp_flash)),
        reactions
    )
end
