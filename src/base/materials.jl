"""
    Reaction(; ν, r, Δhᵣ)

Describe one reaction used by a [`MaterialSource`](@ref).

# Arguments

None.

# Keywords

  - `ν`: Stoichiometric coefficients, ordered consistently with the material source's
    components.
  - `r`: Reaction-rate function `r(p, T, xᵢ)` in mol/s, where `p` is pressure,
    `T` is temperature, and `xᵢ` is the component mole-fraction vector.
  - `Δhᵣ`: Reaction-enthalpy function `Δhᵣ(T)` in J/mol.

# Fields

  - `ν`: Stoichiometric-coefficient vector.
  - `r`: Reaction-rate function.
  - `Δhᵣ`: Reaction-enthalpy function.

# Examples

```jldoctest
julia> reaction = Reaction(;
           ν = [-1.0, 1.0],
           r = (p, T, xᵢ) -> 0.1 * xᵢ[1],
           Δhᵣ = T -> -5.0e4,
       );

julia> reaction.ν
2-element Vector{Float64}:
 -1.0
  1.0
```
"""
struct Reaction
    ν::Vector{Float64}              # Stoichiometry
    r::Function                     # Reaction rate (defined as f(p,T,xᵢ)) in mol/s
    Δhᵣ::Function                   # Reaction enthalpy (defined as f(T)) in J/mol
    Reaction(; ν, r, Δhᵣ) = new(ν, r, Δhᵣ)
end

@public Reaction

"""
    MaterialSource(components; Mw, molar_density, VT_enthalpy,
        VT_internal_energy = ..., pressure = ..., VT_entropy = ...,
        tp_flash = ..., reactions = ..., name = ...)

Provide thermodynamic property functions and reaction data for ProcessSimulator
components. The supplied functions are evaluated symbolically when a component is
constructed, so they must accept symbolic inputs and return expressions compatible
with ModelingToolkit.

# Arguments

  - `components`: A component name or vector of component names. Its order defines
    the order of every composition and stoichiometry vector passed to this source.

# Keywords

  - `Mw`: Required molar masses in kg/mol, either one number for a single component
    or one value per component.
  - `molar_density`: Required function `molar_density(p, T, xᵢ; phase = ...)`
    returning molar density in mol/m^3.
  - `VT_enthalpy`: Required function `VT_enthalpy(ϱ, T, xᵢ)` returning molar
    enthalpy in J/mol.
  - `VT_internal_energy`: Function `VT_internal_energy(ϱ, T, xᵢ)` returning molar
    internal energy in J/mol. It throws when a component needs it and it was omitted.
  - `pressure`: Function `pressure(ϱ, T, xᵢ; phase = ...)` returning pressure in Pa.
  - `VT_entropy`: Function `VT_entropy(ϱ, T, xᵢ)` returning molar entropy in
    J/(mol K).
  - `tp_flash`: Function `tp_flash(p, T, xᵢ; kwargs...)` implementing a
    temperature-pressure flash calculation.
  - `reactions`: Vector of [`Reaction`](@ref) values. Defaults to an empty vector.
  - `name`: Source name. Defaults to the component names joined by `_`.

# Fields

  - `name`: Source name.
  - `components`: Ordered component names.
  - `N_c`: Number of components.
  - `Mw`: Ordered molar masses in kg/mol.
  - `pressure`, `molar_density`: Thermodynamic pressure and molar-density functions.
  - `VT_internal_energy`, `VT_enthalpy`, `VT_entropy`: Molar thermodynamic property
    functions in volume-temperature coordinates.
  - `tp_flash`: Temperature-pressure flash function.
  - `reaction`: Reaction descriptions.

# Examples

```jldoctest
julia> source = MaterialSource("helium";
           Mw = 0.004,
           molar_density = (p, T, xᵢ; kwargs...) -> p / (8.314 * T),
           VT_enthalpy = (ϱ, T, xᵢ) -> 20.0 * T,
           VT_internal_energy = (ϱ, T, xᵢ) -> 12.0 * T,
       );

julia> source.N_c
1
```
"""
struct MaterialSource
    name::String                    # Name of the material source
    components::Vector{String}      # Component names
    N_c::Int                        # Number of components
    Mw::Vector{Float64}             # Molar weight in kg/mol
    pressure::Function              # Pressure function (defined as f(ϱ,T,xᵢ;kwargs...)) in Pa
    molar_density::Function         # Molar density function (defined as f(p,T,xᵢ;kwargs...)) in mol/m³
    VT_internal_energy::Function    # Internal energy function (defined as f(ϱ,T,xᵢ;kwargs...)) in J/mol
    VT_enthalpy::Function           # Enthalpy function (defined as f(ϱ,T,xᵢ;kwargs...)) in J/mol
    VT_entropy::Function            # Entropy function (defined as f(ϱ,T,xᵢ;kwargs...)) in J/(mol K)
    tp_flash::Function              # Flash function (defined as f(p,T,xᵢ;kwargs...))
    reaction::Vector{Reaction}      # Reaction struct
end

"""
    MaterialSource(components; kwargs...)

Construct a [`MaterialSource`](@ref). See [`MaterialSource`](@ref) for the required
thermodynamic-function contract, keyword arguments, fields, and an example.
"""
function MaterialSource(components::Union{String, Vector{String}}; kwargs...)
    components = components isa String ? [components] : components

    # Check for mandatory keyword arguments
    mandatory = [:Mw, :molar_density, :VT_enthalpy]
    [
        haskey(kwargs, k) || throw(ArgumentError("Missing keyword argument $k"))
            for k in mandatory
    ]

    N_c = length(components)
    Mw = kwargs[:Mw] isa Number ? [kwargs[:Mw]] : kwargs[:Mw]
    length(Mw) == N_c ||
        throw(ArgumentError("Length of Mw must be equal to the number of components"))
    name = haskey(kwargs, :name) ? kwargs[:name] : join(components, "_")

    f_NA(field) = error("Function $field not defined in MaterialSource")

    return MaterialSource(
        name,
        components,
        N_c,
        Mw,
        get(kwargs, :pressure, (a, T, n; kws...) -> f_NA(:pressure)),
        kwargs[:molar_density],
        get(kwargs, :VT_internal_energy, (a, T, n; kws...) -> f_NA(:VT_internal_energy)),
        kwargs[:VT_enthalpy],
        get(kwargs, :VT_entropy, (a, T, n; kws...) -> f_NA(:VT_entropy)),
        get(kwargs, :tp_flash, (a, T, n; kws...) -> f_NA(:tp_flash)),
        get(kwargs, :reactions, Reaction[])
    )
end

@public MaterialSource
