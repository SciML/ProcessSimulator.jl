abstract type AbstractMaterialSource end

struct MaterialSource <: AbstractMaterialSource
    name::String                    # Name of the material source
    components::Vector{String}      # Component names
    N_c::Int                        # Number of components
    Mw::Vector{Float64}             # Molar weight in kg/mol
    pressure::Function              # Pressure function (defined as f(ϱ,T,xᵢ;kwargs...))
    molar_density::Function         # Molar density function (defined as f(p,T,xᵢ;kwargs...))
    VT_internal_energy::Function    # Internal energy function (defined as f(ϱ,T,xᵢ;kwargs...))
    VT_enthalpy::Function           # Enthalpy function (defined as f(ϱ,T,xᵢ;kwargs...))
    VT_entropy::Function            # Entropy function (defined as f(ϱ,T,xᵢ;kwargs...))
    tp_flash::Function              # Flash function (defined as f(p,T,xᵢ;kwargs...))
end

function MaterialSource(components::Union{String,Vector{String}}; kwargs...)
    components = components isa String ? [components] : components
    
    # Check for mandatory keyword arguments
    mandatory = [:Mw, :molar_density, :VT_enthalpy]
    [haskey(kwargs, k) || throw(ArgumentError("Missing keyword argument $k")) for k in mandatory]

    N_c = length(components)
    length(kwargs[:Mw]) == N_c || throw(ArgumentError("Length of Mw must be equal to the number of components"))
    name = haskey(kwargs, :name) ? kwargs[:name] : join(components, "_")

    f_NA(field) = error("Function $field not defined in MaterialSource")

    MaterialSource(
        name,
        components,
        N_c,
        kwargs[:Mw] isa Number ? [kwargs[:Mw]] : kwargs[:Mw],
        get(kwargs, :pressure, (a,T,n;kws...) -> f_NA(:pressure)),
        kwargs[:molar_density],
        get(kwargs, :VT_internal_energy, (a,T,n;kws...) -> f_NA(:VT_internal_energy)),
        kwargs[:VT_enthalpy],
        get(kwargs, :VT_entropy, (a,T,n;kws...) -> f_NA(:VT_entropy)),
        get(kwargs, :tp_flash, (a,T,n;kws...) -> f_NA(:tp_flash))
    )
end