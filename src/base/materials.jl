abstract type AbstractMaterialSource end

struct MaterialSource <: AbstractMaterialSource
    name::String                    # Name of the material source
    components::Vector{String}      # Component names
    Mw::Vector{Float64}             # Molar weight in kg/mol
    pressure::Function              # Pressure function (defined as f(ϱ,T,nᵢ;kwargs...))
    molar_density::Function         # Molar density function (defined as f(p,T,nᵢ;kwargs...))
    VT_internal_energy::Function    # Internal energy function (defined as f(ϱ,T,nᵢ;kwargs...))
    VT_enthalpy::Function           # Enthalpy function (defined as f(ϱ,T,nᵢ;kwargs...))
    VT_entropy::Function            # Entropy function (defined as f(ϱ,T,nᵢ;kwargs...))
    tp_flash::Function              # Flash function (defined as f(p,T,nᵢ;kwargs...))
end