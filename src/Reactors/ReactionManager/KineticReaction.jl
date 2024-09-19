# Definindo a estrutura
struct KineticReactionNetwork{T}
    Coef_Cr::Union{Vector{T}, Matrix{T}} # Stoichiometric coefficients of each component in each reaction (-)
    Do_r::Union{Vector{T}, Matrix{T}} # Forward order of the components ()
    substances_user::Vector{String} # Substances in the reaction network
    Nc::Int # Number of components in the reaction network
    Nri::Int # Number of reactions in the reaction network (r)
    Af_r::T # Arrhenius constant of each reaction (s⁻¹)
    Ef_r::T # Activation energy of each reaction (J.mol⁻¹)
    name::String
end

# Construtor para aceitar palavras-chave
function KineticReactionNetwork(; Coef_Cr::Union{Vector{T}, Matrix{T}}, Do_r::Union{Vector{T}, Matrix{T}},
        substances_user::Vector{String}, Af_r::T, Ef_r::T, name::String) where {T}

    Nc = length(substances_user)
    Nri = Coef_Cr isa Matrix ? size(Coef_Cr, 1) : length(Coef_Cr)

    
    return KineticReactionNetwork{T}(Coef_Cr, Do_r, substances_user, Nc, Nri, Af_r, Ef_r, name)
end

export KineticReactionNetwork