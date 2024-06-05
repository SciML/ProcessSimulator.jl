Base.@kwdef struct KineticReactionNetwork
    Coef_Cr::Array{Float64, 2} #Stoichiometric coefficients of each component in each reaction (-)
    Do_r::Array{Float64, 2} # Forward order of the components ()
    substances_user::Array{String, 1} #Substances in the reaction network
    Nc::Int = size(substances_user)[1] #Number of components in the reaction network
    Nri::Int = size(Coef_Cr, 1) #Number of reactions in the reaction network (r)
    Af_r #Arrhenius constant of each reaction (s⁻¹)
    Ef_r #Activation energy of each reaction (J.mol⁻¹)
    name::String
end

my_reaction = KineticReactionNetwork(;substances_user = ["a","b"], 
Af_r = 2, Ef_r = 2, Coef_Cr = [1.0 1.0], Do_r = [1.0 1.0], name = "ReactionNetwork")

my_reaction.Nc

