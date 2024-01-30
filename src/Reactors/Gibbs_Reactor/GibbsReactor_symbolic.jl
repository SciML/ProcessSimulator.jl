using Symbolics, SymbolicUtils, ModelingToolkit
using Reexport
@reexport using ModelingToolkit

using Clapeyron
using GCIdentifier
using DelimitedFiles
using LinearAlgebra

@__DIR__

components_names = ["carbon dioxide", "carbon monoxide", "hydrogen", "water", "methane"]
abspath(joinpath(@__DIR__, "myShomate.csv"))
syngas_prop_model =  ShomateIdeal(components_names; userlocations = [abspath(joinpath(@__DIR__, "myShomate.csv"))])

n_reactions = 2

function readEnthalpyEntropyData(file_path, component_names)
    # Read formation enthalpy from file
    # file_path: File path of the file containing the formation enthalpy
    # component_names: Names of the components
    # Returns: Formation enthalpy in J/mol

    # initialize the formation enthalpy and entropy vector
    H₀_S₀ = zeros(length(component_names), 2)

    # Read the TSV file
    data = readdlm(file_path, '\t', header = true)

    # Get the column of chemical names
    chemical_names = data[1][:, 2]

    # Find the index of the chemical name in the column
    for i in eachindex(component_names)
        index = findfirst(x -> x == component_names[i], chemical_names)

        # If the chemical name was found, return the corresponding line
        if index !== nothing
            H₀_S₀[i, :] =  data[1][index, 3:end]
        end
    end

    return H₀_S₀

end

H₀_S₀ = readEnthalpyEntropyData(abspath(joinpath(@__DIR__, "Hf0_Sf0.tsv")), components_names)

@register Entropy_TP1(T, P, H₀_S₀, ξ, ν, model, f_feed)

function Entropy_TP1(T, P, H₀_S₀, ξ, ν, model, f_feed)
    # Calculate the entropy of the system (J)      

    f_prod = f_feed + ν*ξ 

    ΔSrxn = dot(H₀_S₀[:, 2], ν*ξ) #Standard entropy of reaction times how much was converted from reactants to products

    S = entropy(model, P, T, f_prod) + ΔSrxn 

    return S

end

function GibbsReactor(; name, n_reactions, v)
end

ν = [4  0; -5 -1; 1 -3; -3  1; 1  1]

vars = @variables begin
    (T = 600.0), [description = "Reactor temperature"]
    (ξ[1:n_reactions] = 0.5), [description = "Extend of each independent reaction"]
    (entropy_tp = 0.0), [description = "Entropy of the system"]
end

pars = @parameters begin
    (P = 101325.0), [description = "Reactor pressure"]
    (T_feed = 500.0), [description = "Feed temperature"]
    (f_feed[1:length(components_names)] = [5.0, 20.0, 5.0, 5.0, 5.0]), [description = "Feed molar fraction"]
end

# TODO: structural_simplify optimization system.
eqns = [
entropy_tp ~ Entropy_TP1(T, P, H₀_S₀, ξ, ν, syngas_prop_model, f_feed)
]
