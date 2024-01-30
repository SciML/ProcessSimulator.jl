using Pkg
Pkg.activate("GibbsReactor")

using ModelingToolkit
using Clapeyron
using GCIdentifier
using DelimitedFiles
using LinearAlgebra

#define components
components_names = ["carbon dioxide", "carbon monoxide", "hydrogen", "water", "methane"]
components_smiles = ["C(=O)=O", "[C-]#[O+]", "[HH]", "O=O", "C"]


#Define the model
my_model = ShomateIdeal(components_names;
 userlocations = ["C:/Users/Vinic/OneDrive/Pos-Doc/Flowsheeting/GibbsReactor/GibbsReactor.jl", "myShomate.csv"])

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

H₀_S₀ = readEnthalpyEntropyData("Hf0_Sf0.tsv", components_names)


# Define reaction coefficients:
# Reaction 1: 5 CO + 3 H20 -> 4 CO2 + H2 + CH4
# Reaction 2: CO + 3 H2 -> CH4 + H20 

ν = [4  0; -5 -1; 1 -3; -3  1; 1  1] # Should be done automatically in the future 
ν*[0.0; 0.0] 

function GibbsFreeEnergy_TP(T, P, H₀_S₀, ξ, ν, model, f_feed)
    # Calculate the Gibbs free energy of the system
    # T: Temperature in K
    # P: Pressure in Pa
    # ξ: Extent of reaction
    # ν: Stoichiometric coefficients
    # my_model: Model of the system
    # Returns: Gibbs free energy in J

    f_prod = f_feed + ν*ξ #Mass balance 
    ΔHrxn = dot(H₀_S₀[:, 1], ν*ξ) #Heat released or taken from Reaction
    ΔSrxn = dot(H₀_S₀[:, 2], ν*ξ) #Standard entropy of reaction
    H = enthalpy(model, P, T, f_prod) + ΔHrxn #Outlet stream enthalpy
    S = entropy(model, P, T, f_prod) + ΔSrxn #Outlet stream entropy
    G = H - T*S   
    return G

end



function Entropy_TP(T, P, H₀_S₀, ξ, ν, model, f_feed)
    # Calculate the entropy of the system (J)      

    f_prod = f_feed + ν*ξ 

    ΔSrxn = dot(H₀_S₀[:, 2], ν*ξ) #Standard entropy of reaction times how much was converted from reactants to products

    S = entropy(model, P, T, f_prod) + ΔSrxn 

    return S

end



function EnergyBalance(T, P, H₀_S₀, ξ, ν, model, p)
    # Calculate the energy balance of the system
    # res: Energy balance error
    # ξ: Extent of reaction
    # p: Parameters of the system
    # Returns: Energy balance error (res) in J

    f_feed, T_feed = p[1], p[2] 

    # Calculate the energy balance

    f_prod = f_feed + ν*ξ #mass balance

    Hin = enthalpy(model, P, T_feed, f_feed)
    Hout = enthalpy(model, P, T, f_prod)
    ΔHrxnX = dot(H₀_S₀[:, 1], ν*ξ)

    res = Hin + ΔHrxnX - Hout 

    return res

end

#Testing the functions:
GibbsFreeEnergy_TP(300, 10^5, H₀_S₀, [0.0, 0.0], ν, my_model, [5.0, 20.0, 5.0, 5.0, 5.0])
Entropy_TP(300, 10^5, H₀_S₀, [0.0, 0.0], ν, my_model, [5.0, 20.0, 5.0, 5.0, 5.0])
EnergyBalance(500, 10^5, H₀_S₀, [0.0, 0.0], ν, my_model, ([5.0, 20.0, 5.0, 5.0, 5.0], 500))

#Minimizing S(T,ξ) with respect to ξ constrained to energy balance (adiabatic reactor):
using OptimizationOptimJL, Optimization, OptimizationMOI, Ipopt

f_feed = [5.0, 20.0, 5.0, 5.0, 5.0] #["carbon dioxide", "carbon monoxide", "hydrogen", "water", "methane"]
T_feed = 500.0 #K
P = 101325.0 #Pa

#Target function
Smin(ξ_T) = -Entropy_TP(ξ_T[end], P, H₀_S₀, ξ_T[1:2], ν, my_model, f_feed)


#Energy balance constraint
E(ξ_T) = EnergyBalance(ξ_T[end], P, H₀_S₀, ξ_T[1:2], ν, my_model, [f_feed, T_feed])

#product flow rate positive
function f_prod_cons(ξ_T, f_feed, ν)
    
    res = f_feed + ν*ξ_T[1:2]

    return res

end


#Initial guess
ξ₀ = [0.000; 0.000]
T₀ = 600.0

#Testing the functions:
Smin([ξ₀ ; T₀])
E([ξ₀ ; T₀])
f_prod_cons([ξ₀ ; T₀], f_feed, ν)


#Lower and upper bounds
lb = [0.00; zeros(size(components_names, 1))]
ub = [0.00; Inf*ones(size(components_names, 1))]

optf = OptimizationFunction((x, p) -> Smin(x), Optimization.AutoForwardDiff(), cons = (res, x, p) -> [E(x), f_prod_cons(x, f_feed, ν)])

prob = OptimizationProblem(optf, [ξ₀; T₀], lcons = lb, ucons = ub)

sol = solve(prob, Ipopt.Optimizer())
