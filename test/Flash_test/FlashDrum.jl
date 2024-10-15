## Dynamic Flash Drum Operation - Ali M. Sahlodin, Harry A. J. Watson, and Paul I. Barton (2016) DOI 10.1002/aic.15378

using Revise
using ProcessSimulator
using ModelingToolkit, OrdinaryDiffEq, Clapeyron
#import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs, scalarize
using ModelingToolkit: t_nounits as t, D_nounits as D
using ProcessSimulator: matcon

substances = ["benzene", "toluene", "nitrogen"]

idealmodel = ReidIdeal(substances)

PR_model = PR(substances, idealmodel = idealmodel, translation = RackettTranslation)

@named Flash = ThreePortDrum(; substances = substances, 
    non_condensables = [false, false, true],
    non_volatiles = [false, false, false],
    Ac = 0.1,
    model = PR_model)

@named source = MaterialSource(;substances_user =  substances,
model = PR_model,
P_user = 1.5*101325.0,
T_user = 350.0,
Fₜ_user = 2.0,
zₜ_user = [0.5, 0.5, 0.0], 
guesses = Dict((:zᵂⱼᵢ) => ones(3, 4)*0.25))


for eq in equations(Flash)
    println(eq)
end