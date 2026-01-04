using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs
using ModelingToolkit: t_nounits as t, D_nounits as D
using JSON
using NonlinearSolve
using Test

substances = ["water", "methanol", "propyleneglycol", "methyloxirane"]
properties = Dict(subs => load_component_properties(subs) for subs in substances)
idealmodel = ReidIdeal(substances; userlocations = read_reidcp(properties, substances))
PCSAFT_model = PCPSAFT(substances, idealmodel = idealmodel)

@named source = MaterialSource(;
    substances_user = substances,
    model = PCSAFT_model,
    P_user = 101325, T_user = 297.0,
    Fₜ_user = (36.3 + 453.6 + 45.4) * 1.0e3 / 3600,
    zₜ_user = [0.8473, 1.0 - (0.0678 + 0.8473), 0.0, 0.0678],
    guesses = Dict((:zᵂⱼᵢ) => ones(3, 4) / 0.25)
)

#@named myDisplay = Display(; Nc = 4)

#connections = ODESystem([connect(source.Out, myDisplay.InPort)], t, [], [] ; name = :connection)

#system = compose(connections, [source, myDisplay])

sys = structural_simplify(source)
variables = get_unknowns(sys)
u0 = [x => 0.5 for x in variables]
prob = SteadyStateProblem(sys, u0)
sol = solve(prob, SSRootfind())
sol[source.Hⱼ]
