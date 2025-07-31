## Dynamic Flash Drum Operation - Ali M. Sahlodin, Harry A. J. Watson, and Paul I. Barton (2016) DOI 10.1002/aic.15378

using Revise
using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
#import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs, scalarize
using ModelingToolkit: t_nounits as t, D_nounits as D
using ProcessSimulator: matcon
using NonlinearSolve

substances = ["benzene", "toluene", "nitrogen"]

idealmodel = ReidIdeal(substances)

PR_model = PR(substances, idealmodel = idealmodel, translation = RackettTranslation)

molar_density(PR_model, 1.5*101325.0, 400.0, 1.0*[0.5, 0.5, 0.5], phase = :vapor)

x, n, G = Clapeyron.tp_flash(PR_model, 1.5*101325.0, 50.0, [0.0, 1.0, 1.0], DETPFlash())

bubble_temperature(PR_model, 1.5*101325.0, [0.0, 1.0, 1.0])

@named Drum_101 = ThreePortDrum(; substances = substances,
    Vtot_ = 0.2,
    non_condensables = [false, false, false],
    non_volatiles = [false, false, false],
    Ac = 0.1,
    model = PR_model)

for eq in equations(Drum_101)
    println(eq)
end

@named source = MaterialSource(; substances_user = substances,
    model = PR_model,
    P_user = 1.5*101325.0,
    T_user = 377.0, #Kelvin
    Fₜ_user = 2.0, #mol/s
    zₜ_user = [0.5, 0.5, 0.0],
    guesses = Dict((:zᵂⱼᵢ) => ones(3, 3)*0.33))

for eq in equations(source)
    println(eq)
end

@named Liquid_Valve = LinearValve(;
    Nc = 3,
    CV = 0.141,
    model = PR_model)

@named Vapor_Valve = LinearValve(;
    Nc = 3,
    CV = 0.01,
    model = PR_model)

cons = [connect(source.Out, Drum_101.In), connect(Drum_101.Outᴸ, Liquid_Valve.In),
    connect(Drum_101.Outᵍ, Vapor_Valve.In)]

check_and_block = [(Drum_101.Vᴸ <= 0.01*0.2) => [Liquid_Valve.θ ~ 0.0]
                   (Drum_101.Vᵍ <= 0.02*0.2) => [Vapor_Valve.θ ~ 0.0]
                   (Drum_101.P - Liquid_Valve.P <= 0.0) => [Liquid_Valve.θ ~ 0.0]
                   (Drum_101.P - Vapor_Valve.P <= 0.0) => [Vapor_Valve.θ ~ 0.0]]

flowsheet = ODESystem(
    cons, t; name = :myflowsheet, systems = [Drum_101, source, Liquid_Valve, Vapor_Valve],
    discrete_events = check_and_block)

sistema = structural_simplify(flowsheet)

u0 = [sistema.Drum_101.Nᵢ[1] => 0.0, sistema.Drum_101.Nᵢ[2] => 0.0,
    sistema.Drum_101.Nᵢ[3] => 10.0,
    sistema.Drum_101.T => 377.0, sistema.Liquid_Valve.θ => 0.5, sistema.Vapor_Valve.θ => 0.5]

unknowns(sistema)
alg_equations(sistema)

prob = ODEProblem(sistema, u0, (0.0, 0.01))

dt = 1e-7
sol = solve(prob,
    ImplicitEuler(nlsolve = NLNewton(check_div = false, always_new = true, relax = 4/10, max_iter = 100));
    dt,
    adaptive = false)
sol = solve(prob, FBDF())

sol[sistema.Drum_101.yᵢ]
a = sum.(sol[sistema.Drum_101.xᵢ]) - sum.(sol[sistema.Drum_101.yᵢ])
b = sol[sistema.Drum_101._0_Nᴸ]
c = sol[sistema.Drum_101._0_Nᵍ]

my_med(a, b, c) = a + b + c - min(a, b, c) - max(a, b, c)

my_med(a, b, c)
