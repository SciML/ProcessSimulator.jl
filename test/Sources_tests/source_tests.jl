using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
using NonlinearSolve
using Test
model = PCPSAFT(["propyleneglycol", "methyloxirane", "water", "methanol"], idealmodel = BasicIdeal)
xᵢⱼ, nᵢⱼ, G = tp_flash(model, 40_000, 300.15, [1.0, 1.0, 1.0, 1.0], DETPFlash(; equilibrium = :vle))
@time dew_pressure(model, 300.15, [1.0, 1.0, 1.0, 1.0])

@named source = MaterialSource(; substances_user =  ["propyleneglycol", "methyloxirane", "water", "methanol"],
model = model,
 P_user = 101325, T_user = 298.0,
  Fₜ_user = 100.0, zₜ_user = [0.5, 0.5, 0.5, 0.5]/2.0
  )

equations(source)
for eq in equations(source)
    println(eq)
end
sys = structural_simplify(source, simplify = true)
variables = states(sys)
u0 = [x => 0.5 for x in variables]
prob = SteadyStateProblem(sys, u0, checkbounds=true)
sol = solve(prob, SSRootfind())
sol[states(source)[1]] # Trouble calculating the other states after structural simplification - why?
for ob in observed(sys)
    println(ob)
end

model = PR(["methane", "propane"], idealmodel = ReidIdeal)

enthalpy(model, 1e-10, 298.00, [1.0, 1.0])