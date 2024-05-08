using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
using NonlinearSolve
using Test

@named source = MaterialSource(; substances_user = ["methane", "carbon monoxide"], P_user = 2e6, T_user = 160.0, Fₜ = 100.0, zₜ = [1.0, 0.0])
equations(source)
sys = structural_simplify(source, simplify = true)
variables = states(sys)
u0 = [x => 0.5 for x in variables]
prob = SteadyStateProblem(sys, u0, checkbounds=true)
sol = solve(prob, SSRootfind())
sol[states(source)[1]] # Trouble calculating the other states after structural simplification - why?
