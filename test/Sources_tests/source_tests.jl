using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using NonlinearSolve
using Test

@named source = MaterialSource(; substances_user = ["methane", "carbon monoxide"], P_user = 1.8e6, T_user = 150.0, Fₜ = 100.0, zₜ = [0.7, 0.3])
equations(source)
sys = structural_simplify(source, simplify = true)
variables = states(sys)
u0 = [x => 0.5 for x in variables]
prob = SteadyStateProblem(sys, u0)
sol = solve(prob, SSRootfind())
sol[states(source)[1]] # Trouble calculating the other states after structural simplification - why?