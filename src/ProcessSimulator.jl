module ProcessSimulator


using ModelingToolkit, JSON , DifferentialEquations
import ModelingToolkit: scalarize, equations, get_unknowns
using Clapeyron
using NonlinearSolve

# Write your package code here.
export Gibbs
include("Reactors/Gibbs.jl")


export MaterialSource
include("Sources/MaterialSource.jl")
include("Sources/utils.jl")

end
