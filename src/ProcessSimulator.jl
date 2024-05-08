module ProcessSimulator


using ModelingToolkit, JSON , DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns
using Clapeyron
using NonlinearSolve


export MaterialSource
include("Sources/MaterialSource.jl")
include("Sources/utils.jl")

end
