module ProcessSimulator


using ModelingToolkit, JSON , DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron
using NonlinearSolve


export MaterialSource
include("utils")
include("Sources/MaterialSource.jl")
include("Sources/Sourceutils.jl")


end
