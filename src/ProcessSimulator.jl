module ProcessSimulator


using ModelingToolkit, JSON, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron

include("utils")

include("Sources/MaterialSource.jl")

include("Reactors/ReactionManager/KineticReaction.jl")

include("Reactors/CSTR.jl")

include("HeatExchange/Jacket.jl")

include("Sources/Sourceutils.jl")


end
