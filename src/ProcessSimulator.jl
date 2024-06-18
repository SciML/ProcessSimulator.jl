module ProcessSimulator


using ModelingToolkit, JSON , DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron, Symbolics

include("utils")
export load_component_properties, read_reidcp, my_model, enthalpy_simple, molar_density_simple

include("Sources/MaterialSource.jl")
export MaterialSource

include("Reactors/ReactionManager/KineticReaction.jl")
export KineticReactionNetwork

include("Reactors/SimplifiedCSTR.jl")
export SimpleCSTR

include("Sources/Sourceutils.jl")
export Display, matcon


end
