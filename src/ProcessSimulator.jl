module ProcessSimulator


using ModelingToolkit, JSON , DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron, Symbolics


export load_component_properties, read_reidcp, my_model, enthalpy_simple, molar_density_simple, MaterialSource, SimpleCSTR, my_model, KineticReactionNetwork, Display
include("utils")
include("Reactors/ReactionManager/KineticReaction.jl")
include("Reactors/SimplifiedCSTR.jl")

export MaterialSource
include("Sources/MaterialSource.jl")
include("Sources/Sourceutils.jl")
end
