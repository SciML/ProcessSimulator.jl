module ProcessSimulator

using JSON, JLD2
using Symbolics 
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Surrogates
using Clapeyron
import Clapeyron: enthalpy, molar_density

include("utils")
export load_component_properties, read_reidcp, my_model, enthalpy_simple, molar_density_simple

include("Sources/MaterialSource.jl")
export MaterialSource

include("Reactors/ReactionManager/KineticReaction.jl")
export KineticReactionNetwork

include("Reactors/CSTR.jl")
export CSTR

include("Sources/Sourceutils.jl")
export Display


end
