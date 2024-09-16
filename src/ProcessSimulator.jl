module ProcessSimulator


using ModelingToolkit, JSON, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron

include("utils")
export load_component_properties, read_reidcp

include("Sources/MaterialSource.jl")
export MaterialSource

include("Reactors/ReactionManager/KineticReaction.jl")
export KineticReactionNetwork

include("Reactors/CSTR.jl")
export CSTR

include("Sources/Sourceutils.jl")
export thermal_energy_connector


end
