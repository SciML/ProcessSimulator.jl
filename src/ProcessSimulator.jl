module ProcessSimulator

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize, equations, get_unknowns

# Base
include("base/materials.jl")
include("base/base_components.jl")
include("base/utils.jl")

# Fluid handling
include("fluid_handling/compressors.jl")
include("fluid_handling/heat_exchangers.jl")

# Reactors
include("reactors/types.jl")
include("reactors/reaction_sets.jl")
include("reactors/CSTR.jl")
include("reactors/PFR.jl")
include("reactors/GibbsReactor.jl")

end
