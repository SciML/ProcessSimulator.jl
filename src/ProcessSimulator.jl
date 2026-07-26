module ProcessSimulator

using ModelingToolkit: ModelingToolkit, @component, @connector, @named,
    @parameters, @variables, Equation, Flow, ODESystem
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciMLPublic: @public
using Symbolics: scalarize

# Base
include("base/materials.jl")
include("base/base_components.jl")
include("base/utils.jl")

# Fluid handling
include("fluid_handling/compressors.jl")
include("fluid_handling/heat_exchangers.jl")

# Reactors
include("reactors/CSTR.jl")

end
