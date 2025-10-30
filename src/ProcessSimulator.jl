module ProcessSimulator


using LinearAlgebra
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize, equations, get_unknowns
using Symbolics

# Abstract types for unit operations
abstract type AbstractUnitOperation end
abstract type AbstractReactor <: AbstractUnitOperation end
abstract type AbstractSeparator <: AbstractUnitOperation end

# Base
include("utils/utils.jl")
include("base/basecomponents.jl")
include("pressure_drop/valve.jl")
include("separation/Adsorption.jl")
include("separation/FlashDrum.jl")
include("reactors/CSTR.jl")
include("base/print.jl")
include("base/solution_formatter.jl")
end
