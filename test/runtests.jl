using ProcessSimulator
using Test
using SafeTestsets

@safetestset "Base components" begin
    include("base/simple_steady_state.jl")
end

@safetestset "Reactors" begin
    include("reactors/simple_cstr.jl")
end

@safetestset "Reactor Types" begin
    include("reactors/reactor_types_test.jl")
end
