using ProcessSimulator
using Test
using SafeTestsets

@safetestset "Simple Steady State" begin
     include("base/simple_steady_state.jl")
end
