using ProcessSimulator
using Test
using SafeTestsets

@safetestset "Gibbs reactor" begin
     include("Reactor_tests/gibbs_tests.jl")
end
