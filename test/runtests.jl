using ProcessSimulator
using Test
using SafeTestsets

@safetestset "Interface Compatibility" begin
    include("interface_tests.jl")
end

@safetestset "Base components" begin
    include("base/simple_steady_state.jl")
end

@safetestset "Reactors" begin
    include("reactors/simple_cstr.jl")
end
