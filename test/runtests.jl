using ProcessSimulator
using Test
using SafeTestsets

# New architecture tests - updated for the refactored codebase
@safetestset "Valve" begin
    include("base/valve_test.jl")
end

@safetestset "Steady-State CSTR" begin
    include("reactors/ss_cstr_test.jl")
end

@safetestset "Flash Drum" begin
    include("separation/flash_drum_test.jl")
end
