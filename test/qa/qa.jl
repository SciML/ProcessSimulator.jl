using ProcessSimulator
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(ProcessSimulator)
end

@testset "JET" begin
    JET.test_package(ProcessSimulator; target_defined_modules = true)
end
