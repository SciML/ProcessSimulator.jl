using Pkg
using SafeTestsets
using Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "QA"
    Pkg.activate("qa")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("qa/qa.jl")
end

if GROUP == "All" || GROUP == "Core"
    @safetestset "Base components" begin
        include("base/simple_steady_state.jl")
    end

    @safetestset "Reactors" begin
        include("reactors/simple_cstr.jl")
    end
end
