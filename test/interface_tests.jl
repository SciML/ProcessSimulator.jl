using ProcessSimulator
using Test

const PS = ProcessSimulator

@testset "Interface Compatibility" begin
    @testset "Parameterized Reaction type" begin
        # Test Float64 (default)
        ν_f64 = [-1.0, -1.0, 1.0, 0.0]
        r_f64 = PS.Reaction(
            ν = ν_f64,
            r = (p, T, x) -> 1.0,
            Δhᵣ = (T) -> 1000.0
        )
        @test eltype(r_f64.ν) == Float64
        @test r_f64 isa PS.Reaction{Float64}

        # Test BigFloat
        ν_bf = BigFloat[-1.0, -1.0, 1.0, 0.0]
        r_bf = PS.Reaction(
            ν = ν_bf,
            r = (p, T, x) -> BigFloat("1.0"),
            Δhᵣ = (T) -> BigFloat("1000.0")
        )
        @test eltype(r_bf.ν) == BigFloat
        @test r_bf isa PS.Reaction{BigFloat}

        # Test Float32
        ν_f32 = Float32[-1.0, -1.0, 1.0, 0.0]
        r_f32 = PS.Reaction(
            ν = ν_f32,
            r = (p, T, x) -> Float32(1.0),
            Δhᵣ = (T) -> Float32(1000.0)
        )
        @test eltype(r_f32.ν) == Float32
        @test r_f32 isa PS.Reaction{Float32}
    end

    @testset "Parameterized MaterialSource type" begin
        # Test Float64 (default)
        ms_f64 = PS.MaterialSource("test";
            Mw = 0.004,
            molar_density = (p, T, x; kwargs...) -> p / (8.314 * T),
            VT_enthalpy = (ϱ, T, x) -> 1000.0 * T
        )
        @test eltype(ms_f64.Mw) == Float64
        @test ms_f64 isa PS.MaterialSource{Float64}

        # Test BigFloat
        ms_bf = PS.MaterialSource("test";
            Mw = BigFloat("0.004"),
            molar_density = (p, T, x; kwargs...) -> p / (BigFloat("8.314") * T),
            VT_enthalpy = (ϱ, T, x) -> BigFloat("1000.0") * T
        )
        @test eltype(ms_bf.Mw) == BigFloat
        @test ms_bf isa PS.MaterialSource{BigFloat}

        # Test Float32
        ms_f32 = PS.MaterialSource("test";
            Mw = Float32(0.004),
            molar_density = (p, T, x; kwargs...) -> p / (Float32(8.314) * T),
            VT_enthalpy = (ϱ, T, x) -> Float32(1000.0) * T
        )
        @test eltype(ms_f32.Mw) == Float32
        @test ms_f32 isa PS.MaterialSource{Float32}

        # Test with vector input
        ms_vec = PS.MaterialSource(["A", "B"];
            Mw = [0.018, 0.028],
            molar_density = (p, T, x; kwargs...) -> p / (8.314 * T),
            VT_enthalpy = (ϱ, T, x) -> 1000.0 * T
        )
        @test length(ms_vec.Mw) == 2
        @test ms_vec.N_c == 2
    end

    @testset "MaterialSource with Reactions" begin
        # Test that reactions preserve their type parameter
        ν_bf = BigFloat[-1.0, -1.0, 1.0, 0.0]
        r_bf = PS.Reaction(
            ν = ν_bf,
            r = (p, T, x) -> BigFloat("1.0"),
            Δhᵣ = (T) -> BigFloat("1000.0")
        )

        ms_with_reaction = PS.MaterialSource(
            ["A", "B", "C", "D"];
            Mw = BigFloat[0.05808, 0.01801, 0.07609, 0.03204],
            molar_density = (p, T, x; kwargs...) -> BigFloat("1000"),
            VT_enthalpy = (ϱ, T, x) -> BigFloat("100") * T,
            reactions = [r_bf]
        )
        @test eltype(ms_with_reaction.Mw) == BigFloat
        @test length(ms_with_reaction.reaction) == 1
        @test eltype(ms_with_reaction.reaction[1].ν) == BigFloat
    end

    @testset "Backward compatibility" begin
        # Ensure existing code patterns still work
        Mw = 0.004
        R = 2.1e3 * Mw
        cₚ = 5.2e3 * Mw
        cᵥ = cₚ - R

        matsource = PS.MaterialSource("helium";
            Mw = 0.004,
            molar_density = (p, T, x; kwargs...) -> p / (R * Mw * T),
            VT_internal_energy = (ϱ, T, x) -> cᵥ * T,
            VT_enthalpy = (ϱ, T, x) -> cₚ * T,
            VT_entropy = (ϱ, T, x) -> cᵥ * log(T) + R * log(1 / ϱ)
        )

        @test matsource.name == "helium"
        @test matsource.N_c == 1
        @test matsource.Mw[1] == 0.004
        @test eltype(matsource.Mw) == Float64
    end
end
