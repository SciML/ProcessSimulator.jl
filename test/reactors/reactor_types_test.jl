using ProcessSimulator
using Test

const PS = ProcessSimulator

@testset "Reactor Type Hierarchy" begin
    # Test abstract type hierarchy
    @test PS.AbstractEquilibriumReactor <: PS.AbstractReactor
    @test PS.AbstractKineticReactor <: PS.AbstractReactor
end

@testset "ReactionSet Utilities" begin
    # Define a simple reaction: A -> B
    r1 = PS.Reaction(
        ν = [-1.0, 1.0],
        r = (p, T, x) -> 1e-3 * x[1],  # First order in A
        Δhᵣ = (T) -> -50000.0  # Exothermic, 50 kJ/mol
    )

    rs = PS.ReactionSet([r1], ["A", "B"]; name = "SimpleReaction")

    @test PS.nreactions(rs) == 1
    @test PS.ncomponents(rs) == 2

    ν_matrix = PS.stoichiometry_matrix(rs)
    @test ν_matrix[1, 1] == -1.0  # A is consumed
    @test ν_matrix[1, 2] == 1.0   # B is produced
end

@testset "ArrheniusKinetics" begin
    # Create Arrhenius kinetics for a typical reaction
    # A = 1e10 s^-1, Ea = 50 kJ/mol
    k = PS.ArrheniusKinetics(1e10, 50000.0)

    # Test rate constant at different temperatures
    k_300 = PS.rate_constant(k, 300.0)
    k_400 = PS.rate_constant(k, 400.0)

    # Higher temperature should give higher rate constant
    @test k_400 > k_300

    # k = A * exp(-Ea/(R*T)) where R = 8.314 J/(mol·K)
    # At 300K: k = 1e10 * exp(-50000/(8.314*300)) ≈ 1e10 * exp(-20.05) ≈ 19.7
    expected_k_300 = 1e10 * exp(-50000.0 / (8.314 * 300.0))
    @test isapprox(k_300, expected_k_300, rtol = 1e-6)

    # At 400K: k = 1e10 * exp(-50000/(8.314*400)) ≈ 1e10 * exp(-15.03) ≈ 2.97e3
    expected_k_400 = 1e10 * exp(-50000.0 / (8.314 * 400.0))
    @test isapprox(k_400, expected_k_400, rtol = 1e-6)
end
