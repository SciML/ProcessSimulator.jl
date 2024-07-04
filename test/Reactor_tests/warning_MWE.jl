using ModelingToolkit
using DifferentialEquations
using Clapeyron
using ModelingToolkit: t_nounits as t, D_nounits as D, scalarize

my_function(p, T, N::AbstractArray) = 1.0

@register_symbolic my_function(p, T, N::AbstractArray)

@component function StirredHeatedTank(; model, name)

    Nc = 2

    vars = @variables begin
        N1(t), [description = "molar hold up - component i (mol)", guess = 10000/Nc]
        N2(t), [description = "molar hold up - component i (mol)", guess = 10000/Nc]
        P(t), [description = "Pressure (Pa)", guess = 101325.0]
        T(t), [description = "Temperature (K)", guess = 297.0]
        V(t), [description = "Volume (m³)", guess = 2.0]
        rho(t), [description = "Density (kg/m³)", guess = 1000.0]
    end

    mass_balance = [D(N1) ~ -N1/V, D(N2) ~ -N2/V, P  ~ 501325.0, T ~ 297.0]
    properties = [rho ~ my_function(model, P, T, [N1, N2]), N1 + N2 ~ rho*V] 

    eqs = [mass_balance...; properties...]

    ODESystem([eqs...;], t, collect(Iterators.flatten(vars)), []; name)
end



tank = StirredHeatedTank(; model = PCSAFT_model, name = :tank)

sistema = structural_simplify(tank; check_consistency = true)

u0 = [sistema.N1 => 1.9*57252.65, sistema.N2 => 0.0]
equations(sistema)
prob = ODEProblem(sistema, u0, (0.0, 1.0))
sol = solve(prob, FBDF(autodiff = true), abstol =  1e-8, reltol = 1e-8)

#= Nc = 12

# Include subscript 0 in the range
subscripts = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']

# Adjusted generation of symbols with subscript
# For numbers 10 and above, concatenate subscripts for each digit
Cs = [Symbol("N", i <= 9 ? subscripts[i + 1] : join([subscripts[div(i, 10) + 1], subscripts[mod(i, 10) + 1]])) for i in 0:Nc-1]

# Directly create variables from symbols
variables = collect(Iterators.Flatten([@eval @variables $(Cs[i])(t) for i in 1:Nc])) =#
