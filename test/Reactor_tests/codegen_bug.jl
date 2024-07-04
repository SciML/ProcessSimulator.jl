using ModelingToolkit
using Symbolics
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations

myf(x::Num, sym::Symbol) = -1.0

@register_symbolic myf(x::Num, sym::Symbol)

@component function myc(;name)

    vars = @variables begin
       x(t)
    end
    eqs = [x ~ myf(x, :myfun)]

    ODESystem([eqs...;], t, collect(Iterators.flatten(vars)), []; name)
end

component = myc(; name = :myc)
simple_prob = structural_simplify(component)
u0 = [simple_prob.x => -1.0]
prob = ODEProblem(simple_prob, u0, (0.0, 1.0))
sol = solve(prob, FBDF())