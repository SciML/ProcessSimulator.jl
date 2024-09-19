using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
using NonlinearSolve
import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs, scalarize
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

thermal_fluid_model = IAPWS95()
substances =  ["water"]

@named source = MaterialSource(;substances_user =  substances,
model = thermal_fluid_model,
P_user = 2*101325.0, T_user = 288.7,
Fₜ_user = 126.,
zₜ_user = [1.], 
guesses = Dict((:zᵂⱼᵢ) => ones(3, 1)*0.25))

@named J_101 = Jacket(; substances_user = substances,
    phase = :liquid,
    thermal_fluid_model = thermal_fluid_model,
    heat_transfer_coef = 914.
   )


@named energy_con = thermal_energy_connector()

energy_source = [energy_con.T ~ 350., energy_con.A ~ 9.23]
energy_source_to_J_101 = [connect(energy_con, J_101.EnergyCon)]
source_to_J_101 = [connect(source.Out, J_101.In)]

flowsheet = ODESystem([energy_source...; energy_source_to_J_101...;source_to_J_101...], t; name = :mycon, systems = [source, J_101, energy_con])   
flowsheet_ = structural_simplify(flowsheet)


variables = get_unknowns(flowsheet_)
u0 = [variables[1] => .5, variables[2] => 300., variables[3] => .5, variables[4] => 25.]
prob = SteadyStateProblem(flowsheet_, u0)

sol = solve(prob, SSRootfind())


for eq in unknowns(J_101)
    println(eq)
end

for eq in equations(flowsheet_)
    println(eq)
end