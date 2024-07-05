using ProcessSimulator
import ProcessSimulator: predict_ρ, predict_h
using ModelingToolkit, DifferentialEquations, Clapeyron
import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs, scalarize
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test, JLD2


# Test for CSTR

#---------- Source objct
substances = ["water", "methanol", "propyleneglycol","methyloxirane"]
properties = Dict(subs => load_component_properties(subs) for subs in substances)
read_reidcp(properties, substances)
idealmodel = ReidIdeal(substances; userlocations = read_reidcp(properties, substances))
PCSAFT_model = PCPSAFT(substances, idealmodel = idealmodel)


@named source = MaterialSource(;substances_user =  substances,
model = PCSAFT_model,
P_user = 5*101325, T_user = 297.0,
Fₜ_user = (36.3 + 453.6 + 45.4)*1e3/3600,
zₜ_user = [0.8473, 1.0 - (0.0678 + 0.8473), 1e-15, 0.0678], guesses = Dict((:zᵂⱼᵢ) => ones(3, 4)*0.25))


#----------------- CSTR object
Reaction = KineticReactionNetwork(;substances_user = substances, 
Af_r = 4.71e9, Ef_r = 32400*1055.6/453.6, Coef_Cr = [-1.0 0.0 1.0 -1.0], 
Do_r = [0.0 0.0 0.0 1.0], name = "Propyleneglycol synthesis")


R_101_guess = Dict(:Cᵢ => [57252.65, 1e-10,  1e-10,  1e-10], 
:V => 1.9, :Q_out => (36.3 + 453.6 + 45.4)*1e3/3600, 
:F_out => (36.3 + 453.6 + 45.4)*18/3600,
 :H => enthalpy(PCSAFT_model, 5*101325, 297.0, [1.9*57252.65, 1e-10, 1e-10, 1e-10]))

rho_surrogate = load_object("src/database/surrogates/rho_surrogate.jld");
h_surrogate = load_object("src/database/surrogates/H_surrogate.jld");
surrogate_model = surrogates(rho_surrogate, h_surrogate);

@variables T, N[1:4]
predict_ρ(surrogate_model, [T; N...])

@named R_101 = CSTR(; substances_user = substances, 
    phase = :liquid, 
    model = surrogate_model,
    Reaction = Reaction,
    ninports = 1, 
    Ac = 1.93, #m²
    height_out_port = 0.0, #m
    guesses = R_101_guess 
    )

for eq in equations(R_101)
    println(eq)
end

cons = [connect(source.Out, R_101.InPorts1)]

flowsheet = ODESystem(cons, t; name = :myflowsheet, systems = [R_101, source])   
sistema = structural_simplify(flowsheet; check_consistency = true)
#Inspecting flowsheet simplified equations
equations(R_101)
equations(sistema)
unknowns(sistema)
alg_equations(sistema)

u0 = [sistema.R_101.Nᵢ[1] => 1.9*57252.65, sistema.R_101.Nᵢ[2] => 1e-10,
      sistema.R_101.Nᵢ[3] => 1e-10, sistema.R_101.Nᵢ[4] => 1e-10,
      sistema.R_101.T => 297.0]

prob = ODEProblem(sistema, u0, (0.0, 1.0))

sol = solve(prob, KenCarp47(), abstol = 1e-8, reltol = 1e-8)
using Plots
plotly()
plot(sol.t, sol[sistema.R_101.M], label = "Molar density")
plot(sol.t, sol[sistema.R_101.ρʷ], label = "Molar density")
plot(sol.t, sol[sistema.R_101.T], label = "Temperature")
plot(sol.t, sol[sistema.R_101.V])
plot(sol.t, sol[sistema.R_101.Cᵢ[1]], label = "water")
plot(sol.t, sol[sistema.R_101.Nᵢ[2]], label = "methanol")
#plot(sol.t, sol[sistema.R_101.V], label = "Volume")
plot(sol.t, sol[sistema.R_101.Nᵢ[3]], label = "propyleneglycol")
plot(sol.t, sol[sistema.R_101.Nᵢ[4]], label = "methyloxirane")


solve(prob, Rodas)
#= using ProcessSimulator
using Test

y = Gibbs(2.0)

@test y == 4.0 =#