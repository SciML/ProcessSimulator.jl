using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs, scalarize
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test


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
zₜ_user = [0.8473, 1.0 - (0.0678 + 0.8473), 0.0, 0.0678], guesses = Dict((:zᵂⱼᵢ) => ones(3, 4)*0.25))


#----------------- CSTR object
Reaction = KineticReactionNetwork(;substances_user = substances, 
Af_r = 4.71e9, Ef_r = 32400*1055.6/453.6, Coef_Cr = [-1.0 0.0 1.0 -1.0], 
Do_r = [0.0 0.0 0.0 1.0], name = "Propyleneglycol synthesis")

#= #explicit enthalpy and density model
h_0 = [-43500.64079257049, -37861.75409062767, -66701.57862546679, -28309.177907687146]
Cps = [69.21, 80.66, 192.50, 118.1] # At 298.00 K
pho_coef = Dict("a" => [78821.04345816873, 38133.33588802956, 18453.26055238924, 25211.86290505424], 
"b" => [-114.80261704286386, -85.22435577664774, -26.084451261000222, -73.31971618413455],
"c" => [0.208885275623327, 0.21416405813277503, 0.046784549601356355, 0.18909331028746998],
"d" => [-0.00022293440498512672, -0.0002675908503439664, -4.722426797584051e-5, -0.00022514899466957527])


mymodel = my_model(Cps, pho_coef, h_0) =#

R_101_guess = Dict(:C1 => 57252.65, :C2 => 0.0, :C3 => 0.0, :C4 => 0.0, 
:V => 1.9, :Q_out => (36.3 + 453.6 + 45.4)*1e3/3600, 
:F_out => (36.3 + 453.6 + 45.4)*18/3600,
 :H => enthalpy(PCSAFT_model, 5*101325, 297.0, [1.9*57252.65, 0.0, 0.0, 0.0]))

@named R_101 = SimpleCSTR(; substances_user = substances, 
    phase = :liquid, 
    model = PCSAFT_model,
    Reaction = Reaction,
    ninports = 1, 
    Ac = 1.93, #m²
    height_out_port = 0.0, #m
    guesses = R_101_guess 
    )


#Inspecting reactor equations
for eq in full_equations(R_101)
    println(eq)
end
get_defaults(R_101)

@which Symbolics.connect(source.Out, R_101.InPorts1)

cons = [connect(source.Out, R_101.InPorts1)]

flowsheet = ODESystem(cons, t; name = :myflowsheet, systems = [R_101, source])   
sistema = structural_simplify(flowsheet; check_consistency = true)
#Inspecting flowsheet simplified equations
equations(sistema)
unknowns(sistema)
alg_equations(sistema)

u0 = [sistema.R_101.N1 => 1.9*57252.65, sistema.R_101.N2 => 0.0,
 sistema.R_101.N3 => 0.0, sistema.R_101.N4 => 0.0,
 sistema.R_101.T => 297.0]

prob = ODEProblem(sistema, u0, (0.0, 1.0))
sol = solve(prob, FBDF(autodiff = false), abstol =  1e-8, reltol = 1e-8)
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