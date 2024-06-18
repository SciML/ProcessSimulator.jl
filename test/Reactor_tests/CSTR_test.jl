using ProcessSimulator
using ModelingToolkit, DifferentialEquations, Clapeyron
import ModelingToolkit: get_unknowns, get_observed, get_defaults, get_eqs
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test


# Test for CSTR

#---------- Source objct
substances = ["water", "methanol", "propyleneglycol","methyloxirane"]
properties = Dict(subs => load_component_properties(subs) for subs in substances)
idealmodel = ReidIdeal(substances; userlocations = read_reidcp(properties, substances))
PCSAFT_model = PCPSAFT(substances, idealmodel = idealmodel)

@named source = MaterialSource(;substances_user =  substances,
model = PCSAFT_model,
P_user = 101325, T_user = 297.0,
Fₜ_user = (36.3 + 453.6 + 45.4)*1e3/3600,
zₜ_user = [0.8473, 1.0 - (0.0678 + 0.8473), 0.0, 0.0678], guesses = Dict((:zᵂⱼᵢ) => ones(3, 4)*0.25))


#----------------- CSTR object
Reaction = KineticReactionNetwork(;substances_user = substances, 
Af_r = 4.71e9, Ef_r = 32400*1055.6/453.6, Coef_Cr = [-1.0 0.0 1.0 -1.0], 
Do_r = [1.0 0.0 0.0 1.0], name = "Propyleneglycol synthesis")

#explicit enthalpy and density model
Cps = [69.21, 80.66, 192.50, 118.1] # At 298.00 K
pho_coef = Dict("a" => [78821.04345816873, 38133.33588802956, 18453.26055238924, 25211.86290505424], 
"b" => [-114.80261704286386, -85.22435577664774, -26.084451261000222, -73.31971618413455],
"c" => [0.208885275623327, 0.21416405813277503, 0.046784549601356355, 0.18909331028746998],
"d" => [-0.00022293440498512672, -0.0002675908503439664, -4.722426797584051e-5, -0.00022514899466957527])


mymodel = my_model(Cps, pho_coef)


@named Display_CSTR = Display(; Nc = 4)

@named R_101 = SimpleCSTR(; substances_user = substances, 
    phase = :liquid, 
    model = mymodel,
    Reaction = Reaction,
    ninports = 1, 
    Ac = 1.93, #m²
    height_out_port = 0.0, #m
    guesses = Dict(:Cᵢ => [57252.65, 0.0, 0.0, 0.0], :V => 1.9, :F_out => (36.3 + 453.6 + 45.4)*1e3/3600, :Fʷ_out => (36.3 + 453.6 + 45.4)*18/3600)
    )

#cons = [connect(source.Out, R_101.InPorts1),
#       connect(R_101.Out, Display_CSTR.InPort)]

@which connect(source.Out, R_101.InPorts1)

cons = [connect(source.Out, R_101.InPorts1)]

flowsheet = ODESystem(cons, t; name = :mycon, systems = [R_101, source])   
sistema = structural_simplify(flowsheet)
equations(sistema)
unknowns(sistema)
alg_equations(sistema)

u0 = [sistema.R_101.Nᵢ[1] => 1.9*57252.65, sistema.R_101.Nᵢ[2] => 0.0,
 sistema.R_101.Nᵢ[3] => 0.0, sistema.R_101.Nᵢ[4] => 0.0,
 sistema.R_101.T => 297.0]

prob = ODAEProblem(sistema, u0, (0.0, 1.0))
sol = solve(prob, FBDF(autodiff = false), abstol =  1e-6, reltol = 1e-6)
sol[sistema.R_101.InPorts1.P]
using Plots
plot(sol[sistema.source.Out.z₁])


solve(prob, Rodas)
#= using ProcessSimulator
using Test

y = Gibbs(2.0)

@test y == 4.0 =#