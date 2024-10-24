using ProcessSimulator
using ModelingToolkit, OrdinaryDiffEq, Clapeyron
using ModelingToolkit: t_nounits as t, D_nounits as D, scalarize
using ProcessSimulator: matcon


# Test for CSTR

#---------- Source objct
substances = ["water", "methanol", "propyleneglycol","methyloxirane"]
idealmodel = ReidIdeal(substances; userlocations = read_reidcp(substances))
PCSAFT_model = PCPSAFT(substances, idealmodel = idealmodel)

fugacity_coefficient(PCSAFT_model, 101325.0, 340.15, [0.0, 0.33, 0.33, 0.33])

@named source = MaterialSource(;substances_user =  substances,
model = PCSAFT_model,
P_user = 101325.0*2, T_user = 297.0,
Fₜ_user = (36.3 + 453.6 + 45.4)*1e3/3600,
zₜ_user = [0.8473, 1.0 - (0.0678 + 0.8473), 0.0, 0.0678], 
guesses = Dict((:zᵂⱼᵢ) => ones(3, 4)*0.25))


#----------------- CSTR object
Reaction = KineticReactionNetwork(substances_user = substances, 
Af_r = 4.71e9, Ef_r = 32400*1055.6/453.6, Coef_Cr = [-1.0 0.0 1.0 -1.0], 
Do_r = [0.0 0.0 0.0 1.0], name = "Propyleneglycol synthesis")

@named R_101 = CSTR(; substances_user = substances, 
    phase = :liquid, 
    model = PCSAFT_model,
    Reaction = Reaction,
    Ac = 1.93, #m²
    height_out_port = 0.0, #m
    guesses = Dict(:Cᵢ => [57252.65, 0.0, 0.0, 0.0], :V => 1.9,
    :F_out => (36.3 + 453.6 + 45.4)*1e3/3600, :Fʷ_out => (36.3 + 453.6 + 45.4)*18/3600)
    )


cons = [connect(source.Out, R_101.In)]

flowsheet = ODESystem(cons, t; name = :mycon, systems = [R_101, source])   
sistema = structural_simplify(flowsheet)

u0 = [sistema.R_101.Nᵢ[1] => 1.9*57252.65, sistema.R_101.Nᵢ[2] => 1e-10,
 sistema.R_101.Nᵢ[3] => 1e-10, sistema.R_101.Nᵢ[4] => 1e-10,
 sistema.R_101.T => 297.15]

prob = ODEProblem(sistema, u0, (0.0, 2*3600.0))
@time sol = solve(prob, QNDF(autodiff = true))
using Plots

plot(sol.t, sol[sistema.R_101.R[3]], label = "Rate of reaction")
plot(sol.t, sol[sistema.R_101.fᵢ], label = "fugacity")
plot(sol.t, sol[sistema.R_101.T], label = "Temperature")
plot(sol.t, sol[sistema.R_101.Cᵢ[1]], label = "water")
plot(sol.t, sol[sistema.R_101.V], label = "Volume")
plot(sol.t, sol[sistema.R_101.Cᵢ[3]], label = "propyleneglycol")
plot(sol.t, sol[sistema.R_101.Cᵢ[4]], label = "methyloxirane")
