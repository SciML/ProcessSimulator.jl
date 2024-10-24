using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t

const PS = ProcessSimulator

# Material Parameters
ϱ = [14.8,55.3,13.7,24.7].*1e3          # mol/m^3   
cₚ = [35,18,46,19.5]*4.184              # J/mol/K

# Reaction
ν = [-1.0, -1.0, 1.0, 0.0]
Δhᵣ = 20013*4.184                       # J/mol

matsource = PS.MaterialSource(["propylene oxide (A)","water (B)","propylene glycol (C)","methanol (M)"];
    Mw                  = [0.05808, 0.01801, 0.07609, 0.03204],
    molar_density       = (p, T, x; kwargs...) -> sum([ϱ[i]*x[i] for i in 1:4]),
    VT_enthalpy         = (ϱ, T, x) -> sum([cₚ[i]*x[i] for i in 1:4])*T,
    VT_internal_energy  = (ϱ, T, x) -> sum([cₚ[i]*x[i] for i in 1:4])*T,
    reactions           = [PS.Reaction(
        ν               = ν,
        r               = (p,T,x) -> 2.73e-4*exp(9059*(1/297-1/T))*x[1],
        Δhᵣ             = (T) -> Δhᵣ,
    )],
)

# Parameters and boundary conditions
V_cstr = 1.89                            # m^3
F = [36.3e3, 453.6e3, 0, 45.4e3]./3600  # mol/s

# Initial conditions
T0 = 297. 
n0 = 55.3e3*V_cstr
m0_inlet = F' * matsource.Mw

# Heat exchanger
UA = 7262*4184/3600                     # J/K/s
T1_coolant = 288.71                     # K
cₚ_coolant = 18*4.184                   # J/mol/K   
n_coolant = 453.6e3/3600                # mol/s

# Create flowsheet
@named inlet = PS.Port(matsource)
@named cstr = PS.CSTR(matsource)
@named outlet = PS.Port(matsource)

# Connect the flowsheet
eqs = [
    connect(inlet.c,cstr.cv.c1),
    connect(cstr.cv.c2,outlet.c),
    cstr.Q ~ -n_coolant*cₚ_coolant*(cstr.cv.T-T1_coolant)*(1-exp(-UA/(n_coolant*cₚ_coolant))),
]

@named flowsheet_ = ODESystem(eqs, t, [], [], systems=[inlet,cstr,outlet])

pars = [
    cstr.cv.V => V_cstr,
]

inp = [
    inlet.T => 297.0,
    inlet.p => 1e5,
    inlet.n => sum(F),
    inlet.xᵢ[1] => F[1]/sum(F),
    inlet.xᵢ[2] => F[2]/sum(F),
    inlet.xᵢ[3] => F[3]/sum(F),
    outlet.p => 1e5,
]

u0 = [
    cstr.cv.U => matsource.VT_internal_energy(NaN, T0, [0,1,0,0])*n0,
    cstr.cv.nᵢ[1,1] => 0.0,
    cstr.cv.nᵢ[1,2] => n0,
    cstr.cv.nᵢ[1,3] => 0.0,
    cstr.cv.nᵢ[1,4] => 0.0,
    inlet.m => m0_inlet,
    cstr.cv.T => 297.,
    cstr.cv.c2.n => 0.0,
    outlet.m => 0.,
]

# u0 = [
#     cstr.cv.U => matsource.VT_internal_energy(NaN, T0, [0,1,0,0])*n0,
#     cstr.cv.nᵢ[1,1] => 0.0,
#     cstr.cv.nᵢ[1,2] => n0,
#     cstr.cv.nᵢ[1,3] => 0.0,
#     cstr.cv.nᵢ[1,4] => 0.0,
#     cstr.cv.T => 297.,
#     inlet.m => m0_inlet,
#     outlet.m => 0.,
# ]

flowsheet,idx = structural_simplify(flowsheet_,(first.(inp),[]))

# prob_steady = SteadyStateProblem(flowsheet, vcat(pars,inp), u0)
# sol_steady = solve(prob_steady)

prob = ODEProblem(flowsheet, u0, (0, 4).*3600., vcat(pars,inp))#; guesses = [outlet.m => -m0_inlet])
@time sol = solve(prob, QNDF(), abstol =  1e-6, reltol = 1e-6)


# -------------------------- Plotting and evaluation ------------------------- #
using Plots, PrettyTables
pythonplot()
# gr()
# plotlyjs()
closeall()

t_h = sol.t./3600

ps = [plot(;framestyle=:box,xlabel="t / h",xlims=(0.0,maximum(sol.t)*1.02/3600)) for i in 1:6]
plot!(ps[1],t_h,sol[cstr.cv.T],label="T",ylabel="T / K")
[plot!(ps[2],t_h, sol[cstr.cv.nᵢ[1,i]]./V_cstr./1e3;label="x[$i]",ylabel="Cᵢ / kmol/m³") for i in 1:4]
plot!(ps[3],t_h,sol[cstr.cv.n]./1e3,label="n",xlabel="t / h",ylabel="n / kmol")
plot!(ps[4],t_h,matsource.reaction[1].r.(NaN,sol[cstr.cv.T],sol[cstr.cv.xᵢ]),label="r",xlabel="t / h",ylabel="r / kmol")
plot!(ps[5],t_h,sol[cstr.cv.ΔnR[3]]./1e3,label="Q",xlabel="t / h",ylabel="ΔnR / kmol")
plot!(ps[6],t_h,sol[cstr.cv.c1.h],label="h_in",xlabel="t / h",ylabel="ΔE / Js")
plot!(ps[6],t_h,sol[cstr.cv.c2.h],label="h_out") 
plot!(ps[6],t_h,sol[cstr.cv.ΔHᵣ],label="ΔHᵣ")
plot!(ps[6],t_h,abs.(sol[cstr.cv.q1.Q]),label="|Q|")
fig = plot(ps...;layout=(2,3),size=(1200,800))
display(fig)

# Plot Fogler (Figures E13-3.1+3.2) 
plot()
plts_Fogler = [plot(;framestyle=:box,xlabel="t / h",xlims=(0,4)) for i in 1:2]
plot!(plts_Fogler[1],t_h, sol[cstr.cv.nᵢ[1,1]]./V_cstr./1e3; ylabel="Cₐ / kmol/m³")
plot!(plts_Fogler[2],t_h, sol[cstr.cv.T]; ylabel="T / K")
fig_Fogler = plot(plts_Fogler...;layout=(1,2))
display(fig_Fogler)

# Plots masss and volume
V_sim = sol[cstr.cv.n] ./ sol[cstr.cv.ϱ[1]]
m_sim = [(ni * matsource.Mw)[1] for ni in sol[cstr.cv.nᵢ]]

fig_Vm = plot(
    plot(t_h,V_sim./V_sim[1];xlabel="t / h",ylabel="V / m³",framestyle=:box),
    plot(t_h,m_sim;xlabel="t / h",ylabel="m / kg",framestyle=:box);
    layout=(2,1)
)
display(fig_Vm)

# Table
nᵢ_sol = [sol[cstr.cv.nᵢ[1,i]] for i in 1:4]
header =[
    "Variable",     "Initial Fogler",   "Initial Sim.",                 "Final Fogler",     "Final Sim." 
]
data = [
    "Ca"            0.0                 nᵢ_sol[1][1]/V_cstr/1e3     0.658258            nᵢ_sol[1][end]/V_cstr/1e3;
    "Cb"            55.3                nᵢ_sol[2][1]/V_cstr/1e3     34.06019            nᵢ_sol[2][end]/V_cstr/1e3;
    "Cc"            0.0                 nᵢ_sol[3][1]/V_cstr/1e3     2.247301            nᵢ_sol[3][end]/V_cstr/1e3;
    "Nb"            104.517             nᵢ_sol[2][1]/1e3            64.37375            nᵢ_sol[2][end]/1e3;
    "Nc"            0.0                 nᵢ_sol[3][1]/1e3            4.2474              nᵢ_sol[3][end]/1e3;
    "Nm"            0.0                 nᵢ_sol[4][1]/1e3            6.868166            nᵢ_sol[4][end]/1e3;
    "Qr2"           39920*4184/3600     sol[cstr.cv.q1.Q][1]        205900              sol[cstr.cv.q1.Q][end];
    "T"             297.0               sol[cstr.cv.T][1]           331.4976            sol[cstr.cv.T][end];
]

pretty_table(data; header=header)