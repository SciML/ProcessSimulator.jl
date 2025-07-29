using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t

const PS = ProcessSimulator

# Material Parameters
ϱ = [14.8, 55.3, 13.7, 24.7] .* 1e3          # mol/m^3   
cₚ = [35, 18, 46, 19.5] * 4.184              # J/mol/K

# Reaction
ν = [-1.0, -1.0, 1.0, 0.0]
Δhᵣ = 20013 * 4.184                       # J/mol

matsource = PS.MaterialSource(
    ["propylene oxide (A)", "water (B)", "propylene glycol (C)", "methanol (M)"];
    Mw = [0.05808, 0.01801, 0.07609, 0.03204],
    molar_density = (p, T, x; kwargs...) -> sum([ϱ[i] * x[i] for i in 1:4]),
    VT_enthalpy = (ϱ, T, x) -> sum([cₚ[i] * x[i] for i in 1:4]) * T,
    VT_internal_energy = (ϱ, T, x) -> sum([cₚ[i] * x[i] for i in 1:4]) * T,
    reactions = [PS.Reaction(
        ν = ν,
        r = (p, T, x) -> 2.73e-4 * exp(9059 * (1 / 297 - 1 / T)) * x[1],
        Δhᵣ = (T) -> Δhᵣ
    )]
)

# Parameters and boundary conditions
V_cstr = 1.89                            # m^3
F = [36.3e3, 453.6e3, 0, 45.4e3] ./ 3600  # mol/s

# Initial conditions
T0 = 297.0
n0 = 55.3e3 * V_cstr
m0_inlet = F' * matsource.Mw

# Heat exchanger
UA = 7262 * 4184 / 3600                     # J/K/s
T1_coolant = 288.71                     # K
cₚ_coolant = 18 * 4.184                   # J/mol/K   
n_coolant = 453.6e3 / 3600                # mol/s

# Create flowsheet
@named inlet = PS.Port(matsource)
@named cstr = PS.CSTR(matsource; flowtype = "const. mass")
@named outlet = PS.Port(matsource)

# Connect the flowsheet
eqs = [
    connect(inlet.c, cstr.cv.c1),
    connect(cstr.cv.c2, outlet.c),
    cstr.Q ~ -n_coolant * cₚ_coolant * (cstr.cv.T - T1_coolant) *
             (1 - exp(-UA / (n_coolant * cₚ_coolant)))
]

@named flowsheet_ = ODESystem(eqs, t, [], [], systems = [inlet, cstr, outlet])

pars = []

inp = [
    inlet.T => 297.0,
    inlet.p => 1e5,
    inlet.n => sum(F),
    inlet.xᵢ[1] => F[1] / sum(F),
    inlet.xᵢ[2] => F[2] / sum(F),
    inlet.xᵢ[3] => F[3] / sum(F),
    outlet.p => 1e5
]

u0 = [
    cstr.cv.U => matsource.VT_internal_energy(NaN, T0, [0, 1, 0, 0]) * n0,
    cstr.cv.nᵢ[1, 1] => 0.0,
    cstr.cv.nᵢ[1, 2] => n0,
    cstr.cv.nᵢ[1, 3] => 0.0,
    cstr.cv.nᵢ[1, 4] => 0.0,
    inlet.m => m0_inlet,
    cstr.cv.T => 297.0,
    cstr.cv.c2.n => 0.0,
    outlet.m => 0.0
]

flowsheet, idx = structural_simplify(flowsheet_, (first.(inp), []))

prob = ODEProblem(flowsheet, u0, (0, 2) .* 3600.0, vcat(inp))#; guesses = [outlet.m => -m0_inlet])
sol = solve(prob, QNDF(), abstol = 1e-6, reltol = 1e-6)

(Tmax, iTmax) = findmax(sol[cstr.cv.T])

@test Tmax≈356.149 atol=1e-3
@test sol.t[iTmax]≈2823.24 atol=1e-2

if isinteractive()
    # Plots
    using Plots

    ps = [plot(; framestyle = :box, xlabel = "t / h", xlims = (0.0, 2.0)) for i in 1:2]
    plot!(ps[1], sol.t / 3600, sol[cstr.cv.T], label = "T", ylabel = "T / K")
    [plot!(ps[2], sol.t / 3600, sol[cstr.cv.xᵢ[1, i]];
         label = matsource.components[i], ylabel = "xᵢ / mol/mol") for i in 1:4]
    plot(ps...; layout = (1, 2), size = (800, 400))
end
