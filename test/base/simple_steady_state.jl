using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t

const PS = ProcessSimulator

# Create material sources
Mw = 0.004          # kg/mol
R = 2.1e3*Mw        # J/(mol K)
cₚ = 5.2e3*Mw       # J/(mol K)
cᵥ = cₚ - R

# Perfect gas equations
matsource = PS.MaterialSource("helium";
    Mw = 0.004,
    molar_density = (p, T, x; kwargs...) -> p/(R*Mw*T),
    VT_internal_energy = (ϱ, T, x) -> cᵥ*T,
    VT_enthalpy = (ϱ, T, x) -> cₚ*T,
    VT_entropy = (ϱ, T, x) -> cᵥ*log(T) + R*log(1/ϱ)
)

# Compressor
@named inlet = PS.Port(matsource)
@named comp = PS.SimpleAdiabaticCompressor(matsource)
@named outlet = PS.Port(matsource)

con_eqs_comp = [
    connect(inlet.c, comp.cv.c1),
    connect(comp.cv.c2, outlet.c)
]
@named compressor_ = ODESystem(con_eqs_comp, t, [], []; systems = [inlet, comp, outlet])

inp_comp = [
    inlet.m => 1.0,
    inlet.T => 300.0,
    inlet.p => 1e5,
    outlet.p => 1e6
]
out_comp = [comp.W]
compressor, idx = structural_simplify(compressor_, (first.(inp_comp), out_comp))

u0_comp = [
    (u => 1.0 for u in unknowns(compressor))...,
    comp.ηᴱ => 0.8
]

prob_comp = SteadyStateProblem(compressor, inp_comp, u0_comp)
sol_comp = solve(prob_comp)

@test sol_comp[comp.W] ≈ 2.3933999e6 rtol=1e-5

# Create flowsheet
@named comp_12 = PS.SimpleAdiabaticCompressor(matsource)
@named s2 = PS.MaterialStream(matsource)
@named heat_22⁺ = PS.SimpleIsobaricHeatExchanger(matsource)
@named s2⁺ = PS.MaterialStream(matsource)
@named heat_2⁺3 = PS.SimpleIsobaricHeatExchanger(matsource)
@named s3 = PS.MaterialStream(matsource)
@named turb_34 = PS.SimpleAdiabaticCompressor(matsource)
@named s4 = PS.MaterialStream(matsource)
@named heat_44⁺ = PS.SimpleIsobaricHeatExchanger(matsource)
@named s4⁺ = PS.MaterialStream(matsource)
@named heat_4⁺1 = PS.SimpleIsobaricHeatExchanger(matsource)

systems = [comp_12, heat_22⁺, heat_2⁺3, turb_34, heat_44⁺, heat_4⁺1]
streams = [inlet, s2, s2⁺, s3, s4, s4⁺, outlet]

con_eqs = vcat(
    connect(inlet.c, comp_12.cv.c1),
    connect(comp_12.cv.c2, s2.c1),
    connect(s2.c2, heat_22⁺.cv.c1),
    connect(heat_22⁺.cv.c2, s2⁺.c1),
    connect(s2⁺.c2, heat_2⁺3.cv.c1),
    connect(heat_2⁺3.cv.c2, s3.c1),
    connect(s3.c2, turb_34.cv.c1),
    connect(turb_34.cv.c2, s4.c1),
    connect(s4.c2, heat_44⁺.cv.c1),
    connect(heat_44⁺.cv.c2, s4⁺.c1),
    connect(s4⁺.c2, heat_4⁺1.cv.c1),
    connect(heat_4⁺1.cv.c2, outlet.c)
)
append!(con_eqs, [0.0 ~ heat_2⁺3.Q + heat_4⁺1.Q])

@named flowsheet_ = ODESystem(con_eqs, t, [], []; systems = [systems..., streams...])

inp = [
    inlet.n => 1.0,
    inlet.T => 298.15,
    inlet.p => 1e5,
    s2.p => 1e6,
    s2⁺.T => 298.15,
    s3.T => 253.15,
    outlet.T => 298.15,
    outlet.p => 1e5
]
out = [
    heat_44⁺.Q,
    comp_12.W,
    turb_34.W
]

flowsheet, idx = structural_simplify(flowsheet_, (first.(inp), out))

u0 = [u => 1.0 for u in unknowns(flowsheet)]

prob = SteadyStateProblem(flowsheet, inp, u0)
sol = solve(prob)

ε_KM = abs(sol[heat_44⁺.Q])/abs(sol[turb_34.W] + sol[comp_12.W])

@test round(sol[s2.T], digits = 2) ≈ 755.58
@test round(sol[s4.T], digits = 2) ≈ 99.89
@test round(ε_KM, digits = 3) ≈ 0.504
