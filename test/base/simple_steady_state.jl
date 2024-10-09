using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t

const PS = ProcessSimulator

# Create material sources
Mw = 0.004          # kg/mol
R = 2.1e3*Mw        # J/(mol K)
cₚ = 5.2e3*Mw       # J/(mol K)
cᵥ = cₚ - R

matsource = PS.MaterialSource(
    "helium",
    ["helium"],[0.004],
    (ϱ,T,n) -> ϱ*R*T, 
    (p,T,n;kwargs...) -> p/(R*Mw*T),
    (ϱ,T,n) -> cᵥ*T,
    (ϱ,T,n) -> cₚ*T,
    (ϱ,T,n) -> cᵥ*log(T) + R*log(1/ϱ),
    (p,T,n) -> NaN
)

# Create flowsheet
# Components
comp_12 = PS.SimpleAdiabaticCompressor(matsource; name=:comp_12)
heat_22⁺ = PS.SimpleIsobaricHeatExchanger(matsource; name=:heat_22⁺)  
heat_2⁺3 = PS.SimpleIsobaricHeatExchanger(matsource; name=:heat_2⁺2)
turb_34 = PS.SimpleAdiabaticCompressor(matsource; name=:turb_34)
heat_44⁺ = PS.SimpleIsobaricHeatExchanger(matsource; name=:heat_44⁺)
heat_4⁺1 = PS.SimpleIsobaricHeatExchanger(matsource; name=:heat_4⁺1)
sys = [comp_12,heat_22⁺,heat_2⁺3,turb_34,heat_44⁺,heat_4⁺1]

# Connections
con_eqs = Equation[]
append!(con_eqs,PS.connect_states(comp_12.cv.s2,heat_22⁺.cv.s1; is_stream=true))
append!(con_eqs,PS.connect_states(heat_22⁺.cv.s2,heat_2⁺3.cv.s1; is_stream=true))
append!(con_eqs,PS.connect_states(heat_2⁺3.cv.s2,turb_34.cv.s1; is_stream=true))
append!(con_eqs,PS.connect_states(turb_34.cv.s2,heat_44⁺.cv.s1; is_stream=true))
append!(con_eqs,PS.connect_states(heat_44⁺.cv.s2,heat_4⁺1.cv.s1; is_stream=true))

# Additional connections (internal heat exchanger)
append!(con_eqs,[heat_22⁺.cv.Qs[1] ~ -heat_4⁺1.cv.Qs[1]])

@named flowsheet_ = ODESystem(con_eqs, t, [], []; systems=sys)

giv = [
    comp_12.cv.s1.nᵢ[1] => 1.0,         # T1
    comp_12.cv.s1.T => 298.15,          # T1
    comp_12.cv.s1.p => 1e5,             # p1
    comp_12.cv.s2.p => 1e6,             # p2
    heat_22⁺.cv.s2.T => 298.15,         # T2⁺ 
    turb_34.cv.s1.T => 253.15,          # T3
    turb_34.cv.s2.p => 1e5,             # p4
    heat_44⁺.cv.s2.T => 253.15,         # T4⁺
    comp_12.ηᴱ => 1.0,                  # ηᴱ
    turb_34.ηᴱ => 1.0                   # ηᴱ
]
unk = [
    heat_44⁺.cv.Qs[1],
    comp_12.cv.Ws[1],
    turb_34.cv.Ws[1]
]

flowsheet,idx = structural_simplify(flowsheet_, (first.(giv), first.(unk)))

u0 = [
    comp_12.cv.s2.T => 1.0,
    comp_12.cv_s.s2.T => 1.0,
    (comp_12.cv.Ws)[1] => -1.0,
    turb_34.cv.s2.T => 1.0,
    turb_34.cv_s.s2.T => 1.0,
    (turb_34.cv.Ws)[1] => 1.0,
    heat_4⁺1.cv.s2.T => 1.0,
]

prob = SteadyStateProblem(flowsheet,giv,u0)
sol = solve(prob)

ε_KM = abs(sol[heat_44⁺.cv.Qs[1]])/abs(sol[turb_34.cv.Ws[1]] + sol[comp_12.cv.Ws[1]])

(T₁,T₂,T₂₊,T₃,T₄,T₄₊) = (
    sol.ps[comp_12.cv.s1.T],
    sol[comp_12.cv.s2.T],
    sol.ps[heat_22⁺.cv.s2.T],
    sol.ps[turb_34.cv.s1.T],
    sol[turb_34.cv.s2.T],
    sol.ps[heat_44⁺.cv.s2.T]
)

@test round(T₂,digits=2) ≈ 755.58
@test round(T₄,digits=2) ≈ 99.89
@test round(ε_KM,digits=3) ≈ 0.504