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