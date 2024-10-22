using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t

const PS = ProcessSimulator

# Material Parameters
ϱ = [14.8,55.3,13.7,24.7].*1e3           # mol/m^3   
cₚ = [146.440,75.312,192.464,81.588]    # J/mol/K

# Reaction
ν = [-1.0, -1.0, 1.0, 0.0]
Δhᵣ = -83.734392e6                       # J/mol

matsource = PS.MaterialSource(["propylene oxide (A)","water (B)","propylene glycol (C)","methanol (M)"];
    Mw                  = [0.05808, 0.018015, 0.07609, 0.03204],
    molar_density       = (p, T, x; kwargs...) -> sum([ϱ[i]*x[i] for i in 1:4]),
    VT_enthalpy         = (ϱ, T, x) -> sum([cₚ[i]*x[i] for i in 1:4])*T,
    VT_internal_energy  = (ϱ, T, x) -> sum([cₚ[i]*x[i] for i in 1:4])*T,
    reactions           = [PS.Reaction(
        ν               = ν,
        r               = (p,T,x) -> -2.73e-4*exp(9059*(1/297-1/T))*x[1],
        Δhᵣ             = x -> Δhᵣ,
    )],
)

F = [36.3e3, 453.6e3, 0, 45.4e3]./3600  # mol/s
TF = 297.04                             # K

UA = 8440.06                            # J/K/s
T_coolant = 288.71                      # K
cₚ_coolant = 75.312                     # J/mol/K   


# Create flowsheet
@named inlet = PS.Port(matsource)
@named cstr = PS.CSTR(matsource)
@named outlet = PS.Port(matsource)

# Connect the flowsheet
eqs = [
    connect(inlet.c,cstr.cv.c1),
    connect(cstr.cv.c2,outlet.c),
]

@named flowsheet_ = ODESystem(eqs, t, [], [], systems=[inlet,cstr,outlet])

pars = [
    cstr.cv.V => 1.89,
]

inp = [
    inlet.T => 297.04,
    inlet.p => 1e5,
    inlet.n => sum(F),
    inlet.xᵢ[1] => F[1]/sum(F),
    inlet.xᵢ[2] => F[2]/sum(F),
    inlet.xᵢ[3] => F[3]/sum(F),
    cstr.cv.q1.Q => 0.0,
    outlet.p => 1e5,
]
out = [
    cstr.cv.T,
    cstr.cv.xᵢ[1,1],
    cstr.cv.xᵢ[1,2],
    cstr.cv.xᵢ[1,3],
]

u0 = [
    cstr.cv.nᵢ[1,1] => 0.0,
    cstr.cv.nᵢ[1,2] => 55.3e3*Dict(pars...)[cstr.cv.V],
    cstr.cv.nᵢ[1,3] => 0.0,
    cstr.cv.nᵢ[1,4] => 0.0,
    cstr.cv.T => 297.,
    cstr.cv.ΔnR[1] => 0.0,
    cstr.cv.ΔnR[2] => 0.0,
    cstr.cv.ΔnR[3] => 0.0,
]

flowsheet,idx = structural_simplify(flowsheet_,(first.(inp),out))

prob = ODEProblem(flowsheet, u0, (0.0, 1.4e4), collect(Iterators.flatten([pars,inp])))