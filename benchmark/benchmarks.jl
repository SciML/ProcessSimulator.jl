using ProcessSimulator, BenchmarkTools
using ModelingToolkit
using ModelingToolkit: t_nounits as t

const SUITE = BenchmarkGroup()
const PS = ProcessSimulator

# Material Parameters
ϱ = [14.8, 55.3, 13.7, 24.7] .* 1.0e3          # mol/m^3
cₚ = [35, 18, 46, 19.5] * 4.184              # J/mol/K

ν = [-1.0, -1.0, 1.0, 0.0]
Δhᵣ = 20013 * 4.184                          # J/mol

matsource = PS.MaterialSource(
    ["propylene oxide (A)", "water (B)", "propylene glycol (C)", "methanol (M)"];
    Mw = [0.05808, 0.01801, 0.07609, 0.03204],
    molar_density = (p, T, x; kwargs...) -> sum([ϱ[i] * x[i] for i in 1:4]),
    VT_enthalpy = (ϱ, T, x) -> sum([cₚ[i] * x[i] for i in 1:4]) * T,
    VT_internal_energy = (ϱ, T, x) -> sum([cₚ[i] * x[i] for i in 1:4]) * T,
    reactions = [
        PS.Reaction(
            ν = ν,
            r = (p, T, x) -> 2.73e-4 * exp(9059 * (1 / 297 - 1 / T)) * x[1],
            Δhᵣ = (T) -> Δhᵣ,
        ),
    ]
)

UA = 7262 * 4184 / 3600
T1_coolant = 288.71
cₚ_coolant = 18 * 4.184
n_coolant = 453.6e3 / 3600

# =============================================================================
# Component construction
# =============================================================================

SUITE["components"] = BenchmarkGroup()

SUITE["components"]["port"] = @benchmarkable PS.Port($matsource; name = :inlet)
SUITE["components"]["cstr"] = @benchmarkable PS.CSTR(
    $matsource; flowtype = "const. mass", name = :reactor
)

# =============================================================================
# Flowsheet assembly + system compilation (the representative workload:
# MTK structural_simplify on a connected reactor flowsheet)
# =============================================================================

F = [36.3e3, 453.6e3, 0, 45.4e3] ./ 3600  # mol/s
x_frac = F ./ sum(F)

function build_flowsheet()
    @named inlet = PS.Port(matsource)
    @named cstr = PS.CSTR(matsource; flowtype = "const. mass")
    @named outlet = PS.Port(matsource)
    eqs = [
        connect(inlet.c, cstr.cv.c1),
        connect(cstr.cv.c2, outlet.c),
        cstr.Q ~
            -n_coolant * cₚ_coolant * (cstr.cv.T - T1_coolant) *
            (1 - exp(-UA / (n_coolant * cₚ_coolant))),
    ]
    return @named flowsheet = ODESystem(
        eqs, t, [], []; systems = [inlet, cstr, outlet]
    )
end

@named inlet = PS.Port(matsource)
@named cstr = PS.CSTR(matsource; flowtype = "const. mass")
@named outlet = PS.Port(matsource)
eqs = [
    connect(inlet.c, cstr.cv.c1),
    connect(cstr.cv.c2, outlet.c),
    cstr.Q ~
        -n_coolant * cₚ_coolant * (cstr.cv.T - T1_coolant) *
        (1 - exp(-UA / (n_coolant * cₚ_coolant))),
]
@named flowsheet_ = ODESystem(eqs, t, [], []; systems = [inlet, cstr, outlet])

inp = [
    inlet.T => 297.0,
    inlet.p => 1.0e5,
    inlet.n => sum(F),
    inlet.xᵢ[1] => x_frac[1],
    inlet.xᵢ[2] => x_frac[2],
    inlet.xᵢ[3] => x_frac[3],
    inlet.xᵢ[4] => x_frac[4],
    outlet.p => 1.0e5,
]

SUITE["flowsheet"] = BenchmarkGroup()

SUITE["flowsheet"]["build"] = @benchmarkable build_flowsheet()
SUITE["flowsheet"]["simplify"] = @benchmarkable structural_simplify(
    $flowsheet_, ($([first(i) for i in inp]), [])
)
