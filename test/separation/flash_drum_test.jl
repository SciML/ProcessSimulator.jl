using ProcessSimulator
using Clapeyron
using ModelingToolkit
using OrdinaryDiffEq
using NonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

# ==================== Steady-State Flash Drum Test ====================
# Based on EMSO manual example (Section 3.2.4)
# Feed: 496.3 kmol/h at 338 K and 507.1 kPa
# Flash conditions: 2.5 atm, 315.87 K

components = ["1,3-butadiene", "isobutene", "n-pentane", "1-pentene", "1-hexene", "benzene"]
model = PR(components, idealmodel = ReidIdeal)
medium = EoSBased(components = components, eosmodel = model)

@named S1 = Boundary_pTzn(
    medium = medium,
    p = 507.1e3,
    T = 338.0,
    z = [0.2379, 0.3082, 0.09959, 0.1373, 0.08872, 0.1283],
    flowrate = 496.3/3600,
    flowbasis = :molar
)

flash_state = pTNVState(2.5*101325.0, 315.87, [0.2379, 0.3082, 0.09959, 0.1373, 0.08872, 0.1283], base = :Pressure)

@named FL1 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = flash_state,
    pressure = 2.5*101325.0,
    Q = 0.0
)

@named LiquidPort = ConnHouse(medium = medium)
@named VaporPort = ConnHouse(medium = medium)

connection_set = [
    connect(S1.odesystem.OutPort, FL1.odesystem.InPort),
    connect(FL1.odesystem.LiquidOutPort, LiquidPort.InPort),
    connect(FL1.odesystem.VaporOutPort, VaporPort.InPort)
]

@named sys = System(connection_set, t, [], []; 
    systems = [S1.odesystem, FL1.odesystem, LiquidPort, VaporPort])

SteadyStateFlash = mtkcompile(sys)

default_guesses = guesses(SteadyStateFlash)
guesses_Flash = Dict(
    FL1.odesystem.VaporOutPort.ṅ[1] => 396.3/3600,
    FL1.odesystem.LiquidOutPort.ṅ[1] => 100.0/3600,
)
merged_guesses = merge(default_guesses, guesses_Flash)

prob = NonlinearProblem(SteadyStateFlash, merged_guesses)
@time sol = solve(prob, abstol = 1e-8, reltol = 1e-8)

