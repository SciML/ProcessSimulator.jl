using ProcessSimulator
using Clapeyron
using ModelingToolkit
using OrdinaryDiffEq
using NonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

# ==================== Splitter Test ====================
# Test that a splitter correctly divides flow according to split ratio
# and maintains energy and mass balance

components = ["methane", "ethane", "propane"]
model = SRK(components, idealmodel = ReidIdeal)
medium = EoSBased(components = components, eosmodel = model)

# Inlet boundary - vapor phase feed
@named FEED1 = Boundary_pTzn(
    medium = medium,
    p = 30e5,
    T = 300.0,
    z = [0.5, 0.3, 0.2],
    flowrate = 100.0,  # 1 mol/s
    flowbasis = :molar
)

# Splitter state guess
splitter_state = pTNVState(30e5, 310.0, [0.6, 0.3, 0.1], base = :Pressure)

# Split ratio: 0.6 means 60% goes to OutPort, 40% goes to OutPort2
@named SP1 = Splitter(
    medium = medium,
    state = splitter_state,
    split_ratio = 0.7)

# Fixed pressure outlets
@named P1 = ConnHouse(medium = medium)
@named P2 = ConnHouse(medium = medium)

connection_set = [
    connect(FEED1.odesystem.OutPort, SP1.odesystem.InPort),
    connect(SP1.odesystem.OutPort, P1.InPort),
    connect(SP1.odesystem.OutPort2, P2.InPort)
]

@named sys = System(connection_set, t, [], []; 
    systems = [FEED1.odesystem, SP1.odesystem, P1, P2])
ModelingToolkit.flatten_equations(equations(expand_connections(sys)))
equations(expand_connections(sys))
for eq in ModelingToolkit.flatten_equations(equations(expand_connections(sys)))
    println(eq)
end

SplitterSystem = mtkcompile(sys)

default_guesses = guesses(SplitterSystem)
guesses_Flash = Dict(
    SP1.odesystem.OutPort.z[1, 1] => 0.3,
    SP1.odesystem.OutPort.z[2, 1] => 0.3,
   SP1.odesystem.OutPort.z[3, 1] => 0.4,
   SP1.odesystem.OutPort2.z[3, 1] => 0.4
)
guess = merge(default_guesses, guesses_Flash)

prob = NonlinearProblem(SplitterSystem, guess, use_scc = true)
@time sol = solve(prob, NewtonRaphson(autodiff = AutoFiniteDiff()), abstol = 1e-6, reltol = 1e-6)

# Verify results
inlet_flow = sol[SplitterSystem.S1.n_dot]
outlet1_flow = abs(sol[SplitterSystem.P1.Port.ṅ[1]])
outlet2_flow = abs(sol[SplitterSystem.P2.Port.ṅ[1]])

split_ratio = 0.6

# Test flow split (within 1% tolerance)
@test abs(outlet1_flow - inlet_flow * split_ratio) / inlet_flow < 0.01
@test abs(outlet2_flow - inlet_flow * (1 - split_ratio)) / inlet_flow < 0.01

# Test overall mass balance
@test abs(inlet_flow - outlet1_flow - outlet2_flow) / inlet_flow < 0.01

# Test that compositions remain the same
inlet_z = [0.5, 0.3, 0.2]
for i in 1:length(components)
    outlet1_z = sol[SplitterSystem.P1.Port.z[i, 1]]
    outlet2_z = sol[SplitterSystem.P2.Port.z[i, 1]]
    @test abs(outlet1_z - inlet_z[i]) < 0.01
    @test abs(outlet2_z - inlet_z[i]) < 0.01
end

println("Splitter Test Passed!")
println("Inlet flow: ", inlet_flow, " mol/s")
println("Outlet 1 flow: ", outlet1_flow, " mol/s")
println("Outlet 2 flow: ", outlet2_flow, " mol/s")
println("Split ratio: ", split_ratio)
