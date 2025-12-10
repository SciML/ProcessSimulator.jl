# Reactors

ProcessSimulator.jl provides reactor models for chemical process simulation.

## Continuous Stirred Tank Reactor (CSTR)

The CSTR model represents a well-mixed reactor with chemical reactions.

### Steady-State CSTR

```julia
using ProcessSimulator

# Define reaction system
@named rxn_system = ReactionSystem(
    reactions = [
        Reaction(k_forward, [A], [B])
    ],
    components = ["A", "B"],
    medium = medium
)

# Create CSTR
@named reactor = SteadyStateCSTR(
    medium = medium,
    state = initial_state,
    V = 1.0,              # Volume (m³)
    Q = 0.0,              # Heat input (W)
    reaction_system = rxn_system
)
```

### Features

- **Material Balance**: Component mole balances with reaction terms
- **Energy Balance**: Enthalpy balance with heat input and reaction heat
- **VLE**: Automatic vapor-liquid equilibrium calculations
- **Flexible Reactions**: Support for multiple reactions with various kinetics

## Reaction Systems

Define chemical reactions:

```julia
@named rxn = ReactionSystem(
    reactions = [
        Reaction(1e5, ["A", "B"], ["C"]),        # Forward reaction
        Reaction(1e3, ["C"], ["A", "B"])         # Reverse reaction
    ],
    components = ["A", "B", "C"],
    medium = medium
)
```

### Reaction Kinetics

- Power law kinetics
- Arrhenius temperature dependence
- Custom rate expressions

## Dynamic CSTR

For transient simulations:

```julia
@named reactor = DynamicCSTR(
    medium = medium,
    state = initial_state,
    V = 1.0,
    Q_func = t -> 1000.0 * sin(t),  # Time-varying heat input
    reaction_system = rxn_system
)

# Solve with ODE solver
using OrdinaryDiffEq
prob = ODEProblem(compiled_sys, u0, tspan)
sol = solve(prob, Rodas5())
```

## Example: Exothermic Reaction

```julia
# A → B (exothermic)
components = ["A", "B"]
eos = PR(components)
medium = EoSBased(components = components, eosmodel = eos)

# Feed stream (pure A)
@named feed = Boundary_pTzn(
    medium = medium,
    p = 10e5,
    T = 350.0,
    z = [1.0, 0.0],
    flowrate = 10.0,
    flowbasis = :molar
)

# Reaction: A → B with rate k = 1e5 s⁻¹
@named rxn = ReactionSystem(
    reactions = [Reaction(1e5, ["A"], ["B"])],
    components = components,
    medium = medium
)

# Adiabatic CSTR
initial_state = pTNVState(10e5, 350.0, [0.5, 0.5], base = :Pressure)
@named cstr = SteadyStateCSTR(
    medium = medium,
    state = initial_state,
    V = 0.1,
    Q = 0.0,  # Adiabatic
    reaction_system = rxn
)

# Connect and solve
@named product = ConnHouse(medium = medium)
connections = [
    connect(feed.odesystem.OutPort, cstr.odesystem.InPort),
    connect(cstr.odesystem.OutPort, product.InPort)
]

@named sys = System(connections, t, [], [];
    systems = [feed.odesystem, cstr.odesystem, product])

compiled = mtkcompile(sys)
prob = NonlinearProblem(compiled, guesses(compiled))
sol = solve(prob, NewtonRaphson())

# Results
conversion = (1 - sol[cstr.ControlVolumeState.z[1, 1]]) * 100
T_reactor = sol[cstr.ControlVolumeState.T]
println("Conversion: ", conversion, "%")
println("Reactor Temperature: ", T_reactor, " K")
```
