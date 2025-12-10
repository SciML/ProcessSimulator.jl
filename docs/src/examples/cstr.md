# CSTR Example

This example demonstrates a continuous stirred tank reactor (CSTR) with chemical reaction.

## Problem Statement

Produce ethylene glycol by reacting ethylene oxide with water in an adiabatic CSTR:

**Reaction:**
```
C₂H₄O + H₂O → C₂H₆O₂
(ethylene oxide + water → ethylene glycol)
```

**Feed Conditions:**
- Flow rate: 100 mol/s
- Pressure: 5 atm (505 kPa)
- Temperature: 350.15 K
- Composition: 40% ethylene oxide, 60% water

**Reactor:**
- Volume: 1.0 m³
- Operation: Adiabatic (Q = 0)
- Limiting reactant: Ethylene oxide

## Setup

Import required packages:

```julia
using ProcessSimulator
using Clapeyron
using ModelingToolkit
using NonlinearSolve
using ModelingToolkit: t_nounits as t
```

## Define Thermodynamics

This system requires a composite model for accurate liquid-phase activity:

```julia
components = ["ethylene oxide", "water", "ethylene glycol"]

# Composite ideal model for reference states
idealmodel = CompositeModel(components,
    liquid = RackettLiquid,
    gas = ReidIdeal(components, reference_state = :formation),
    saturation = DIPPR101Sat,
    hvap = DIPPR106HVap
)

# NRTL activity coefficient model for liquid phase
activity = NRTL(components)

# Complete composite model
model = CompositeModel(components, fluid = idealmodel, liquid = activity)

# Create medium
medium = EoSBased(components = components, eosmodel = model)
```

**Why Composite Model?**
- **Ideal Gas**: Reid ideal gas for vapor phase
- **Liquid Activity**: NRTL model for non-ideal liquid interactions
- **Properties**: DIPPR correlations for saturation and heat of vaporization

## Verify Thermodynamics

Test the model before simulation:

```julia
# Check bubble pressure calculation
p_bubble = bubble_pressure(model, 350.15, [0.4, 0.6, 0.0])
println("Bubble pressure at 350.15 K: ", p_bubble[1]/1e5, " bar")
```

## Define Feed Stream

Create the feed boundary:

```julia
@named S1 = Boundary_pTzn(
    medium = medium,
    p = 5 * 101325.0,    # 5 atm
    T = 350.15,          # 350.15 K
    z = [0.4, 0.6, 0.0], # No glycol in feed
    flowrate = 100.0,
    flowbasis = :molar
)
```

## Define Reaction

Create the power law reaction:

```julia
rxn1 = PowerLawReaction(
    components = components,
    stoichiometry = Dict(
        "ethylene oxide" => -1.0,
        "water" => -1.0,
        "ethylene glycol" => 1.0
    ),
    order = Dict(
        "ethylene oxide" => 1.0,
        "water" => 0.0
    ),
    A = 0.5,    # Pre-exponential factor
    Eₐ = 0.0    # Activation energy (J/mol)
)
```

**Reaction Kinetics:**
- Rate = A × [ethylene oxide]¹ × [water]⁰
- First order in ethylene oxide, zero order in water
- No temperature dependence (Eₐ = 0)

## Create Reactor

Define initial state and reactor:

```julia
# Initial state guess
reactor_state = pTNVState(
    5 * 101325.0,           # Pressure (Pa)
    350.15,                 # Temperature (K)
    ones(3) / 3.0,          # Equal composition guess
    base = :Pressure
)

# CSTR
@named R1 = FixedVolumeSteadyStateCSTR(
    medium = medium,
    reactionset = rxn1,
    limiting_reactant = "ethylene oxide",
    state = reactor_state,
    volume = 1.0,    # 1 m³
    W = 0.0,         # No shaft work
    Q = nothing      # Adiabatic
)
```

**Note**: `limiting_reactant` helps the solver by identifying which reactant determines conversion.

## Define Product Stream

```julia
@named sink = ConnHouse(medium = medium)
```

## Connect Components

```julia
connection_set = [
    connect(S1.odesystem.OutPort, R1.odesystem.InPort),
    connect(R1.odesystem.OutPort, sink.InPort)
]
```

Flowsheet diagram:
```
Feed (S1) → CSTR (R1) → Product (sink)
```

## Build System

```julia
@named sys = System(connection_set, t, [], [];
    systems = [S1.odesystem, R1.odesystem, sink])

AdiabaticVolumeReactor = mtkcompile(sys)
```

## Set Initial Guesses

Customize guesses for better convergence:

```julia
default_guesses = guesses(AdiabaticVolumeReactor)
guesses_Reactor = Dict(
    R1.odesystem.V[3] => 1e-8,          # Small vapor volume
    R1.odesystem.OutPort.ṅ[1] => 100.0  # Outlet flow rate
)
merged_guesses = merge(default_guesses, guesses_Reactor)
```

## Solve

```julia
prob = NonlinearProblem(AdiabaticVolumeReactor, merged_guesses)
@time sol = solve(prob, abstol = 1e-8, reltol = 1e-8)
```

## Extract Results

### Reactor Conditions

```julia
T_reactor = sol[AdiabaticVolumeReactor.R1.ControlVolumeState.T]
p_reactor = sol[AdiabaticVolumeReactor.R1.ControlVolumeState.p]
Q_reactor = sol[AdiabaticVolumeReactor.R1.Q]

println("Reactor Temperature: ", T_reactor, " K")
println("Reactor Pressure: ", p_reactor/1e5, " bar")
println("Heat Duty: ", Q_reactor)
```

### Conversion

```julia
# Feed and product moles of ethylene oxide

conversion = sol[AdiabaticVolumeReactor.R1.X]*100.0

println("\nEthylene Oxide Conversion: ", round(conversion, digits = 3), "%")
```

### Reaction Rate

```julia
r_reaction = sol[AdiabaticVolumeReactor.R1.r[1]]  # Reaction rate (mol/s)
println("Reaction Rate: ", r_reaction, " mol/s")
```

## Print Summary

Use the built-in summary function:

```julia
# Method 1: Using unit_ops dictionary
unit_ops = Dict(:R1 => :CSTR, :S1 => :Feed)
print_flowsheet_summary(sol, AdiabaticVolumeReactor, unit_ops, components)

# Method 2: Direct component specification
print_flowsheet_summary(sol, AdiabaticVolumeReactor, components, R1, S1)
```

This provides a formatted table with:
- Stream compositions
- Flow rates
- Temperatures and pressures
- Phase information

### CSTR Series

Chain multiple reactors:

```julia
@named R1 = FixedVolumeSteadyStateCSTR(...)
@named R2 = FixedVolumeSteadyStateCSTR(...)

connections = [
    connect(feed.odesystem.OutPort, R1.odesystem.InPort),
    connect(R1.odesystem.OutPort, R2.odesystem.InPort),
    connect(R2.odesystem.OutPort, product.InPort)
]
```

## Troubleshooting

### Convergence Issues

1. **Adjust Initial Guesses**: Try different temperature/composition guesses
2. **Check Reaction Rate**: Ensure kinetics are reasonable (not too fast/slow)
3. **Verify Thermodynamics**: Test `bubble_pressure` and `flash` separately
4. **Use Different Solver**: Try `FastShortcutNonlinearPolyalg()`

### Physical Constraints

- **Positive Compositions**: Ensure all z[i] ≥ 0
- **Temperature Limits**: Check for unreasonable temperature rise
- **Phase Stability**: Verify reactor operates in single phase if assumed

## Complete Code

The full code is available in `test/reactors/ss_cstr_test.jl`.
