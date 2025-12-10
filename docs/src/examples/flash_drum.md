# Flash Drum Example

This example demonstrates a complete vapor-liquid flash separation based on the EMSO manual example (Section 3.2.4).

## Problem Statement

Separate a hydrocarbon mixture in a flash drum:

- **Feed**: 496.3 kmol/h at 338 K and 507.1 kPa
- **Flash Conditions**: 2.5 atm (253 kPa), adiabatic
- **Components**: 1,3-butadiene, isobutene, n-pentane, 1-pentene, 1-hexene, benzene

The goal is to calculate the vapor and liquid product compositions, flow rates, and flash temperature.

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

Create the equation of state model:

```julia
components = [
    "1,3-butadiene",
    "isobutene",
    "n-pentane",
    "1-pentene",
    "1-hexene",
    "benzene"
]

# Soave-Redlich-Kwong EoS with Reid ideal gas
model = SRK(components, idealmodel = ReidIdeal)
medium = EoSBased(components = components, eosmodel = model)
```

**Why SRK?** The Soave-Redlich-Kwong equation of state provides good accuracy for light hydrocarbons at moderate pressures.

## Define Feed Stream

Create the feed boundary condition:

```julia
@named S1 = Boundary_pTzn(
    medium = medium,
    p = 507.1e3,    # 507.1 kPa
    T = 338.0,      # 338 K
    z = [0.2379, 0.3082, 0.09959, 0.1373, 0.08872, 0.1283],
    flowrate = 496.3/3600,  # Convert kmol/h to kmol/s
    flowbasis = :molar
)
```

The feed composition breakdown:
- 1,3-butadiene: 23.79%
- Isobutene: 30.82%
- n-pentane: 9.96%
- 1-pentene: 13.73%
- 1-hexene: 8.87%
- Benzene: 12.83%

## Create Flash Drum

Define the initial state guess:

```julia
flash_state = pTNVState(
    2.5e5,      # Flash pressure (Pa)
    315.87,     # Initial temperature guess (K)
    [0.2379, 0.3082, 0.09959, 0.1373, 0.08872, 0.1283],
    base = :Pressure
)
```

Create the flash drum:

```julia
@named FL1 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = flash_state,
    pressure = flash_state.p,
    Q = nothing  # Adiabatic operation
)
```

**Note**: The composition guess uses the feed composition. In a larger flowsheet, this would depend on upstream conditions.

## Define Product Streams

Create connection houses for products:

```julia
@named LiquidPort = ConnHouse(medium = medium)
@named VaporPort = ConnHouse(medium = medium)
```

## Connect Components

Build the flowsheet connections:

```julia
connection_set = [
    connect(S1.odesystem.OutPort, FL1.odesystem.InPort),
    connect(FL1.odesystem.LiquidOutPort, LiquidPort.InPort),
    connect(FL1.odesystem.VaporOutPort, VaporPort.InPort)
]
```

Flowsheet diagram:
```
Feed (S1) → Flash Drum (FL1) → Vapor (VaporPort)
                              ↓
                           Liquid (LiquidPort)
```

## Build System

Create the complete system:

```julia
@named sys = System(connection_set, t, [], [];
    systems = [S1.odesystem, FL1.odesystem, LiquidPort, VaporPort])

SteadyStateFlash = mtkcompile(sys)
```

## Set Initial Guesses

Customize initial guesses for better convergence:

```julia
default_guesses = guesses(SteadyStateFlash)
guesses_Flash = Dict(
    FL1.odesystem.V[2] => 1.0/FL1.medium.Guesses.ρ[2],
    FL1.odesystem.V[3] => 1.0/FL1.medium.Guesses.ρ[3],
    FL1.odesystem.nᴸⱽ[2] => 1.0,  # Vapor holdup guess
)
merged_guesses = merge(default_guesses, guesses_Flash)
```

**Tip**: Volume guesses are based on estimated densities from the thermodynamic flash.

## Solve

Create and solve the nonlinear problem:

```julia
prob = NonlinearProblem(SteadyStateFlash, merged_guesses, use_scc = true)
sol = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()))
```

The `use_scc = true` option enables strongly connected components analysis for better performance on large systems.

## Extract Results

### Flash Drum Conditions

```julia
T_flash = sol[SteadyStateFlash.FL1.ControlVolumeState.T]
p_flash = sol[SteadyStateFlash.FL1.ControlVolumeState.p]
Q_flash = sol[SteadyStateFlash.FL1.Q]

println("Flash Temperature: ", T_flash, " K")
println("Flash Pressure: ", p_flash/1e5, " bar")
println("Heat Duty: ", Q_flash, " W")
```

### Product Compositions

```julia
# Vapor composition
z_vapor = [sol[SteadyStateFlash.VaporPort.p_T_z_n.z[i]] for i in 1:6]
n_vapor = sol[SteadyStateFlash.VaporPort.p_T_z_n.n]

# Liquid composition
z_liquid = [sol[SteadyStateFlash.LiquidPort.p_T_z_n.z[i]] for i in 1:6]
n_liquid = sol[SteadyStateFlash.LiquidPort.p_T_z_n.n]

# Vapor fraction
β = n_vapor / (n_vapor + n_liquid)

println("\nVapor Fraction: ", β)
println("\nVapor Composition:")
for (i, comp) in enumerate(components)
    println("  ", comp, ": ", z_vapor[i])
end

println("\nLiquid Composition:")
for (i, comp) in enumerate(components)
    println("  ", comp, ": ", z_liquid[i])
end
```

### Flow Rates

```julia
println("\nFlow Rates:")
println("  Vapor: ", n_vapor * 3600, " kmol/h")
println("  Liquid: ", n_liquid * 3600, " kmol/h")
```

## Expected Results

For this system, you should obtain:

- **Flash Temperature**: ~315-316 K (adiabatic cooling from feed)
- **Vapor Fraction**: Higher for lighter components (butadiene, isobutene)
- **Liquid Fraction**: Higher for heavier components (hexene, benzene)

The light components (butadiene, isobutene) will be enriched in the vapor phase, while heavier components (hexene, benzene) will concentrate in the liquid phase.

## Validation

Compare with EMSO manual results (Section 3.2.4):

```julia
# Expected values from EMSO
emso_T_flash = 315.87  # K
emso_vapor_fraction = 0.xx  # Replace with actual value

error_T = abs(T_flash - emso_T_flash) / emso_T_flash * 100
println("\nTemperature Error: ", error_T, "%")
```

## Troubleshooting

If convergence fails:

1. **Adjust Temperature Guess**: Try temperatures between feed and flash conditions
2. **Change Solver**: Use `FastShortcutNonlinearPolyalg()` for automatic solver selection
3. **Check EoS Parameters**: Ensure Clapeyron has parameters for all components
4. **Verify Feed Conditions**: Ensure feed is subcooled or two-phase capable at flash pressure

## Extensions

### Add Heat Exchanger

Pre-cool feed before flash:

```julia
@named cooler = HeatExchanger(
    medium = medium,
    Q = -50000.0  # 50 kW cooling
)
```

### Multi-Stage Flash

Create a flash train:

```julia
# First flash at high pressure
@named FL1 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state1,
    pressure = 5e5
)

# Second flash at lower pressure
@named FL2 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state2,
    pressure = 2.5e5
)

# Connect liquid from FL1 to FL2
connections = [
    connect(feed.odesystem.OutPort, FL1.odesystem.InPort),
    connect(FL1.odesystem.LiquidOutPort, FL2.odesystem.InPort)
]
```

## Complete Code

The full code is available in `test/separation/flash_drum_test.jl`.
