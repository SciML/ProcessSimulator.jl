# Separation Units

ProcessSimulator.jl provides models for vapor-liquid separation equipment.

## Flash Drums

Flash drums perform vapor-liquid separation at specified conditions.

### Fixed Pressure Steady-State Flash Drum

The most common flash drum configuration maintains constant pressure:

```julia
using ProcessSimulator
using Clapeyron

# Define components and thermodynamics
components = ["methane", "ethane", "propane"]
eos = PR(components)
medium = EoSBased(components = components, eosmodel = eos)

# Feed stream
@named feed = Boundary_pTzn(
    medium = medium,
    p = 30e5,           # 30 bar
    T = 300.0,          # 300 K
    z = [0.3, 0.3, 0.4],
    flowrate = 100.0,
    flowbasis = :molar
)

# Flash drum at 10 bar
initial_state = pTNVState(10e5, 280.0, [0.3, 0.3, 0.4], base = :Pressure)
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = initial_state,
    p = 10e5,
    Q = 0.0  # Adiabatic
)

# Product streams
@named vapor = ConnHouse(medium = medium)
@named liquid = ConnHouse(medium = medium)

# Connect
connections = [
    connect(feed.odesystem.OutPort, flash.odesystem.InPort),
    connect(flash.odesystem.VaporOutPort, vapor.InPort),
    connect(flash.odesystem.LiquidOutPort, liquid.InPort)
]

# Build system
@named sys = System(connections, t, [], [];
    systems = [feed.odesystem, flash.odesystem, vapor, liquid])

# Solve
compiled = mtkcompile(sys)
prob = NonlinearProblem(compiled, guesses(compiled))
sol = solve(prob, NewtonRaphson())
```

### Flash Drum Features

- **Three-Port Control Volume**: Feed, vapor product, liquid product
- **VLE Calculations**: Automatic phase equilibrium using Clapeyron.jl
- **Energy Balance**: Adiabatic or with heat input/removal
- **Material Balance**: Component mole balances
- **Flexible Initial Guesses**: Multiple state specification options

## Accessing Results

Extract flash drum results from the solution:

```julia
# Temperatures
T_flash = sol[flash.ControlVolumeState.T]
T_vapor = sol[vapor.p_T_z_n.T]
T_liquid = sol[liquid.p_T_z_n.T]

# Pressures
p_flash = sol[flash.ControlVolumeState.p]

# Compositions
z_vapor = [sol[vapor.p_T_z_n.z[i]] for i in 1:3]
z_liquid = [sol[liquid.p_T_z_n.z[i]] for i in 1:3]

# Flow rates
n_vapor = sol[vapor.p_T_z_n.n]
n_liquid = sol[liquid.p_T_z_n.n]

# Vapor fraction
vapor_fraction = n_vapor / (n_vapor + n_liquid)

println("Flash Temperature: ", T_flash, " K")
println("Vapor Fraction: ", vapor_fraction)
println("Vapor Composition: ", z_vapor)
println("Liquid Composition: ", z_liquid)
```

## Initial State Specifications

Flash drums require initial state guesses. Multiple options are available:

### Pressure-Temperature-Composition-Volume State
```julia
state = pTNVState(
    10e5,                    # Pressure (Pa)
    280.0,                   # Temperature (K)
    [0.3, 0.3, 0.4],        # Composition
    base = :Pressure         # Use pressure-based flash
)
```

### Pressure-Enthalpy State
```julia
state = pHNVState(
    10e5,                    # Pressure (Pa)
    -50000.0,               # Enthalpy (J/mol)
    [0.3, 0.3, 0.4],        # Composition
    base = :Pressure
)
```

### Temperature-Volume State
```julia
state = TVNState(
    280.0,                   # Temperature (K)
    0.001,                   # Molar volume (m³/mol)
    [0.3, 0.3, 0.4]         # Composition
)
```

## Operating Modes

### Adiabatic Flash
No heat exchange with surroundings:
```julia
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = initial_state,
    p = 10e5,
    Q = 0.0  # Adiabatic
)
```

### Flash with Heat Input
Add or remove heat:
```julia
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = initial_state,
    p = 10e5,
    Q = -50000.0  # Remove 50 kW
)
```

## Complete Example

See the [Flash Drum Example](../examples/flash_drum.md) for a detailed walkthrough of:
- Setting up components and thermodynamics
- Defining feed conditions
- Creating the flash drum
- Solving the flowsheet
- Analyzing results
- Validating against reference data

## Other Separation Units

### Adsorbers (Experimental)

Adsorption units for gas purification:
```julia
@named adsorber = Adsorber(
    medium = medium,
    adsorbent = activated_carbon,
    bed_volume = 1.0
)
```

*Note: Adsorption models are under development.*

## Troubleshooting

### Convergence Issues

If the flash drum doesn't converge:

1. **Check Initial Guesses**: Ensure temperature and pressure are in reasonable range
2. **Verify EoS Model**: Some EoS models work better for certain mixtures
3. **Adjust Solver**: Try different nonlinear solvers:
   ```julia
   sol = solve(prob, FastShortcutNonlinearPolyalg())
   ```
4. **Check Phase Stability**: Ensure feed conditions allow VLE

### Common Errors

- **DimensionMismatch**: Usually from incorrect composition array size
- **Flash Failed**: Initial guess too far from solution, adjust T or p guess
- **Negative Mole Numbers**: Unphysical solution, check model equations

## Advanced Topics

### Custom Flash Specifications

Create custom flash configurations by extending base components:
- Fixed temperature flash drums
- Multi-stage flash trains
- Flash with chemical reactions
