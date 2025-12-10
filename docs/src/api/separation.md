# Separation Units API

Reference documentation for separation equipment.

## Flash Drums

### FixedPressureSteadyStateFlashDrum

Vapor-liquid flash drum operating at fixed pressure.

```@docs
FixedPressureSteadyStateFlashDrum
```

**Parameters:**
- `medium::EoSBased`: Thermodynamic medium
- `state`: Initial state specification (pTNVState, pHNVState, TVNState)
- `pressure::Real`: Operating pressure (Pa)
- `Q`: Heat duty (W), `nothing` for calculated/adiabatic

**Connectors:**
- `InPort`: Feed stream inlet
- `VaporOutPort`: Vapor product outlet (port 2)
- `LiquidOutPort`: Liquid product outlet (port 3)

**Internal Components:**
- `ControlVolume`: ThreePortControlVolume_SteadyState
- `ControlVolumeState`: Thermodynamic state inside drum

**Variables:**
- `p`: Operating pressure (Pa) - fixed
- `T`: Flash temperature (K) - calculated
- `V[2]`: Vapor volume (m³)
- `V[3]`: Liquid volume (m³)
- `nᴸⱽ[2]`: Vapor holdup (mol)
- `nᴸⱽ[3]`: Liquid holdup (mol)
- `Q`: Heat duty (W)
- `z[i]`: Overall composition
- `ϕ[i,j]`: Fugacity coefficients (j=2:vapor, j=3:liquid)

**Equations:**

Material balance:
```
ṅ_in[i] = ṅ_vapor[i] + ṅ_liquid[i]
```

Energy balance:
```
Ḣ_in + Q = Ḣ_vapor + Ḣ_liquid
```

VLE:
```
K[i] = y[i] / x[i] = ϕ_liquid[i] / ϕ_vapor[i]
```

Pressure constraint:
```
p = pressure (fixed)
```

**Example:**
```julia
# Initial guess
flash_state = pTNVState(
    10e5,           # Pressure (Pa)
    300.0,          # Temperature guess (K)
    [0.3, 0.4, 0.3], # Composition
    base = :Pressure
)

# Flash drum
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = flash_state,
    pressure = 10e5,  # 10 bar
    Q = nothing       # Adiabatic
)

# Connect
@named vapor = ConnHouse(medium = medium)
@named liquid = ConnHouse(medium = medium)

connections = [
    connect(feed.odesystem.OutPort, flash.odesystem.InPort),
    connect(flash.odesystem.VaporOutPort, vapor.InPort),
    connect(flash.odesystem.LiquidOutPort, liquid.InPort)
]
```

**Accessing Results:**
```julia
# Flash conditions
T_flash = sol[flash.ControlVolumeState.T]
p_flash = sol[flash.ControlVolumeState.p]

# Vapor product
T_vapor = sol[vapor.p_T_z_n.T]
y = [sol[vapor.p_T_z_n.z[i]] for i in 1:n]
n_vapor = sol[vapor.p_T_z_n.n]

# Liquid product
T_liquid = sol[liquid.p_T_z_n.T]
x = [sol[liquid.p_T_z_n.z[i]] for i in 1:n]
n_liquid = sol[liquid.p_T_z_n.n]

# Heat duty (if Q was nothing)
Q = sol[flash.Q]

# Vapor fraction
β = n_vapor / (n_vapor + n_liquid)
```

## Flash Calculations

### Vapor Fraction

```julia
function vapor_fraction(sol, flash_drum, vapor_port, liquid_port)
    n_V = sol[vapor_port.p_T_z_n.n]
    n_L = sol[liquid_port.p_T_z_n.n]
    return n_V / (n_V + n_L)
end
```

### K-Values

Calculate equilibrium K-values:

```julia
function k_values(sol, flash_drum, vapor_port, liquid_port, n_components)
    y = [sol[vapor_port.p_T_z_n.z[i]] for i in 1:n_components]
    x = [sol[liquid_port.p_T_z_n.z[i]] for i in 1:n_components]
    return y ./ x
end
```

### Flash Quality

For phase quality (0 = saturated liquid, 1 = saturated vapor):

```julia
function flash_quality(sol, vapor_port, liquid_port)
    n_V = sol[vapor_port.p_T_z_n.n]
    n_L = sol[liquid_port.p_T_z_n.n]
    
    # Mass-based quality (requires MW)
    # q = m_V / (m_V + m_L)
    
    # Molar-based quality
    q = n_V / (n_V + n_L)
    
    return q
end
```

## Initial State Specifications

Different flash specifications:

### Isothermal Flash
```julia
# Known temperature, calculate Q
state = pTNVState(p_flash, T_known, z_feed, base = :Pressure)
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state,
    pressure = p_flash,
    Q = nothing  # Will be calculated
)
```

### Adiabatic Flash
```julia
# Unknown temperature, Q = 0
state = pTNVState(p_flash, T_guess, z_feed, base = :Pressure)
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state,
    pressure = p_flash,
    Q = 0.0  # Adiabatic
)
```

### Flash with Heat Input
```julia
# Known Q, calculate T
state = pTNVState(p_flash, T_guess, z_feed, base = :Pressure)
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state,
    pressure = p_flash,
    Q = -100000.0  # 100 kW cooling
)
```

## Customization

### Initial Guesses

Improve convergence with better guesses:

```julia
default_guesses = guesses(compiled_system)

custom_guesses = Dict(
    flash.odesystem.V[2] => 1.0 / flash.medium.Guesses.ρ[2],  # Vapor volume
    flash.odesystem.V[3] => 1.0 / flash.medium.Guesses.ρ[3],  # Liquid volume
    flash.odesystem.nᴸⱽ[2] => 1.0,  # Vapor holdup
    flash.odesystem.nᴸⱽ[3] => 1.0   # Liquid holdup
)

merged = merge(default_guesses, custom_guesses)
prob = NonlinearProblem(compiled_system, merged)
```

### Solver Options

Different solvers for difficult cases:

```julia
# Newton-Raphson with finite differences
sol = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()))

# Trust region method
sol = solve(prob, TrustRegion())

# Automatic solver selection
sol = solve(prob, FastShortcutNonlinearPolyalg())

# With tolerances
sol = solve(prob, NewtonRaphson(), abstol=1e-8, reltol=1e-8)
```

## Validation

### Energy Balance Check

Verify energy conservation:

```julia
function check_energy_balance(sol, flash, feed, vapor, liquid)
    # Inlet enthalpy flow
    H_in = sol[feed.odesystem.OutPort.Ḣ]
    
    # Outlet enthalpy flows
    H_vapor = sol[flash.odesystem.VaporOutPort.Ḣ]
    H_liquid = sol[flash.odesystem.LiquidOutPort.Ḣ]
    
    # Heat duty
    Q = sol[flash.Q]
    
    # Check balance
    error = abs(H_in + Q - H_vapor - H_liquid)
    relative_error = error / abs(H_in) * 100
    
    println("Energy Balance Error: ", relative_error, "%")
    return relative_error < 0.1  # Less than 0.1%
end
```

### Material Balance Check

Verify component conservation:

```julia
function check_material_balance(sol, flash, feed, vapor, liquid, n_components)
    for i in 1:n_components
        n_in = sol[feed.odesystem.OutPort.ṅ[i]]
        n_vapor = sol[flash.odesystem.VaporOutPort.ṅ[i]]
        n_liquid = sol[flash.odesystem.LiquidOutPort.ṅ[i]]
        
        error = abs(n_in - n_vapor - n_liquid) / n_in * 100
        println("Component ", i, " balance error: ", error, "%")
    end
end
```

## Advanced Topics

### Multi-Stage Flash

Cascade flash drums:

```julia
@named flash1 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state1,
    pressure = 20e5,  # High pressure
    Q = 0.0
)

@named flash2 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state2,
    pressure = 10e5,  # Lower pressure
    Q = 0.0
)

@named flash3 = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state3,
    pressure = 5e5,   # Lowest pressure
    Q = 0.0
)

connections = [
    # Feed to first flash
    connect(feed.odesystem.OutPort, flash1.odesystem.InPort),
    
    # Liquid from flash1 to flash2
    connect(flash1.odesystem.LiquidOutPort, flash2.odesystem.InPort),
    
    # Liquid from flash2 to flash3
    connect(flash2.odesystem.LiquidOutPort, flash3.odesystem.InPort),
    
    # Collect products
    connect(flash1.odesystem.VaporOutPort, vapor1.InPort),
    connect(flash2.odesystem.VaporOutPort, vapor2.InPort),
    connect(flash3.odesystem.VaporOutPort, vapor3.InPort),
    connect(flash3.odesystem.LiquidOutPort, liquid_final.InPort)
]
```

### Flash Train Analysis

Calculate overall recovery:

```julia
function overall_recovery(sol, feed, product_ports, component_idx)
    n_feed = sol[feed.odesystem.OutPort.ṅ[component_idx]]
    
    total_product = sum([
        sol[port.p_T_z_n.n] * sol[port.p_T_z_n.z[component_idx]]
        for port in product_ports
    ])
    
    return total_product / n_feed * 100
end
```

## Troubleshooting

### Common Issues

**Convergence Failure:**
- Adjust initial temperature guess
- Check if feed is single-phase at flash pressure
- Try different solver
- Verify EoS parameters exist for all components

**Negative Mole Numbers:**
- Improve initial guesses
- Check thermodynamic consistency
- Verify feed composition sums to 1.0

**DimensionMismatch:**
- Ensure composition arrays match number of components
- Check connector dimensions

### Debugging Tips

```julia
# Check initial flash
using Clapeyron
(x, y, β), state = pt_flash(medium.eosmodel, p_flash, T_guess, z_feed)
println("Initial flash state: ", state)
println("Vapor fraction: ", β)
println("Vapor comp: ", y)
println("Liquid comp: ", x)

# Check if feed is flashable
T_bubble, _ = bubble_temperature(medium.eosmodel, p_flash, z_feed)
T_dew, _ = dew_temperature(medium.eosmodel, p_flash, z_feed)
println("Bubble T: ", T_bubble, " K")
println("Dew T: ", T_dew, " K")
println("Flash T must be between these values for two-phase")
```

## See Also

- [Separation Guide](../guide/separation.md) - Conceptual overview
- [Flash Drum Example](../examples/flash_drum.md) - Complete example
- [Base Components API](base_components.md) - Control volume details
- [Thermodynamics API](thermodynamics.md) - VLE calculations
