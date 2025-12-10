# Base Components API

Reference documentation for base components, connectors, and control volumes.

## Connectors

### StreamConnector

Two-phase stream connector for material and energy transfer.

```@docs
StreamConnector
```

**Variables:**
- `ṅ[i]`: Component molar flow rates (mol/s)
- `Ḣ`: Enthalpy flow rate (W)

**Example:**
```julia
@connector StreamConnector begin
    ṅ(t)[1:n], [description = "Molar flow rate"]
    Ḣ(t), [description = "Enthalpy flow rate"]
end
```

**Usage in Components:**
```julia
@named InPort = StreamConnector(n = length(components))
@named OutPort = StreamConnector(n = length(components))
```

## Boundaries

### Boundary_pTzn

Fixed pressure, temperature, composition, and flow rate boundary.

```@docs
Boundary_pTzn
```

**Parameters:**
- `medium::EoSBased`: Thermodynamic medium
- `p::Real`: Pressure (Pa)
- `T::Real`: Temperature (K)
- `z::Vector{Real}`: Mole fractions
- `flowrate::Real`: Flow rate
- `flowbasis::Symbol`: `:molar` or `:mass`

**Example:**
```julia
@named feed = Boundary_pTzn(
    medium = medium,
    p = 30e5,
    T = 350.0,
    z = [0.4, 0.6],
    flowrate = 100.0,
    flowbasis = :molar
)
```

**Connector:**
- `OutPort`: Stream connector providing the specified conditions

### ConnHouse

Connection house for collecting stream properties.

```@docs
ConnHouse
```

Collects properties from a stream connector into a named tuple state.

**Parameters:**
- `medium::EoSBased`: Thermodynamic medium

**Example:**
```julia
@named product = ConnHouse(medium = medium)

# After solving, access properties
p = sol[product.p_T_z_n.p]
T = sol[product.p_T_z_n.T]
z = [sol[product.p_T_z_n.z[i]] for i in 1:n]
n = sol[product.p_T_z_n.n]
```

**Connector:**
- `InPort`: Stream connector accepting inlet stream

**Properties Collected:**
- `p`: Pressure
- `T`: Temperature
- `z[i]`: Mole fractions
- `n`: Total molar flow rate

## Control Volumes

### ThreePortControlVolume_SteadyState

Three-port control volume for steady-state separation equipment.

```@docs
ThreePortControlVolume_SteadyState
```

**Ports:**
- `InPort`: Feed inlet
- `VaporOutPort`: Vapor outlet (port 2)
- `LiquidOutPort`: Liquid outlet (port 3)

**Parameters:**
- `medium::EoSBased`: Thermodynamic medium
- `state`: Initial state specification
- `Q`: Heat duty (W), `nothing` for calculated

**Variables:**
- `V[i]`: Phase volumes (m³)
- `nᴸⱽ[i]`: Phase holdups (mol)
- `U`: Internal energy (J)
- `Nᵢ[i]`: Component inventories (mol)

**Equations:**
1. Material balance: `ṅ_in[i] = ṅ_vapor[i] + ṅ_liquid[i]`
2. Energy balance: `Ḣ_in + Q = Ḣ_vapor + Ḣ_liquid`
3. VLE: Phase equilibrium for two-phase systems
4. Volume constraint: `V_total = V_vapor + V_liquid`

**Example:**
```julia
initial_state = pTNVState(10e5, 300.0, [0.5, 0.5], base = :Pressure)

@named cv = ThreePortControlVolume_SteadyState(
    medium = medium,
    state = initial_state
)
```

**Used in:**
- `FixedPressureSteadyStateFlashDrum`
- Custom separation equipment

### TwoPortControlVolume_SteadyState

Two-port control volume for steady-state equipment with single outlet.

**Ports:**
- `InPort`: Feed inlet
- `OutPort`: Product outlet

**Parameters:**
- `medium::EoSBased`: Thermodynamic medium
- `state`: Initial state specification
- `Q`: Heat duty (W), `nothing` for calculated
- `W`: Shaft work (W)

**Variables:**
- `V`: Volume (m³)
- `n`: Holdup (mol)
- `U`: Internal energy (J)

**Equations:**
1. Material balance: `ṅ_in[i] = ṅ_out[i]`
2. Energy balance: `Ḣ_in + Q + W = Ḣ_out`

**Example:**
```julia
@named cv = TwoPortControlVolume_SteadyState(
    medium = medium,
    state = initial_state,
    Q = 0.0,
    W = 0.0
)
```

**Used in:**
- `FixedVolumeSteadyStateCSTR`
- Heat exchangers
- Valves

## State Specifications

Components use state specifications for initial guesses:

### pTNVState
```julia
state = pTNVState(p, T, z; base = :Pressure)
```

### pHNVState
```julia
state = pHNVState(p, H, z; base = :Pressure)
```

### TVNState
```julia
state = TVNState(T, V, z)
```

See [Thermodynamics API](thermodynamics.md) for details.

## Helper Functions

### connect

ModelingToolkit connection function:

```julia
connections = [
    connect(source.OutPort, destination.InPort)
]
```

### System

Create MTK system from connections:

```julia
@named sys = System(
    connections,
    t,
    [],
    [];
    systems = [comp1.odesystem, comp2.odesystem, ...]
)
```

### mtkcompile

Compile MTK system for solving:

```julia
compiled = mtkcompile(sys)
compiled = mtkcompile(sys, use_scc = true)  # Use strongly connected components
```

## Solving

### Initial Guesses

Get default guesses:
```julia
guess_dict = guesses(compiled_system)
```

Merge with custom guesses:
```julia
custom = Dict(var => value, ...)
merged = merge(guesses(compiled_system), custom)
```

### NonlinearProblem

Create nonlinear problem:
```julia
prob = NonlinearProblem(compiled_system, guess_dict)
prob = NonlinearProblem(compiled_system, guess_dict, use_scc = true)
```

### Solve

```julia
sol = solve(prob, NewtonRaphson())
sol = solve(prob, FastShortcutNonlinearPolyalg())
sol = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()))
```

### Access Results

```julia
# Variable value
value = sol[compiled_system.component.variable]

# Array variable
values = [sol[compiled_system.component.array[i]] for i in 1:n]
```

## Output Formatting

### print_flowsheet_summary

Print formatted flowsheet results:

```julia
# Using unit operations dictionary
unit_ops = Dict(:R1 => :CSTR, :F1 => :Feed)
print_flowsheet_summary(sol, compiled_sys, unit_ops, components)

# Direct component specification
print_flowsheet_summary(sol, compiled_sys, components, comp1, comp2, ...)
```

Outputs:
- Component flow rates and compositions
- Temperatures and pressures
- Phase information
- Heat duties

## See Also

- [Components Guide](../guide/components.md) - Conceptual overview
- [Media Guide](../guide/media.md) - Thermodynamics integration
- [Flash Drum Example](../examples/flash_drum.md) - Using control volumes
