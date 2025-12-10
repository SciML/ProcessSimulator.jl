# Thermodynamics API

Reference documentation for thermodynamic functions and state specifications.

## Medium Definition

```@docs
EoSBased
```

## State Specifications

### Pressure-Temperature State

```@docs
pTNVState
```

Specify state by pressure, temperature, composition, and optionally volume.

**Parameters:**
- `p::Real`: Pressure (Pa)
- `T::Real`: Temperature (K)
- `z::Vector{Real}`: Mole fractions
- `base::Symbol`: Flash type (`:Pressure` or `:Temperature`)

**Example:**
```julia
state = pTNVState(10e5, 300.0, [0.5, 0.5], base = :Pressure)
```

### Pressure-Enthalpy State

```@docs
pHNVState
```

Specify state by pressure, enthalpy, and composition.

**Parameters:**
- `p::Real`: Pressure (Pa)
- `H::Real`: Molar enthalpy (J/mol)
- `z::Vector{Real}`: Mole fractions
- `base::Symbol`: Flash type

**Example:**
```julia
state = pHNVState(10e5, -50000.0, [0.5, 0.5], base = :Pressure)
```

### Temperature-Volume State

```@docs
TVNState
```

Specify state by temperature, volume, and composition.

**Parameters:**
- `T::Real`: Temperature (K)
- `V::Real`: Molar volume (m³/mol)
- `z::Vector{Real}`: Mole fractions

**Example:**
```julia
state = TVNState(300.0, 0.001, [0.5, 0.5])
```

## Thermodynamic Property Functions

### VLE Calculations

Calculate vapor-liquid equilibrium:

```julia
# Bubble point
T_bubble, (x, y, β) = bubble_temperature(eos_model, p, z)
p_bubble, (x, y, β) = bubble_pressure(eos_model, T, z)

# Dew point
T_dew, (x, y, β) = dew_temperature(eos_model, p, z)
p_dew, (x, y, β) = dew_pressure(eos_model, T, z)

# PT flash
(x, y, β), state = pt_flash(eos_model, p, T, z)
```

**Returns:**
- `x`: Liquid composition
- `y`: Vapor composition
- `β`: Vapor fraction
- `state`: Flash state (`:liquid`, `:vapor`, `:two-phase`)

### Enthalpy Calculations

```julia
# Total enthalpy
H = enthalpy(eos_model, p, T, z)

# Phase-specific enthalpy
H_liquid = enthalpy(eos_model, V_liquid, T, z, phase = :liquid)
H_vapor = enthalpy(eos_model, V_vapor, T, z, phase = :vapor)
```

### Entropy Calculations

```julia
# Total entropy
S = entropy(eos_model, p, T, z)

# Phase-specific entropy
S_liquid = entropy(eos_model, V_liquid, T, z, phase = :liquid)
```

### Volume Calculations

```julia
# Molar volume
V = volume(eos_model, p, T, z)

# Phase-specific volume
V_liquid = volume(eos_model, p, T, z, phase = :liquid)
V_vapor = volume(eos_model, p, T, z, phase = :vapor)
```

### Density Calculations

```julia
# Mass density
ρ_mass = mass_density(eos_model, p, T, z)

# Molar density
ρ_molar = molar_density(eos_model, p, T, z)
```

## Equation of State Models

### Cubic EoS

Supported cubic equations of state from Clapeyron.jl:

**Peng-Robinson:**
```julia
eos = PR(components)
eos = PR(components, alpha = PRAlpha, mixing = vdW1fRule)
```

**Soave-Redlich-Kwong:**
```julia
eos = SRK(components)
eos = SRK(components, idealmodel = ReidIdeal)
```

**Redlich-Kwong:**
```julia
eos = RK(components)
```

### SAFT Models

Statistical Associating Fluid Theory models:

**PC-SAFT:**
```julia
eos = PCSAFT(components)
```

**SAFT-γ Mie:**
```julia
eos = SAFTgammaMie(components)
```

### Composite Models

Combine different models for different phases:

```julia
idealmodel = CompositeModel(components,
    liquid = RackettLiquid,
    gas = ReidIdeal(components, reference_state = :formation),
    saturation = DIPPR101Sat,
    hvap = DIPPR106HVap
)

activity = NRTL(components)
model = CompositeModel(components, fluid = idealmodel, liquid = activity)
```

### Ideal Models

**Ideal Gas:**
```julia
eos = IdealModel(components)
eos = ReidIdeal(components, reference_state = :formation)
```

## Transport Properties

### Viscosity

```julia
# Dynamic viscosity (Pa·s)
μ = viscosity(eos_model, p, T, z)

# Phase-specific viscosity
μ_liquid = viscosity(eos_model, p, T, z, phase = :liquid)
μ_vapor = viscosity(eos_model, p, T, z, phase = :vapor)
```

### Thermal Conductivity

```julia
# Thermal conductivity (W/m/K)
k = thermal_conductivity(eos_model, p, T, z)

# Phase-specific thermal conductivity
k_liquid = thermal_conductivity(eos_model, p, T, z, phase = :liquid)
```

### Surface Tension

```julia
# Surface tension (N/m)
σ = surface_tension(eos_model, T, z)
```

## Initial Guess Management

### Resolve Guesses

Automatically update thermodynamic guesses:

```julia
resolve_guess!(medium, state)
```

This function:
1. Performs flash calculation based on state specification
2. Updates `medium.Guesses` with VLE results
3. Stores phase compositions, densities, and enthalpies

**Called automatically by**: Boundary conditions and components during initialization

### Manual Guess Creation

Create custom guesses:

```julia
guesses = EosBasedGuesses(
    eos_model,
    p,
    T,
    z,
    Val(:Pressure)  # or Val(:Temperature)
)
```

## Utility Functions

### Component Properties

```julia
# Molecular weight
MW = molar_mass(eos_model, z)

# Critical properties
Tc = critical_temperature(eos_model)
pc = critical_pressure(eos_model)
Vc = critical_volume(eos_model)

# Acentric factor
ω = acentric_factor(eos_model)
```

### Unit Conversions

```julia
# Temperature
T_C = celsius(T_K)
T_K = kelvin(T_C)

# Pressure
p_bar = bar(p_Pa)
p_Pa = pascal(p_bar)

# Flow basis
n_molar = molar_flow(m_mass, MW)
m_mass = mass_flow(n_molar, MW)
```

## See Also

- [Media Guide](../guide/media.md) - Conceptual overview of thermodynamics
- [Components Guide](../guide/components.md) - Using thermodynamics in components
- [Flash Drum Example](../examples/flash_drum.md) - VLE flash example
