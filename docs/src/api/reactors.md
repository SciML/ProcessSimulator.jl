# Reactors API

Reference documentation for reactor components.

## Continuous Stirred Tank Reactor (CSTR)

### FixedVolumeSteadyStateCSTR

Steady-state CSTR with fixed volume.

```@docs
FixedVolumeSteadyStateCSTR
```

**Parameters:**
- `medium::EoSBased`: Thermodynamic medium
- `reactionset`: Single reaction or array of reactions
- `limiting_reactant::String`: Name of limiting reactant component (only relevant for single reactions but still required)
- `state`: Initial state specification (pTNVState, pHNVState, etc.)
- `volume::Real`: Reactor volume (m³)
- `W::Real`: Shaft work input (W), typically 0 for CSTR
- `Q`: Heat duty (W), `nothing` for adiabatic

**Connectors:**
- `InPort`: Feed stream inlet
- `OutPort`: Product stream outlet

**Variables:**
- `ControlVolumeState`: Thermodynamic state inside reactor
  - `p`: Pressure (Pa)
  - `T`: Temperature (K)
  - `z[i]`: Mole fractions
  - `ϕ[i]`: Vaporized Fraction
- `r[j]`: Reaction rates (mol/s) for each reaction
- `V`: Reactor volume (m³)
- `n`: Total molar holdup (mol)
- `U`: Internal energy (J)
- `Q`: Heat duty (W)

**Equations:**

Material balance:
```
ṅ_in[i] + Σ(νᵢⱼ × r[j]) = ṅ_out[i]
```

Energy balance:
```
Ḣ_in + Q + W + Σ(ΔH_rxn[j] × r[j]) = Ḣ_out
```

VLE (if two-phase):
```
K[i] = y[i] / x[i] = ϕ_liquid[i] / ϕ_vapor[i]
```

**Example:**
```julia
# Define reaction
rxn = PowerLawReaction(
    components = ["A", "B", "C"],
    stoichiometry = Dict("A" => -1.0, "B" => -1.0, "C" => 1.0),
    order = Dict("A" => 1.0, "B" => 1.0),
    A = 1e5,
    Eₐ = 50000.0
)

# Create CSTR
initial_state = pTNVState(10e5, 350.0, [0.33, 0.33, 0.34], base = :Pressure)

@named reactor = FixedVolumeSteadyStateCSTR(
    medium = medium,
    reactionset = rxn,
    limiting_reactant = "A",
    state = initial_state,
    volume = 1.0,
    W = 0.0,
    Q = nothing  # Adiabatic
)
```

**Accessing Results:**
```julia
# Reactor conditions
T = sol[reactor.ControlVolumeState.T]
p = sol[reactor.ControlVolumeState.p]
z = [sol[reactor.ControlVolumeState.z[i]] for i in 1:n]

# Reaction rates
rates = [sol[reactor.r[j]] for j in 1:n_reactions]

# Heat duty (if Q was nothing)
Q = sol[reactor.Q]

# Product flow
n_out = sol[reactor.OutPort.ṅ[1]]  # Total outlet flow
```

## Reaction Definitions

### PowerLawReaction

Power law kinetic model.

```@docs
PowerLawReaction
```

**Parameters:**
- `components::Vector{String}`: Component names
- `stoichiometry::Dict{String, Real}`: Stoichiometric coefficients (negative for reactants)
- `order::Dict{String, Real}`: Reaction orders for each component
- `A::Real`: Pre-exponential factor (units depend on overall order)
- `Eₐ::Real`: Activation energy (J/mol)

**Rate Expression:**
```
r = A × exp(-Eₐ/RT) × Π(Cᵢ^orderᵢ)
```

**Example:**
```julia
# A + 2B → C
rxn = PowerLawReaction(
    components = ["A", "B", "C"],
    stoichiometry = Dict(
        "A" => -1.0,
        "B" => -2.0,
        "C" => 1.0
    ),
    order = Dict(
        "A" => 1.0,
        "B" => 2.0,
        "C" => 0.0
    ),
    A = 1e6,      # mol^(-2) m^6 s^-1
    Eₐ = 75000.0  # J/mol
)
```

### ArrheniusReaction

Arrhenius kinetic model (alias for PowerLawReaction).

**Example:**
```julia
rxn = ArrheniusReaction(
    components = ["A", "B"],
    stoichiometry = Dict("A" => -1.0, "B" => 1.0),
    order = Dict("A" => 1.0),
    A = 5e4,
    Eₐ = 60000.0
)
```

### Custom Reactions

Extend the reaction framework:

```julia
struct CustomReaction
    components::Vector{String}
    stoichiometry::Dict{String, Float64}
    # Custom parameters
end

function rate(rxn::CustomReaction, T, p, C)
    # Implement custom rate law
    # Return reaction rate (mol/s)
end
```

## Reaction Systems

### ReactionSystem

Container for multiple reactions.

**Creating:**
```julia
# Multiple reactions
rxns = [rxn1, rxn2, rxn3]

# Use in CSTR
@named reactor = FixedVolumeSteadyStateCSTR(
    medium = medium,
    reactionset = rxns,  # Array of reactions
    limiting_reactant = "A",
    ...
)
```

**Accessing Individual Rates:**
```julia
r1 = sol[reactor.r[1]]  # Rate of first reaction
r2 = sol[reactor.r[2]]  # Rate of second reaction
```

## Conversion Calculations

Calculate conversion from CSTR results:

```julia
# Single reactant
function conversion(sol, reactor, feed, component_index)
    n_in = sol[feed.odesystem.OutPort.ṅ[component_index]]
    n_out = sol[reactor.OutPort.ṅ[component_index]]
    return (n_in - n_out) / n_in * 100  # Percentage
end

# Using compositions
function conversion_from_composition(z_feed, z_product, n_feed, n_product, i)
    moles_in = z_feed[i] * n_feed
    moles_out = z_product[i] * n_product
    return (moles_in - moles_out) / moles_in * 100
end
```

## Selectivity and Yield

For multiple products:

```julia
# Selectivity (moles of desired product / moles of limiting reactant consumed)
function selectivity(sol, reactor, feed, product_idx, reactant_idx)
    reactant_consumed = sol[feed.odesystem.OutPort.ṅ[reactant_idx]] - 
                       sol[reactor.OutPort.ṅ[reactant_idx]]
    product_formed = sol[reactor.OutPort.ṅ[product_idx]] - 
                    sol[feed.odesystem.OutPort.ṅ[product_idx]]
    
    return product_formed / reactant_consumed
end

# Yield (moles of product / moles of reactant fed)
function yield(sol, reactor, feed, product_idx, reactant_idx)
    product_formed = sol[reactor.OutPort.ṅ[product_idx]] - 
                    sol[feed.odesystem.OutPort.ṅ[product_idx]]
    reactant_fed = sol[feed.odesystem.OutPort.ṅ[reactant_idx]]
    
    return product_formed / reactant_fed * 100  # Percentage
end
```

## Heat Duty Calculation

For non-adiabatic operation:

```julia
# If Q is specified
@named reactor = FixedVolumeSteadyStateCSTR(..., Q = -50000.0)  # 50 kW cooling

# If Q is calculated (Q = nothing)
@named reactor = FixedVolumeSteadyStateCSTR(..., Q = nothing)
Q_required = sol[reactor.Q]  # Heat duty needed for specified T
```

## Residence Time

Calculate mean residence time:

```julia
function residence_time(sol, reactor, component_index = 1)
    V = sol[reactor.V]  # Reactor volume
    n = sol[reactor.n]  # Total holdup
    ṅ = sol[reactor.OutPort.ṅ[component_index]]
    
    # Molar residence time
    τ_molar = n / ṅ  # seconds
    
    # Volumetric residence time (approximate)
    F = ṅ / sum([sol[reactor.OutPort.ṅ[i]] for i in 1:length(components)])
    τ_vol = V / F  # m³/(mol/s) - needs conversion to time
    
    return τ_molar
end
```

## Advanced Topics

### Temperature Control

Implement temperature control via heat duty:

```julia
# Target temperature
T_target = 360.0  # K

# Solve with temperature constraint
# (Requires adding constraint to system)
```

### Pressure Drop

Add pressure drop correlation:

```julia
# Custom component with pressure drop
@component function CSTRWithPressureDrop(...)
    @extend cv = TwoPortControlVolume_SteadyState(...)
    
    @equations begin
        # Pressure drop equation
        OutPort.p ~ InPort.p - ΔP
    end
end
```

### Multiple Phases

CSTR automatically handles VLE if conditions are appropriate:

```julia
# Check phase state
ϕ_vapor = [sol[reactor.ControlVolumeState.ϕ[i, 2]] for i in 1:n]
ϕ_liquid = [sol[reactor.ControlVolumeState.ϕ[i, 3]] for i in 1:n]
```

## See Also

- [Reactors Guide](../guide/reactors.md) - Conceptual overview
- [CSTR Example](../examples/cstr.md) - Complete example
- [Base Components API](base_components.md) - Control volume details
- [Thermodynamics API](thermodynamics.md) - Property calculations
