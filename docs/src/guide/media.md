# Media & Thermodynamics

ProcessSimulator.jl uses Clapeyron.jl for thermodynamic property calculations. This page explains how to define and use fluid media.

## Defining a Medium

A medium represents the thermodynamic model for your fluids:

```julia
using Clapeyron

# Define components
components = ["methane", "ethane", "propane", "n-butane"]

# Choose an equation of state
eos_model = PR(components)  # Peng-Robinson

# Create the medium
medium = EoSBased(components = components, eosmodel = eos_model)
```

## Available Equations of State

ProcessSimulator.jl supports all Clapeyron.jl EoS models:

- **Cubic EoS**: `PR`, `SRK`, `RK`, `vdW`
- **SAFT**: `PCSAFT`, `sPCSAFT`, `SAFTVRMie`
- **Activity Models**: `UNIFAC`, `NRTL`, `Wilson`

```julia
# Peng-Robinson with Reid ideal gas
model = PR(components, idealmodel = ReidIdeal)

# Soave-Redlich-Kwong
model = SRK(components)

# PC-SAFT
model = PCSAFT(components)
```

## Thermodynamic States

Define the state of a stream or unit:

### Pressure-Temperature-Composition (pTz)

```julia
state = pTNVState(
    5e5,              # Pressure (Pa)
    350.0,            # Temperature (K)
    [0.3, 0.4, 0.3],  # Mole fractions
    base = :Pressure   # Calculate from P,T,z
)
```

### Volume-Temperature-Composition (VTN)

```julia
state = pTNVState(
    nothing,          # Pressure to be calculated
    350.0,            # Temperature (K)
    [10.0, 15.0, 5.0], # Molar amounts (mol)
    base = :Volume    # Calculate from V,T,N
)
state.V = 0.1  # Set volume (m³)
```

## Property Calculations

The medium provides automatic property calculations:

```julia
# After creating a flash calculation
medium, state, phase = resolve_guess!(medium, state)

# Access properties from medium.Guesses
p = medium.Guesses.p          # Pressure
T = medium.Guesses.T          # Temperature
ρ = medium.Guesses.ρ          # Densities [overall, liquid, vapor]
h = medium.Guesses.h          # Enthalpies [overall, liquid, vapor]
x = medium.Guesses.x          # Compositions [overall, liquid, vapor]
ϕ = medium.Guesses.ϕ          # Phase fractions [liquid, vapor]
```

## Multi-Phase Equilibrium

ProcessSimulator.jl automatically handles VLE calculations through Clapeyron.jl (At the moment):

```julia
# Flash calculation at given P,T
sol = TP_flash(medium.EoSModel, p, T, z)
ϕ = sol[1]        # Phase fractions
x = sol[2]        # Phase compositions

# Get phase properties
x_liquid = flash_mol_fractions_liquid(medium.EoSModel, p, T, z)
x_vapor = flash_mol_fractions_vapor(medium.EoSModel, p, T, z)
vapor_frac = flash_vaporized_fraction(medium.EoSModel, p, T, z)
```

## Transport Properties

Define mass and heat transfer coefficients:

```julia
medium = EoSBased(
    components = components,
    eosmodel = eos_model,
    transportmodel = TransportModel(
        mass_transfer = ConstantMassTransferCoeff([0.5, 0.5, 0.5]),
        heat_transfer = ConstantHeatTransferCoeff(10.0)
    )
)
```
