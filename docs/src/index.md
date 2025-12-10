```@meta
CurrentModule = ProcessSimulator
```

# ProcessSimulator.jl

ProcessSimulator.jl is a Julia package for process modeling and simulation using the ModelingToolkit.jl ecosystem. It provides a framework for building equation-oriented models of chemical processes with advanced thermodynamic property calculations.

## Features

- **Thermodynamic Models**: Integration with Clapeyron.jl for equation of state (EoS) based property calculations
- **Component-Based Modeling**: Modular components for reactors, separators, and unit operations
- **Symbolic Equations**: Built on ModelingToolkit.jl for symbolic equation manipulation
- **Flexible Connections**: Stream-based connections between unit operations
- **Multiple Phases**: Support for multi-phase systems (liquid-vapor equilibrium)

## Package Components

### Core Modules

- **Media & Thermodynamics**: EoS-based fluid property calculations
- **Base Components**: Control volumes, connectors, and boundary conditions
- **Reactors**: CSTR and other reactor models
- **Separation Units**: Flash drums, distillation columns (in development)

## Quick Example

```julia
using ProcessSimulator
using Clapeyron
using ModelingToolkit

# Define components and EoS model
components = ["methane", "ethane", "propane"]
model = PR(components)
medium = EoSBased(components = components, eosmodel = model)

# Create a flash drum
state = pTNVState(5e5, 300.0, [0.3, 0.4, 0.3], base = :Pressure)
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = state,
    pressure = 3e5,
    Q = nothing
)
```

## Installation

```julia
using Pkg
Pkg.add("ProcessSimulator")
```

## Getting Help

- **Documentation**: This documentation site
- **Issues**: [GitHub Issues](https://github.com/SciML/ProcessSimulator.jl/issues)
- **Discussions**: [GitHub Discussions](https://github.com/SciML/ProcessSimulator.jl/discussions)
