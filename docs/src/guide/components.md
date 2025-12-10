# Components

ProcessSimulator.jl provides various component models for building process flowsheets.

## Base Components

### Connectors

Connectors link unit operations via material streams:

#### PhZConnector

Standard connector for pressure-enthalpy-composition streams:

```julia
@connector function PhZConnector_(;medium, name)
    # Carries: pressure, enthalpy, composition, flow rates
end
```

### Boundary Conditions

#### Fixed Stream Boundary

Define inlet/outlet streams with fixed conditions:

```julia
@named feed = Boundary_pTzn(
    medium = medium,
    p = 10e5,                    # Pressure (Pa)
    T = 350.0,                   # Temperature (K)
    z = [0.3, 0.4, 0.3],        # Mole fractions
    flowrate = 100.0,            # Flow rate
    flowbasis = :molar           # :molar, :mass, or :volumetric
)
```

#### Connection House

Simple connector for stream endpoints:

```julia
@named outlet = ConnHouse(medium = medium)
```

## Control Volumes

Control volumes represent the material and energy balances within a unit.

### Two-Port Control Volume

For single-inlet, single-outlet units:

```julia
@named CV = TwoPortControlVolume_(medium = medium)
# Has: InPort, OutPort, ControlVolumeState
```

### Three-Port Control Volume

For flash drums with separate liquid and vapor outlets:

```julia
@named CV = ThreePortControlVolume_SteadyState(medium = medium)
# Has: InPort, LiquidOutPort, VaporOutPort, ControlVolumeState
```

## Connecting Components

Use ModelingToolkit's `connect` function to link components:

```julia
connections = [
    connect(feed.odesystem.OutPort, unit.odesystem.InPort),
    connect(unit.odesystem.OutPort, product.InPort)
]

@named sys = System(connections, t, [], [];
    systems = [feed.odesystem, unit.odesystem, product])
```

## Building Custom Components

Create custom components using the `@component` macro:

```julia
@component function MyUnit(; medium, name, parameter1)
    # Define systems (control volumes, connectors)
    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        CV = TwoPortControlVolume_(medium = medium)
    end
    
    # Define variables
    vars = @variables begin
        custom_var(t), [description = "Custom variable"]
    end
    
    # Define parameters
    pars = @parameters begin
        param1 = parameter1
    end
    
    # Define equations
    eqs = [
        # Your model equations here
        OutPort.p ~ InPort.p - 1000.0  # Pressure drop example
    ]
    
    return ODESystem(eqs, t, vars, pars; name, systems = [systems...])
end
```

## State Variables

All components have access to thermodynamic state through `ControlVolumeState`:

- `p`: Pressure
- `T`: Temperature  
- `z[:, j]`: Mole fractions (component i, phase j)
- `ϕ[j]`: Phase fractions
- `ρ[j]`: Molar densities
- `h[j]`: Molar enthalpies
