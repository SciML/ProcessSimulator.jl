# Getting Started

This guide will help you get started with ProcessSimulator.jl for process modeling and simulation.

## Installation

ProcessSimulator.jl requires Julia 1.9 or later. Install it using the Julia package manager:

```julia
using Pkg
Pkg.add("ProcessSimulator")
```

## Basic Concepts

ProcessSimulator.jl uses a component-based modeling approach built on ModelingToolkit.jl:

1. **Media**: Define thermodynamic models for fluids
2. **Components**: Build unit operations (reactors, separators)
3. **Connections**: Connect components via streams
4. **System**: Assemble and solve the complete flowsheet

## Your First Simulation

Let's create a simple steady-state flash drum:

```julia
using ProcessSimulator
using Clapeyron
using ModelingToolkit
using NonlinearSolve
using ModelingToolkit: t_nounits as t

# Step 1: Define the fluid medium
components = ["methane", "ethane"]
eos_model = PR(components)
medium = EoSBased(components = components, eosmodel = eos_model)

# Step 2: Create a feed stream
@named feed = Boundary_pTzn(
    medium = medium,
    p = 50e5,           # 50 bar
    T = 300.0,          # 300 K
    z = [0.5, 0.5],     # Equal molar composition
    flowrate = 100.0,   # mol/s
    flowbasis = :molar
)

# Step 3: Create a flash drum
flash_state = pTNVState(50e5, 300.0, [0.5, 0.5], base = :Pressure)
@named flash = FixedPressureSteadyStateFlashDrum(
    medium = medium,
    state = flash_state,
    pressure = 30e5,    # Flash at 30 bar
    Q = nothing         # Adiabatic
)

# Step 4: Create outlet streams
@named liquid_out = ConnHouse(medium = medium)
@named vapor_out = ConnHouse(medium = medium)

# Step 5: Connect the flowsheet
connections = [
    connect(feed.odesystem.OutPort, flash.odesystem.InPort),
    connect(flash.odesystem.LiquidOutPort, liquid_out.InPort),
    connect(flash.odesystem.VaporOutPort, vapor_out.InPort)
]

# Step 6: Build and solve the system
@named sys = System(connections, t, [], [];
    systems = [feed.odesystem, flash.odesystem, liquid_out, vapor_out])

compiled_sys = mtkcompile(sys)
prob = NonlinearProblem(compiled_sys, guesses(compiled_sys))
sol = solve(prob, NewtonRaphson())

# Step 7: Extract results
T_flash = sol[flash.ControlVolumeState.T]
liquid_flow = sol[flash.LiquidOutPort.ṅ[1]]
vapor_flow = sol[flash.VaporOutPort.ṅ[1]]

println("Flash Temperature: ", T_flash, " K")
println("Liquid Flow: ", liquid_flow, " mol/s")
println("Vapor Flow: ", vapor_flow, " mol/s")
```

## Next Steps

- Learn about [Media & Thermodynamics](guide/media.md)
- Explore [Component Models](guide/components.md)
- See more [Examples](examples/flash_drum.md)
