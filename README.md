# ProcessSimulator.jl

[![Build Status](https://github.com/SciML/ProcessSimulator.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SciML/ProcessSimulator.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/SciML/ProcessSimulator.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/SciML/ProcessSimulator.jl)

A Julia package for process modeling and simulation built on the [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) ecosystem. It provides component-based, equation-oriented modeling of reactors, separators, and other unit operations, with thermodynamic support via [Clapeyron.jl](https://github.com/viniviena/Clapeyron.jl).

## Installation

Install Julia (≥ 1.10), then add the package from the Julia REPL:

```julia
using Pkg
Pkg.add("ProcessSimulator")
```

Thermodynamic calculations require Clapeyron.jl. Install the fork used by this package:

```julia
Pkg.add(url = "https://github.com/viniviena/Clapeyron.jl")
```

Install the remaining solver dependencies:

```julia
Pkg.add(["ModelingToolkit", "NonlinearSolve"])
```

## Example: Adiabatic CSTR — Ethylene Glycol Production

This example models the liquid-phase hydration of ethylene oxide to ethylene glycol in a steady-state adiabatic CSTR.

**Reaction:**
```
C₂H₄O + H₂O → C₂H₆O₂
(ethylene oxide + water → ethylene glycol)
```

**Feed:** 100 mol/s, 5 atm, 350.15 K — 40 mol% ethylene oxide, 60 mol% water
**Reactor:** 1.0 m³, adiabatic (Q = 0), limiting reactant: ethylene oxide

### 1. Imports

```julia
using ProcessSimulator
using Clapeyron
using ModelingToolkit
using NonlinearSolve
using ModelingToolkit: t_nounits as t
```

### 2. Thermodynamic Model

A composite model combines an ideal reference state with NRTL liquid-phase activity:

```julia
components = ["ethylene oxide", "water", "ethylene glycol"]

idealmodel = CompositeModel(components,
    liquid     = RackettLiquid,
    gas        = ReidIdeal(components, reference_state = :formation),
    saturation = DIPPR101Sat,
    hvap       = DIPPR106HVap
)

activity = NRTL(components)
model    = CompositeModel(components, fluid = idealmodel, liquid = activity)

# Quick sanity check
bubble_pressure(model, 350.15, [0.4, 0.6, 0.0])

medium = EoSBased(components = components, eosmodel = model)
```

### 3. Feed Stream

```julia
@named S1 = Boundary_pTzn(
    medium    = medium,
    p         = 5 * 101325.0,    # Pa
    T         = 350.15,          # K
    z         = [0.4, 0.6, 0.0],
    flowrate  = 100.0,
    flowbasis = :molar
)
```

### 4. Reaction Kinetics

Power-law rate: r = 0.5 × [C₂H₄O]¹ × [H₂O]⁰ (first order in EO, no temperature dependence)

```julia
rxn1 = PowerLawReaction(
    components   = components,
    stoichiometry = Dict(
        "ethylene oxide"  => -1.0,
        "water"           => -1.0,
        "ethylene glycol" =>  1.0
    ),
    order = Dict(
        "ethylene oxide" => 1.0,
        "water"          => 0.0
    ),
    A  = 0.5,
    Eₐ = 0.0
)
```

### 5. Reactor

```julia
reactor_state = pTNVState(5 * 101325.0, 350.15, ones(3) / 3.0, base = :Pressure)

@named R1 = FixedVolumeSteadyStateCSTR(
    medium            = medium,
    reactionset       = rxn1,
    limiting_reactant = "ethylene oxide",
    state             = reactor_state,
    volume            = 1.0,   # m³
    W                 = 0.0,
    Q                 = nothing  # adiabatic
)
```

### 6. Flowsheet Assembly

```julia
@named sink = ConnHouse(medium = medium)

connection_set = [
    connect(S1.odesystem.OutPort, R1.odesystem.InPort),
    connect(R1.odesystem.OutPort, sink.InPort)
]

@named sys = System(connection_set, t, [], [];
    systems = [S1.odesystem, R1.odesystem, sink])

AdiabaticVolumeReactor = mtkcompile(sys, use_scc = true)
```

### 7. Solve

```julia
default_guesses = guesses(AdiabaticVolumeReactor)
guesses_Reactor = Dict(
    R1.odesystem.V[3]           => 1e-8,
    R1.odesystem.OutPort.ṅ[1]  => 100.0,
)
merged_guesses = merge(default_guesses, guesses_Reactor)

prob = NonlinearProblem(AdiabaticVolumeReactor, merged_guesses)
sol  = solve(prob, abstol = 1e-8, reltol = 1e-8)
```

### 8. Results

```julia
T_reactor  = sol[AdiabaticVolumeReactor.R1.ControlVolumeState.T]
conversion = sol[AdiabaticVolumeReactor.R1.X] * 100.0
r_reaction = sol[AdiabaticVolumeReactor.R1.r[1]]

println("Reactor Temperature : ", T_reactor, " K")
println("EO Conversion       : ", round(conversion, digits = 3), " %")
println("Reaction Rate       : ", r_reaction, " mol/s")
```

Print a full flowsheet table:

```julia
unit_ops = Dict(:R1 => :CSTR, :S1 => :Feed)
print_flowsheet_summary(sol, AdiabaticVolumeReactor, unit_ops, components)
```

The complete runnable script is in [test/reactors/ss_cstr_test.jl](test/reactors/ss_cstr_test.jl).

## Features

- **Equation-oriented modeling** via ModelingToolkit.jl — symbolic manipulation and automatic differentiation
- **Thermodynamics** — equation-of-state and activity-coefficient models through Clapeyron.jl
- **Reactors** — steady-state and dynamic CSTR with power-law kinetics
- **Separation** — flash drums and adsorption columns
- **Stream connections** — declarative `connect` API for building flowsheets
- **Summary utilities** — formatted tables for stream and unit-operation results

## License

MIT — see [LICENSE](LICENSE).
