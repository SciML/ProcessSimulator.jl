"""
    Abstract type hierarchy for reactors in ProcessSimulator.jl

The reactor type hierarchy provides a structured way to categorize different
reactor models based on their fundamental modeling approach:

- `AbstractReactor`: Base type for all reactor models
- `AbstractEquilibriumReactor`: Reactors based on thermodynamic equilibrium (Gibbs minimization)
- `AbstractKineticReactor`: Reactors based on reaction kinetics (rate equations)
"""

"""
    AbstractReactor

Base abstract type for all reactor components. Reactors transform input streams
through chemical reactions to produce output streams.
"""
abstract type AbstractReactor end

"""
    AbstractEquilibriumReactor <: AbstractReactor

Abstract type for equilibrium-based reactors. These reactors assume chemical
equilibrium is reached and determine the outlet composition by minimizing
Gibbs free energy subject to elemental balance constraints.

Examples include:
- GibbsReactor: General equilibrium reactor using Gibbs minimization
- EquilibriumReactor: Simplified equilibrium reactor with specified reactions
"""
abstract type AbstractEquilibriumReactor <: AbstractReactor end

"""
    AbstractKineticReactor <: AbstractReactor

Abstract type for kinetics-based reactors. These reactors model the rate
of chemical reactions using rate expressions (e.g., Arrhenius kinetics).

Examples include:
- CSTR (Continuous Stirred Tank Reactor): Perfect mixing assumption
- PFR (Plug Flow Reactor): No axial mixing, 1D spatial variation
- BatchReactor: Time-varying composition in closed vessel
"""
abstract type AbstractKineticReactor <: AbstractReactor end

export AbstractReactor, AbstractEquilibriumReactor, AbstractKineticReactor
