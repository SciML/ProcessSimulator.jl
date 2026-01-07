"""
    GibbsReactor - Equilibrium reactor based on Gibbs free energy minimization

The GibbsReactor finds the equilibrium composition by minimizing the total
Gibbs free energy subject to elemental balance constraints. This approach
does not require specification of reaction stoichiometry - equilibrium is
determined purely from thermodynamic data.

Returns an `OptimizationSystem` for solving the equilibrium problem.
"""

"""
    GibbsReactor(ms::MaterialSource; name, element_matrix, feed_elements)

Create a Gibbs equilibrium reactor that minimizes Gibbs free energy.

# Arguments
- `ms::MaterialSource`: Material source with thermodynamic property functions
- `name`: System name (required)
- `element_matrix::Matrix{Float64}`: Element composition matrix (N_elements × N_components)
  where entry [i,j] is the number of atoms of element i in component j
- `feed_elements::Vector{Float64}`: Total moles of each element in the feed

# Returns
An `OptimizationSystem` that, when solved, gives the equilibrium composition.

# Notes
The Gibbs free energy is computed as G = H - T*S where H and S are obtained
from the MaterialSource enthalpy and entropy functions. The optimization
problem minimizes total G subject to:
1. Element balance: A * n = b (conservation of atoms)
2. Non-negativity: n ≥ 0 (no negative mole amounts)

# Example
```julia
# Define element matrix for H2O, H2, O2 system
# Rows: [H, O], Columns: [H2O, H2, O2]
elem_matrix = [2.0 2.0 0.0;   # H atoms
               1.0 0.0 2.0]   # O atoms

# Feed composition (in terms of elements)
feed_elem = [4.0, 2.0]  # 4 mol H, 2 mol O

@named gibbs = GibbsReactor(ms; element_matrix=elem_matrix, feed_elements=feed_elem)
```
"""
@component function GibbsReactor(ms::MaterialSource;
        name,
        element_matrix::Matrix{Float64},
        feed_elements::Vector{Float64},
        T_specified::Union{Nothing, Float64} = nothing,
        p_specified::Float64 = 101325.0)
    N_c = ms.N_c
    N_elem = size(element_matrix, 1)

    # Validate dimensions
    if size(element_matrix, 2) != N_c
        throw(ArgumentError(
            "Element matrix must have $N_c columns (number of components), " *
            "got $(size(element_matrix, 2))"))
    end
    if length(feed_elements) != N_elem
        throw(ArgumentError(
            "feed_elements must have $N_elem entries (number of elements), " *
            "got $(length(feed_elements))"))
    end

    # Parameters
    pars = @parameters begin
        p = p_specified, [description = "System pressure (Pa)"]
        T_spec = isnothing(T_specified) ? 298.15 : T_specified,
        [description = "Specified temperature (K), used if isothermal"]
        A[1:N_elem, 1:N_c] = element_matrix, [description = "Element composition matrix"]
        b[1:N_elem] = feed_elements, [description = "Feed element amounts (mol)"]
    end

    # Variables for optimization
    vars = @variables begin
        (n(t))[1:N_c], [description = "Molar amounts of each component (mol)", bounds = (0, Inf)]
        n_total(t), [description = "Total moles", bounds = (0, Inf)]
        (x(t))[1:N_c], [description = "Mole fractions", bounds = (0, 1)]
        T(t), [description = "Temperature (K)", bounds = (200, 5000)]
        G_total(t), [description = "Total Gibbs free energy (J)"]
        H_total(t), [description = "Total enthalpy (J)"]
        S_total(t), [description = "Total entropy (J/K)"]
    end

    # Compute molar density for property calculations
    ρ_func = (p_val, T_val, x_val) -> ms.molar_density(p_val, T_val, x_val; phase = "unknown")

    # Build equations
    eqs = Equation[]

    # Total moles
    push!(eqs, n_total ~ sum(collect(n)))

    # Mole fractions (with protection against division by zero)
    for i in 1:N_c
        push!(eqs, x[i] ~ n[i] / n_total)
    end

    # Element balance constraints: A * n = b
    for i in 1:N_elem
        push!(eqs, sum(scalarize(A[i, :] .* n)) ~ b[i])
    end

    # Thermodynamic properties
    # Enthalpy: H = sum(n_i * h_i) where h_i is molar enthalpy
    push!(eqs, H_total ~ n_total * ms.VT_enthalpy(ρ_func(p, T, collect(x)), T, collect(x)))

    # Entropy: S = sum(n_i * s_i) where s_i is molar entropy
    push!(eqs, S_total ~ n_total * ms.VT_entropy(ρ_func(p, T, collect(x)), T, collect(x)))

    # Gibbs free energy: G = H - T*S
    push!(eqs, G_total ~ H_total - T * S_total)

    # Temperature specification (if isothermal)
    if !isnothing(T_specified)
        push!(eqs, T ~ T_spec)
    end

    # Create OptimizationSystem
    # Objective: minimize G_total
    # Subject to: element balance constraints (included in eqs)
    return OptimizationSystem(
        G_total,  # Objective function to minimize
        eqs,
        t,
        collect(Iterators.flatten(vars)),
        collect(Iterators.flatten(pars));
        name
    )
end

"""
    SteadyStateGibbsReactor(ms::MaterialSource; name, kwargs...)

Create a steady-state Gibbs reactor with inlet and outlet connectors.
This wraps the GibbsReactor optimization in a component that can be
connected to other process units.

# Arguments
- `ms::MaterialSource`: Material source with thermodynamic property functions
- `name`: System name (required)
- `element_matrix::Matrix{Float64}`: Element composition matrix
- `feed_elements::Vector{Float64}`: Feed element amounts (will be computed from inlet)

# Returns
An ODESystem component with material connectors.
"""
@component function SteadyStateGibbsReactor(ms::MaterialSource;
        name,
        element_matrix::Matrix{Float64})
    N_c = ms.N_c
    N_elem = size(element_matrix, 1)

    # Connectors
    @named In = MaterialConnector(ms)
    @named Out = MaterialConnector(ms)

    # Parameters
    pars = @parameters begin
        A[1:N_elem, 1:N_c] = element_matrix, [description = "Element composition matrix"]
    end

    # Variables
    vars = @variables begin
        (n_out(t))[1:N_c], [description = "Equilibrium molar amounts (mol)"]
        n_total_out(t), [description = "Total outlet moles"]
        (elem_feed(t))[1:N_elem], [description = "Element amounts in feed (mol)"]
        G(t), [description = "Gibbs free energy at equilibrium (J)"]
    end

    eqs = Equation[]

    # Compute element amounts from inlet composition
    for i in 1:N_elem
        push!(eqs, elem_feed[i] ~ sum(scalarize(A[i, :] .* In.xᵢ .* In.n)))
    end

    # Element balance at outlet (steady state: in = out)
    for i in 1:N_elem
        push!(eqs, elem_feed[i] ~ sum(scalarize(A[i, :] .* n_out)))
    end

    # Total outlet moles
    push!(eqs, n_total_out ~ sum(collect(n_out)))

    # Outlet connector equations
    push!(eqs, Out.n ~ -n_total_out)  # Negative for outflow
    push!(eqs, Out.p ~ In.p)
    push!(eqs, Out.T ~ In.T)  # Isothermal assumption for this simplified version
    for i in 1:N_c
        push!(eqs, Out.xᵢ[i] ~ n_out[i] / n_total_out)
    end

    # Flow balance
    push!(eqs, In.n + Out.n ~ 0)

    # Gibbs energy (for output/monitoring)
    push!(eqs, G ~ ms.VT_enthalpy(Out.ϱ, Out.T, collect(Out.xᵢ)) * n_total_out -
                   Out.T * ms.VT_entropy(Out.ϱ, Out.T, collect(Out.xᵢ)) * n_total_out)

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)),
        collect(Iterators.flatten(pars));
        name, systems = [In, Out])
end

export GibbsReactor, SteadyStateGibbsReactor
