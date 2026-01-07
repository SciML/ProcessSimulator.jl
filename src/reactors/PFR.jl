"""
    PFR - Plug Flow Reactor

The Plug Flow Reactor (PFR) models a tubular reactor with no axial mixing.
The composition varies along the reactor length as reactions proceed.
This implementation uses spatial discretization to convert the PDEs into
a system of ODEs.

Key assumptions:
- No radial gradients (1D model)
- No axial dispersion (plug flow)
- Steady-state or dynamic operation
- Constant cross-sectional area (can be extended)
"""

"""
    PFR(ms::MaterialSource; name, L, A, N_nodes, i_reacts, adiabatic)

Create a Plug Flow Reactor component.

# Arguments
- `ms::MaterialSource`: Material source with thermodynamic properties and reactions
- `name`: System name
- `L::Float64`: Reactor length (m)
- `A::Float64`: Cross-sectional area (m²)
- `N_nodes::Int`: Number of discretization nodes along length
- `i_reacts`: Indices of reactions to include (default: all)
- `adiabatic::Bool`: If true, no heat transfer; if false, Q must be specified (default: true)

# Returns
An `ODESystem` representing the discretized PFR.

# Notes
The PFR is discretized into N_nodes spatial nodes. At each node, mass and
energy balances are solved. The spatial coordinate z varies from 0 to L.

For steady-state operation, use with `structural_simplify` and set
time derivatives to zero.

# Example
```julia
@named pfr = PFR(ms; L=2.0, A=0.01, N_nodes=20)
```
"""
@component function PFR(ms::MaterialSource;
        name,
        L::Float64,
        A::Float64,
        N_nodes::Int = 20,
        i_reacts = 1:length(ms.reaction),
        adiabatic::Bool = true)
    # Validate
    N_nodes >= 2 || throw(ArgumentError("N_nodes must be at least 2"))
    L > 0 || throw(ArgumentError("Reactor length L must be positive"))
    A > 0 || throw(ArgumentError("Cross-sectional area A must be positive"))

    N_c = ms.N_c
    N_r = length(ms.reaction[i_reacts])
    Δz = L / (N_nodes - 1)  # Node spacing
    V_node = A * Δz  # Volume per node

    # Connectors
    @named In = MaterialConnector(ms)
    @named Out = MaterialConnector(ms)

    # Parameters
    pars = @parameters begin
        L_reactor = L, [description = "Reactor length (m)"]
        A_cross = A, [description = "Cross-sectional area (m²)"]
        dz = Δz, [description = "Node spacing (m)"]
        V_seg = V_node, [description = "Volume per segment (m³)"]
    end

    # Variables at each node
    vars = @variables begin
        # State variables at each node
        (T_node(t))[1:N_nodes], [description = "Temperature at each node (K)"]
        (p_node(t))[1:N_nodes], [description = "Pressure at each node (Pa)"]
        (xᵢ_node(t))[1:N_nodes, 1:N_c], [description = "Mole fractions at each node"]
        (ϱ_node(t))[1:N_nodes], [description = "Molar density at each node (mol/m³)"]
        (n_node(t))[1:N_nodes], [description = "Molar flow at each node (mol/s)"]
        (h_node(t))[1:N_nodes], [description = "Molar enthalpy at each node (J/mol)"]

        # Reaction rates at each node
        (r_node(t))[1:N_nodes, 1:N_r], [description = "Reaction rates at each node (mol/s)"]
        (R_node(t))[1:N_nodes, 1:N_c], [description = "Species production rates (mol/m³/s)"]
        (ΔHᵣ_node(t))[1:N_nodes], [description = "Heat of reaction at each node (J/s)"]

        # Outlet variables
        F_out(t), [description = "Outlet molar flow (mol/s)"]
        T_out(t), [description = "Outlet temperature (K)"]

        # Heat transfer (if not adiabatic)
        (Q_node(t))[1:N_nodes], [description = "Heat transfer at each node (J/s)"]
    end

    eqs = Equation[]

    # Boundary conditions at inlet (node 1)
    push!(eqs, T_node[1] ~ In.T)
    push!(eqs, p_node[1] ~ In.p)
    push!(eqs, n_node[1] ~ In.n)
    for i in 1:N_c
        push!(eqs, xᵢ_node[1, i] ~ In.xᵢ[i])
    end

    # Equations at each node
    for j in 1:N_nodes
        # Density from EOS
        push!(eqs, ϱ_node[j] ~ ms.molar_density(p_node[j], T_node[j],
            collect(xᵢ_node[j, :]); phase = "unknown"))

        # Enthalpy from EOS
        push!(eqs, h_node[j] ~ ms.VT_enthalpy(ϱ_node[j], T_node[j], collect(xᵢ_node[j, :])))

        # Reaction rates at this node
        for (ir, reac) in enumerate(ms.reaction[i_reacts])
            push!(eqs, r_node[j, ir] ~ reac.r(p_node[j], T_node[j], collect(xᵢ_node[j, :])))
        end

        # Species production rates
        for i in 1:N_c
            push!(eqs, R_node[j, i] ~ sum([r_node[j, ir] * ms.reaction[i_reacts][ir].ν[i]
                                           for ir in 1:N_r]))
        end

        # Heat of reaction
        push!(eqs, ΔHᵣ_node[j] ~ sum([r_node[j, ir] * ms.reaction[i_reacts][ir].Δhᵣ(T_node[j])
                                      for ir in 1:N_r]) * V_seg)

        # Heat transfer (adiabatic or specified)
        if adiabatic
            push!(eqs, Q_node[j] ~ 0.0)
        end
    end

    # Spatial derivatives using finite differences (upwind scheme)
    # dF_i/dz = R_i * A  (component molar flow)
    # dT/dz depends on energy balance
    for j in 2:N_nodes
        # Component balances: d(n*x_i)/dz = R_i * A
        # Using upwind finite difference: (n*x_i)[j] - (n*x_i)[j-1] = R_i * A * dz
        for i in 1:N_c
            push!(eqs, n_node[j] * xᵢ_node[j, i] - n_node[j - 1] * xᵢ_node[j - 1, i] ~
                       R_node[j, i] * V_seg)
        end

        # Total molar flow (sum of component flows)
        push!(eqs, n_node[j] ~ sum([n_node[j] * xᵢ_node[j, i] for i in 1:N_c]))

        # Pressure drop (simplified: assume small or specify separately)
        # For now, assume isobaric
        push!(eqs, p_node[j] ~ p_node[j - 1])

        # Energy balance: d(n*h)/dz = ΔHᵣ * A + Q
        push!(eqs, n_node[j] * h_node[j] - n_node[j - 1] * h_node[j - 1] ~
                   ΔHᵣ_node[j] + Q_node[j])
    end

    # Normalize mole fractions (ensure sum = 1)
    for j in 2:N_nodes
        push!(eqs, 1.0 ~ sum([xᵢ_node[j, i] for i in 1:N_c]))
    end

    # Outlet conditions
    push!(eqs, F_out ~ n_node[N_nodes])
    push!(eqs, T_out ~ T_node[N_nodes])

    # Outlet connector
    push!(eqs, Out.n ~ -F_out)
    push!(eqs, Out.T ~ T_out)
    push!(eqs, Out.p ~ p_node[N_nodes])
    push!(eqs, Out.ϱ ~ ϱ_node[N_nodes])
    push!(eqs, Out.h ~ h_node[N_nodes])
    for i in 1:N_c
        push!(eqs, Out.xᵢ[i] ~ xᵢ_node[N_nodes, i])
    end

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)),
        collect(Iterators.flatten(pars));
        name, systems = [In, Out])
end

"""
    SteadyStatePFR(ms::MaterialSource; name, L, A, N_nodes, i_reacts, adiabatic)

Create a steady-state Plug Flow Reactor. This is a convenience wrapper that
sets up the PFR for steady-state analysis (d/dt = 0).

All arguments are the same as PFR.
"""
@component function SteadyStatePFR(ms::MaterialSource;
        name,
        L::Float64,
        A::Float64,
        N_nodes::Int = 20,
        i_reacts = 1:length(ms.reaction),
        adiabatic::Bool = true)
    # For steady state, we simply create the PFR - the time derivatives
    # will be zero in the discretized spatial equations
    return PFR(ms; name, L, A, N_nodes, i_reacts, adiabatic)
end

"""
    SimplePFR(ms::MaterialSource; name, L, A, Q_total)

A simplified PFR model without spatial discretization. Uses a single
control volume with residence time based on volume and flow rate.
Suitable for preliminary calculations.

# Arguments
- `ms::MaterialSource`: Material source
- `name`: System name
- `L::Float64`: Reactor length (m)
- `A::Float64`: Cross-sectional area (m²)
- `Q_total::Float64`: Total volumetric flow rate (m³/s) - optional, can be from inlet

# Returns
An `ODESystem` for the simplified PFR.
"""
@component function SimplePFR(ms::MaterialSource;
        name,
        L::Float64,
        A::Float64,
        i_reacts = 1:length(ms.reaction))
    N_c = ms.N_c
    N_r = length(ms.reaction[i_reacts])
    V_reactor = L * A

    # Connectors
    @named In = MaterialConnector(ms)
    @named Out = MaterialConnector(ms)

    # Parameters
    pars = @parameters begin
        V = V_reactor, [description = "Reactor volume (m³)"]
        L_reactor = L, [description = "Reactor length (m)"]
    end

    # Variables
    vars = @variables begin
        τ(t), [description = "Residence time (s)"]
        T(t), [description = "Outlet temperature (K)"]
        p(t), [description = "Pressure (Pa)"]
        (xᵢ(t))[1:N_c], [description = "Outlet mole fractions"]
        ϱ(t), [description = "Molar density (mol/m³)"]
        h(t), [description = "Molar enthalpy (J/mol)"]
        n_out(t), [description = "Outlet molar flow (mol/s)"]
        (r(t))[1:N_r], [description = "Reaction rates (mol/s/m³)"]
        (R(t))[1:N_c], [description = "Species production rates (mol/s/m³)"]
        ΔHᵣ(t), [description = "Heat of reaction (J/s)"]
    end

    eqs = Equation[]

    # Residence time
    push!(eqs, τ ~ V * ϱ / In.n)

    # Pressure (isobaric)
    push!(eqs, p ~ In.p)

    # Density
    push!(eqs, ϱ ~ ms.molar_density(p, T, collect(xᵢ); phase = "unknown"))

    # Enthalpy
    push!(eqs, h ~ ms.VT_enthalpy(ϱ, T, collect(xᵢ)))

    # Reaction rates (using average composition)
    for (ir, reac) in enumerate(ms.reaction[i_reacts])
        # Use inlet composition for rate calculation (first-order approximation)
        push!(eqs, r[ir] ~ reac.r(p, T, collect(In.xᵢ)))
    end

    # Species production rates
    for i in 1:N_c
        push!(eqs, R[i] ~ sum([r[ir] * ms.reaction[i_reacts][ir].ν[i] for ir in 1:N_r]))
    end

    # Component balance (steady state): F_in * x_in + R * V = F_out * x_out
    for i in 1:N_c
        push!(eqs, In.n * In.xᵢ[i] + R[i] * V ~ n_out * xᵢ[i])
    end

    # Total molar balance
    push!(eqs, n_out ~ In.n + sum([R[i] * V for i in 1:N_c]))

    # Heat of reaction
    push!(eqs, ΔHᵣ ~ sum([r[ir] * ms.reaction[i_reacts][ir].Δhᵣ(T) * V for ir in 1:N_r]))

    # Energy balance (adiabatic): F_in * h_in = F_out * h_out + ΔHᵣ
    push!(eqs, In.n * In.h ~ n_out * h + ΔHᵣ)

    # Mole fraction normalization
    push!(eqs, 1.0 ~ sum([xᵢ[i] for i in 1:N_c]))

    # Outlet connector
    push!(eqs, Out.n ~ -n_out)
    push!(eqs, Out.T ~ T)
    push!(eqs, Out.p ~ p)
    push!(eqs, Out.ϱ ~ ϱ)
    push!(eqs, Out.h ~ h)
    for i in 1:N_c
        push!(eqs, Out.xᵢ[i] ~ xᵢ[i])
    end

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)),
        collect(Iterators.flatten(pars));
        name, systems = [In, Out])
end

export PFR, SteadyStatePFR, SimplePFR
