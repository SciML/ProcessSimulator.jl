# Base components

@connector function State(ms::MaterialSource;N_c=1,n_bnds=(-Inf,Inf),phase="unknown",name)
    vars = @variables begin
        T(t),               [description="temperature", unit="K", output=true]
        p(t),               [description="pressure", unit="Pa", output=true]
        ϱ(t),               [description="density", unit="mol/m³", output=true]
        (nᵢ(t))[1:N_c],     [description="molar flow", unit="mol/s", bounds=n_bnds, output=true]
        n(t),               [description="total molar flow", unit="mol/s", output=true]
    end

    eqs = [
        # Density
        ϱ ~ ms.molar_density(p,T,nᵢ;phase=phase)
        n ~ sum(nᵢ)
    ]

    pars = @parameters begin
    end

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name)
end

"""
    MaterialStream(ms::MaterialSource;N_c=1, phase="unknown", flowdir=:generic, name)

Create a material stream connector (between points `A` and `B`). Constant pressure and temperature are assumed. 
...
"""
@connector function MaterialStream(ms::MaterialSource;N_c=1, phase="unknown", flowdir=:AB, name)
    if flowdir == :bidir
        boundsA = boundsB = (-Inf,Inf)
    elseif flowdir == :AB
        boundsA = (-Inf,0)
        boundsB = (0,Inf)
    elseif flowdir == :BA
        boundsA = (0,Inf)
        boundsB = (-Inf,0)
    else
        error("Invalid type")
    end

    @named A = State(ms;N_c=N_c,n_bnds=boundsA,phase=phase)
    @named B = State(ms;N_c=N_c,n_bnds=boundsB,phase=phase)

    eqs = [
        scalarize(0.0 .~ A.nᵢ .+ B.nᵢ)...   # Mass balance
        A.T ~ B.T                           # Thermal equilibrium
        A.p ~ B.p                           # Mechanical equilibrium
    ]
    
    ODESystem(eqs, t; systems=[A,B], name)
end

"""
    EnergyStream(;name)

Create an energy stream connector (between poinnts `A` and `B`).
...
"""
@connector function HeatStream(; name)
    vars = @variables begin
        QA(t),              [description="heat A", unit="J/s", output=true]
        QB(t),              [description="heat B", unit="J/s", output=true]
    end

    eqs = [
        0.0 ~ QA + QB
    ]

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), []; name)
end


@connector function WorkStream(; name)
    vars = @variables begin
        WA(t),              [description="work A", unit="J/s", output=true]
        WB(t),              [description="work B", unit="J/s", output=true]
    end

    eqs = [
        0.0 ~ WA + WB
    ]

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), []; name)
end

"""
    ControlVolume(ms::MaterialSource; ... name)

Create a control volume with constant pressure and temperature and chemical equilibrium.
...
""",
@component function TPControlVolume(ms::MaterialSource;
                                    N_states,
                                    N_heats=0,
                                    N_works=0,
                                    reactions=Vector{ODESystem}[],
                                    N_ph=1, phases=repeat(["unknown"],N_ph),
                                    name)
    # Init 
    N_c = length(ms.components)
    states = [State(ms;N_c=N_c) for i in 1:N_states]

    vars = @variables begin
        T(t),                   [description="temperature", unit="K", output=true, bounds=(0,Inf)]
        p(t),                   [description="pressure", unit="Pa", output=true, bounds=(0,Inf)]
        ϱ(t)[1:N_ph],           [description="density", unit="mol/m³", output=true, bounds=(0,Inf)]
        (nᵢ(t))[1:N_ph,1:N_c],  [description="molar holdup", unit="mol", output=true, bounds=(0,Inf)]
        n(t),                   [description="total molar holdup", unit="mol", output=true, bounds=(0,Inf)]
        U(t),                   [description="internal energy", unit="J", output=true]
        ΔH(t),                  [description="enthalpy difference inlets/outlets", unit="J/s", output=true]
        ΔHᵣ(t),                 [description="enthalpy difference reactions", unit="J/s", output=true]
        ΔE(t),                  [description="added/removed heat or work", unit="J/s", output=true]
    end
    if N_heats > 0
        @variables (Qs(t))[1:N_heats], [description="Heats added/removed", unit="J/s", output=true]
        vars = vcat(vars, Qs)
    else
        Qs = 0.0
    end
    if N_works > 0
        @variables (Ws(t))[1:N_works], [description="Work added/removed", unit="J/s", output=true]
        vars = vcat(vars, Ws)
    else
        Ws = 0.0
    end

    eqs = [
        ΔH ~ sum([ms.VT_enthalpy(st.ϱ,st.T,st.nᵢ) for st in states])
        ΔHᵣ ~ 0.0 # TODO: Reactions -> sum([reaction.ΔH for reaction in reactions])
        ΔE ~ sum(Qs) + sum(Ws)
        # Energy balance
        D(U) ~ ΔH + ΔE + ΔHᵣ
        # Mole balance
        [D(sum(nᵢ[:,j])) ~ sum([st.nᵢ[j] for st in states]) for j in 1:N_c][:]
        n ~ sum(nᵢ)
        # TODO: Definition of chemical equilibrium and reactions outside of the control volume obejct?
        # Thermodynamic system properties
        [ϱ[i] ~ ms.molar_density(p,T,nᵢ[i,:];phase=phases[i]) for i in 1:N_ph]
        U ~ sum([ms.VT_internal_energy(ϱ[i],T,nᵢ[i,:]) for i in 1:N_ph])
    ]

    ODESystem([eqs...], t, collect(Iterators.flatten(vars)), []; name, systems=sys)
end

@component function SimpleControlVolume(ms::MaterialSource;
                                        N_states,
                                        N_heats=0,
                                        N_works=0,
                                        N_ph=1, phases=repeat(["unknown"],N_ph),
                                        name)
    # Init 
    N_c = length(ms.components)
    states = [State(ms;N_c=N_c,name=Symbol("s$i")) for i in 1:N_states]

    vars = @variables begin
        ΔH(t), [description="Enthalpy difference inlets/outlets", unit="J/s", output=true]
        ΔE(t), [description="Added/removed heat or work", unit="J/s", output=true]
    end
    if N_heats > 0
        @variables (Qs(t))[1:N_heats], [description="Heats added/removed", unit="J/s", output=true]
        vars = vcat(vars, Qs)
    else
        Qs = zeros(1)
    end
    if N_works > 0
        @variables (Ws(t))[1:N_works], [description="Work added/removed", unit="J/s", output=true]
        vars = vcat(vars, Ws)
    else
        Ws = zeros(1)
    end

    eqs = [
        ΔH ~ sum([ms.VT_enthalpy(st.ϱ,st.T,st.nᵢ)*sign(st.nᵢ[1]) for st in states])
        ΔE ~ sum([Q for Q in Qs]) + sum([W for W in Ws])
        # Energy balance
        0.0 ~ ΔH + ΔE
        # Mole balance
        [0.0 ~ sum([st.nᵢ[j] for st in states]) for j in 1:N_c][:]
    ]

    ODESystem([eqs...], t, collect(Iterators.flatten(vars)), []; name, systems=states)
end