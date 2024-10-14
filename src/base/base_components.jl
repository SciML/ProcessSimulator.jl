@component function MaterialFlow(ms::MaterialSource;N_c=1,phase="unknown",name)

    @named s = TϱState(N_c=N_c)
    @named c = PHxConnector(ms;N_c=N_c,phase=phase)
    sub = [s,c]

    eqs = [
        c.p ~ ms.pressure(s.ϱ,s.T,c.xᵢ)
        c.h ~ ms.VT_enthalpy(s.ϱ,s.T,c.xᵢ)
        1.0 ~ sum([xᵢ for xᵢ in c.xᵢ])
        scalarize(c.xᵢ .~ s.nᵢ ./ c.n)
    ]

    ODESystem(eqs, t, [], []; name, systems=sub)
end

@component function TϱState(;N_c=1, name)
    vars = @variables begin
        T(t),               [description="temperature", output=true]    #, unit=u"K"]
        ϱ(t),               [description="density", output=true]        #, unit=u"mol m^-3"]
        nᵢ(t)[1:N_c],       [description="mole fractions", output=true] #, unit=u"mol s^-1"]
    end

    return ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)
end

@connector function PHxConnector(ms::MaterialSource;N_c=1,phase="unknown",name)
    vars = @variables begin
        p(t),               [description="pressure"]                        #, unit=u"Pa"]
        h(t),               [description="molar enthalpy", connect=Stream]  #, unit=u"J mol^-1"]
        xᵢ(t)[1:N_c]=1/N_c, [description="mole fractions", connect=Stream]  #, unit=u"mol mol^-1"]
        n(t),               [description="total molar flow", connect=Flow]  #, unit=u"mol s^-1"]
    end

    return ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)
end

@connector function HeatConnector(;name)
    vars = @variables begin
        Q(t)#,              [description="heat flux", connect=Flow]  #, unit=u"J s^-1"]
    end

    return ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)
end

@connector function WorkConnector(;name)
    vars = @variables begin
        W(t)#,              [description="power", connect=Flow]      #, unit=u"J s^1"]
    end

    return ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)
end

@component function TPControlVolume(ms::MaterialSource;
                                    N_states,
                                    N_heats=0,
                                    N_works=0,
                                    reactions=Vector{ODESystem}[],
                                    N_ph=1, phases=repeat(["unknown"],N_ph),
                                    name)
    # Init 
    N_c = length(ms.components)
    flows = [MaterialFlow(ms;N_c=N_c,name=Symbol("s$i")) for i in 1:N_states]
    works = [WorkConnector(name=Symbol("w$i")) for i in 1:N_works]
    heats = [HeatConnector(name=Symbol("q$i")) for i in 1:N_heats]

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

    eqs = [
        ΔH ~ sum([ms.VT_enthalpy(st.ϱ,st.T,st.nᵢ) for st in states])
        ΔHᵣ ~ 0.0 # TODO: Reactions -> sum([reaction.ΔH for reaction in reactions])
        ΔE ~ (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) + (isempty(works) ? 0.0 : sum([w.W for w in works]))
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
    flows = [MaterialFlow(ms;N_c=N_c,name=Symbol("s$i")) for i in 1:N_states]
    works = [WorkConnector(name=Symbol("w$i")) for i in 1:N_works]
    heats = [HeatConnector(name=Symbol("q$i")) for i in 1:N_heats]

    vars = @variables begin
        ΔH(t), [description="Enthalpy difference inlets/outlets", unit="J/s", output=true]
        ΔE(t), [description="Added/removed heat or work", unit="J/s", output=true]
    end

    eqs = [
        ΔH ~ sum([f.c.h*f.c.n for f in flows])
        ΔE ~ (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) + (isempty(works) ? 0.0 : sum([w.W for w in works]))
        # Energy balance
        0.0 ~ ΔH + ΔE
        # Mole balance
        [0.0 ~ sum([f.s.nᵢ[j] for f in flows]) for j in 1:N_c][:]
    ]

    ODESystem([eqs...], t, collect(Iterators.flatten(vars)), []; name, systems=states)
end