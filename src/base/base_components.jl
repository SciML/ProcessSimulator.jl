
@component function MaterialStream(ms::MaterialSource;phase="unknown",name)
    @named c1 = MaterialConnector(ms)
    @named c2 = MaterialConnector(ms)
    mcons = [c1,c2]
    
    vars = @variables begin
        T(t),               [description="temperature"]         #, unit=u"K"]
        ϱ(t),               [description="density"]             #, unit=u"mol m^-3"]
        p(t),               [description="pressure"]            #, unit=u"Pa"]
        h(t),               [description="molar enthalpy"]      #, unit=u"J mol^-1"]
        s(t),               [description="molar entropy"]       #, unit=u"J mol^-1 K^-1"]
        n(t),               [description="molar flow "]         #, unit=u"mol s^-1"]
        xᵢ(t)[1:ms.N_c],    [description="mole fractions"]      #, unit=u"mol mol^-1"]
    end

    eqs = [
        # EOS 
        ϱ ~ ms.molar_density(p,T,xᵢ;phase=phase),
        h ~ ms.VT_enthalpy(ϱ,T,xᵢ),
        s ~ ms.VT_entropy(ϱ,T,xᵢ),
        1.0 ~ sum([xᵢ for xᵢ in xᵢ]),
        # Connector "1"
        T ~ c1.T,
        ϱ ~ c1.ϱ,
        p ~ c1.p,
        h ~ c1.h,
        s ~ c1.s,
        n ~ c1.n,
        scalarize(xᵢ .~ c1.xᵢ)...,
        # Connector "2"
        T ~ c2.T,
        ϱ ~ c2.ϱ,
        p ~ c2.p,
        h ~ c2.h,
        s ~ c2.s,
        n ~ -c2.n,
        scalarize(xᵢ .~ c2.xᵢ)...,   
    ]

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)), []; name, systems=mcons)
end

@connector function MaterialConnector(ms::MaterialSource;name)
    vars = @variables begin
        T(t),               [description="temperature", output=true]        #, unit=u"K"]
        ϱ(t),               [description="density", output=true]            #, unit=u"mol m^-3"]
        p(t),               [description="pressure"]                        #, unit=u"Pa"]
        h(t),               [description="molar enthalpy", output=true]     #, unit=u"J mol^-1"]
        s(t),               [description="molar entropy", output=true]      #, unit=u"J mol^-1 K^-1"]
        xᵢ(t)[1:ms.N_c],    [description="mole fractions", output=true]     #, unit=u"mol mol^-1"]
        n(t),               [description="total molar flow", connect=Flow]  #, unit=u"mol s^-1"]
    end

    return ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)
end

@component function HeatConnector(;name)
    vars = @variables begin
        Q(t),              [description="heat flux"]  #, unit=u"J s^-1"]
    end

    return ODESystem(Equation[], t, vars, []; name)
end

@component function WorkConnector(;name)
    vars = @variables begin
        W(t),              [description="power"]      #, unit=u"J s^1"]
    end

    return ODESystem(Equation[], t, vars, []; name)
end

@component function SimpleControlVolume(ms::MaterialSource;
                                        N_mcons,
                                        N_heats=0,
                                        N_works=0,
                                        name)
    # Init 
    mcons = [MaterialConnector(ms;name=Symbol("c$i")) for i in 1:N_mcons]
    works = [WorkConnector(name=Symbol("w$i")) for i in 1:N_works]
    heats = [HeatConnector(name=Symbol("q$i")) for i in 1:N_heats]

    vars = @variables begin
        ΔH(t), [description="Enthalpy difference inlets/outlets"]  #, unit=u"J/s"]
        ΔE(t), [description="Added/removed heat or work"]          #, unit=u"J/s"]
    end

    eqs = Equation[
        ΔH ~ sum([c.h*c.n for c in mcons])
        ΔE ~ (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) + (isempty(works) ? 0.0 : sum([w.W for w in works]))
        # Energy balance
        0.0 ~ ΔH + ΔE
        # Mole balance
        [0.0 ~ sum([c.xᵢ[j]*c.n for c in mcons]) for j in 1:ms.N_c][:]
    ]

    return ODESystem(eqs, t, vars, []; name, systems=[mcons...,works...,heats...])
end

@component function TPControlVolume(ms::MaterialSource;
                                    N_mcons,
                                    N_heats=0,
                                    N_works=0,
                                    N_ph=1, phases=repeat(["unknown"],N_ph),
                                    name)
    # Init 
    mcons = [MaterialConnector(ms;name=Symbol("c$i")) for i in 1:N_mcons]
    works = [WorkConnector(name=Symbol("w$i")) for i in 1:N_works]
    heats = [HeatConnector(name=Symbol("q$i")) for i in 1:N_heats]

    vars = @variables begin
        T(t),                       [description="temperature", bounds=(0,Inf)]        #, unit=u"K"]
        p(t),                       [description="pressure", bounds=(0,Inf)]           #, unit=u"Pa"]
        ϱ(t)[1:N_ph],               [description="density", bounds=(0,Inf)]            #, unit=u"mol m^-3"]
        (nᵢ(t))[1:N_ph,1:ms.N_c],   [description="molar holdup", bounds=(0,Inf)]       #, unit=u"mol"]
        n(t),                       [description="total molar holdup", bounds=(0,Inf)] #, unit=u"mol"]
        U(t),                       [description="internal energy"]                    #, unit=u"J"]
        ΔH(t),                      [description="enthalpy difference inlets/outlets"] #, unit=u"J/s"]
        ΔE(t),                      [description="added/removed heat or work"]         #, unit=u"J/s"]
    end

    pars = @parameters begin
        V,                          [description="volume", bounds=(0,Inf)]             #, unit=u"m^3"]
    end

    eqs = [
        ΔH ~ sum([c.h*c.n for c in mcons])
        ΔE ~ (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) + (isempty(works) ? 0.0 : sum([w.W for w in works]))
        # Energy balance
        D(U) ~ ΔH + ΔE + ΔHᵣ
        # Mole balance
        [D(sum(nᵢ[:,j])) ~ sum([c.xᵢ[j]*c.n for c in mcons]) for j in 1:ms.N_c][:]
        n ~ sum(nᵢ)
        # Thermodynamic system properties
        [ϱ[i] ~ ms.molar_density(p,T,nᵢ[i,:];phase=phases[i]) for i in 1:N_ph]
        U ~ sum([ms.VT_internal_energy(ϱ[i],T,nᵢ[i,:]) for i in 1:N_ph])
    ]

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems=[mcons...,works...,heats...])
end