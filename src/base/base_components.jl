@component function MaterialStream(ms::MaterialSource; phase = "unknown", name)
    @named c1 = MaterialConnector(ms)
    @named c2 = MaterialConnector(ms)
    mcons = [c1, c2]

    vars = @variables begin
        T(t), [description = "temperature"]         #, unit=u"K"]
        ϱ(t), [description = "density"]             #, unit=u"mol m^-3"]
        p(t), [description = "pressure"]            #, unit=u"Pa"]
        h(t), [description = "molar enthalpy"]      #, unit=u"J mol^-1"]
        n(t), [description = "molar flow "]         #, unit=u"mol s^-1"]
        xᵢ(t)[1:ms.N_c], [description = "mole fractions"]      #, unit=u"mol mol^-1"]
    end

    eqs = [
        # EOS
        ϱ ~ ms.molar_density(p, T, xᵢ; phase = phase),
        h ~ ms.VT_enthalpy(ϱ, T, xᵢ),
        1.0 ~ sum(collect(xᵢ)),
        # Connector "1"
        T ~ c1.T,
        ϱ ~ c1.ϱ,
        p ~ c1.p,
        h ~ c1.h,
        n ~ c1.n,
        scalarize(xᵢ .~ c1.xᵢ)...,
        # Connector "2"
        T ~ c2.T,
        ϱ ~ c2.ϱ,
        p ~ c2.p,
        h ~ c2.h,
        n ~ -c2.n,
        scalarize(xᵢ .~ c2.xᵢ)...,
    ]

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)), []; name, systems = mcons)
end

@connector function MaterialConnector(ms::MaterialSource; name)
    vars = @variables begin
        T(t), [description = "temperature", output = true]        #, unit=u"K"]
        ϱ(t), [description = "density", output = true]            #, unit=u"mol m^-3"]
        p(t), [description = "pressure"]                        #, unit=u"Pa"]
        h(t), [description = "molar enthalpy", output = true]     #, unit=u"J mol^-1"]
        xᵢ(t)[1:ms.N_c], [description = "mole fractions", output = true]     #, unit=u"mol mol^-1"]
        n(t), [description = "total molar flow", connect = Flow]  #, unit=u"mol s^-1"]
    end

    return ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)
end

@component function HeatConnector(; name)
    vars = @variables begin
        Q(t), [description = "heat flux"]  #, unit=u"J s^-1"]
    end

    return ODESystem(Equation[], t, vars, []; name)
end

@component function WorkConnector(; name)
    vars = @variables begin
        W(t), [description = "power"]      #, unit=u"J s^1"]
    end

    return ODESystem(Equation[], t, vars, []; name)
end

@component function SimpleControlVolume(
        ms::MaterialSource;
        N_mcons,
        N_heats = 0,
        N_works = 0,
        name
    )
    # Init
    mcons = [MaterialConnector(ms; name = Symbol("c$i")) for i in 1:N_mcons]
    works = [WorkConnector(name = Symbol("w$i")) for i in 1:N_works]
    heats = [HeatConnector(name = Symbol("q$i")) for i in 1:N_heats]

    vars = @variables begin
        ΔH(t), [description = "Enthalpy difference inlets/outlets"]  #, unit=u"J/s"]
        ΔE(t), [description = "Added/removed heat or work"]          #, unit=u"J/s"]
    end

    eqs = Equation[
        ΔH ~ sum([c.h * c.n for c in mcons])
        ΔE ~
            (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) +
            (isempty(works) ? 0.0 : sum([w.W for w in works]))
        # Energy balance
        0.0 ~ ΔH + ΔE
        # Mole balance
        [0.0 ~ sum([c.xᵢ[j] * c.n for c in mcons]) for j in 1:ms.N_c][:]
    ]

    return ODESystem(eqs, t, vars, []; name, systems = [mcons..., works..., heats...])
end

@component function TPControlVolume(
        ms::MaterialSource;
        N_mcons,
        N_heats = 0,
        N_works = 0,
        phases = ["unknown"],
        reactive = false,
        name
    )
    # Init
    N_ph = length(phases)
    mcons = [MaterialConnector(ms; name = Symbol("c$i")) for i in 1:N_mcons]
    works = [WorkConnector(name = Symbol("w$i")) for i in 1:N_works]
    heats = [HeatConnector(name = Symbol("q$i")) for i in 1:N_heats]

    vars = @variables begin
        T(t), [description = "temperature", bounds = (0, Inf)]         #, unit=u"K"]
        p(t), [description = "pressure", bounds = (0, Inf)]            #, unit=u"Pa"]
        ϱ(t)[1:N_ph], [description = "density", bounds = (0, Inf)]             #, unit=u"mol m^-3"]
        (nᵢ(t))[1:N_ph, 1:ms.N_c], [description = "molar holdup", bounds = (0, Inf)]        #, unit=u"mol"]
        (xᵢ(t))[1:N_ph, 1:ms.N_c], [description = "mole fractions", bounds = (0, 1)]        #, unit=u"mol mol^-1"]
        n(t), [description = "total molar holdup", bounds = (0, Inf)]  #, unit=u"mol"]
        U(t), [description = "internal energy"]                     #, unit=u"J"]
        ΔH(t), [description = "enthalpy difference inlets/outlets"]  #, unit=u"J/s"]
        ΔE(t), [description = "added/removed heat or work"]          #, unit=u"J/s"]
        V(t), [description = "volume", bounds = (0, Inf)]             #, unit=u"m^3"]
    end
    reactive ?
        append!(
            vars, @variables begin
                ΔnR(t)[1:ms.N_c], [description = "molar holdup change by reaction"]     #, unit=u"mol"])
                ΔHᵣ(t), [description = "enthalpy of reaction"]                #, unit=u"J/s"]
            end
        ) : nothing

    eqs = [
        ΔH ~ sum([c.h * c.n for c in mcons]),
        ΔE ~
            (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) +
            (isempty(works) ? 0.0 : sum([w.W for w in works])),
        # Energy balance
        D(U) ~ ΔH + ΔE + (reactive ? ΔHᵣ : 0.0),
        # Mole balance
        [
            D(sum(collect(nᵢ[:, i]))) ~
                sum([c.xᵢ[i] * c.n for c in mcons]) +
                (reactive ? ΔnR[i] : 0.0) for i in 1:ms.N_c
        ]...,
        n ~ sum(collect([nᵢ...])),
        [xᵢ[j, i] ~ nᵢ[j, i] / sum(collect(nᵢ[j, :])) for i in 1:ms.N_c, j in 1:N_ph]...,
        # Thermodynamic system properties
        [
            ϱ[j] ~ ms.molar_density(p, T, collect(xᵢ[j, :]); phase = phases[j])
                for j in 1:N_ph
        ]...,
        U ~ sum(
            [
                ms.VT_internal_energy(ϱ[j], T, collect(xᵢ[j, :])) * sum(collect(nᵢ[j, :]))
                    for j in 1:N_ph
            ]
        ),
        V ~ n / sum(collect(ϱ)),
    ]

    return ODESystem(
        eqs, t, collect(Iterators.flatten(vars)), [];
        name, systems = [mcons..., works..., heats...]
    )
end
