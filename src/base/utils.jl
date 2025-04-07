#= @component function Port(ms;phase="unknown", name)
    @named c = PhXConnector(ms)

    vars = @variables begin
        T(t),               [description="inlet temperature"]   #, unit=u"K"]
        ϱ(t),               [description="inlet density"]       #, unit=u"mol m^-3"]
        p(t),               [description="inlet pressure"]      #, unit=u"Pa"]
        n(t),               [description="inlet molar flow"]    #, unit=u"mol s^-1"]
        m(t),               [description="inlet mass flow"]     #, unit=u"kg s^-1"]
        xᵢ(t)[1:ms.N_c],    [description="inlet mole fractions"]#, unit=u"mol mol^-1"] 
    end

    eqs = Equation[
        # EOS 
        ϱ ~ ms.molar_density(p,T,xᵢ;phase=phase),
        1.0 ~ sum(collect(xᵢ)),
        c.h ~ ms.VT_enthalpy(ϱ,T,xᵢ),
        # Connector
        T ~ c.T,
        ϱ ~ c.ϱ,
        p ~ c.p,
        scalarize(xᵢ .~ c.xᵢ)...,
        n ~ m / sum(ms.Mw[i] * xᵢ[i] for i in 1:ms.N_c),
        0.0 ~ n + c.n,
     ]

    return ODESystem(eqs, t, collect(Iterators.flatten(vars)), []; name, systems=[c])
end
 =#