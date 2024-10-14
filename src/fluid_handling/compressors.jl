@component function SimpleAdiabaticCompressor(ms::MaterialSource; name)

    # Subsystems
    @named cv = SimpleControlVolume(ms;N_states=2,N_works=1)
    sys = [cv]

    # Variables
    vars = @variables begin
        T2_s(t),            [description="temperature", output=true]
        ϱ2_s(t),            [description="density", output=true]
    end

    pars = @parameters begin
        ηᴱ,                 [description="efficiency", output=true]
    end

    # Equations
    eqs = [
        ms.VT_entropy(cv.f1.s.ϱ,cv.f1.s.T,cv.f1.c.xᵢ) ~ ms.VT_entropy(cv.f2.s.ϱ,cv.f2.s.T,cv.f2.c.xᵢ)
        ϱ2_s ~ ms.molar_density(cv.f2.c.p,T2_s,cv.f2.c.xᵢ)
        ηᴱ ~ cv.Ws[1] / (ms.VT_enthalpy(ϱ2_s,T2_s,cv.f2.c.xᵢ) - ms.VT_enthalpy(cv.f1.s.ϱ,cv.f1.s.T,cv.f1.c.xᵢ))
    ]

    return ODESystem(eqs, t, vars, pars; name, systems=sys)
end