@component function SimpleAdiabaticCompressor(ms::MaterialSource; name)

    # Subsystems
    @named cv = SimpleControlVolume(ms; N_mcons = 2, N_works = 1)

    # Variables
    vars = @variables begin
        W(t), [description="power"]       #, unit=u"J s^1"]
        T2_s(t), [description="temperature"] #, unit=u"K"]
        ϱ2_s(t), [description="density"]     #, unit=u"mol m^-3"]
    end

    pars = @parameters begin
        ηᴱ = 1.0, [description="efficiency"]
    end

    # Equations
    eqs = [
        cv.w1.W ~ W,
        ms.VT_entropy(cv.c1.ϱ, cv.c1.T, cv.c1.xᵢ) ~
        ms.VT_entropy(cv.c2.ϱ, cv.c2.T, cv.c2.xᵢ),
        ϱ2_s ~ ms.molar_density(cv.c2.p, T2_s, cv.c2.xᵢ),
        ηᴱ ~ cv.w1.W / ((ms.VT_enthalpy(ϱ2_s, T2_s, cv.c2.xᵢ) - cv.c1.h)*cv.c1.n)
    ]

    return ODESystem(eqs, t, vars, pars; name, systems = [cv])
end
