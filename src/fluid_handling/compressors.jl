@component function SimpleAdiabaticCompressor(ms::MaterialSource; name)

    # Subsystems
    @named cv = SimpleControlVolume(ms;N_states=2,N_works=1)
    @named cv_s = SimpleControlVolume(ms;N_states=2,N_works=1)
    sys = [cv,cv_s]

    # Variables
    vars = @variables begin
        ηᴱ(t),              [description="efficiency", unit="-", output=true]
    end

    pars = []

    # Equations
    eqs = [
        cv.s1.p ~ cv_s.s1.p,
        cv.s1.T ~ cv_s.s1.T,
        cv.s2.p ~ cv_s.s2.p,
        scalarize(cv.s1.nᵢ .~ cv_s.s1.nᵢ)...,
        ms.VT_entropy(cv.s1.ϱ,cv.s1.T,cv.s1.nᵢ) ~ ms.VT_entropy(cv_s.s2.ϱ,cv_s.s2.T,cv_s.s2.nᵢ),
        ηᴱ ~ cv.Ws[1] / cv_s.Ws[1]
    ]
    eqs = Equation[eqs...]

    return ODESystem([eqs...], t, vars, pars; name, systems=sys)
end