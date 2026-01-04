@component function SimpleIsobaricHeatExchanger(ms::MaterialSource; name)

    # Subsystems
    @named cv = SimpleControlVolume(ms; N_mcons = 2, N_heats = 1)

    # Variables
    vars = @variables begin
        Q(t), [description = "heat flux"] #, unit=u"J s^-1"]
    end

    # Equations
    eqs = [
        cv.q1.Q ~ Q,
        cv.c1.p ~ cv.c2.p,
    ]

    return ODESystem(eqs, t, vars, []; systems = [cv], name)
end
