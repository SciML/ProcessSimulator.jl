@component function SimpleIsobaricHeatExchanger(ms::MaterialSource; name)

    # Subsystems
    @named cv = SimpleControlVolume(ms;N_states=2,N_heats=1)
    sys = [cv]

    # Variables
    vars = @variables begin
    end

    # Equations
    eqs = [
        cv.s1.p ~ cv.s2.p,
    ]

    return ODESystem(eqs, t, vars, []; systems=sys, name)
end