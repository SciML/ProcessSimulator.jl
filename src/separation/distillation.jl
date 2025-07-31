@component function DistillationColumn(ms::MaterialSource; N_stages, i_feeds, name)
    # Subsystems
    stages = [SimpleStage(ms; name = "stage$i", add_flows = sum(i .==
                                                                [
                  1, N_stages])+sum(i .== i_feeds)) for i in 1:N_stages]

    # Connect stages
    # ...

    return ODESystem(; name)
end

@component function SimpleStage(ms::MaterialSource; name, add_flows)
    # Subsystem
    @named cv = TPControlVolume(ms; N_states = 4)

    # VLE 
    # ...

    return ODESystem(; name)
end
