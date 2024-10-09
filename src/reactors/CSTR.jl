@component function CSTR(ms:MaterialSource,reac::Reaction;name)
    # Subsystems
    @named cv = SimpleControlVolume(ms;)

    # Stoichiometry

    return ODESystem()
end