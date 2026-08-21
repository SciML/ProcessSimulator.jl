using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    reaction_rate = (p, T, xᵢ) -> 0.1 * xᵢ[1]
    reaction_enthalpy = T -> -5.0e4
    molar_density = (p, T, xᵢ; kwargs...) -> p / (8.314 * T)
    enthalpy = (ϱ, T, xᵢ) -> 20.0 * T
    internal_energy = (ϱ, T, xᵢ) -> 12.0 * T
    entropy = (ϱ, T, xᵢ) -> log(T)

    @compile_workload begin
        reaction = Reaction(
            ν = [-1.0, 1.0], r = reaction_rate, Δhᵣ = reaction_enthalpy
        )
        material = MaterialSource(
            ["reactant", "product"];
            Mw = [0.01, 0.02],
            molar_density,
            VT_enthalpy = enthalpy,
            VT_internal_energy = internal_energy,
            VT_entropy = entropy,
            reactions = [reaction],
        )

        Port(material; name = :inlet)
        MaterialStream(material; name = :stream)
        SimpleAdiabaticCompressor(material; name = :compressor)
        SimpleIsobaricHeatExchanger(material; name = :exchanger)
        CSTR(material; name = :reactor, flowtype = "const. volume")
    end
end
