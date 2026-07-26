"""
    SimpleIsobaricHeatExchanger(ms::MaterialSource; name)

Create a two-sided heat-exchanger component with an isobaric material path and one
heat connector.

# Arguments

  - `ms`: Material source used for the component's thermodynamic relations.

# Keywords

  - `name`: Required ModelingToolkit component name; use
    `@named exchanger = SimpleIsobaricHeatExchanger(ms)`.

# Examples

```jldoctest
julia> source = MaterialSource("helium"; Mw = 0.004,
           molar_density = (p, T, xᵢ; kwargs...) -> p / (8.314 * T),
           VT_enthalpy = (ϱ, T, xᵢ) -> 20.0 * T);

julia> exchanger = SimpleIsobaricHeatExchanger(source; name = :exchanger);

julia> nameof(exchanger)
:exchanger
```
"""
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

@public SimpleIsobaricHeatExchanger
