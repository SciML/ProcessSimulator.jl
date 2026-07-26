"""
    CSTR(ms::MaterialSource; name, i_reacts = 1:length(ms.reaction), flowtype = "")

Create a continuous stirred-tank reactor using the reactions and thermodynamic
functions stored in `ms`.

# Arguments

  - `ms`: Material source that owns component data and [`Reaction`](@ref)s.

# Keywords

  - `name`: Required ModelingToolkit component name; `@named reactor = CSTR(ms)`
    supplies it automatically.
  - `i_reacts`: Indices of `ms.reaction` active in the reactor. Defaults to every
    reaction in the material source.
  - `flowtype`: Optional holdup constraint. Use `"const. mass"` for constant mass,
    `"const. volume"` for constant volume, or `""` for no additional constraint.

# Examples

```jldoctest
julia> source = MaterialSource("helium"; Mw = 0.004,
           molar_density = (p, T, xᵢ; kwargs...) -> p / (8.314 * T),
           VT_enthalpy = (ϱ, T, xᵢ) -> 20.0 * T,
           VT_internal_energy = (ϱ, T, xᵢ) -> 12.0 * T);

julia> reactor = CSTR(source; name = :reactor, flowtype = "const. volume");

julia> nameof(reactor)
:reactor
```
"""
@component function CSTR(ms::MaterialSource; name, i_reacts = 1:length(ms.reaction), flowtype = "")
    # Subsystems
    @named cv = TPControlVolume(
        ms; N_mcons = 2, N_heats = 1, N_works = 1, phases = ["liquid"], reactive = true
    )

    # Variables
    vars = @variables begin
        Q(t), [description = "heat flux"] #, unit=u"J s^-1"]
    end

    # Parameters
    pars = @parameters begin
        W = 0.0, [description = "work"]      #, unit=u"J s^1"]
    end

    # Equations
    eqs = Equation[
        cv.w1.W ~ W,
        cv.q1.Q ~ Q,
        [
            cv.c2.xᵢ[i] ~ cv.xᵢ[1, i]
                for i in 1:(ms.N_c - 1)
        ]...,
        cv.c2.p ~ cv.p,
        cv.c2.T ~ cv.T,     # Change of moles by reaction
        [
            cv.ΔnR[i] ~ reac.r(
                    cv.p,
                    cv.T,
                    collect(
                        cv.xᵢ[
                            1, :,
                        ]
                    )
                ) * reac.ν[i] * cv.n
                for i in 1:ms.N_c,
                reac in ms.reaction[i_reacts]
        ]...,     # Enthalpy of reaction
        [
            cv.ΔHᵣ ~ reac.r(
                    cv.p,
                    cv.T,
                    collect(
                        cv.xᵢ[
                            1, :,
                        ]
                    )
                ) * reac.Δhᵣ(cv.T) * cv.n
                for reac in ms.reaction[i_reacts]
        ]...,
    ]
    if flowtype == "const. mass"
        push!(
            eqs,
            0.0 ~
                cv.c1.n * sum(ms.Mw[i] * cv.c1.xᵢ[i] for i in 1:ms.N_c) +
                cv.c2.n * sum(ms.Mw[i] * cv.c2.xᵢ[i] for i in 1:ms.N_c)
        )
    elseif flowtype == "const. volume"
        push!(eqs, D(cv.V) ~ 0.0)
    end

    return ODESystem(eqs, t, vars, pars; name, systems = [cv])
end

@public CSTR
