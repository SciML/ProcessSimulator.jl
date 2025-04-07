@component function CSTR(ms::FluidMedium, Reactions::Union{AbstractReaction, Vector{AbstractReaction}}; name, N_heats::Int, N_works::Int, flowtype = "const. volume")
    
    # Subsystems - I think that CSTR would be an extension of the TPControlVolume instead of a composition

    N_Reactions = length(Reactions)

    @named cv = TwoPortControlVolume(ms;
                    N_heats = N_heats, 
                    N_works = N_works, 
                    phases = ["liquid"],
                    N_SinkSource = N_Reactions)

    # Equations
    eqs = Equation[
        ## Connector equations
        [c2.xᵢ[i] ~ xᵢ[1, i] for i in 1:ms.N_c - 1]...,
        c2.p ~ p,
        c2.T ~ T,

        # Change of moles by each reaction
        [ifelse(isa(Reactions, Vector), [rᵥ[j, :] .~ Rate(Reactions[i], V*xᵢ[1, :]..., T, V) for j in 1:N_reaction]..., rᵥ[1, :] .~ Rate(Reactions, V*xᵢ[1, :]..., T, V))]...,


    ]

    if flowtype == "const. mass"

        push!(eqs, 0.0 ~ cv.c1.n * sum(ms.Mw[i]*cv.c1.xᵢ[i] for i in 1:ms.N_c) + cv.c2.n * sum(ms.Mw[i]*cv.c2.xᵢ[i] for i in 1:ms.N_c))

    elseif flowtype == "const. volume"

        push!(eqs, D(V) ~ 0.0)

    end

    return extend(ODESystem(eqs, t, vars, pars; name), cv)
end

