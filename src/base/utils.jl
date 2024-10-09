function connect_states(state_1,state_2; is_stream=false)
    ndir = is_stream ? -1 : 1
    con_eqs = [
        state_1.T ~ state_2.T
        state_1.p ~ state_2.p
        # state_1.ϱ ~ state_2.ϱ
        scalarize(state_1.nᵢ .~ ndir * state_2.nᵢ)...
    ]
    return con_eqs
end