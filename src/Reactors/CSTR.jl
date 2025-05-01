@component function AdiabaticCSTR(;medium, reactions, P, W = 0.0, Q_rate = 0.0, n_out, flowbasis = :molar, phase = "liquid", name)

    @named CV = TwoPortControlVolume_(medium = medium)
    @unpack Nᵢ, V, InPort, OutPort, ControlVolumeState, rₐ, rᵥ, Q, p, Wₛ = CV

    vars = @variables begin
        cᵢ(t)[1:medium.FluidConstants.Nc],    [description = "bulk concentrations"]   #, unit=u"mol m^-3"]
        X(t)[1:medium.FluidConstants.Nc],     [description = "conversion"]            #, unit=u"mol mol^-1"]
    end

    #reactions = ifelse(length(reactions) == 1, [reactions], reactions)
    
    eqs = [
        scalarize(cᵢ .~ Nᵢ/V[1])...
        Wₛ ~ W
        scalarize(rₐ[:, 2:end] .~ 0.0)...
        Q ~ Q_rate
        ControlVolumeState.p ~ p
        p ~ P
        OutPort.ṅ[1] ~ OutPort.ṅ[2] + OutPort.ṅ[3]
    ]

    if phase == "liquid"
        eq_reaction = [
            scalarize(ControlVolumeState.z[:, 3] .~ ones(medium.FluidConstants.Nc)/medium.FluidConstants.Nc)...
            ControlVolumeState.ϕ[1] ~ 1.0
            scalarize(rᵥ[:, 2] .~ sum(Rate.(reactions, cᵢ, ControlVolumeState.T)))...
            scalarize(rᵥ[:, 3] .~ 0.0)...
            OutPort.ṅ[2] ~ -XToMolar(n_out)
            #D(V[2]) ~ 0.0
            OutPort.ṅ[3] ~ -0.0
            scalarize(X .~ (InPort.ṅ[2].*InPort.z₂ .+ OutPort.ṅ[2].*OutPort.z₂)./(InPort.ṅ[2].*InPort.z₂ .+ 1e-8))...
        ]
    elseif phase == "vapor"
        eq_reaction = [
            ControlVolumeState.z[:, 2] ~ ones(medium.FluidConstants.Nc)/medium.FluidConstants.Nc
            ControlVolumeState.ϕ[2] ~ 1.0
            scalarize(rᵥ[:, 2] .~ 0.0)...
            scalarize(rᵥ[:, 3] .~ sum(Rate.(reactions, cᵢ, ControlVolumeState.T)))...
            OutPort.ṅ[2] ~ -0.0
            OutPort.ṅ[3] ~ -XToMolar(n_out)
            #V[3] ~ v
            scalarize(X .~ (InPort.ṅ[3].*InPort.z₃ .+ OutPort.ṅ[3].*OutPort.z₃)./(InPort.ṅ[3].*InPort.z₃ .+ 1e-8))...
        ]
    end

    pars = []

    return extend(ODESystem([eqs...;eq_reaction...], t, collect(Iterators.flatten(vars)), pars; name), CV)

end



reaction1 = PowerLawReaction(["water", "methanol"], [-1.0, 0.0], [2.0, 0.0], 1e-10, 10_000.0)

@component function Reactor(;medium, reactions, P, n_in, n_out, W = 0.0, Q_rate = 0.0, phase = "liquid", name)

    systems = @named begin
        inlet_stream = FixedBoundary_pTzn_(medium = medium, p = P, T = 300.15, z = [0.8, 0.2], ṅ = n_in)
        tank = CSTR(medium = medium, reactions = reactions, P = P, n_out = n_out, W = W, Q_rate = Q_rate, phase = phase)
    end

    vars = []

    pars = []

    connections = [
        connect(inlet_stream.OutPort, tank.InPort)
        ]

    return ODESystem(connections, t, vars, pars; name, systems = [systems...])

end

@named R_101 = Reactor(medium = medium, reactions = [reaction1], P = 101325.0, n_in = 10.0, n_out = 10.0, W = 0.0, Q_rate = 0.0, phase = "liquid")

sistem = structural_simplify(R_101)

#unknowns_ = unknowns(sistem)

u0 = [sistem.tank.Nᵢ[1] => 60.0, sistem.tank.Nᵢ[2] => 60.0,
 sistem.tank.ControlVolumeState.T => 300.0]

guesses_ = [
sistem.tank.V[2] => 0.2,
sistem.tank.V[3] => 0.0001,
sistem.tank.nᴸⱽ[2] => 0.0]

prob = ODEProblem(sistem, u0, (0.0, 100.0), guesses = guesses_);

sol = solve(prob, FBDF(autodiff = false))

plot(sol.t, sol[sistem.tank.X[1]])


