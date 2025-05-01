using ModelingToolkit
using Clapeyron
using DynamicQuantities
using LinearAlgebra, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize, equations, get_unknowns
using NonlinearSolve

#Building media
ideal = JobackIdeal(["water", "methanol"])
ideal.params.reference_state

model = Clapeyron.PR(["water", "methanol"], idealmodel = ReidIdeal)
bubble_pressure(model, 350.15, [0.8, 0.2])
dew_pressure(model, 350.15, [0.8, 0.2])
#flash_cl = Clapeyron.tp_flash(model, 1.01325e5, 350.15, [0.33, 0.33, 0.34])
#enthalpy(model, 1.01325e5, 350.15, [0.8, 0.2], phase = "liquid")

guess = EosBasedGuesses(model, 1.01325e5, 350.15, [0.8, 0.2])
medium = EoSBased(BasicFluidConstants([0.01801528, 0.1801528]), model, guess)
medium.Guesses.ρ

### ------ Reservoir test
@named stream = FixedBoundary_pTzn_(medium = medium, p = 1.01325e5, T = 350.15, z = [.8, .2], ṅ = 10.0)

simple_stream = structural_simplify(stream)

prob = SteadyStateProblem(simple_stream, [])
@time sol = solve(prob, SSRootfind())
sol[stream.OutPort.ṅ]

### ------ ControlVolume

@component function HeatedTank_(;medium, Q̇, pressure, ṅ_out, name)

        @named CV = TwoPortControlVolume_(medium = medium)
        @unpack ControlVolumeState, OutPort, rₐ, rᵥ, Q, p, Wₛ, nᴸⱽ  = CV

        eqs = [
            Wₛ ~ 0.0
            scalarize(rᵥ[:, 2:end] .~ 0.0)...
            scalarize(rₐ[:, 2:end] .~ 0.0)...
            Q ~ Q̇
            p ~ pressure
            OutPort.ṅ[1] ~ -ṅ_out
            OutPort.ṅ[2] ~ -ṅ_out
            OutPort.ṅ[3] ~ -1e-8
            ControlVolumeState.p ~ p
            scalarize(ControlVolumeState.z[:, 3] ~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...
            scalarize(nᴸⱽ[1]/sum(nᴸⱽ) ~ flash_vaporized_fraction(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1]))[1])

        ]

        pars = []

        vars = []


    return extend(ODESystem(eqs, t, vars, pars; name), CV)
end

@named tank = HeatedTank_(medium = medium, Q̇ = 100.0, pressure = 1.01325e5, ṅ_out = 10.0)

@component function PerfectFlowHeatedTank_(; medium, ṅ_in, ṅ_out, p, Q, name)

    systems = @named begin
        pT_Boundary = FixedBoundary_pTzn_(medium = medium, p = p, T = 300.15, z = [0.8, 0.2], ṅ = ṅ_in)
        tank = HeatedTank_(medium = medium, Q̇ = Q, pressure = p, ṅ_out = ṅ_out)
    end

    vars = []

    pars = []

    connections = [
                connect(pT_Boundary.OutPort, tank.InPort)
                ]
    
    return ODESystem(connections, t, vars, pars; name = name, systems = [systems...])

end


@named WaterTank = PerfectFlowHeatedTank_(medium = medium, ṅ_in = 10.0, ṅ_out = 10.0, p = 1.01325e5, Q = 100.0)

sistem = structural_simplify(WaterTank)

length(alg_equations(sistem))

equations(sistem)
defaults(sistem)
guesses_ = [
 sistem.tank.ControlVolumeState.z[1, 2] => 0.2, 
 sistem.tank.ControlVolumeState.z[2, 2] => 0.2, 
 sistem.tank.V[2] => 50.0,
  sistem.tank.V[3] => 100.0,
  sistem.tank.nᴸⱽ[2] => 50.0]

u0 = [sistem.tank.Nᵢ[1] => 80.0, sistem.tank.Nᵢ[2] => 80.0,
 sistem.tank.ControlVolumeState.T => 300.0]

prob = ODEProblem(sistem, u0, (0.0, 100.0), guesses = guesses_);
ssprob = SteadyStateProblem(sistem, [guesses_...; u0])

sol = solve(prob, Rodas42(autodiff = false))

sol_ss = solve(ssprob, SSRootfind())

plot(sol.t, sol[sistem.tank.ControlVolumeState.T])

initialization_equations(sistem)
equations(prob.f.initializationprob.f.sys)



















D(U) ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1])

scalarize(D(nᵢ[:, 1]) .~ InPort.ṅ[1].*InPort.z₁ + (OutPort.ṅ[2].*OutPort.z₂ + OutPort.ṅ[3].*OutPort.z₃) .+ collect(rᵥ[:, 2:end]*V[2:end]))

scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))

scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))

scalarize(nᵢ[:, 1] .~ sum(collect(nᵢ[:, 2:end]), dims = 2))

scalarize(sum(collect(ControlVolumeState.ϕ)) ~ 1.0)

scalarize(ControlVolumeState.z .~ nᵢ ./ sum(collect(nᵢ), dims = 1))

ControlVolumeState.p ~ p

[ControlVolumeState.ϕ[j - 1] ~ sum(collect(nᵢ[:, j]), dims = 1)./sum(collect(nᵢ[:, 1]), dims = 1) for j in 2:medium.FluidConstants.nphases]

U ~ (OutPort.h[1] - ControlVolumeState.p/ControlVolumeState.ρ[1])*sum(collect(nᵢ[:, 1])) 

V[1] ~ sum(collect(V[2:end]))
        
V[2]*ControlVolumeState.ρ[2] ~ sum(collect(nᵢ[:, 2]))

V[3]*ControlVolumeState.ρ[3] ~ sum(collect(nᵢ[:, 3]))
        
        
# Outlet port properties

[OutPort.h[j] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.FluidConstants.nphases]
        
OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))
         
OutPort.p ~ ControlVolumeState.p

scalarize(OutPort.z₁ .~ ControlVolumeState.z[:, 1])
scalarize(OutPort.z₂ .~ ControlVolumeState.z[:, 2])
scalarize(OutPort.z₃ .~ ControlVolumeState.z[:, 3])