using ModelingToolkit
using Clapeyron
using DynamicQuantities
using LinearAlgebra, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize, equations, get_unknowns
using NonlinearSolve

#Building media
#model = Clapeyron.NRTL(["water", "methanol", "ethanol"], puremodel = PR(["water", "methanol", "ethanol"], idealmodel = ReidIdeal))
#model = Clapeyron.PR(["water", "methanol"])

model = Clapeyron.PR(["water", "methanol"], idealmodel = ReidIdeal)
bubble_pressure(model, 298.15, [0.5, 0.5])
dew_pressure(model, 298.15, [0.55, 0.5])
#flash_cl = Clapeyron.tp_flash(model, 1.01325e5, 350.15, [0.33, 0.33, 0.34])
#enthalpy(model, 1.01325e5, 350.15, [0.8, 0.2], phase = "liquid")

h(p) = TP_flash(model, p, 298.15, [0.5, 0.5])
sol = TP_flash(model, 101325.0, 350.15, [10.0, 5])[2]
rho = molar_density(model, 101325.0, 350.15, sol[:, 3])
ρT_enthalpy(model, rho, 350.15, sol[:, 3])
enthalpy(model, 101325.0, 350.15, sol[:, 3])

FiniteDifferences.central_fdm(3, 1)(h, 150_000.0)

guess = EosBasedGuesses(model, 1.01325e5, 298.15, [0.8, 0.2])
medium = EoSBased(BasicFluidConstants([0.01801528, 0.1801528]), model, guess)
medium.Guesses


#= @named InPort = PhZConnector_(medium = medium)
@named OutPort = PhZConnector_(medium = medium)
@named ControlVolumeState = ρTz_ThermodynamicState_(medium = medium) =#

### ------ Reservoir test
@named stream = FixedBoundary_pTzn_(medium = medium, p = 1.01325e5, T = 350.15, z = [.8, .2], ṅ = 10.0)

keys(guesses(stream)) |> collect

simple_stream = structural_simplify(stream)
#= u0 = [simple_stream.ControlVolumeState.z[1, 1] => 0.8, simple_stream.ControlVolumeState.z[2, 1] => 0.2
,simple_stream.ControlVolumeState.z[1, 3] => 0.5, simple_stream.ControlVolumeState.z[2, 3] => 0.5] =#

prob = SteadyStateProblem(simple_stream, [])
@time sol = solve(prob, SSRootfind())
sol[stream.ControlVolumeState.ϕ]

### ------ ControlVolume


@component function HeatedTank_(;medium, Q̇, pressure, ṅ_out, name)

        @named CV = TwoPortControlVolume_(medium = medium)
        @unpack OutPort, rₐ, rᵥ, Q, p, Wₛ  = CV

        eqs = [
            Wₛ ~ 0.0
            scalarize(rᵥ[:, 2:end] .~ 0.0)...
            scalarize(rₐ[:, 2:end] .~ 0.0)...
            Q ~ Q̇
            p ~ pressure
            OutPort.ṅ[1] ~ -ṅ_out
            OutPort.ṅ[2] ~ -ṅ_out
            OutPort.ṅ[3] ~ -1e-8

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

u0 = [sistem.tank.Nᵢ[1] => 150.0, sistem.tank.Nᵢ[2] => 80.0,
 sistem.tank.ControlVolumeState.T => 300.0]

prob = ODEProblem(sistem, u0, (0.0, 10.0));
prob.u0

sol = solve(prob, FBDF(autodiff = false))

plot(sol.t, sol[sistem.tank.ControlVolumeState.z[2, 2]])

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