using ProcessSimulator
using Clapeyron
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test



#= @testset "isenthalpic valve" begin =#
#Building media
components = ["carbon dioxide", "methane"]

model = cPR(components, idealmodel = ReidIdeal)

p__ = 50.0*101325.0 # Pa
T__ = 273.15 + 25.0 # K
z__ = [0.5, 0.5] # Mole fractions

medium = EoSBased(components = components, eosmodel = model)

@named S1 = FixedBoundary_pTz_(medium = medium, p = p__, T = T__, z = z__)

@named V1 = Valve(medium = medium, 
state = pTNVState(0.5*p__, T__, z__),
Cv = 4.4e-6,
f = x -> x/sqrt(abs(x) + 1e-8),
flowrate_guess = 1e-3,
flowbasis = :volume)


@named P1 = ConstantPressure(medium = medium, p = 0.5*p__)

connection_set = [
    connect(S1.OutPort, V1.odesystem.InPort),
    connect(V1.odesystem.OutPort, P1.Port)
]

@named sys = System(connection_set, t, [], []; systems = [S1, V1.odesystem, P1])


simple_sys = mtkcompile(sys)

u0 = [V1.odesystem.opening => 0.5]
valve_guesses = [V1.odesystem.InPort.ṅ[1] => V1.molar_flowrate_guess]
prob = ODEProblem(simple_sys, u0, (0.0, 0.001), guesses = valve_guesses, use_scc = false)
sol = solve(prob, FBDF(autodiff = false), abstol = 1e-6, reltol = 1e-6)


T_valve_out = first(sol[V1.odesystem.ControlVolumeState.T]) 

reference_T_valve_out = 290.39815

@test abs(reference_T_valve_out - T_valve_out)/reference_T_valve_out < 0.05

#= end =#