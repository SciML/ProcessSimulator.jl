using ProcessSimulator
using Clapeyron
using ModelingToolkit
using OrdinaryDiffEq
using NonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test


#Clapeyron model definition
components = ["ethylene oxide", "water", "ethylene glycol"]
idealmodel = CompositeModel(components, 
    liquid = RackettLiquid, 
    gas = ReidIdeal(components, reference_state = :formation), 
    saturation = DIPPR101Sat, 
    hvap = DIPPR106HVap)
activity = NRTL(components)

model = CompositeModel(components, fluid = idealmodel, liquid = activity)

bubble_pressure(model, 350.15, [0.4, 0.6, 0.0]) #Just to test if model is working
#Standard medium building (Will be overwritten inside each equipment with state or guesses)
medium = EoSBased(components = components, eosmodel = model)


@named S1 = Boundary_pTzn(medium = medium, p = 5*101325.0, T = 350.15, z = [0.4, 0.6, 0.0], flowrate = 100.0, flowbasis = :molar)

rxn1 = PowerLawReaction(components = components,
    stoichiometry = Dict("ethylene oxide" => -1.0, "water" => -1.0, "ethylene glycol" => 1.0),
    order = Dict("ethylene oxide" => 1.0, "water" => 0.0),
    A = 0.5, Eₐ = 0.0
)

reactor_state = pTNVState(5*101325.0, 350.15, ones(3)/3.0, base = :Pressure) 

@named R1 = FixedVolumeSteadyStateCSTR(medium = medium, reactionset = rxn1,
                        limiting_reactant = "ethylene oxide",
                         state = reactor_state, volume = 1.0, W = 0.0, Q = nothing)



@named sink = ConnHouse(medium = medium)

connection_set = [
    connect(S1.odesystem.OutPort, R1.odesystem.InPort),
    connect(R1.odesystem.OutPort, sink.InPort)
]

@named sys = System(connection_set, t, [], []; systems = [S1.odesystem, R1.odesystem, sink])

AdiabaticVolumeReactor = mtkcompile(sys, use_scc = true)

default_guesses = guesses(AdiabaticVolumeReactor)
guesses_Reactor = Dict(
    R1.odesystem.V[3] => 1e-8,
    R1.odesystem.OutPort.ṅ[1] => 100.0,
) #These guesses are problem specific

merged_guesses = merge(default_guesses, guesses_Reactor)

prob = NonlinearProblem(AdiabaticVolumeReactor, merged_guesses)

@time sol = solve(prob, abstol = 1e-8, reltol = 1e-8)

unit_ops = Dict(:R1 => :CSTR, :S1 => :Feed)

print_flowsheet_summary(sol, AdiabaticVolumeReactor, unit_ops, components)
print_flowsheet_summary(sol, AdiabaticVolumeReactor, components, R1, S1)