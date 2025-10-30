using ModelingToolkit
using Clapeyron
using LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize, equations, get_unknowns
using NonlinearSolve
using OrdinaryDiffEq

#Building media
components = ["carbon dioxide", "methane"]

model = SRK(components, idealmodel = ReidIdeal)

#model = ReidIdeal(components)

p__ = 1.0*101325.0 # Pa
T__ = 273.15 + 25.0 # K
z__ = [0.5, 0.5] # Mole fractions

masstransfermodel = ConstantMassTransferCoeff([5e-1, 3e-1]); #This assumes the same correlation for all interfaces you might have with the control volume.
heat_transfer_model = ConstantHeatTransferCoeff(10.0);
viscosity_model = ChapmanEnskogModel(components, ref = "Poling et al. (2001)")
fluid_transport_model = TransportModel(masstransfermodel, heat_transfer_model, viscosity_model)
medium = EoSBased(components = components, eosmodel = model, transportmodel = fluid_transport_model, state = pTzState(p__, T__, z__))


### ------ Reservoir test

@named reservoir = FixedBoundary_pTzn_(medium = medium, p = p__, T = T__, z = z__, flowrate = -5e-4, flowbasis = :volume)
@named sink = ConnHouse(medium = medium)

connections = [connect(reservoir.OutPort, sink.port)]

@named sys = System(connections, t; systems = [reservoir, sink])

prob = NonlinearProblem(simple_stream, guesses(simple_stream))
@time sol = solve(prob, RobustMultiNewton())
sol[reservoir.ControlVolumeState.ϕ]




#Testing adsorption interface
## Solid Eos
adsorbenteos = SolidEoSModel(750.0, 273.15, [0.0, 0.0, 0.0, 0.0], 0.0, 273.15, [935.0, 0.0, 0.0, 0.0]) #J/kg/K

# Isotherm
iso_c1 = LangmuirS1(1e-8, 1e-8, -30_000.1456)
iso_c2 = LangmuirS1(10.0, 1e-9, -30_000.1456)
all_iso = ExtendedLangmuir(iso_c1, iso_c2)


# Mass transfer in film
masstransfermodel = HomogeneousDiffusivityCoeff([0.0, 1e-8])
heattransfermodel = ConstantHeatTransferCoeff(0.0)
transportmodel = TransportModel(masstransfermodel, heattransfermodel)


#Defining adsorbent model
adsorbent = Adsorbent(adsorbent_name = "XYZ",
                     particle_size = 2e-3,
                     components = components,
                     isotherm = all_iso,
                     EoSModel = adsorbenteos,
                     transport_model = transportmodel) 


#Set constants
porosity = 0.5
V = 5e-3*10.0
#p = 101325.0
phase = "vapor"
solidmedium = adsorbent
fluidmedium = medium
mass_of_adsorbent = V * solidmedium.EoSModel.ρ_T0 * porosity
A  = V*(1.0 - porosity)*area_per_volume(solidmedium) #Interfacial area


_0 = 1e-8
Ntot = 5.0
guess_adsorber = EosBasedGuesses(model, Clapeyron.pressure(model, V, T__, [Ntot, _0]), T__, [Ntot, _0]/(Ntot + _0))
medium_adsorber = EoSBased(constants, model, fluid_transport_model, guess_adsorber)
#Well mixed adsorber test in vapor phase with valve
@named boundary = FixedBoundary_pTzn_(medium = fluidmedium, p = p__, T = T__, z = z__, flowrate = -0.10, flowbasis = :molar)
@named tank =  WellMixedAdsorber(fluidmedium = medium_adsorber, solidmedium = solidmedium, porosity = 0.5, p = p__, V = V, phase = "vapor")
@named valve = Valve_(medium = medium_adsorber, Cv = 8.4e-6, ΔP_f = x -> √(x) , phase = "vapor")
@named sink = ConstantPressure(medium = medium_adsorber, p = 0.2*p__)
@named flowsink = ConstantFlowRate(medium = fluidmedium, flowrate = 0.3, flowbasis = :molar)


topography = [connect(boundary.OutPort, tank.mobilephase.InPort),
              connect(tank.mobilephase.OutPort, valve.InPort),
              connect(valve.OutPort, sink.port)
]

#= perfect_flow = [connect(boundary.OutPort, tank.mobilephase.InPort),
                connect(tank.mobilephase.OutPort, flowsink.port)] =#


@named flowsheet = System(topography, t, [], [], systems = [boundary, tank, valve, sink])

#@named perfect_flow_system = System(perfect_flow, t, [], [], systems = [tank, boundary, flowsink])

ModelingToolkit.flatten_equations(equations(expand_connections(flowsheet)))

simple_flowsheet = mtkcompile(flowsheet)

alg_equations(simple_flowsheet)
unknowns(simple_flowsheet)

see_g = guesses(simple_flowsheet)


u0_flowsheet = [tank.mobilephase.Nᵢ[1] => Ntot,
                tank.mobilephase.Nᵢ[2] => _0,
                tank.mobilephase.ControlVolumeState.T => T__,
                tank.stationaryphase.Nᵢ[1] => Ntot,
                tank.stationaryphase.Nᵢ[2] => _0,
                tank.stationaryphase.ControlVolumeState.T => T__,
                valve.opening => 0.5]


missing_guesses = [
                   tank.mobilephase.ControlVolumeState.p => guess_adsorber.p,
                   tank.mobilephase.V[2] => _0,
                   valve.InPort.ṅ[1] => 0.1
]


prob = ODEProblem(simple_flowsheet, u0_flowsheet, (0.0, 15500.0), guesses = missing_guesses, use_scc = false)

@time sol = solve(prob, FBDF(autodiff = false), abstol = 1e-10, reltol = 1e-10)




sol[tank.mobilephase.Nᵢ]
sol[tank.mobilephase.InPort.ṅ[1]]
sol[tank.stationaryphase.Nᵢ[2]]
plot(sol.t, sol[valve.InPort.ṅ[1]], label = "total molar flow rate")
sol[tank.interface.SolidSurface.InPort.ϕₘ[1]]
sol[tank.mobilephase.ControlVolumeState.p]
sol[tank.mobilephase.V[1]]


sol[tank.interface.SolidSurface.InPort.ϕₘ[2]] .+ sol[tank.interface.FluidSurface.OutPort.ϕₘ[2]]


plot(sol.t, sol[tank.interface.FluidSurface.InPort.ϕₕ])
plot(sol.t, sol[tank.interface.SolidSurface.InPort.ϕₘ[1]])
plot(sol.t, sol[tank.interface.FluidSurface.InPort.ϕₘ[2]])
plot(sol.t, sol[tank.stationaryphase.Q], label = "methane")
plot(sol.t, sol[tank.stationaryphase.Nᵢ[2]], label = "methane")
plot(sol.t, sol[tank.mobilephase.U], label = "mobile phase")
plot(sol.t, sol[tank.mobilephase.ControlVolumeState.T], label = "mobile phase")
plot(sol.t, sol[tank.mobilephase.Nᵢ[1]])
plot(sol.t, sol[tank.mobilephase.ControlVolumeState.z[1, 1]], label = "mobile phase")

Clapeyron.pressure(model, V, 300.15, 0.22*[0.0, 0.0, 1.0])
Clapeyron.volume(model, 101325.0, 300.15, 0.22*[0.0, 0.0, 1.0], phase = :stable)















### ---- 
#= @named InPort = PhZConnector_(medium = medium)
@named OutPort = PhZConnector_(medium = medium)
@named ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
vars = @variables begin
    rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "mass source or sink - volumetric  basis"]
    rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                                      [description = "molar source or sink - through surface in contact with other phases"]      
    Nᵢ(t)[1:medium.Constants.Nc],                                                                  [description = "molar holdup"]
    nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                        [description = "molar holdup in each phase, excluding the overall phase"]              
    U(t),                                                                                          [description = "internal energy holdup"]                                               
    p(t),                                                                                          [description = "pressure"]  
    V(t)[1:medium.Constants.nphases],                                                              [description = "volume"]    
    A(t),                                                                                          [description = "control volume area"]
    Q(t),                                                                                          [description = "heat flux"]
    Wₛ(t),                                                                                         [description = "shaft work"]
end

D(U) ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) + Q + Wₛ
[D(Nᵢ[i]) ~ InPort.ṅ[1]*InPort.z[i, 1] + sum(dot(collect(OutPort.ṅ[2:end]), collect(OutPort.z[i:i, 2:end]))) + sum(collect(rᵥ[i, 2:end].*V[2:end])) + rₐ[i, 1]*A for i in 1:medium.Constants.Nc]


scalarize(sum(OutPort.ṅ[2:end].*OutPort.z[:, 2:end], dims = 2))

sum(dot(OutPort.ṅ[2:end], OutPort.z[i:i, 2:end]))

scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2)) =#

