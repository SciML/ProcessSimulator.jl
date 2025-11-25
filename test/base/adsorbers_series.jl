using ModelingToolkit
using Clapeyron
using LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: scalarize, equations, get_unknowns
using NonlinearSolve
using OrdinaryDiffEq
using EntropyScaling

#Building media
components = ["carbon dioxide", "methane"]

#model = cPR(components, idealmodel = ReidIdeal, translation = RackettTranslation(components))

model = ReidIdeal(components)

p__ = 2.0*101325.0 # Pa
T__ = 273.15 + 35.0 # K
z__ = [0.2, 0.8] # Mole fractions

#Fluid medium
guess = EosBasedGuesses(model, p__/2.0, T__, z__)
constants = BasicFluidConstants(components)
masstransfermodel = ConstantMassTransferCoeff([5e-1, 3e-1]);
heat_transfer_model = ConstantHeatTransferCoeff(10.0);
viscosity_model = ChapmanEnskogModel(components, ref = "Poling et al. (2001)")
fluid_transport_model = TransportModel(masstransfermodel, heat_transfer_model, viscosity_model)
medium = EoSBased(constants, model, fluid_transport_model, guess)


#Solid medium
## Solid Eos
adsorbenteos = SolidEoSModel(750.0, 273.15, [0.0, 0.0, 0.0, 0.0], 0.0, 273.15, [935.0, 0.0, 0.0, 0.0]) #J/kg/K

# Isotherm
iso_c1 = LangmuirS1(1e-8, 1e-8, -30_000.1456)
iso_c2 = LangmuirS1(10.0, 1e-9, -30_000.1456)
all_iso = ExtendedLangmuir(iso_c1, iso_c2)


# Mass transfer in film
masstransfermodel = HomogeneousDiffusivityCoeff([0.0, 1e-8])
heattransfermodel = ConstantHeatTransferCoeff(10.0)
transportmodel = TransportModel(masstransfermodel, heattransfermodel)


#Defining adsorbent model
adsorbent = Adsorbent(adsorbent_name = "XYZ",
                     particle_size = 2e-3,
                     components = components,
                     isotherm = all_iso,
                     EoSModel = adsorbenteos,
                     transport_model = transportmodel) 


#Tank info
phase = "vapor"
tank = CylindricalTank(D = 0.2, L = 2.0, porosity = 0.4)
discretize!(tank, 4)
V = volume_(tank)


#Two well mixed adsorbers connected by ergun pressure drop
@named reservoir = FixedBoundary_pTz_(medium = medium, p = p__, T = T__, z = z__)
@named adsorber1 = WellMixedAdsorber(fluidmedium = medium, solidmedium = adsorbent, porosity = tank.porosity, p = p__, V = V, phase = phase)
@named pdrop2 = ErgunDrop(medium = medium, solidmedium = adsorbent, tank = tank, phase = phase)
@named adsorber2 = WellMixedAdsorber(fluidmedium = medium, solidmedium = adsorbent, porosity = tank.porosity, p = p__, V = V, phase = phase)
@named outletvalve = Valve_(medium = medium, Cv = 8.4e-6, ΔP_f = x -> √(x) , phase = phase)
@named pressure_sink = ConstantPressure(medium = medium, p = p__/2.0)
