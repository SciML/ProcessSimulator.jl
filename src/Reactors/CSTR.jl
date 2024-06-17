@component function CSTR(; substances_user = ["methanol", "propylene oxide", "water"], 
    Nc = length(substances_user), 
    phase = :liquid, 
    model = PR(substances_user),
    properties = Dict(subs => load_component_properties(subs) for subs in substances_user),
    Reaction = KineticReactionNetwork(Af_r = [1.0], Ef_r = [1.0], Do_r = [1.0 1.0], name = :DefaultReaction),
    ninports = 1,  
    InPorts = [matcon(; Nc = Nc, name = Symbol("InPorts$i")) for i in 1:ninports],
    Ac = 1.0, 
    height_out_port = 0.0
   )
    
    #Numerical variables
    
    #Constants
    GravitationalConst = 9.81 # m²/s
    gramsToKilograms = 10^(-3)
    Rᵍ = 8.314 # J/(mol K)
    
    #Properties of individual substances
    MWs = [properties[subs]["MW"] for subs in substances_user]
    
    #Connection and reaction constants
    Nri = defaults(Reaction)[@nonamespace Reaction.Nr]

    pars = @parameters begin 
        ## How to inherent parameters from Reaction?
        height_out = height_out_port, [description = "Height of the outlet stream port with reference from the bottom of the tank (m)"]
        N_InPorts = ninports
        Nr = Nri, [description = "Number of reactions"]
        Af_r[1:Nri] = defaults(Reaction)[@nonamespace Reaction.Af_r], [description = "Arrhenius constant of each reaction at given temperature ()"]
        Coef_Cr[1:Nri, 1:Nc] = defaults(Reaction)[@nonamespace Reaction.Coef_Cr], [description = "Stoichiometric coefficients of each component in each reaction (-)"]
        Do_r[1:Nri, 1:Nc] = defaults(Reaction)[@nonamespace Reaction.Do_r], [description = "Forward order of the components (-)"]
        Ef_r[1:Nri] = defaults(Reaction)[@nonamespace Reaction.Ef_r], [description = "Activation energy of each reaction at given temperature ()"]
        N = Nc, [description = "Number of components"]
        A = Ac, [description = "Cross sectional area of the tank (m²)"]
    end   
    
    
    OutPorts = @named begin
        Out = matcon(; Nc = Nc) 
    end   

    vars = @variables begin

    T(t), [description = "Temperature of vessel contents (K)"]
    M(t), [description = "Mass holdup in the tank (kg)"]
    N(t), [description = "Total molar holdup in the tank (kmol)"]
    V(t), [description = "Volume holdup in the tank (m³)"]
    (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component in the tank (mol)"]
    (Cᵢ(t))[1:Nc], [description = "Concentration of each component in the tank (mol/m³)"]
    (ρ(t)), [description = "Molar Density of the fluid in the tank (mol/m³)"]
    (ρʷ(t)), [description = "Mass Density of the fluid in the tank (kg/m³)"]
    #h_out(t), [description = "Outlet specific enthalpy (J/mol)"]
    MW(t), [description = "Molecular weight of fluid in the tank (kg/kmol)"]    
    P_out(t), [description = "Pressure at the outlet stream level (Pa)"]
    H(t), [description = "Enthalpy holdup of the fluid in the tank (J)"] 
    S(t), [description = "Entropy holdup of the fluid in the tank (J/K)"]
    F_out(t), [description = "Outlet molar flow rate (mol/s)"] 
    Fʷ_out(t), [description = "Outlet molar flow rate (mol/s)"]
    Q_out(t), [description = "Outlet volumetric flow rate (m³/s)"] # DoF
    height(t), [description = "Liquid level in vessel measured from bottom of the tank (m)"]

    (Cᵢ_in(t))[1:Nc, 1:Ni_InPorts], [description = "Inlet concentration of each component (mol/m³)"] # DoF through inlet stream
    (F_in(t))[1:Ni_InPorts], [description = "Inlet molar flow rate (mol/s)"] # DoF through inlet stream
    (Q_in(t))[1:Ni_InPorts], [description = "Inlet volumetric flow rate(s) (m³/s)"]
    (T_in(t))[1:Ni_InPorts], [description = "Inlet temperature (K)"] # DoF through inlet stream
    (h_in(t))[1:Ni_InPorts], [description = "Inlet specific enthalpy (J/mol)"]
    (ρ_in(t))[1:Ni_InPorts], [description = "Inlet density (mol/m³)"]
    (ρʷ_in(t))[1:Ni_InPorts], [description = "Inlet density (mol/m³)"]
    P_atm(t), [description = "Tank pressure (Pa)"] # Pontetial DoF and equal to inlet pressures.
    Q̇(t), [description = "Heat transfer rate (J/s)"] # Potential DoF
    (r(t))[1:Nri, 1:Nc], [description = "Rate of each reaction for each component (mol/s/m³)"]
    (R(t))[1:Nc], [description = "Overall reaction rate (mol/s/m³)"]
    end

    

    #Reaction equations
    reaction_rate = [r[i, j] ~ Af_r[i]*exp(-Ef_r[i]/(Rᵍ*T))*(Cᵢ[:].^(Do_r[j, :])) for i in 1:Nri for j in 1:Nc] # If there's an inert, the order is just zero, but it has to be written
    overall_reaction_rate = [R[i] ~ sum(r[:, i].*Coef_Cr[:, i]) for i in 1:Nc]  # If there's an inert, the coefficient is just zero, but it has to be written

   
    #Inlet connector variables's equations
    atm_pressure = [P_atm ~ InPorts[1].P]
    mass_density_eqs = [ρʷ_in[j] ~ InPorts[j].ρʷ for j in 1:Ni_InPorts]
    molar_density_eqs = [ρ_in[j] ~ InPorts[j].ρ for j in 1:Ni_InPorts]
    inletenthalpy = [h_in[j] ~ InPorts[j].H for j in 1:Ni_InPorts]
    inletconcentrations = [Cᵢ_in[i, j] ~ InPorts[j].z₁[i]*ρ_in[j] for j in 1:Ni_InPorts for i in 1:Nc]
    inlettemperature_eqs = [T_in[j] ~ InPorts[j].T for j in 1:Ni_InPorts]
    inletmolarflow_eqs = [F_in[j] ~ InPorts[j].F for j in 1:Ni_InPorts]
    volumetricflow_eqs = [Q_in[j] ~ F_in[j] / ρ_in[j] for j in 1:Ni_InPorts]

    eqs = [reaction_rate; overall_reaction_rate; atm_pressure; mass_density_eqs; molar_density_eqs; inletenthalpy; inletconcentrations; inlettemperature_eqs; inletmolarflow_eqs; volumetricflow_eqs]

    #Outlet connector equations:
    out_conn = [Out.P ~ P_out
                Out.T ~ T
                Out.F ~ F_out
                Out.Fʷ ~ Fʷ_out
                Out.H ~ H/N
                Out.S ~ S/N
                Out.ρʷ ~ ρʷ
                Out.ρ ~ ρ
                scalarize(Out.z₁ .~ Nᵢ/N)...
                Out.MW[1] ~ MW
    ]

    out_conn_phases = [if phase == :liquid
                            scalarize(Out.z₂ .~ 0.0)...
                            scalarize(Out.z₃ .~ z₁)...
                            Out.MW[2] ~ 0.0
                            Out.MW[3] ~ Out.MW[1]
                            Out.α_g ~ 0.0
                    elseif phase == :vapor
                            scalarize(Out.z₂ .~ z₁)...
                            scalarize(Out.z₃ .~ 0.0)...
                            Out.MW[2] ~ Out.MW[1]
                            Out.MW[3] ~ 0.0
                            Out.α_g ~ 1.0
                    end
            ]



    #balances
    mass_balance = [D(M) ~ sum(Q_in.*ρʷ_in) - Q_out*ρʷ]
    component_balance = [d(Nᵢ[i]) ~ sum(Q_in.[:]*Cᵢ_in[i, :]) - Q_out*Cᵢ[i] + R[i]*V for i in 1:Nc]
    energy_balance = [D(H(t)) ~ Q_in*ρ_in*h_in - Q_out*h_out*ρ + Q̇]
    mass_volume_eq = [ρʷ*V ~ M]
    mol_holdup = [N ~ sum(Nᵢ)]
    mol_to_concentration = [Nᵢ .~ Cᵢ*V]
    height_to_volume = [phase = :vapor ? height ~ 0.0 : height*A ~ V]
    volumetricflow_to_molarflow = [Q_out ~ F_out/ρ]
    volumetricflow_to_massflow = [Q_out ~ Fʷ_out/ρʷ]
    

    #Thermodynamic properties (outlet)
    pressure_out = [phase == :liquid ? P_out ~ P_atm + ρʷ*GravitationalConst*(height - height_out) : P_out ~ P_atm] #Estimation considering static pressure for liquids (May be off as tank is stirred and not static)
    density_eqs = [ρ ~ molar_density(model, P_out, T, Nᵢ; phase = phase), ρʷ ~ mass_density(model, P, T, Nᵢ; phase = phase)]
    globalEnthalpy_eq = [H ~ enthalpy(model, P_out, T, Nᵢ; phase = phase)]
    molar_mass = [MW ~ sum(MWs[i]*Nᵢ[i]/N)*gramsToKilograms]
    entropy_eq = [S ~ entropy(model, P_atm, T, Nᵢ; phase = phase)]

    eqs = [reaction_rate...; overall_reaction_rate...; atm_pressure...; mass_density_eqs...; molar_density_eqs...; inletenthalpy...; inletconcentrations...; inlettemperature_eqs...; inletmolarflow_eqs...; volumetricflow_eqs...; out_conn...;
    out_conn_phases...; mass_balance...; component_balance...; energy_balance...; mass_volume_eq...; mol_holdup...; mol_to_concentration...; height_to_volume...; volumetricflow_to_molarflow...; volumetricflow_to_massflow...;
    pressure_out...; density_eqs...; globalEnthalpy_eq...; molar_mass...; entropy_eq...]

end

#= 
using ModelingToolkit, JSON, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron, GCIdentifier
using NonlinearSolve
using JSON
@parameters t
D = Differential(t)
include("C:/Users/Vinic/OneDrive/Pos-Doc/Flowsheeting/ProcessModeling/ProcessSimulator.jl/src/utils")
include("C:/Users/Vinic/OneDrive/Pos-Doc/Flowsheeting/ProcessModeling/ProcessSimulator.jl/src/Sources/Sourceutils.jl")
substances_user = ["water", "methanol", "propyleneglycol","methyloxirane"]
Nc = size(substances_user, 1)
#substances_user = ["water"]
properties = Dict(subs => load_component_properties(subs) for subs in substances_user)
# Function to extract parameters for ReidIdeal model
read_reidcp(properties, substances_user)


cp_params = (a = [36.54206320678348, 39.19437678197436, 25.7415, 34.91747774761048], b = [-0.03480434051958945, -0.05808483585041852, 0.2355, -0.014935581577635826], c = [0.000116818199785053, 0.0003501220208504329, 0.0001578, 0.000756101594841365], d = [-1.3003819534791665e-7, -3.6941157412454843e-7, -4.0939e-7, -1.0894144551347726e-6], e = [5.2547403746728466e-11, 1.276270011886522e-10, 2.1166e-10, 4.896983427747592e-10])
idealmodel = ReidIdeal(["water", "methanol", "propyleneglycol","methyloxirane"]; userlocations = cp_params)
pcpsaft = PCPSAFT(["water", "methanol", "propyleneglycol","methyloxirane"], idealmodel = idealmodel)
phase = :liquid
model = pcpsaft
ΔT = 10.
Ts = Base._linspace(298.00, 370.0, 20) |> collect
println(Ts)
z = [1.0, 1.0, 1.0, 1.0]
isobaric_heat_capacity(pcpsaft, 101325., 298.00, z, phase = :liquid)
bubble_temperature(pcpsaft, 5*101325., z)
molar_density(pcpsaft, 5*101325, 350.15, z, phase = :liquid)
rhos = [molar_density(pcpsaft, 5*101325, T, z, phase = :liquid) for T in Ts]
enthalpy(pcpsaft, eps(1.), 298.00, z)

isobaric_heat_capacity(IAPWS95(), 101325, 298.0, 1.)

cp_w = (a = [36.54206320678348], b = [-0.03480434051958945], c = [0.000116818199785053], d = [-1.3003819534791665e-7], e = [5.2547403746728466e-11])
ideal_water = ReidIdeal(["water"]; userlocations = cp_w)
pcpwater = PCPSAFT(["water"], idealmodel = ideal_water)
isobaric_heat_capacity(pcpwater, 101325, 298.00, 1., phase = :liquid)
enthalpy(pcpwater, 101325, 298.00, 1., phase = :liquid)


Reaction = KineticReactionNetwork(;substances_user = substances_user, 
Af_r = 4.71e9, Ef_r = 32400*1055.6/453.6, Coef_Cr = [-1.0 0.0 1.0 -1.0], 
Do_r = [1.0 0.0 0.0 1.0], name = "Propyleneglycol synthesis")
#["water", "methanol", "propyleneglycol","methyloxirane"]
ninports = 1
InPorts = [matcon(; Nc = Nc, name = Symbol("InPorts$i")) for i in 1:ninports]

Ac = 1.0
height_out_port = 0.0
Nri =  Reaction.Nri

#Constants
GravitationalConst = 9.81 # m²/s
gramsToKilograms = 10^(-3)
Rᵍ = 8.314 # J/(mol K)

#Properties of individual substances
MWs = [properties[subs]["MW"] for subs in substances_user]
ΔH₀f = [properties[subs]["IGHF"]/10^3 for subs in substances_user] # (IG formation enthalpy) J/mol


pars = @parameters begin 
height_out = height_out_port, [description = "Height of the outlet stream port with reference from the bottom of the tank (m)"]
N_InPorts = ninports
Nr = Nri, [description = "Number of reactions"]
Af_r[1:Nri] = Reaction.Af_r, [description = "Arrhenius constant of each reaction at given temperature ()"]
Coef_Cr[1:Nri, 1:Nc] = Reaction.Coef_Cr, [description = "Stoichiometric coefficients of each component in each reaction (-)"]
Do_r[1:Nri, 1:Nc] = Reaction.Do_r, [description = "Forward order of the components (-) "]
Ef_r[1:Nri] = Reaction.Ef_r, [description = "Activation energy of each reaction at given temperature ()"]
N = Nc, [description = "Number of components"]
A = Ac, [description = "Cross sectional area of the tank (m²)"]
end

@variables v
    
OutPorts = @named begin
    Out = matcon(; Nc = Nc) 
end   

vars = @variables begin

    M(t), [description = "Mass holdup in the tank (kg)"]
    N(t), [description = "Total molar holdup in the tank (kmol)"]
    V(t), [description = "Volume holdup in the tank (m³)"]
    (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component in the tank (mol)"]
    (Cᵢ(t))[1:Nc], [description = "Concentration of each component in the tank (mol/m³)"]
    (ρ(t)), [description = "Molar Density of the fluid in the tank (mol/m³)"]
    (ρʷ(t)), [description = "Mass Density of the fluid in the tank (kg/m³)"]
    #h_out(t), [description = "Outlet specific enthalpy (J/mol)"]
    MW(t), [description = "Molecular weight of fluid in the tank (kg/kmol)"]  
    T(t), [description = "Temperature of vessel contents (K)"]  
    P_out(t), [description = "Pressure at the outlet stream level (Pa)"]
    H(t), [description = "Enthalpy holdup of the fluid in the tank (J)"] 
    S(t), [description = "Entropy holdup of the fluid in the tank (J/K)"]
    F_out(t), [description = "Outlet molar flow rate (mol/s)"] 
    Fʷ_out(t), [description = "Outlet molar flow rate (mol/s)"]
    Q_out(t), [description = "Outlet volumetric flow rate (m³/s)"] # DoF
    height(t), [description = "Liquid level in vessel measured from bottom of the tank (m)"]

    (Cᵢ_in(t))[1:Nc, 1:Ni_InPorts], [description = "Inlet concentration of each component (mol/m³)"] # DoF through inlet stream
    (F_in(t))[1:Ni_InPorts], [description = "Inlet molar flow rate (mol/s)"] # DoF through inlet stream
    (Q_in(t))[1:Ni_InPorts], [description = "Inlet volumetric flow rate(s) (m³/s)"]
    (T_in(t))[1:Ni_InPorts], [description = "Inlet temperature (K)"] # DoF through inlet stream
    (h_in(t))[1:Ni_InPorts], [description = "Inlet specific enthalpy (J/mol)"]
    (ρ_in(t))[1:Ni_InPorts], [description = "Inlet density (mol/m³)"]
    (ρʷ_in(t))[1:Ni_InPorts], [description = "Inlet density (mol/m³)"]
    P_atm(t), [description = "Tank pressure (Pa)"] # Equal to inlet pressures.
    Q̇(t), [description = "Heat transfer rate (J/s)"] # Potential DoF
    (r(t))[1:Nri], [description = "Rate of each reaction for each component (mol/s/m³)"]
    (R(t))[1:Nc], [description = "Overall reaction rate (mol/s/m³)"]
end

    

#Reaction equations
reaction_rate = [r[i] ~ Af_r[i]*exp(-Ef_r[i]/(Rᵍ*T))*prod(scalarize((Cᵢ[:].^Do_r[i, :]))) for i in 1:Nri] # If there's an inert, the order is just zero, but it has to be written
overall_reaction_rate = [R[i] ~ sum(scalarize(r[:].*Coef_Cr[:, i])) for i in 1:Nc]  # If there's an inert, the coefficient is just zero, but it has to be written

   
#Inlet connector variables's equations
atm_pressure = [P_atm ~ InPorts[1].P]
mass_density_eqs = [ρʷ_in[j] ~ InPorts[j].ρʷ for j in 1:Ni_InPorts]
molar_density_eqs = [ρ_in[j] ~ InPorts[j].ρ for j in 1:Ni_InPorts]
inletenthalpy = [h_in[j] ~ InPorts[j].H for j in 1:Ni_InPorts]
inletconcentrations = [Cᵢ_in[i, j] ~ InPorts[j].z₁[i]*ρ_in[j] for j in 1:Ni_InPorts for i in 1:Nc]
inlettemperature_eqs = [T_in[j] ~ InPorts[j].T for j in 1:Ni_InPorts]
inletmolarflow_eqs = [F_in[j] ~ InPorts[j].F for j in 1:Ni_InPorts]
volumetricflow_eqs = [Q_in[j] ~ F_in[j] / ρ_in[j] for j in 1:Ni_InPorts]

eqs = [reaction_rate; overall_reaction_rate; atm_pressure; mass_density_eqs; molar_density_eqs; inletenthalpy; inletconcentrations; inlettemperature_eqs; inletmolarflow_eqs; volumetricflow_eqs]

#Outlet connector equations:
out_conn = [Out.P ~ P_out
            Out.T ~ T
            Out.F ~ F_out
            Out.Fʷ ~ Fʷ_out
            Out.H ~ H/N
            Out.S ~ S/N
            Out.ρʷ ~ ρʷ
            Out.ρ ~ ρ
            scalarize(Out.z₁ .~ Nᵢ/N)...
            Out.MW[1] ~ MW
]

    if phase == :liquid
    out_conn_phases = [
                    scalarize(Out.z₂ .~ 0.0)...
                    scalarize(Out.z₃ .~ Out.z₁)...
                    Out.MW[2] ~ 0.0
                    Out.MW[3] ~ Out.MW[1]
                    Out.α_g ~ 0.0]

    elseif phase == :vapor
        out_conn_phases = [
        scalarize(Out.z₂ .~ Out.z₁)...
        scalarize(Out.z₃ .~ 0.0)...
        Out.MW[2] ~ Out.MW[1]
        Out.MW[3] ~ 0.0
        Out.α_g ~ 1.0]
    end



    #balances
    mass_balance = [D(M) ~ sum(scalarize(Q_in.*ρʷ_in)) - Q_out*ρʷ]
    component_balance = [D(Nᵢ[i]) ~ sum(scalarize(Q_in.*Cᵢ_in[i, :])) - Q_out*Cᵢ[i] + R[i]*V for i in 1:Nc] #Neglectable loss to vapor phase head space
    energy_balance = [D(H) ~ sum(scalarize(Q_in.*ρ_in.*h_in)) - Q_out*H/N*ρ + Q̇]
    mass_volume_eq = [ρʷ*V ~ M]
    mol_holdup = [N ~ sum(scalarize(Nᵢ))]
    mol_to_concentration = [scalarize(Nᵢ .~ Cᵢ*V)...]
    height_to_volume = [height*A ~ V]
    volumetricflow_to_molarflow = [Q_out ~ F_out/ρ]
    volumetricflow_to_massflow = [Q_out ~ Fʷ_out/ρʷ]
  
    #Thermodynamic properties (outlet)
    pressure_out = [phase == :liquid ? P_out ~ P_atm + ρʷ*GravitationalConst*(height - height_out) : P_out ~ P_atm] #Estimation considering static pressure (May be off as tank is agitated and not static)
    density_eqs = [ρ ~ molar_density(model, P_out, T, scalarize(Nᵢ); phase = phase), ρʷ ~ mass_density(model, P_out, T, scalarize(Nᵢ); phase = phase)]
    globalEnthalpy_eq = [H ~ enthalpy(model, P_out, T, Nᵢ; phase = phase)]
    molar_mass = [MW ~ sum(MWs[i]*Nᵢ[i]/N)*gramsToKilograms]
    entropy_eq = [S ~ entropy(model, P_atm, T, Nᵢ; phase = phase)]

    eqs = [reaction_rate...; overall_reaction_rate...; atm_pressure...; mass_density_eqs...; molar_density_eqs...; inletenthalpy...; inletconcentrations...; inlettemperature_eqs...; inletmolarflow_eqs...; volumetricflow_eqs...; out_conn...;
    out_conn_phases...; mass_balance...; component_balance...; energy_balance...; mass_volume_eq...; mol_holdup...; mol_to_concentration...; height_to_volume...; volumetricflow_to_molarflow...; volumetricflow_to_massflow...;
    pressure_out...; density_eqs...; globalEnthalpy_eq...; molar_mass...; entropy_eq...] =#