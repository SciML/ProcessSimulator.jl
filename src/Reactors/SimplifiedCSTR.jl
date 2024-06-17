@component function SimpleCSTR(; substances_user, 
    Nc = length(substances_user), 
    phase, 
    model,
    Reaction,
    ninports,  
    Ac, 
    height_out_port,
    name
   )
    
#Numerical variables
    
#Constants
#GravitationalConst = 9.81 # m²/s
gramsToKilograms = 10^(-3)
Rᵍ = 8.314 # J/(mol K)
Nri = Reaction.Nri

#Properties of individual substances
properties = Dict(subs => load_component_properties(subs) for subs in substances_user)
MWs = [properties[subs]["MW"] for subs in substances_user]
ΔH₀f = [properties[subs]["IGHF"]/10^3 for subs in substances_user] # (IG formation enthalpy) J/mol


pars = @parameters begin 
Af_r[1:Nri] = Reaction.Af_r, [description = "Arrhenius constant of each reaction at given temperature ()"]
Coef_Cr[1:Nri, 1:Nc] = Reaction.Coef_Cr, [description = "Stoichiometric coefficients of each component in each reaction (-)"]
Do_r[1:Nri, 1:Nc] = Reaction.Do_r, [description = "Forward order of the components (-) "]
Ef_r[1:Nri] = Reaction.Ef_r, [description = "Activation energy of each reaction at given temperature ()"]
A = Ac, [description = "Cross sectional area of the tank (m²)"]
end

#Ports creation
 
InPorts = [matcon(; Nc = Nc, name = Symbol("InPorts$i")) for i in 1:ninports]

OutPorts = @named begin
    Out = matcon(; Nc = Nc) 
end   

vars = @variables begin
    M(t), [description = "Mass holdup in the tank (kg)"]
    N(t), [description = "Total molar holdup in the tank (kmol)"]
    V(t), [description = "Volume holdup in the tank (m³)"]
    (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component in the tank (mol)"]
    (Cᵢ(t))[1:Nc], [description = "Concentration of each component in the tank (mol/m³)"]
    ρ(t), [description = "Molar Density of the fluid in the tank (mol/m³)"]
    ρʷ(t), [description = "Mass Density of the fluid in the tank (kg/m³)"]
    MW(t), [description = "Molecular weight of fluid in the tank (kg/kmol)"]  
    T(t), [description = "Temperature of vessel contents (K)"]  
    P_out(t), [description = "Pressure at the outlet stream level (Pa)"]
    H(t), [description = "Enthalpy holdup of the fluid in the tank (J)"] 
    #S(t), [description = "Entropy holdup of the fluid in the tank (J/K)"]
    F_out(t), [description = "Outlet molar flow rate (mol/s)"] 
    Fʷ_out(t), [description = "Outlet molar flow rate (mol/s)"]
    Q_out(t), [description = "Outlet volumetric flow rate (m³/s)"] # DoF
    height(t), [description = "Liquid level in vessel measured from bottom of the tank (m)"]

    (Cᵢ_in(t))[1:Nc, 1:ninports], [description = "Inlet concentration of each component (mol/m³)"] # DoF through inlet stream
    (F_in(t))[1:ninports], [description = "Inlet molar flow rate (mol/s)"] # DoF through inlet stream
    (Q_in(t))[1:ninports], [description = "Inlet volumetric flow rate(s) (m³/s)"]
    (T_in(t))[1:ninports], [description = "Inlet temperature (K)"] # DoF through inlet stream
    (h_in(t))[1:ninports], [description = "Inlet specific enthalpy (J/mol)"]
    (ρ_in(t))[1:ninports], [description = "Inlet density (mol/m³)"]
    (ρʷ_in(t))[1:ninports], [description = "Inlet density (mol/m³)"]
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
mass_density_eqs = [ρʷ_in[j] ~ InPorts[j].ρʷ for j in 1:ninports]
molar_density_eqs = [ρ_in[j] ~ InPorts[j].ρ for j in 1:ninports]
inletenthalpy = [h_in[j] ~ InPorts[j].H for j in 1:ninports]
inletconcentrations = [Cᵢ_in[i, j] ~ InPorts[j].z₁[i]*ρ_in[j] for j in 1:ninports for i in 1:Nc]
inlettemperature_eqs = [T_in[j] ~ InPorts[j].T for j in 1:ninports]
inletmolarflow_eqs = [F_in[j] ~ InPorts[j].F for j in 1:ninports]
volumetricflow_eqs = [Q_in[j] ~ F_in[j] / ρ_in[j] for j in 1:ninports]


#Outlet connector equations:
out_conn = [Out.P ~ - P_out
            Out.T ~ - T
            Out.F ~ - F_out
            Out.Fʷ ~ - Fʷ_out
            Out.H ~ - H/N
            Out.S ~  0.0
            Out.ρʷ ~ - ρʷ
            Out.ρ ~ - ρ
            scalarize(Out.z₁ .~ - Nᵢ/N)...
            Out.MW[1] ~ - MW
]

if phase == :liquid
out_conn_phases = [
                scalarize(Out.z₂ .~ 0.0)...
                scalarize(Out.z₃ .~ - Nᵢ/N)...
                Out.MW[2] ~  0.0
                Out.MW[3] ~ -MW
                Out.α_g ~  0.0]

elseif phase == :vapor
    out_conn_phases = [
    scalarize(Out.z₂ .~ - Nᵢ/N)...
    scalarize(Out.z₃ .~  0.0)...
    Out.MW[2] ~ -MW
    Out.MW[3] ~ 0.0
    Out.α_g ~ -1.0]
end


#balances
#mass_balance = [D(M) ~ sum(scalarize(Q_in.*ρʷ_in)) - Q_out*ρʷ]
component_balance = [D(Nᵢ[i]) ~ sum(scalarize(Q_in[:].*Cᵢ_in[i, :])) - Q_out*Cᵢ[i] + R[i]*V for i in 1:Nc] #Neglectable loss to vapor phase head space
energy_balance = [D(H) ~ sum(scalarize(Q_in.*ρ_in.*h_in)) - Q_out*ρ*H/N + Q̇]
jacket_energy_balance = [Q̇ ~ -2.27*4184.0*(T - 288.7)*(1.0 - exp(-8440.26/(2.27*4184)))] 
mass_volume_eq = [ρʷ*V ~ M, ρ*V ~ N] #Modified to calculate volume from N
mol_holdup = [N ~ sum(scalarize(Nᵢ))]
mol_to_concentration = [scalarize(Nᵢ .~ Cᵢ*V)...]
height_to_volume = [height*A ~ V]
volumetricflow_to_molarflow = [Q_out ~ F_out/ρ Q_out ~ sum(scalarize(Q_in))] # Modified
volumetricflow_to_massflow = [Q_out ~ Fʷ_out/ρʷ]
  
#Thermodynamic properties (outlet)
pressure_out = [phase == :liquid ? P_out ~ P_atm : P_out ~ P_atm] #Estimation considering static pressure (May be off as tank is agitated and not static)
density_eqs = [ρ ~ molar_density_simple(model, P_out, T, Nᵢ)] 
mass_density = [ρʷ ~ ρ*MW]
globalEnthalpy_eq = [H ~ enthalpy_simple(model, P_out, T, Nᵢ) + sum(scalarize(ΔH₀f.*Nᵢ))]
molar_mass = [MW ~ sum(scalarize(MWs.*Nᵢ))/N*gramsToKilograms]


eqs = [reaction_rate...; overall_reaction_rate...; atm_pressure...; mass_density_eqs...; molar_density_eqs...; inletenthalpy...; inletconcentrations...; inlettemperature_eqs...; inletmolarflow_eqs...; volumetricflow_eqs...; out_conn...;
out_conn_phases...; component_balance...; energy_balance...; jacket_energy_balance...; mass_volume_eq...; mol_holdup...; mol_to_concentration...; height_to_volume...; volumetricflow_to_molarflow...; volumetricflow_to_massflow...;
pressure_out...; density_eqs...; mass_density...; globalEnthalpy_eq...; molar_mass...]

unfold_pars = []
for par in pars
    unfold_pars = [unfold_pars...; par...]
end

unfold_vars = []
for var in vars
    unfold_vars = [unfold_vars...; var...]
end

ODESystem([eqs...;], t, unfold_vars, unfold_pars; name, systems = [InPorts...; OutPorts])

end

#= 



struct my_model
    Cp
    ρ_coefs
end

function enthalpy_simple(m::my_model, P, T, N)
    sum(m.Cp[i]*N[i]*(T - 298.15) for i in eachindex(N))
end

function molar_density_simple(m::my_model, P, T, N)
        sum(N)/sum(N[i]/(m.ρ_coefs["a"][i] + m.ρ_coefs["b"][i]*T + m.ρ_coefs["c"][i]*T^2  + m.ρ_coefs["d"][i]*T^3) for i in eachindex(N))  
end

@register_symbolic enthalpy_simple(m::my_model, P::Num, T::Num, N::Vector{Num})
@register_symbolic molar_density_simple(m::my_model, P::Num, T::Num, N::Vector{Num}) =#


#= using ModelingToolkit, JSON, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics
import ModelingToolkit: scalarize, equations, get_unknowns, defaults
using Clapeyron, GCIdentifier
using NonlinearSolve
using JSON
@parameters t
D = Differential(t)
include("C:/Users/Vinic/OneDrive/Pos-Doc/Flowsheeting/ProcessModeling/ProcessSimulator.jl/src/utils")
include("C:/Users/Vinic/OneDrive/Pos-Doc/Flowsheeting/ProcessModeling/ProcessSimulator.jl/src/Sources/Sourceutils.jl")
include("C:/Users/Vinic/OneDrive/Pos-Doc/Flowsheeting/ProcessModeling/ProcessSimulator.jl/src/Reactors/ReactionManager/KineticReaction.jl")
substances_user = ["water", "methanol", "propyleneglycol","methyloxirane"]
Nc = size(substances_user, 1)
properties = Dict(subs => load_component_properties(subs) for subs in substances_user)
# Function to extract parameters for ReidIdeal model
Cps = [69.21, 80.66, 192.50, 118.1] # At 298.00 K
pho_coef = Dict("a" => [78821.04345816873, 38133.33588802956, 18453.26055238924, 25211.86290505424], 
"b" => [-114.80261704286386, -85.22435577664774, -26.084451261000222, -73.31971618413455],
"c" => [0.208885275623327, 0.21416405813277503, 0.046784549601356355, 0.18909331028746998],
"d" => [-0.00022293440498512672, -0.0002675908503439664, -4.722426797584051e-5, -0.00022514899466957527])

phase = :liquid

struct my_model
    Cp
    ρ_coefs
end

function enthalpy_simple(m::my_model, P, T, N)
    sum(m.Cp[i]*N[i]*(T - 298.15) for i in eachindex(N))
end

function molar_density_simple(m::my_model, P, T, N)
        sum(N)/sum(N[i]/(m.ρ_coefs["a"][i] + m.ρ_coefs["b"][i]*T + m.ρ_coefs["c"][i]*T^2  + m.ρ_coefs["d"][i]*T^3) for i in eachindex(N))  
end


mymodel = my_model(Cps, pho_coef)
enthalpy_simple(mymodel, 101325, 298.15, [1.0, 1.0, 1.0, 1.0])
molar_density_simple(mymodel, 101325, 350.15, [1.0, 1.0, 1.0, 1.0])

@register_symbolic enthalpy_simple(m::my_model, P::Num, T::Num, N::Vector{Num})
@register_symbolic molar_density_simple(m::my_model, P::Num, T::Num, N::Vector{Num})


Reaction = KineticReactionNetwork(;substances_user = substances_user, 
Af_r = 4.71e9, Ef_r = 32400*1055.6/453.6, Coef_Cr = [-1.0 0.0 1.0 -1.0], 
Do_r = [1.0 0.0 0.0 1.0], name = "Propyleneglycol synthesis")
#["water", "methanol", "propyleneglycol","methyloxirane"]
ninports = 1
InPorts = [matcon(; Nc = Nc, name = Symbol("InPorts$i")) for i in 1:ninports]

Ac = 26.51 # m²
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
    #S(t), [description = "Entropy holdup of the fluid in the tank (J/K)"]
    F_out(t), [description = "Outlet molar flow rate (mol/s)"] 
    Fʷ_out(t), [description = "Outlet molar flow rate (mol/s)"]
    Q_out(t), [description = "Outlet volumetric flow rate (m³/s)"] # DoF
    height(t), [description = "Liquid level in vessel measured from bottom of the tank (m)"]

    (Cᵢ_in(t))[1:Nc, 1:ninports], [description = "Inlet concentration of each component (mol/m³)"] # DoF through inlet stream
    (F_in(t))[1:ninports], [description = "Inlet molar flow rate (mol/s)"] # DoF through inlet stream
    (Q_in(t))[1:ninports], [description = "Inlet volumetric flow rate(s) (m³/s)"]
    (T_in(t))[1:ninports], [description = "Inlet temperature (K)"] # DoF through inlet stream
    (h_in(t))[1:ninports], [description = "Inlet specific enthalpy (J/mol)"]
    (ρ_in(t))[1:ninports], [description = "Inlet density (mol/m³)"]
    (ρʷ_in(t))[1:ninports], [description = "Inlet density (mol/m³)"]
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
mass_density_eqs = [ρʷ_in[j] ~ InPorts[j].ρʷ for j in 1:ninports]
molar_density_eqs = [ρ_in[j] ~ InPorts[j].ρ for j in 1:ninports]
inletenthalpy = [h_in[j] ~ InPorts[j].H for j in 1:ninports]
inletconcentrations = [Cᵢ_in[i, j] ~ InPorts[j].z₁[i]*ρ_in[j] for j in 1:ninports for i in 1:Nc]
inlettemperature_eqs = [T_in[j] ~ InPorts[j].T for j in 1:ninports]
inletmolarflow_eqs = [F_in[j] ~ InPorts[j].F for j in 1:ninports]
volumetricflow_eqs = [Q_in[j] ~ F_in[j] / ρ_in[j] for j in 1:ninports]


#Outlet connector equations:
out_conn = [Out.P ~ P_out
            Out.T ~ T
            Out.F ~ - F_out
            Out.Fʷ ~ - Fʷ_out
            Out.H ~ H/N
            Out.S ~ 0.0
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
energy_balance = [D(H) ~ sum(scalarize(Q_in.*ρ_in.*h_in)) - Q_out*ρ*H/N + Q̇]
jacket_energy_balance = [Q̇ ~ -2.27*4184.0*(T - 288.7)*(1.0 - exp(-8440.26/(2.27*4184)))] 
mass_volume_eq = [ρ*V ~ N]
mol_holdup = [N ~ sum(scalarize(Nᵢ))]
mol_to_concentration = [scalarize(Nᵢ .~ Cᵢ*V)...]
height_to_volume = [height*A ~ V]
volumetricflow_to_molarflow = [Q_out ~ F_out/ρ]
volumetricflow_to_massflow = [Q_out ~ Fʷ_out/ρʷ]
  
#Thermodynamic properties (outlet)
pressure_out = [phase == :liquid ? P_out ~ P_atm + ρʷ*GravitationalConst*(height - height_out) : P_out ~ P_atm] #Estimation considering static pressure (May be off as tank is agitated and not static)
density_eqs = [ρ ~ molar_density_simple(mymodel, P_out, T, Nᵢ)] 
mass_density = [ρʷ ~ ρ*MW]
globalEnthalpy_eq = [H ~ enthalpy_simple(mymodel, P_out, T, Nᵢ) + sum(scalarize(ΔH₀f.*Nᵢ))]
molar_mass = [MW ~ sum(scalarize(MWs.*Nᵢ)/N)*gramsToKilograms]


eqs = [reaction_rate...; overall_reaction_rate...; atm_pressure...; mass_density_eqs...; molar_density_eqs...; inletenthalpy...; inletconcentrations...; inlettemperature_eqs...; inletmolarflow_eqs...; volumetricflow_eqs...; out_conn...;
out_conn_phases...; mass_balance...; component_balance...; energy_balance...; jacket_energy_balance...; mass_volume_eq...; mol_holdup...; mol_to_concentration...; height_to_volume...; volumetricflow_to_molarflow...; volumetricflow_to_massflow...;
pressure_out...; density_eqs...; mass_density...; globalEnthalpy_eq...; molar_mass...]

eqs_cstr = [reaction_rate...; overall_reaction_rate...; atm_pressure...; mass_density_eqs...; molar_density_eqs...; inletenthalpy...; inletconcentrations...; inlettemperature_eqs...; 
inletmolarflow_eqs...; volumetricflow_eqs...; mass_balance...; component_balance...; energy_balance...; jacket_energy_balance...; mass_volume_eq...; mol_holdup...; mol_to_concentration...; height_to_volume...; volumetricflow_to_molarflow...; volumetricflow_to_massflow...;
pressure_out...; density_eqs...; mass_density...; globalEnthalpy_eq...; molar_mass...]

unfold_pars = []
for par in pars
    unfold_pars = [unfold_pars...; par...]
end

unfold_vars = []
for var in vars
    unfold_vars = [unfold_vars...; var...]
end =#