@component function FlashDrum(; substances_user, 
    Nc = length(substances_user), 
    model,
    ninports,  
    Volume, 
    name,
    guesses 
   )
  

    
#Constants
gramsToKilograms = 10^(-3)
Rᵍ = 8.314 # J/(mol K)
Nri = Reaction.Nri

#Properties of individual substances
properties = Dict(subs => load_component_properties(subs) for subs in substances_user)
MWs = [properties[subs]["MW"] for subs in substances_user]


pars = @parameters begin 
V = Volume, [description = "Drum volume (m3)"]
end

#Ports creation
 
InPorts = [matcon(; Nc = Nc, name = Symbol("InPorts$i")) for i in 1:ninports]

OutPorts = @named begin
    Out = matcon(; Nc = Nc) 
end   

vars = @variables begin
    M(t)[1:3], [description = "Mass holdup in the tank (kg)"]
    N(t)[1:3], [description = "Total molar holdup in the tank (kmol)"]
    V(t)[1:3], [description = "Volume holdup in the tank (m³)", guess = guesses[:V]]
    (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component in the tank (mol)"]
    (Cᵢ(t))[1:Nc], [description = "Concentration of each component in the tank (mol/m³)", guess = guesses[:Cᵢ]]
    ρ(t), [description = "Molar Density of the fluid in the tank (mol/m³)"]
    ρʷ(t), [description = "Mass Density of the fluid in the tank (kg/m³)"]
    MW(t), [description = "Molecular weight of fluid in the tank (kg/kmol)"]  
    T(t), [description = "Temperature of vessel contents (K)"]  
    P_out(t), [description = "Pressure at the outlet stream level (Pa)"]
    H(t), [description = "Enthalpy holdup of the fluid in the tank (J)"] 
    #S(t), [description = "Entropy holdup of the fluid in the tank (J/K)"]
    F_out(t), [description = "Outlet molar flow rate (mol/s)", guess = guesses[:F_out]] 
    Fʷ_out(t), [description = "Outlet molar flow rate (mol/s)", guess = guesses[:Fʷ_out]]
    Q_out(t), [description = "Outlet volumetric flow rate (m³/s)"] # DoF
    height(t), [description = "Liquid level in vessel measured from bottom of the tank (m)"]
    P_atm(t), [description = "Tank pressure (Pa)"] # Equal to inlet pressures.
    Q̇(t), [description = "Heat transfer rate (J/s)"] # Potential DoF
    (r(t))[1:Nri], [description = "Rate of each reaction for each component (mol/s/m³)"]
    (R(t))[1:Nc], [description = "Overall reaction rate (mol/s/m³)"]

    (Cᵢ_in(t))[1:Nc, 1:ninports], [description = "Inlet concentration of each component (mol/m³)"] # DoF through inlet stream
    (F_in(t))[1:ninports], [description = "Inlet molar flow rate (mol/s)"] # DoF through inlet stream
    (Q_in(t))[1:ninports], [description = "Inlet volumetric flow rate(s) (m³/s)"]
    (T_in(t))[1:ninports], [description = "Inlet temperature (K)"] # DoF through inlet stream
    (h_in(t))[1:ninports], [description = "Inlet specific enthalpy (J/mol)"]
    (ρ_in(t))[1:ninports], [description = "Inlet density (mol/m³)"]
    (ρʷ_in(t))[1:ninports], [description = "Inlet density (mol/m³)"]

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
            Out.F ~ F_out
            Out.Fʷ ~ Fʷ_out
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
                scalarize(Out.z₃ .~ Nᵢ/N)...
                Out.MW[2] ~  0.0
                Out.MW[3] ~ MW
                Out.α_g ~ 0.0]

elseif phase == :vapor
    out_conn_phases = [
    scalarize(Out.z₂ .~ Nᵢ/N)...
    scalarize(Out.z₃ .~ 0.0)...
    Out.MW[2] ~ MW
    Out.MW[3] ~ 0.0
    Out.α_g ~ 1.0]
end


#balances
#mass_balance = [D(M) ~ sum(scalarize(Q_in.*ρʷ_in)) - Q_out*ρʷ]
component_balance = [D(Nᵢ[i]) ~ sum(scalarize(Q_in[:].*Cᵢ_in[i, :])) - Q_out*Cᵢ[i] + R[i]*V for i in 1:Nc] #Neglectable loss to vapor phase head space
energy_balance = [D(H) ~ sum(scalarize(Q_in.*ρ_in.*h_in)) - Q_out*ρ*H/N + Q̇]
jacket_energy_balance = [Q̇ ~ -2.27*4184.0*(T - 288.7)*(1.0 - exp(-8440.26/(2.27*4184)))] 
mass_volume_eq = [ρʷ*V ~ M, ρ*V ~ N] #Modified to calculate volume from N
mol_holdup = [N ~ sum(collect(Nᵢ))]
mol_to_concentration = [scalarize(Nᵢ .~ Cᵢ*V)...]
height_to_volume = [height*A ~ V]
volumetricflow_to_molarflow = [Q_out ~ F_out/ρ, Q_out ~ sum(scalarize(Q_in))] # Modified
volumetricflow_to_massflow = [Q_out ~ Fʷ_out/ρʷ]
  
#Thermodynamic properties (outlet)
pressure_out = [phase == :liquid ? P_out ~ P_atm : P_out ~ P_atm] #Estimation considering static pressure (May be off as tank is agitated and not static)
density_eqs = [ρ ~ molar_density(model, P_out, T, Nᵢ, phase = "unknown")] 
mass_density = [ρʷ ~ ρ*MW]
globalEnthalpy_eq = [H ~ enthalpy(model, P_out, T, Nᵢ, phase = "unknown") + sum(scalarize(ΔH₀f.*Nᵢ))]
molar_mass = [MW ~ sum(scalarize(MWs.*Nᵢ))/N*gramsToKilograms]


eqs = [reaction_rate...; overall_reaction_rate...; atm_pressure...; mass_density_eqs...; molar_density_eqs...; inletenthalpy...; inletconcentrations...; inlettemperature_eqs...; inletmolarflow_eqs...; volumetricflow_eqs...; out_conn...;
out_conn_phases...; component_balance...; energy_balance...; jacket_energy_balance...; mass_volume_eq...; mol_holdup...; mol_to_concentration...; height_to_volume...; volumetricflow_to_molarflow...; volumetricflow_to_massflow...;
pressure_out...; density_eqs...; mass_density...; globalEnthalpy_eq...; molar_mass...]


ODESystem([eqs...;], t, collect(Iterators.flatten(vars)),
 collect(Iterators.flatten(pars)); name, systems = [InPorts...; OutPorts])

end

