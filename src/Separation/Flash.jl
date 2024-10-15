@component function ThreePortDrum(; substances,
    non_volatiles::Vector{Bool} = falses(length(substances)),
    non_condensables::Vector{Bool} = falses(length(substances)),
    Ac,
    Nc = length(substances), 
    model,
    name,
   )
  

    
#Constants
g = 9.81 # m/s²


pars = []


Ports = @named begin
    In = matcon(; Nc = Nc)
    Outᴸ = matcon(; Nc = Nc)
    Outᵍ = matcon(; Nc = Nc)
    #EnergyCon = thermal_energy_connector()
end    

vars = @variables begin
    (Nᵢᴸ(t))[1:Nc], [description = "Molar holdup of each component in liquid phase (mol)"]
    (Nᵢᵍ(t))[1:Nc], [description = "Molar holdup of each component in gas phase (mol)"]
    (xᵢ(t))[1:Nc], [description = "Molar fraction in liquid phase (-)"]
    (yᵢ(t))[1:Nc], [description = "Molar fraction in gas phase (-)"]
    (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component (mol)"]
    (ϕᵢᴸ(t))[1:Nc], [description = "Fugacity coefficient in liquid phase"]
    (ϕᵢᵍ(t))[1:Nc], [description = "Fugacity coefficient in liquid phase"]
    Vᴸ(t), [description = "Volume of the liquid phase (m³)"]
    Vᵍ(t), [description = "Volume of the gas phase (m³)"]
    ρᴸ(t), [description = "Molar Density of gas phase (mol/m³)"]
    ρᵍ(t), [description = "Molar Density of liquid phase (mol/m³)"]
    Nᴸ(t), [description = "Total molar holdup in the liquid phase (mol)"] 
    Nᵍ(t), [description = "Total molar holdup in the gas phase (mol)"] 

    T(t), [description = "Drum temperature (K)"]  
    P(t), [description = "Drum pressue (Pa)"]
    U(t), [description = "Total internal energy holdup in the tank (L + G) (J)"]
    V(t), [description = "Total volume in the tank (L + G) (m³)"]

    hᴸ_outflow(t), [description = "Outflow Enthalpy of liquid phase (J/mol)"]
    hᵍ_outflow(t), [description = "Outflow Enthalpy of liquid phase (J/mol)"]
    Fᴸ_outflow(t), [description = "Outlet molar flow rate of liquid phase (mol/s)"] 
    Fᵍ_outflow(t), [description = "Outlet molar flow rate of gas phase (mol/s)"] 
    height(t), [description = "Liquid level in vessel measured from bottom of the tank (m)"]
    Q̇(t), [description = "Heat transfer rate (J/s)"]


    _0_Nᴸ(t), [description = "Condition for all gas phase"]
    _0_Nᵍ(t), [description = "Condition for all liquid phase"]

end

    
#Conservation equations

liquid_gas_holdups = [Nᴸ ~ sum(collect(Nᵢᴸ)), Nᵍ ~ sum(collect(Nᵢᵍ))]

xy_fractions = [scalarize(xᵢ.*Nᴸ .~ Nᵢᴸ)...; scalarize(yᵢ.*Nᵍ .~ Nᵢᵍ)]

component_holdups = [Nᵢ[i] ~ Nᵢᴸ[i] + Nᵢᵍ[i] for i in 1:Nc]


component_balance = [D(Nᵢ[i]) ~ In.F*In.z₁[i] - Fᴸ_outflow*xᵢ[i] - Fᵍ_outflow*yᵢ[i] for i in 1:Nc]

energy_balance = [D(U) ~ In.F*In.H - Fᴸ_outflow*hᴸ_outflow - Fᵍ_outflow*hᵍ_outflow + Q̇]

heat_source_equation = [Q̇ ~ 0.0]

#Thermodynamic constraints

internal_energy = [U + P*V ~ Nᴸ*hᴸ_outflow + Nᵍ*hᵍ_outflow] 

total_volume = [V ~ Vᴸ + Vᵍ] #That should not be fixed a priori. One may want to have a perfect pressure control vessel. 

liquid_gas_enthalpy = [hᴸ_outflow ~ enthalpy(model, P, T, xᵢ, phase = "liquid")
 hᵍ_outflow ~ enthalpy(model, P, T, yᵢ, phase = "vapor")]

liquid_gas_densities = [ρᴸ ~ molar_density(model, P, T, Nᵢᴸ, phase = "liquid")
 ρᵍ ~ molar_density(model, P, T, Nᵢᵍ, phase = "vapor")]

liquid_volume = [ρᴸ*Vᴸ ~ Nᴸ]

gas_volume = [ρᵍ*Vᵍ ~ Nᵍ]

fugacities_g = [ϕᵢᵍ[i] ~ fugacity_coefficient(model, P, T, Nᵢᵍ, phase = "vapor")[i] for i in 1:Nc if ~non_condensables[i] & ~non_volatiles[i]]

fugacities_l = [ϕᵢᴸ[i] ~ fugacity_coefficient(model, P, T, Nᵢᴸ, phase = "liquid")[i] for i in 1:Nc if ~non_condensables[i] & ~non_volatiles[i]]

non_smoothness = [_0_Nᵍ ~ Nᵍ/(Nᵍ + Nᴸ)
 _0_Nᴸ ~ Nᵍ/(Nᵍ + Nᴸ) - 1.0
 0.0 ~ _0_Nᵍ + _0_Nᴸ + sum(collect(xᵢ)) - sum(collect(yᵢ)) - min(_0_Nᵍ, _0_Nᴸ, sum(collect(xᵢ)) - sum(collect(yᵢ))) - max(_0_Nᵍ, _0_Nᴸ, sum(collect(xᵢ)) - sum(collect(yᵢ))) 
 ] #Non-smooth formulation 

 equilibrium = [ϕᵢᵍ[i].*yᵢ[i] .~ ϕᵢᴸ[i].*xᵢ[i] for i in 1:Nc if ~non_condensables[i] & ~non_volatiles[i]]

 # Out connectors

 out_conn_L = [Outᴸ.P ~ P + ρᴸ*g*Ac/Vᴸ
 Outᴸ.T ~ T
 Outᴸ.F ~ Fᴸ_outflow
 Outᴸ.H ~ hᴸ_outflow
 scalarize(Outᴸ.z₁ .~ xᵢ)...
 scalarize(Outᴸ.z₂ .~ 0.0)...
 scalarize(Outᴸ.z₃ .~ xᵢ)...
 Outᴸ.α_g ~ 0.0
]

out_conn_g = [Outᵍ.P ~ P
Outᵍ.T ~ T
Outᵍ.F ~ Fᵍ_outflow
Outᵍ.H ~ hᵍ_outflow
scalarize(Outᵍ.z₁ .~ yᵢ)...
scalarize(Outᵍ.z₂ .~ yᵢ)...
scalarize(Outᵍ.z₃ .~ 0.0)...
Outᵍ.α_g ~ 1.0
]



eqs = [liquid_gas_holdups...; xy_fractions...; component_holdups...; component_balance...; energy_balance...; 
heat_source_equation...; internal_energy...; total_volume...; liquid_gas_enthalpy...; liquid_gas_densities...; liquid_volume...;
gas_volume...; fugacities_g...;  fugacities_l...; non_smoothness...; equilibrium...; out_conn_L...; out_conn_g...]


ODESystem([eqs...;], t, collect(Iterators.flatten(vars)), collect(Iterators.flatten(pars)); name, systems = [Ports...])

end

export ThreePortDrum