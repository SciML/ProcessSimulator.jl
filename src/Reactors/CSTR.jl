using ModelingToolkit, Clapeyron
@parameters t
D = Differential(t)

@component function CSTR(; InPorts = [matcon(; Nc = Nc, name = Symbol("InPorts$i")) for i in 1:1], Ac = 1.0, P = 101325.0, phase = :liquid, substances_user = ["methane", "carbon monoxide"], Nc = size(substances_user, 1),
    model = PR(substances_user), Reaction = KineticReactionNetwork(Af_r = [1.0], Ef_r = [1.0], Do_r = [1.0 1.0], name = :DefaultReaction))

    Nri = defaults(Reaction)[@nonamespace Reaction.Nr]
    Ni_InPorts = size(InPorts, 1)

    pars = @parameters begin 
        ## How to inherent parameters from Reaction?
        N_InPorts = Ni_InPorts
        Nr = Nri, [description = "Number of reactions"]
        Af_r[1:Nri] = defaults(Reaction)[@nonamespace Reaction.Af_r], [description = "Arrhenius constant of each reaction at given temperature ()"]
        Coef_Cr[1:Nri, 1:Nc] = defaults(Reaction)[@nonamespace Reaction.Coef_Cr], [description = "Stoichiometric coefficients of each component in each reaction (-)"]
        Do_r[1:Nri, 1:Nc] = defaults(Reaction)[@nonamespace Reaction.Do_r], [description = "Forward order of the components (-) "]
        Ef_r[1:Nri] = defaults(Reaction)[@nonamespace Reaction.Ef_r], [description = "Activation energy of each reaction at given temperature ()"]
        N = Nc, [description = "Number of components"]
        A = Ac, [description = "Cross sectional area of the tank (m²)"]
    end   
    
    @extend Reaction
    
    OutPorts = @named begin
        Out = matcon(; Nc = Nc) 
    end   

    vars = @variables begin
    T(t), [description = "Temperature of vessel contents (K)"]
    M(t), [description = "Mass of fluid in the tank (kg)"]
    N(t), [description = "Mols of fluid in the tank (kmol)"]
    MW(t), [description = "Molecular weight of fluid in the tank (kg/kmol)"]    
    V(t), [description = "Volume of fluid in the tank (m³)"]
    (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component in the tank(mol)"]
    (Cᵢ(t))[1:Nc], [description = "Concentration of each component in the tank (mol/m³)"]
    (ρ(t)), [description = "Molar Density of the fluid in the tank (mol/m³)"]
    (ρʷ(t)), [description = "Mass Density of the fluid in the tank (kg/m³)"]
    (H(t)), [description = "Enthalpy holdup of the fluid in the tank (J)"] 
    (Cᵢ_in(t))[1:Nc, 1:Ni_InPorts], [description = "Inlet concentration of each component (mol/m³)"]
    (Q_in(t))[1:Ni_InPorts], [description = "Inlet volumetric flow rate(s) (m³/s)"]
    Q_out(t), [description = "Inlet volumetric flow rate (m³/s)"]
    Q̃(t), [description = "Heat transfer from fluid in jacket (J/s)"]
    (T_in(t))[1:Ni_InPorts], [description = "Inlet temperature (K)"]
    (h_in(t))[1:Ni_InPorts], [description = "Inlet specific enthalpy (J/mol)"]
    (ρ_in(t))[1:Ni_InPorts], [description = "Inlet density (mol/m³)"]
    (ρʷ_in(t))[1:Ni_InPorts], [description = "Inlet density (mol/m³)"]
    h_out(t), [description = "Outlet specific enthalpy (J/mol)"]
    P_bottom(t), [description = "Bottom pressure of the tank (Pa)"]
    P_atm(t), [description = "Tank pressure (Pa)"]
    (r(t))[1:Nri, 1:Nc], [description = "Rate of each reaction for each component (mol/s/m³)"]
    (R(t))[1:Nc], [description = "Overall reaction rate (mol/s/m³)"]
    h(t), [description = "Liquid level in vessel (m)"]
    end


    #Reaction equations
    reaction_rate = [r[i, j] ~ Af_r[i]*exp(-Ef_r[i]/(R*T))*(Cᵢ[:].^(Do_r[j, :])) for i in 1:Nri for j in 1:Nc] # If there's an inert, the order is just zero, but it has to be written

    overall_reaction_rate = [R[i] ~ sum(r[:, i].*Coef_Cr[:, i]) for i in 1:Nc]  # If there's an inert, the coefficient is just zero, but it has to be written

    
    #Inlet variables

    volumetricflow = [Q_in[j] ~ InPorts[j].F * InPorts[j].ρ for j in 1:Ni_InPorts]

    #balances
    mass_balance = [D(M) ~ Q_in*ρʷ_in - Q_out*ρʷ]

    component_balance = [d(Nᵢ[i]) ~ Q_in*Cᵢ_in[i] - Q_out*Cᵢ[i] + R[i]*V for i in 1:Nc]

    energy_balance = [D(H(t)) ~ Q_in*ρ_in*h_in - Q_out*h_out*ρ + Q̃]

    mass_volume_eq = [ρʷ*V ~ M]

    mol_holdup = [N ~ sum(Nᵢ)]

    mol_to_concentration = [Nᵢ .~ Cᵢ*V]

    #Thermodynamic properties

    globalEnthalpy = [H ~ isobaric_heat_capacity(model, P_atm, T, Nᵢ; phase = phase)]

    molarEnthalpy = [h_out*N ~ H]

end
