@component function MaterialSource(; substances_user = ["methane", "carbon monoxide"],
     Nc = size(substances_user, 1),
     model = PR(substances_user),
    properties = Dict(subs => load_component_properties(subs) for subs in substances_user),
    P_user = 101325.0, T_user = 25.0 + 273.15,
     Fₜ_user = 100.0, zₜ_user = [0.5, 0.5], name)
    #phase 1 is total, 2 is vapor, 3 is liquid

    
    Tcs = [properties[subs]["Tc"] for subs in substances_user] #Critical temperature (K)
    MWs = [properties[subs]["MW"] for subs in substances_user] #Molecular weight (g/mol)
    ΔH₀f = [properties[subs]["IGHF"]/10^3 for subs in substances_user] # (IG formation enthalpy) J/mol
    gramsToKilograms = 10^(-3)

    pars = @parameters begin
        N = Nc, [description = "Number of components"]
        substances[1:Nc] = substances_user, [description = "Component names"]
        P = P_user, [description = "Pressure (Pa)"]
        T = T_user, [description = "Temperature (K)"]
        zₜ[1:Nc] = zₜ_user, [description = "Components molar fraction (-)"]
        Fₜ = Fₜ_user, [description = "Total molar flow rate (mol/s)"]
    end 

    vars = @variables begin
        Tc(t), [description = "Critical temperature (K)"]
        Pc(t), [description = "Critical pressure (Pa)"]
        P_buble(t), [description = "Bubble point pressure (Pa)"]
        P_dew(t), [description = "Dew point pressure (Pa)"]
        α_g(t), [description = "Vapor phase molar fraction"]
        α_l(t), [description = "Liquid phase molar fraction"]
        αᵂ_g(t), [description = "Vapor phase mass fraction"]
        αᵂ_l(t), [description = "Liquid phase mass fraction"]
        (Fⱼ(t))[1:3], [description = "Molar flow rate in each phase j (mol/s)"]
        (Fᵂⱼ(t))[1:3], [description = "Mass flow rate in each phase j (kg/s)"]
        (zⱼᵢ(t))[1:3, 1:Nc], [description = "Component molar fraction in each phase j component i (-)"] 
        (zᵂⱼᵢ(t))[1:3, 1:Nc], [description = "Component mass fraction in each phase j component i (-)"] 
        (Fⱼᵢ(t))[1:3, 1:Nc], [description = "Molar flow rate in each phase j and component i (mol/s)"] 
        (Fᵂⱼᵢ(t))[1:3, 1:Nc], [description = "Mass flow rate in each phase j and component i (mol/s)"]
        (Hⱼ(t))[1:3], [description = "Enthalpy in each phase j at T and P (J/mol)"] 
        (Sⱼ(t))[1:3], [description = "Entropy in each phase j at T and P (J/mol.K)"]
        (ρ(t))[1:3], [description = "Molar density in each phase j (mol/m³)"]
        (ρʷ(t))[1:3], [description = "Mass density in each phase j(mol/m³)"]
        (MWⱼ(t))[1:3], [description = "Molar mass of each phase j (kg/mol)"] 
    end

    systems = @named begin 
        Out = matcon(; Nc = Nc)
    end

    #Connector equations
    eqs_conn = [
        Out.P ~ P 
        Out.T ~ T  
        Out.F ~ - Fₜ # F is negative as it is leaving the component 
        Out.Fʷ ~ - Fᵂⱼ[1]
        Out.H ~ Hⱼ[1] 
        Out.S ~ Sⱼ[1] 
        scalarize(Out.z₁ .- zⱼᵢ[1, :] .~ 0.0)...
        scalarize(Out.z₂ .- zⱼᵢ[2, :] .~ 0.0)...
        scalarize(Out.z₃ .- zⱼᵢ[3, :] .~ 0.0)...
        Out.α_g ~ α_g 
        Out.ρ ~ ρ[1]
        Out.ρʷ ~ ρʷ[1] 
        scalarize(Out.MW .~ MWⱼ)...
    ]

    #Global Mass and Molar balances
    global_mol_balance = [       
        Fⱼ[1] ~ Fₜ  
        Fⱼ[1] - Fⱼ[2] - Fⱼ[3] ~ 0.0 
        #Fᵂⱼ[1] - Fᵂⱼ[2] - Fᵂⱼ[3] ~ 0.0
        α_l + α_g ~ 1.0
        Fⱼ[2] ~ α_g*Fⱼ[1] 

        αᵂ_g + αᵂ_l ~ 1.0 
        αᵂ_g ~ Fᵂⱼ[2]/Fᵂⱼ[1]   
        
    ]

    # Mass based Equations 

    molar_to_mass = [	
        scalarize(Fᵂⱼ[:] .~ sum(Fᵂⱼᵢ[:, :], dims = 2))...
    ]

    molar_to_mass_2 = [scalarize(Fᵂⱼᵢ[:, i] .~ (MWs[i]*gramsToKilograms)*Fⱼᵢ[:, i] ) for i in 1:Nc]
    
    molar_to_mass = [molar_to_mass...; molar_to_mass_2...]
    
    #Component Molar and Mass Balance
    component_balance = [scalarize(Fⱼᵢ[:, i] .~ Fⱼ[:].*zⱼᵢ[:, i]) for i in 1:Nc]


    #Phase check
    #Tci, Pci, Vc = crit_mix(model, zₜ_user) #Critical point of the mixture
    Tci = minimum(Tcs)
    Pci = NaN
    Pdew, Vᵍdew, Vˡdew, xdew = dew_pressure(model, T_user, zₜ_user)
    Pbuble, Vᵍbuble, Vˡbuble, xbuble = bubble_pressure(model, T_user, zₜ_user)

    if T_user ≥ Tci
        pc = [α_g ~ 1.0
            Tc ~ Tci
            Pc ~ Pci
            P_buble ~ Pbuble
            P_dew ~ Pdew
            Hⱼ[1] ~ Hⱼ[2]
            Hⱼ[3] ~ 0.0
            Hⱼ[2] ~ enthalpy(model, P_user, T_user, zₜ_user, phase = :vapor)
            Sⱼ[1] ~ Sⱼ[2]
            Sⱼ[3] ~ 0.0
            Sⱼ[2] ~ entropy(model, P_user, T_user, zₜ_user, phase = :vapor)
            scalarize(zⱼᵢ[2, :] .~ zₜ)...
            scalarize(zⱼᵢ[3, :] .~ 0.0)...
            scalarize(zⱼᵢ[1, :] .~ zₜ)...
            scalarize(zᵂⱼᵢ[3, :] .~ 0.0)... # Mass base 
            scalarize(zᵂⱼᵢ[1, :] .~ zᵂⱼᵢ[2, :])... # Mass base
            scalarize(Fᵂⱼᵢ[2, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[2, :])... # Mass base
            ρ[1] ~ ρ[2]
            ρ[2] ~ molar_density(model, P_user, T_user, zₜ_user, phase = :vapor)
            ρ[3] ~ 0.0
            ρʷ[1] ~ ρʷ[2]
            ρʷ[2] ~ mass_density(model, P_user, T_user, zₜ_user, phase = :vapor)
            ρʷ[3] ~ 0.0
            MWⱼ[1] ~ MWⱼ[2]
            MWⱼ[2] ~ sum(MWs.*zⱼᵢ[2, :]*gramsToKilograms)
            MWⱼ[3] ~ 0.0         
              ]
    else

        if P_user < Pdew
            pc = [α_g ~ 1.0
                Tc ~ Tci
                Pc ~ Pci
                P_buble ~ Pbuble
                P_dew ~ Pdew
                Hⱼ[1] ~ Hⱼ[2]
                Hⱼ[2] ~ enthalpy(model, P_user, T_user, zₜ_user, phase = :vapor) + sum(zₜ_user.*ΔH₀f)
                Hⱼ[3] ~ 0.0
                Sⱼ[1] ~ Sⱼ[2]
                Sⱼ[2] ~ entropy(model, P_user, T_user, zₜ_user, phase = :vapor) 
                Sⱼ[3] ~ 0.0
                scalarize(zⱼᵢ[2, :] .~ zₜ)...
                scalarize(zⱼᵢ[3, :] .~ 0.0)...
                scalarize(zⱼᵢ[1, :] .~ zₜ)...
                scalarize(zᵂⱼᵢ[3, :] .~ 0.0)... # Mass base 
                scalarize(zᵂⱼᵢ[1, :] .~ zᵂⱼᵢ[2, :])... # Mass base
                scalarize(Fᵂⱼᵢ[2, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[2, :])... # Mass base
                ρ[1] ~ ρ[2]
                ρ[2] ~ molar_density(model, P_user, T_user, zₜ_user, phase = :vapor)
                ρ[3] ~ 0.0
                ρʷ[1] ~ ρʷ[2]
                ρʷ[2] ~ mass_density(model, P_user, T_user, zₜ_user, phase = :vapor)
                ρʷ[3] ~ 0.0
                MWⱼ[1] ~ MWⱼ[2]
                MWⱼ[2] ~ sum(MWs.*zⱼᵢ[2, :]*gramsToKilograms)
                MWⱼ[3] ~ 0.0  

            ] 

        elseif P_user > Pbuble
            pc = [α_g  ~ 0.0
                Tc ~ Tci
                Pc ~ Pci
                P_buble ~ Pbuble
                P_dew ~ Pdew
                Hⱼ[1] ~ Hⱼ[3]
                Hⱼ[3] ~ enthalpy(model, P_user, T_user, zₜ_user, phase = :liquid) + sum(zₜ_user.*ΔH₀f)
                Hⱼ[2] ~ 0.0
                Sⱼ[1] ~ Sⱼ[3]
                Sⱼ[3] ~ entropy(model, P_user, T_user, zₜ_user, phase = :liquid)
                Sⱼ[2] ~ 0.0
                scalarize(zⱼᵢ[3, :] .~ zₜ)...
                scalarize(zⱼᵢ[2, :] .~ 0.0)...
                scalarize(zⱼᵢ[1, :] .- zⱼᵢ[3, :] .~ 0.0)...
                scalarize(zᵂⱼᵢ[2, :] .~ 0.0)... # Mass base 
                scalarize(zᵂⱼᵢ[1, :] .- zᵂⱼᵢ[3, :] .~ 0.0)... # Mass base
                scalarize(Fᵂⱼᵢ[3, :] .~ Fᵂⱼ[3].*zᵂⱼᵢ[3, :])... # Mass base
                ρ[1] ~ ρ[3]
                ρ[3] ~ molar_density(model, P_user, T_user, zₜ_user, phase = :liquid)
                ρ[2] ~ 0.0
                ρʷ[1] ~ ρʷ[3]
                ρʷ[3] ~ mass_density(model, P_user, T_user, zₜ_user, phase = :liquid)
                ρʷ[2] ~ 0.0
                MWⱼ[1] ~ MWⱼ[3]
                MWⱼ[3] ~ sum(MWs.*zⱼᵢ[3, :]*gramsToKilograms)
                MWⱼ[2] ~ 0.0

            ]
        else
            #flash calculation of the two phases (assuming liquid and gas only)
            xᵢⱼ, nᵢⱼ, G = tp_flash(model, P_user, T_user, zₜ_user, DETPFlash(; equilibrium = :vle)) #Phase i, component j
            H_l = enthalpy(model, P_user, T_user, xᵢⱼ[1, :], phase = :liquid) 
            H_g = enthalpy(model, P_user, T_user, xᵢⱼ[2, :], phase = :vapor)
            S_l = entropy(model, P_user, T_user, xᵢⱼ[1, :], phase = :liquid) 
            S_g = entropy(model, P_user, T_user, xᵢⱼ[2, :], phase = :vapor)
            V_l = volume(model, P_user, T_user, nᵢⱼ[1, :], phase = :liquid)
            V_g = volume(model, P_user, T_user, nᵢⱼ[2, :], phase = :vapor)

            pc = [α_g ~ sum(nᵢⱼ[2, :])/(sum(nᵢⱼ[1, :]) + sum(nᵢⱼ[2, :])) # Vapor phase is the second entry
                Tc ~ Tci
                Pc ~ Pci
                P_buble ~ Pbuble
                P_dew ~ Pdew
                Hⱼ[2] ~ H_g + sum(xᵢⱼ[2, :].*ΔH₀f)
                Hⱼ[3] ~ H_l + sum(xᵢⱼ[1, :].*ΔH₀f)
                Hⱼ[1] ~ sum(nᵢⱼ[2, :])*Hⱼ[2] + sum(nᵢⱼ[1, :])*Hⱼ[3]
                Sⱼ[2] ~ S_g
                Sⱼ[3] ~ S_l
                Sⱼ[1] ~ sum(nᵢⱼ[2, :])*Sⱼ[2] + sum(nᵢⱼ[1, :])*Sⱼ[3]
                scalarize(zⱼᵢ[1, :] .~ zₜ)... #Global phase
                scalarize(zⱼᵢ[3, :] .~ xᵢⱼ[1, :])... #Liquid phase
                scalarize(zⱼᵢ[2, :] .~ xᵢⱼ[2, :])... #Vapor phase
                scalarize(Fᵂⱼᵢ[1, :] .~ Fᵂⱼ[1].*zᵂⱼᵢ[1, :])... # Mass base
                scalarize(Fᵂⱼᵢ[2, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[2, :])... # Mass base
                scalarize(Fᵂⱼᵢ[3, :] .~ Fᵂⱼ[3].*zᵂⱼᵢ[3, :])... # Mass base
                ρ[1] ~ (sum(nᵢⱼ[2, :]) + sum(nᵢⱼ[1, :]))/(V_l + V_g)
                ρ[2] ~ sum(nᵢⱼ[2, :])/V_g
                ρ[3] ~ sum(nᵢⱼ[1, :])/V_l
                ρʷ[1] ~ (sum(nᵢⱼ[2, :].*MWs)*gramsToKilograms + sum(nᵢⱼ[1, :].*MWs)*gramsToKilograms)/(V_l + V_g)
                ρʷ[2] ~ sum(nᵢⱼ[2, :].*MWs)*gramsToKilograms/(V_g)
                ρʷ[3] ~ sum(nᵢⱼ[1, :].*MWs)*gramsToKilograms/(V_l)
                MWⱼ[1] ~ sum(MWs.*zⱼᵢ[1, :]*gramsToKilograms)
                MWⱼ[2] ~ sum(MWs.*zⱼᵢ[2, :]*gramsToKilograms)
                MWⱼ[3] ~ sum(MWs.*zⱼᵢ[3, :]*gramsToKilograms)

            ]
        end
        
    end

    #phase 1 is total, 2 is vapor, 3 is liquid
    eqs = [eqs_conn..., global_mol_balance..., molar_to_mass..., component_balance..., pc...]
    
    unfold_pars = []
    for par in pars
        unfold_pars = [unfold_pars...; par...]
    end

    unfold_vars = []
    for var in vars
        unfold_vars = [unfold_vars...; var...]
    end

    ODESystem([eqs...;], t, unfold_vars, unfold_pars; name, systems) 
end

