@component function MaterialSource(; substances_user = ["methane", "carbon monoxide"], model = PR(substances_user), properties = Dict(subs => load_component_properties(subs) for subs in substances_user),
    P_user = 101325.0, T_user = 25.0 + 273.15, Fₜ = 100.0, zₜ = [0.5, 0.5], name)
    #phase 1 is total, 2 is vapor, 3 is liquid

    
    MWs = Dict(subs => properties[subs]["MW"] for subs in substances_user) # TODO: Also read parameters from the properties dictionary. So far it's comming from clapeyron.jl
    Nc = size(substances_user, 1) #Numeric number of components
    
    pars = @parameters begin
        N = size(substances_user, 1), [description = "Number of components", tunable = false]
        substances[1:Nc] = substances_user, [description = "Component names", tunable = false]
        P = P_user, [description = "Pressure (Pa)", tunable = false]
        T = T_user, [description = "Temperature (K)",tunable = false]
        MWⱼ[1:Nc] = [MWs[subs] for subs in substances_user], [description = "Molar mass of each phase j (kg/mol)"] # 14
    end 

    vars = @variables begin
        P_buble(t), [description = "Bubble point pressure (Pa)"]
        P_dew(t), [description = "Dew point pressure (Pa)"]
        α_g(t), [description = "Vapor phase molar fraction"]
        α_l(t), [description = "Liquid phase molar fraction"]
        αᵂ_g(t), [description = "Vapor phase mass fraction"]
        αᵂ_l(t), [description = "Liquid phase mass fraction"]
        (Fⱼ(t))[1:3], [description = "Molar flow rate in each phase j (kmol/s)"]
        (Fᵂⱼ(t))[1:3], [description = "Mass flow rate in each phase j (kg/s)"]
        (zⱼᵢ(t))[1:3, 1:Nc], [description = "Component molar fraction in each phase (-)"] # 6
        (zᵂⱼᵢ(t))[1:3, 1:Nc], [description = "Component mass fraction in each phase (-)"] # 6
        (Fⱼᵢ(t))[1:3, 1:Nc], [description = "Molar flow rate in each phase j and component i (kmol/s)"] # 6
        (Fᵂⱼᵢ(t))[1:3, 1:Nc], [description = "Mass flow rate in each phase j and component i (kmol/s)"] #6
        #Cpⱼ[1:3], [description = "Heat capacity in each phase j at T and P (J/mol.K)"]
        (Hⱼ(t))[1:3], [description = "Enthalpy in each phase j at T and P (J/mol)"] #3
        (Sⱼ(t))[1:3], [description = "Entropy in each phase j at T and P (J/mol.K)"] #3

    end

    systems = @named begin 
        Out = matcon(; Nc = Nc)
    end

    #Connector equations
    eqs_conn = [
        Out.P ~ P_user 
        Out.T ~ T_user  
        Out.F ~ - Fₜ # F is negative as it is leaving the component 
        Out.H ~ Hⱼ[1] 
        Out.S ~ Sⱼ[1] 
        scalarize(Out.z₁ .- zⱼᵢ[1, :] .~ 0.0)...
        scalarize(Out.z₂ .- zⱼᵢ[2, :] .~ 0.0)...
        scalarize(Out.z₃ .- zⱼᵢ[3, :] .~ 0.0)...
        Out.α_g ~ α_g 
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

    molar_to_mass_2 = [scalarize(Fᵂⱼᵢ[:, i] .~ MWⱼ[i]*Fⱼᵢ[:, i] ) for i in 1:Nc]
    
    molar_to_mass = [molar_to_mass...; molar_to_mass_2...]
    
    #Component Molar and Mass Balance
    component_balance = [scalarize(Fⱼᵢ[:, i] .~ Fⱼ[:].*zⱼᵢ[:, i]) for i in 1:Nc]


    #Phase check
    Tc, Pc, Vc = crit_mix(model, zₜ) #Critical point of the mixture
    if T_user > Tc
        pc = [α_g ~ 1.0
            P_buble ~ 0.0
            P_dew ~ 0.0
            Hⱼ[1] ~ Hⱼ[2]
            Hⱼ[3] ~ 0.0
            Hⱼ[2] ~ enthalpy(model, P_user, T_user, zₜ, phase = :vapor)
            Sⱼ[1] ~ Sⱼ[2]
            Sⱼ[3] ~ 0.0
            Sⱼ[2] ~ entropy(model, P_user, T_user, zₜ, phase = :vapor)
            scalarize(zⱼᵢ[2, :] .~ zₜ)...
            scalarize(zⱼᵢ[3, :] .~ 0.0)...
            scalarize(zⱼᵢ[1, :] .~ zₜ)...

            scalarize(zᵂⱼᵢ[3, :] .~ 0.0)... # Mass base 
            scalarize(zᵂⱼᵢ[1, :] .~ zᵂⱼᵢ[2, :])... # Mass base
            scalarize(Fᵂⱼᵢ[2, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[2, :])... # Mass base
              ]
    else

        Pdew, Vᵍdew, Vˡdew, xdew = dew_pressure(model, T_user, zₜ)
        Pbuble, Vᵍbuble, Vˡbuble, xbuble = bubble_pressure(model, T_user, zₜ)

        if P_user < Pdew
            pc = [α_g ~ 1.0
                P_buble ~ Pbuble
                P_dew ~ Pdew
                Hⱼ[1] ~ Hⱼ[2]
                Hⱼ[2] ~ enthalpy(model, P_user, T_user, zₜ, phase = :vapor)
                Hⱼ[3] ~ 0.0
                Sⱼ[1] ~ Sⱼ[2]
                Sⱼ[2] ~ entropy(model, P_user, T_user, zₜ, phase = :vapor)
                Sⱼ[3] ~ 0.0
                scalarize(zⱼᵢ[2, :] .~ zₜ)...
                scalarize(zⱼᵢ[3, :] .~ 0.0)...
                scalarize(zⱼᵢ[1, :] .~ zₜ)...
                scalarize(zᵂⱼᵢ[3, :] .~ 0.0)... # Mass base 
                scalarize(zᵂⱼᵢ[1, :] .~ zᵂⱼᵢ[2, :])... # Mass base
                scalarize(Fᵂⱼᵢ[2, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[2, :])... # Mass base
            ] 
        elseif P_user > Pbuble
            pc = [α_g  ~ 0.0
                P_buble ~ Pbuble
                P_dew ~ Pdew
                Hⱼ[1] ~ Hⱼ[3]
                Hⱼ[3] ~ enthalpy(model, P_user, T_user, zₜ, phase = :liquid)
                Hⱼ[2] ~ 0.0
                Sⱼ[1] ~ Sⱼ[3]
                Sⱼ[3] ~ entropy(model, P_user, T_user, zₜ, phase = :liquid)
                Sⱼ[2] ~ 0.0
                scalarize(zⱼᵢ[3, :] .~ zₜ)...
                scalarize(zⱼᵢ[2, :] .~ 0.0)...
                scalarize(zⱼᵢ[1, :] .- zⱼᵢ[3, :] .~ 0.0)...
                scalarize(zᵂⱼᵢ[2, :] .~ 0.0)... # Mass base 
                scalarize(zᵂⱼᵢ[1, :] .- zᵂⱼᵢ[3, :] .~ 0.0)... # Mass base
                scalarize(Fᵂⱼᵢ[3, :] .~ Fᵂⱼ[3].*zᵂⱼᵢ[3, :])... # Mass base

            ]
        else
            #flash calculation of the two phases (assuming liquid and gas only)
            xᵢⱼ, nᵢⱼ, G = tp_flash(model, P_user, T_user, zₜ, RRTPFlash(; equilibrium = :vle)) #Phase i, component j
            H_l = enthalpy(model, P_user, T_user, xᵢⱼ[1, :], phase = :liquid) 
            H_g = enthalpy(model, P_user, T_user, xᵢⱼ[2, :], phase = :vapor)
            S_l = entropy(model, P_user, T_user, xᵢⱼ[1, :], phase = :liquid) 
            S_g = entropy(model, P_user, T_user, xᵢⱼ[2, :], phase = :vapor)

            pc = [α_g ~ sum(nᵢⱼ[2, :])/(sum(nᵢⱼ[1, :]) + sum(nᵢⱼ[2, :])) # Vapor phase is the second entry
                P_buble ~ Pbuble
                P_dew ~ Pdew
                Hⱼ[2] ~ H_g
                Hⱼ[3] ~ H_l
                Hⱼ[1] ~ sum(nᵢⱼ[2, :])*Hⱼ[2] + sum(nᵢⱼ[1, :])*Hⱼ[3]
                Sⱼ[2] ~ S_g
                Sⱼ[3] ~ S_l
                Sⱼ[1] ~ sum(nᵢⱼ[2, :])*Sⱼ[2] + sum(nᵢⱼ[1, :])*Sⱼ[3]
                scalarize(zⱼᵢ[1, :] .~ zₜ)... #Global phase
                scalarize(zⱼᵢ[3, :] .~ xᵢⱼ[1, :])... #Liquid phase
                scalarize(zⱼᵢ[2, :] .~ xᵢⱼ[2, :])... #Vapor phase
                scalarize(Fᵂⱼᵢ[1, :] .~ Fᵂⱼ[1].*zᵂⱼᵢ[1, :])... # Mass base
                scalarize(Fᵂⱼᵢ[2, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[2, :])... # Mass base
                scalarize(Fᵂⱼᵢ[3, :] .~ Fᵂⱼ[2].*zᵂⱼᵢ[3, :])... # Mass base
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

