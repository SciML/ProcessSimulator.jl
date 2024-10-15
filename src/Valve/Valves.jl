@component function LinearValve(; Nc, CV, θ, model, name, phase)

    pars = @parameters begin
        Cᵥ = CV, [description = "Valve coefficient"]
    end

    vars = @variables begin
        #(xᵢ(t))[1:Nc], [description = "Molar fraction in liquid phase (-)"]
        #(yᵢ(t))[1:Nc], [description = "Molar fraction in gas phase (-)"]
        #(ϕᵢᴸ(t))[1:Nc], [description = "Fugacity coefficient in liquid phase"]
        #(ϕᵢᵍ(t))[1:Nc], [description = "Fugacity coefficient in gas phase"]
        ρᴸ(t), [description = "Molar Density of gas phase (mol/m³)"]
        ρᵍ(t), [description = "Molar Density of liquid phase (mol/m³)"]
        θ(t), [description = "Valve opening"]
        T(t), [description = "Exit temperature (K)"]  
        P(t), [description = "Exit pressure (Pa)"]
    
        hᴸ_outflow(t), [description = "Outflow Enthalpy of liquid phase (J/mol)"]
        hᵍ_outflow(t), [description = "Outflow Enthalpy of liquid phase (J/mol)"]
        Fᴸ_outflow(t), [description = "Outlet molar flow rate of liquid phase (mol/s)"] 
        Fᵍ_outflow(t), [description = "Outlet molar flow rate of gas phase (mol/s)"] 
        Q̇(t), [description = "Heat transfer rate (J/s)"]
    
    
        #_0_Nᴸ(t), [description = "Condition for all gas phase"]
        #_0_Nᵍ(t), [description = "Condition for all liquid phase"]
    end


    Ports = @named begin
        In = matcon(; Nc = Nc)
        Out = matcon(; Nc = Nc)
        #EnergyCon = thermal_energy_connector()
    end

    connector_eqs = [Out.P ~ P
                    Out.T ~ T
                    scalarize(Out.z₁ .~ In.z₁)...
                    scalarize(Out.z₂ .~ In.z₂)...
                    scalarize(Out.z₃ .~ In.z₃)...
                    Out.α_g ~ In.α_g]

    heat = [Q̇ ~ 0.0]

    if phase == :liquid
        SS_mass_balance = [Fᴸ_outflow ~ ρᴸ*Cᵥ*θ*max(0.0, (In.P - P)/√(abs(In.P - P)) + eps(1.))
        0.0 ~ In.F*In.H - hᴸ_outflow*Fᴸ_outflow + Q̇
        Fᴸ_outflow ~ In.F
        Fᵍ_outflow ~ 0.0
        Out.F ~ Fᴸ_outflow
        Out.H ~ hᴸ_outflow]

    elseif phase == :vapor  
        SS_mass_balance = [Fᵍ_outflow ~ ρᵍ*Cᵥ*θ*max(0.0, (In.P - P)/√(abs(In.P - P)) + eps(1.)) 
        0.0 ~ In.F*In.H - hᵍ_outflow*Fᵍ_outflow + Q̇
        Fᵍ_outflow ~ In.F
        Fᴸ_outflow ~ 0.0
        Out.F ~ Fᵍ_outflow
        Out.H ~ hᵍ_outflow]
    end

    density = [ρᴸ ~ molar_density(model, P, T, In.z₃, phase = "liquid")
    ρᵍ ~ molar_density(model, P, T, In.z₂, phase = "vapor")]

    enthalpy = [hᴸ_outflow ~ molar_enthalpy(model, P, T, In.z₃, phase = "liquid")
    hᵍ_outflow ~ molar_enthalpy(model, P, T, In.z₂, phase = "vapor")]

    eqs = [SS_mass_balance...; density...; enthalpy...; connector_eqs...; heat...]

    ODESystem([eqs...;], t, collect(Iterators.flatten(vars)), collect(Iterators.flatten(pars)); name, systems = [Ports...])

end

export LinearValve