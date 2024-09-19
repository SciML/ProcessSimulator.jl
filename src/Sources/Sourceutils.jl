@connector function matcon(; Nc, name)

    vars = @variables begin
    P(t), [description = "Pressure (Pa)", output = true] 
    T(t), [description = "Temperature (K)", output = true] 
    F(t), [description = "Molar Flow rate (mol/s)", output = true]  
    Fʷ(t), [description = "Mass Flow rate (kg/s)", output = true] 
    H(t), [description = "Enthalpy (J/mol)", output = true]  
    S(t), [description = "Entropy (J/mol.K)", output = true] # 
    (z₁(t))[1:Nc], [description = "component molar fraction global (mol/mol)", output = true]    
    (z₂(t))[1:Nc], [description = "component molar fraction in vapor phase (mol/mol)", output = true] # 
    (z₃(t))[1:Nc], [description = "component molar fraction in liquid phase (mol/mol)", output = true] #  
    α_g(t), [description = "gas phase fraction (mol/mol)", output = true] # 
    ρ(t), [description = "Molar density (mol/m³)", output = true] # 
    ρʷ(t), [description = "Mass density (kg/m³)", output = true] # 
    (MW(t))[1:3], [description = "Molar mass (g/mol)", output = true] # 
    end

    ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)

end



@connector function thermal_energy_connector(; name)

        vars = @variables begin
        ϕᴱ(t), [description = "Energy flux at the interface (W/m²)", output = true]   
        T(t), [description = "Interface temperature (T)", output = true] 
        A(t), [description = "Actual heat transfer area", output = true]  
        end

    ODESystem(Equation[], t, collect(Iterators.flatten(vars)), []; name)

end

