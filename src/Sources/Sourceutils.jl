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

    unfold_vars = []
    for var in vars
        unfold_vars = [unfold_vars...; var...]
    end

    ODESystem(Equation[], t, unfold_vars, []; name)

end

@component function Display(; Nc, name)
    
    vars = @variables begin
    P(t) # Pressure (Pa) 
    T(t) # Temperature (K)
    F(t) # Molar Flow rate (mol/s)
    Fʷ(t) # Mass Flow rate (kg/s)
    H(t) # Enthalpy (J/mol)
    S(t) # Entropy (J/mol.K)
    (z₁(t))[1:Nc] # component molar fraction global (mol/mol)  
    (z₂(t))[1:Nc] # component molar fraction in vapor phase (mol/mol)
    (z₃(t))[1:Nc] # component molar fraction in liquid phase (mol/mol) 
    α_g(t) # gas phase fraction (mol/mol) 
    ρ(t) # Molar density (mol/m³)
    ρʷ(t) # Mass density (kg/m³)
    (MW(t))[1:3] # Molar mass (g/mol)
    end

    systems = @named begin 
        InPort = matcon(; Nc = Nc)
    end

    eqs_conn = [     
        #InPort data
        InPort.P + P ~ 0.0 
        InPort.T + T ~ 0.0  
        InPort.F + F ~ 0.0 
        InPort.Fʷ + Fʷ ~ 0.0
        InPort.H + H ~ 0.0
        InPort.S + S ~ 0.0
        scalarize(InPort.z₁ .+ z₁ .~ 0.0)...
        scalarize(InPort.z₂ .+ z₂ .~ 0.0)...
        scalarize(InPort.z₃ .+ z₃ .~ 0.0)...
        InPort.α_g + α_g ~ 0.0 
        InPort.ρ + ρ ~ 0.0
        InPort.ρʷ + ρʷ ~ 0.0
        scalarize(InPort.MW + MW .~ 0.0)...
    ]

    unfold_vars = []
    for var in vars
        unfold_vars = [unfold_vars...; var...]
    end

    ODESystem(eqs_conn, t, unfold_vars, []; name, systems)

end



function ModelingToolkit.connect(::Type{matcon}, con...)
        eqs = []
        for i in length(con) - 1
            push!(eqs, con[i].P + con[i+1].P ~ 0.0)
            push!(eqs, con[i].T + con[i+1].T ~ 0.0)  
            push!(eqs, con[i].F + con[i+1].F ~ 0.0) 
            push!(eqs, con[i].Fʷ + con[i+1].Fʷ ~ 0.0)
            push!(eqs, con[i].H + con[i+1].H ~ 0.0) 
            push!(eqs, con[i].S + con[i+1].S ~ 0.0)
            push!(eqs, scalarize(con[i].z₁ .+ con[i+1].z₁ .~ 0.0)...)
            push!(eqs, scalarize(con[i].z₂ .+ con[i+1].z₂ .~ 0.0)...)
            push!(eqs, scalarize(con[i].z₃ .+ con[i+1].z₃ .~ 0.0)...) 
            push!(eqs, con[i].α_g + con[i+1].α_g ~ 0.0)
            push!(eqs, con[i].ρ + con[i+1].ρ ~ 0.0)
            push!(eqs, con[i].ρʷ + con[i+1].ρʷ ~ 0.0)
            push!(eqs, scalarize(con[i].MW .+ con[i+1].MW .~ 0.0)...)
        end
        return eqs
end
