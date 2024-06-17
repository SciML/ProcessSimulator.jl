@connector function matcon(; Nc, name)

    vars = @variables begin
    P(t), [connect = Flow] # Pressure (Pa) 
    T(t), [connect = Flow] # Temperature (K)
    F(t), [connect = Flow] # Molar Flow rate (mol/s)
    Fʷ(t), [connect = Flow] # Mass Flow rate (kg/s)
    H(t), [connect = Flow] # Enthalpy (J/mol)
    S(t), [connect = Flow] # Entropy (J/mol.K)
    (z₁(t))[1:Nc], [connect = Flow] # component molar fraction global (mol/mol)  
    (z₂(t))[1:Nc], [connect = Flow] # component molar fraction in vapor phase (mol/mol)
    (z₃(t))[1:Nc], [connect = Flow] # component molar fraction in liquid phase (mol/mol) 
    α_g(t), [connect = Flow] # gas phase fraction (mol/mol) 
    ρ(t), [connect = Flow] # Molar density (mol/m³)
    ρʷ(t), [connect = Flow] # Mass density (kg/m³)
    (MW(t))[1:3], [connect = Flow] # Molar mass (g/mol)
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
        InPort.P ~ + P 
        InPort.T ~ + T  
        InPort.F ~ + F 
        InPort.Fʷ ~ + Fʷ
        InPort.H ~ + H
        InPort.S ~ + S
        scalarize(InPort.z₁ .- z₁ .~ 0.0)...
        scalarize(InPort.z₂ .- z₂ .~ 0.0)...
        scalarize(InPort.z₃ .- z₃ .~ 0.0)...
        InPort.α_g ~ + α_g 
        InPort.ρ ~ + ρ
        InPort.ρʷ ~ + ρʷ
        scalarize(InPort.MW .~ + MW)...
    ]

    unfold_vars = []
    for var in vars
        unfold_vars = [unfold_vars...; var...]
    end

    ODESystem(eqs_conn, t, unfold_vars, []; name, systems)

end









#= function ModelingToolkit.connect(::Type{matcon}, con...)
        eqs = []
        for i in length(con) - 1
            push!(eqs, con[i].P ~ con[i+1].P )
            push!(eqs, con[i].T ~ con[i+1].T )  
            push!(eqs, con[i].F ~ -con[i+1].F ) 
            push!(eqs, con[i].Fʷ ~ -con[i+1].Fʷ )
            push!(eqs, con[i].H ~ con[i+1].H ) 
            push!(eqs, con[i].S ~ con[i+1].S )
            push!(eqs, scalarize(con[i].z₁ .~ con[i+1].z₁)...)
            push!(eqs, scalarize(con[i].z₂ .~ con[i+1].z₂)...)
            push!(eqs, scalarize(con[i].z₃ .~ con[i+1].z₃)...) 
            push!(eqs, con[i].α_g ~ con[i+1].α_g )
            push!(eqs, con[i].ρ ~ con[i+1].ρ )
            push!(eqs, con[i].ρʷ ~ con[i+1].ρʷ )
            push!(eqs, scalarize(con[i].MW .~ con[i+1].MW)...)
        end
        return eqs
end =#
