@parameters t
D = Differential(t)

@connector function matcon(; Nc = 2, name)
    
    pars = @parameters begin
    N = Nc
    end

    vars = @variables begin
    P(t) #Pressure (Pa) 
    T(t) #Temperature (K)
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

    ODESystem(Equation[], t, [vars...;], pars; name)

end

function ModelingToolkit.connect(::typeof(matcon), con...)
        eqs = []
        for i in length(con) - 1
            push!(eqs, con[i].P ~ con[i+1].P)
            push!(eqs, con[i].T ~ con[i+1].T)  
            push!(eqs, con[i].F ~ -con[i+1].F) 
            push!(eqs, con[i].Fʷ ~ -con[i+1].Fʷ)
            push!(eqs, con[i].H ~ con[i+1].H)
            push!(eqs, con[i].S ~ con[i+1].S)
            push!(eqs, con[i].z ~ con[i+1].z₁)
            push!(eqs, con[i].z ~ con[i+1].z₂)
            push!(eqs, con[i].z ~ con[i+1].z₃) 
            push!(eqs, con[i].α_g ~ con[i+1].α_g)
            push!(eqs, con[i].ρ ~ con[i+1].ρ)
            push!(eqs, con[i].ρʷ ~ con[i+1].ρʷ)
            push!(eqs, con[i].MW ~ con[i+1].MW)
        end
end
