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
    H(t) # Enthalpy (J/mol)
    S(t) # Entropy (J/mol.K)
    z₁(t)[1:Nc] # component molar fraction global (mol/mol)  
    z₂(t)[1:Nc] # component molar fraction in vapor phase (mol/mol)
    z₃(t)[1:Nc] # component molar fraction in liquid phase (mol/mol) 
    α_g(t)
    end

    ODESystem(Equation[], t, [vars...;], pars; name)

end

matcon(Nc = 2, name = "test")



function ModelingToolkit.connect(::typeof(matcon), con...)
        eqs = []
        for i in length(con) - 1
            push!(eqs, con.P[i] ~ con[i+1].P)
            push!(eqs, con.T[i] ~ con[i+1].T)  
            push!(eqs, con.F[i] ~ -con[i+1].F) #Should it be mass or molar flow rate?
            push!(eqs, con.H[i] ~ con[i+1].H)
            push!(eqs, con.S[i] ~ con[i+1].S)
            push!(eqs, con.z[i] ~ con[i+1].z) 
        end
end


"""
    load_component_properties(filepath::String)

Load component properties from a file.

# Arguments
- `filepath::String`: Path to the file containing the component properties.

# Returns
- `properties::Dict`: A dictionary containing the component properties.
"""
function load_component_properties(component_name::String)
    file_path = abspath(joinpath(@__DIR__, "OMChemSimDataBase/$(component_name).json"))
    println(file_path)
    if isfile(file_path)
        return JSON.parsefile(file_path)
    else
        error("Component file does not exist: $file_path")
    end
end