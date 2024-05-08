@parameters t
@component function KineticReactionNetwork(;substances_user = ["methane", "carbon monoxide"], Coef_Cr::Array = [-1.0 1.0], Af_r::Array, Ef_r::Array,
     Do_r::Array, Nc::Int = size(substances_user, 1), Nri::Int = size(Coef_Cr, 1), name)
    
    pars = @parameters begin
        Nr = Nri
        Coef_Cr[1:Nri, 1:Nc] = Coef_Cr, [description = "Stoichiometric coefficients of each component in each reaction (-)"]
        Af_r[1:Nri] = Af_r, [description = "Arrhenius constant of each reaction at given temperature ()"]
        Ef_r[1:Nri] = Ef_r, [description = "Activation energy of each reaction at given temperature ()"]
        Do_r[1:Nri, 1:Nc] = Do_r, [description = "Forward order of the components (-) "]
    end 


    vars = []

    unfold_pars = []
    for par in pars
        unfold_pars = [unfold_pars...; par...]
    end

    
    unfold_vars = []
    for var in vars
        unfold_vars = [unfold_vars...; var...]
    end

    ODESystem(Equation[], t, unfold_vars, unfold_pars; name)
end


a = KineticReactionNetwork(Af_r = [1.0], Ef_r = [1.0], Do_r = [1.0 1.0], name = :Reaction)
parameters(a)[1]

