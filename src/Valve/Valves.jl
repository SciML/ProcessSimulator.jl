@component function LinearValve(; Nc, CV, model, name)

    pars = @parameters begin
        Cᵥ = CV, [description = "Valve coefficient"]
    end

    vars = @variables begin
        θ(t), [description = "Valve opening (-)", guess = 1.0, irreducible = true]
        T(t), [description = "Exit temperature (K)"]  
        P(t), [description = "Exit pressure (Pa)", guess = 101325.0]
        F_outflow(t), [description = "Outlet molar flow rate of liquid phase (mol/s)"] 
        Q̇(t), [description = "Heat transfer rate (J/s)"]
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
                    Out.α_g ~ In.α_g
                    P ~ 101325.0]

    heat = [Q̇ ~ 0.0]


    SS_mass_balance = [F_outflow ~ Cᵥ*θ*(In.P - P)/√(abs(In.P - P) + eps(10.))
    D(θ) ~ 0.0
    In.T ~ Out.T
    F_outflow ~ In.F
    Out.F ~ F_outflow
    Out.H ~ In.H]

    eqs = [SS_mass_balance...; connector_eqs...; heat...]

    ODESystem([eqs...;], t, collect(Iterators.flatten(vars)), collect(Iterators.flatten(pars)); name, systems = [Ports...])

end

export LinearValve