@component function Jacket(;
        substances_user,
        Nc = length(substances_user),
        phase,
        thermal_fluid_model,
        heat_transfer_coef,
        name
    )
    properties = Dict(subs => load_component_properties(subs) for subs in substances_user)
    MWs = [properties[subs]["MW"] for subs in substances_user]

    pars = @parameters begin
        U = heat_transfer_coef
    end

    systems = @named begin
        Out = matcon(; Nc = Nc)
        In = matcon(; Nc = Nc)
        EnergyCon = thermal_energy_connector()
    end

    vars = @variables begin
        LMTD(t), [description = "log mean temperature difference (-)"]
        Tⱼ(t), [description = "Thermal fluid exiting temperature (K)"]
        H(t), [description = "Thermal fluid exiting enthalpy (J/mol)"]
        S(t), [description = "Thermal fluid exiting entropy (J/K/mol)"]
        ρ(t), [description = "Thermal fluid exiting molar density (mol/m³)"]
        ρʷ(t), [description = "Thermal fluid exiting mass density (kg/m³)"]
        P_out(t), [description = "Thermal fluid exiting pressure (Pa)"]
        F_out(t), [description = "Thermal fluid exiting molar flow rate (mol/s)"]
        Fʷ_out(t), [description = "Thermal fluid exiting mass flow rate (mol/s)"]
        Q_out(t), [description = "Thermal fluid exiting volumetric flow rate (m³/s)"]
        Q̇(t), [description = "Heat transfer rate from/to the jacket"]
    end

    enthalpy_eq = [H ~ enthalpy(thermal_fluid_model, P_out, Tⱼ, Out.z₁, phase = "unknown")]

    entropy_eq = [S ~ 0.0]

    densities_eq = [
        ρ ~ molar_density(
            thermal_fluid_model, P_out, Tⱼ, Out.z₁, phase = "unknown"
        )
        ρʷ ~
            mass_density(thermal_fluid_model, P_out, Tⱼ, Out.z₁, phase = "unknown")
    ]

    pressure_drop = [P_out ~ In.P]

    mass_balance = [
        In.Fʷ - Fʷ_out ~ 0.0
        Q_out ~ Fʷ_out / ρʷ
        Q_out ~ F_out / ρ
    ]

    lmtd = [LMTD ~ log((EnergyCon.T - In.T) / (EnergyCon.T - Tⱼ))]

    heat_flux = [Q̇ ~ U * EnergyCon.A * (In.T - Tⱼ) / LMTD]

    energy_balance = [In.H * In.F - H * F_out - Q̇ ~ 0.0] # Quasi steady state assumption (Temperature change is much faster than reactor change)

    #Outlet connector equations:
    out_conn = [
        Out.P ~ P_out
        Out.T ~ Tⱼ
        Out.F ~ F_out
        Out.Fʷ ~ Fʷ_out
        Out.H ~ H
        Out.S ~ S
        Out.ρʷ ~ ρʷ
        Out.ρ ~ ρ
        scalarize(Out.z₁ .~ In.z₁)...
        Out.MW[1] ~ MWs
        EnergyCon.ϕᴱ ~ Q̇
    ]

    if phase == :liquid
        out_conn_phases = [
            scalarize(Out.z₂ .~ 0.0)...
            scalarize(Out.z₃ .~ In.z₁)...
            Out.MW[2] ~ 0.0
            Out.MW[3] ~ MWs
            Out.α_g ~ 0.0
        ]

    elseif phase == :vapor
        out_conn_phases = [
            scalarize(Out.z₂ .~ In.z₁)...
            scalarize(Out.z₃ .~ 0.0)...
            Out.MW[2] ~ MWs
            Out.MW[3] ~ 0.0
            Out.α_g ~ 1.0
        ]
    end

    eqs = [
        enthalpy_eq...; entropy_eq...; densities_eq...; pressure_drop...;
        mass_balance...; lmtd...; heat_flux...; energy_balance...; out_conn...;
        out_conn_phases...
    ]

    ODESystem(
        [eqs...;], t, collect(Iterators.flatten(vars)),
        collect(Iterators.flatten(pars)); name, systems = [Out, In, EnergyCon]
    )
end

export Jacket
