@component function CSTR(; substances_user,
        Nc = length(substances_user),
        phase,
        model,
        Reaction,
        Ac,
        height_out_port,
        name,
        guesses
)

    #Constants
    gramsToKilograms = 10^(-3)
    Rᵍ = 8.314 # J/(mol K)
    Nri = Reaction.Nri

    #Properties of individual substances
    properties = Dict(subs => load_component_properties(subs) for subs in substances_user)
    MWs = [properties[subs]["MW"] for subs in substances_user]
    ΔH₀f = [properties[subs]["IGHF"] / 10^3 for subs in substances_user] # (IG formation enthalpy) J/mol

    pars = @parameters begin
        Af_r[1:Nri] = Reaction.Af_r,
        [description = "Arrhenius constant of each reaction at given temperature ()"]
        Coef_Cr[1:Nri, 1:Nc] = Reaction.Coef_Cr,
        [description = "Stoichiometric coefficients of each component in each reaction (-)"]
        Do_r[1:Nri, 1:Nc] = Reaction.Do_r,
        [description = "Forward order of the components (-) "]
        Ef_r[1:Nri] = Reaction.Ef_r,
        [description = "Activation energy of each reaction at given temperature ()"]
        A = Ac, [description = "Cross sectional area of the tank (m²)"]
    end

    Ports = @named begin
        In = matcon(; Nc = Nc)
        Out = matcon(; Nc = Nc)
        #EnergyCon = thermal_energy_connector()
    end

    vars = @variables begin
        M(t), [description = "Mass holdup in the tank (kg)"]
        N(t), [description = "Total molar holdup in the tank (mol)"]
        V(t), [description = "Volume holdup in the tank (m³)", guess = guesses[:V]]
        (Nᵢ(t))[1:Nc], [description = "Molar holdup of each component in the tank (mol)"]
        (Cᵢ(t))[1:Nc],
        [description = "Concentration of each component in the tank (mol/m³)",
            guess = guesses[:Cᵢ]]

        ρ(t), [description = "Molar Density of the fluid in the tank (mol/m³)"]
        ρʷ(t), [description = "Mass Density of the fluid in the tank (kg/m³)"]
        MW(t), [description = "Molecular weight of fluid in the tank (kg/kmol)"]
        T(t), [description = "Temperature of vessel contents (K)"]
        P_atm(t), [description = "Tank pressure (Pa)"] # Equal to inlet pressures.
        H(t), [description = "Enthalpy holdup of the fluid in the tank (J)"]

        F_out(t), [description = "Outlet molar flow rate (mol/s)", guess = guesses[:F_out]]
        Fʷ_out(t),
        [description = "Outlet molar flow rate (mol/s)", guess = guesses[:Fʷ_out]]
        Q_out(t), [description = "Outlet volumetric flow rate (m³/s)"] # DoF

        Q̇(t), [description = "Heat transfer rate (J/s)"] # Potential DoF
        height(t),
        [description = "Liquid level in vessel measured from bottom of the tank (m)"]
        P_out(t),
        [description = "Pressure at the outlet stream level (Pa)", guess = 101325.0]

        (r(t))[1:Nri], [description = "Rate of each reaction for each component (mol/s/m³)"]
        (R(t))[1:Nc], [description = "Overall reaction rate (mol/s/m³)"]

        (F_in(t)), [description = "Inlet molar flow rate (mol/s)"] # DoF through inlet stream
        (T_in(t)), [description = "Inlet temperature (K)"] # DoF through inlet stream
        (h_in(t)), [description = "Inlet specific enthalpy (J/mol)"]

        fᵢ(t), [description = "fugacity"]
    end

    #Reaction equations
    reaction_rate = [r[i] ~ Af_r[i] * exp(-Ef_r[i] / (Rᵍ * T)) *
                            prod(scalarize((Cᵢ[:] .^ Do_r[i, :]))) for i in 1:Nri] # If there's an inert, the order is just zero, but it has to be written
    overall_reaction_rate = [R[i] ~ sum(scalarize(r[:] .* Coef_Cr[:, i])) for i in 1:Nc]  # If there's an inert, the coefficient is just zero, but it has to be written

    #Inlet connector variables's equations
    atm_pressure = [P_atm ~ In.P]
    inletenthalpy = [h_in ~ In.H]
    inlettemperature_eqs = [T_in ~ In.T]
    inletmolarflow_eqs = [F_in ~ In.F]

    #Outlet connector equations:
    out_conn = [Out.P ~ P_out
                Out.T ~ T
                Out.F ~ F_out
                Out.H ~ H / N
                scalarize(Out.z₁ .~ Nᵢ / N)...]

    if phase == :liquid
        out_conn_phases = [scalarize(Out.z₂ .~ 0.0)...
                           scalarize(Out.z₃ .~ Nᵢ / N)...
                           Out.α_g ~ 0.0]

    elseif phase == :vapor
        out_conn_phases = [scalarize(Out.z₂ .~ Nᵢ / N)...
                           scalarize(Out.z₃ .~ 0.0)...
                           Out.α_g ~ 1.0]
    end

    #balances

    component_balance = [D(Nᵢ[i]) ~ F_in * In.z₁[i] - F_out * Nᵢ[i] / N + R[i] * V
                         for i in 1:Nc] #Neglectable loss to vapor phase head space
    energy_balance = [D(H) ~ F_in * h_in - F_out * H / N + Q̇]
    jacket_energy_balance = [Q̇ ~ -2.27 * 4184.0 * (T - 288.7) *
                                  (1.0 - exp(-8440.26 / (2.27 * 4184)))]
    mass_volume_eq = [ρʷ * V ~ M, ρ * V ~ N]
    mol_holdup = [N ~ sum(collect(Nᵢ))]
    mol_to_concentration = [scalarize(Nᵢ .~ Cᵢ * V)...]
    height_to_volume = [height * A ~ V]
    volumetricflow_to_molarflow = [Q_out ~ F_out / ρ, F_out ~ F_in]
    volumetricflow_to_massflow = [Q_out ~ Fʷ_out / ρʷ]

    #Thermodynamic properties (outlet)
    pressure_out = [phase == :liquid ? P_out ~ P_atm : P_out ~ P_atm]   # Estimation considering static pressure (May be off as tank is agitated and not static)
    density_eqs = [ρ ~ molar_density(model, P_out, T, Nᵢ, phase = "unknown")]
    mass_density = [ρʷ ~ ρ * MW]
    globalEnthalpy_eq = [
        H ~ enthalpy(model, P_out, T, Nᵢ, phase = "unknown") + sum(scalarize(ΔH₀f .* Nᵢ)),
        fᵢ ~ fugacity_coefficient(model, P_out, T, Nᵢ, phase = "liquid")[1]] #fᵢ ~ fugacity_coefficient(model, P_out, T, Nᵢ, "liquid")[1]
    molar_mass = [MW ~ sum(scalarize(MWs .* Nᵢ)) / N * gramsToKilograms]

    eqs = [reaction_rate...; overall_reaction_rate...; atm_pressure...; inletenthalpy...;
           inlettemperature_eqs...; inletmolarflow_eqs...; out_conn...;
           out_conn_phases...; component_balance...; energy_balance...;
           jacket_energy_balance...; mass_volume_eq...; mol_holdup...;
           mol_to_concentration...;
           height_to_volume...; volumetricflow_to_molarflow...;
           volumetricflow_to_massflow...;
           pressure_out...; density_eqs...; mass_density...; globalEnthalpy_eq...;
           molar_mass...]

    ODESystem([eqs...;], t, collect(Iterators.flatten(vars)),
        collect(Iterators.flatten(pars)); name, systems = [Ports...])
end

export CSTR
