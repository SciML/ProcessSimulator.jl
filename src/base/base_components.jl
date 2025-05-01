@component function ρTz_ThermodynamicState_(;medium, name)

    pars = []

    vars = @variables begin
        ϕ(t)[1:medium.Constants.nphases - 1],                                             [description = "phase fraction"]
        ρ(t)[1:medium.Constants.nphases],                                                 [description = "molar density"] 
        z(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                     [description = "mole fraction", irreducible = true] 
        p(t),                                                                                  [description = "pressure"]
        T(t),                                                                                  [description = "Temperature"]
    end

    eq = [
        [ρ[j] ~ PT_molar_density(medium.EoSModel, p, T, collect(z[:, j]), phase = medium.Constants.phaseNames[j]) for j in 2:medium.Constants.nphases]...
        ρ[1] .~ 1.0/(sum([ϕ[j - 1]/ρ[j] for j in 2:medium.Constants.nphases]))...
    ]

    guesses = [z[i, j] => medium.Guesses.x[i, j] for i in 1:medium.Constants.Nc, j in 1:medium.Constants.nphases]

    ODESystem(eq, t, collect(Iterators.flatten(vars)), pars; guesses = guesses, name)

end


@connector function PhZConnector_(;medium, name)

    vars = @variables begin
        p(t),                                                                                      [description = "pressure"]                      
        h(t)[1:medium.Constants.nphases],                                                     [description = "molar enthalpy"]  
        z₁(t)[1:medium.Constants.Nc],                                                         [description = "overall mole fraction"]
        z₂(t)[1:medium.Constants.Nc],                                                         [description = "dense state mole fraction"] 
        z₃(t)[1:medium.Constants.Nc],                                                         [description = "dense mole fraction"]   
        ṅ(t)[1:medium.Constants.nphases],                                                     [description = "molar flow"]              
    end

    pars = []

    eqs = []  # This avoids "BoundsError: attempt to access 0-element Vector{Vector{Any}} at index [0]"
	eq = eqs==[] ? Equation[] : eqs

    ODESystem(eq, t, collect(Iterators.flatten(vars)), pars; name)

end


@connector function SurfaceConnector_(;medium, name)

    vars = @variables begin
        T(t),                                                                                      [description = "Temperature"]                      
        z(t)[1:medium.Constants.Nc],                                                          [description = "overall mole fraction"]  
        ϕₘ(t)[1:medium.Constants.Nc],                                                         [description = "Molar Flow Rate", connect = Flow]
        ϕₕ(t),                                                                                     [description = "Heat Flow Rate", connect = Flow]       
    end

    pars = []

    eqs = []  # This avoids "BoundsError: attempt to access 0-element Vector{Vector{Any}} at index [0]"
	eq = eqs==[] ? Equation[] : eqs

    ODESystem(eq, t, collect(Iterators.flatten(vars)), pars; name)

end

@component function Surface(;medium, name)

        systems = @named begin
            InPort = SurfaceConnector_(medium = medium)
            OutPort = SurfaceConnector_(medium = medium)
        end

        vars = @variables begin
            A(t),                                                                                   [description = "Area"]
            kₘ(t),                                                                                  [description = "mass transfer coefficient"]
            kₕ(t),                                                                                  [description = "heat transfer coefficient"]
        end

        eqs = [
                kₘ ~ mass_transfer_coefficient(medium)
                kₕ ~ heat_transfer_coefficient(medium)
                ϕₕ ~ A*kₕ*(InPort.T - OutPort.T)
                scalarize(InPort.ϕₘ .~ A*kₘ*(InPort.z .- OutPort.z))...
                InPort.ϕₘ + OutPort.ϕₘ ~ 0
                InPort.ϕₕ + OutPort.ϕₕ ~ 0
        ]

        ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])

end



@component function AdsorptionInterface(;fluidmedium, solidmedium, name)

    systems = @named begin
        SolidSurface = Surface(medium = solidmedium)
        FluidSurface = Surface(medium = fluidmedium)
    end

    eqs = [domain_connect(FluidSurface.OutPort, SolidSurface.InPort)
            SolidSurface.InPort.z .~ loading(solidmedium.IsothermModel, FluidSurface.OutPort.z) 
        ]

    ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])

end



#### ------ ControlVolumes ------

@component function HeatConnector(;name)
    vars = @variables begin
        Q(t),              [description = "heat flux"]  #, unit=u"J s^-1"]
    end

    return ODESystem(Equation[], t, vars, []; name)
end

@component function WorkConnector(;name)
    vars = @variables begin
        W(t),              [description = "power"]      #, unit=u"J s^1"]
    end

    return ODESystem(Equation[], t, vars, []; name)
end


@component function TwoPortControlVolume_(;medium, name)


    systems = @named begin
        InPort = PhZConnector_(medium = medium)
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end

    vars = @variables begin
        rᵥ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                            [description = "mass source or sink - volumetric  basis"]
        rₐ(t)[1:medium.Constants.Nc, 1:medium.Constants.nphases],                            [description = "molar source or sink - through surface in contact with other phases"]      
        Nᵢ(t)[1:medium.Constants.Nc],                                                             [description = "molar holdup"]
        nᴸⱽ(t)[1:medium.Constants.nphases - 1],                                                   [description = "molar holdup"]              
        U(t),                                                                                          [description = "internal energy holdup"]                                               
        p(t),                                                                                          [description = "pressure"]  
        V(t)[1:medium.Constants.nphases],                                                         [description = "volume"]    
        Q(t),                                                                                          [description = "heat flux"]
        Wₛ(t),                                                                                         [description = "shaft work"]
    end

    pars = []

 
    eqs = [

        # Energy balance

        D(U) ~ InPort.h[1]*InPort.ṅ[1] + OutPort.h[1]*(OutPort.ṅ[1]) + Q + Wₛ

        # Mole balance

        scalarize(D(Nᵢ) .~ InPort.ṅ[1].*InPort.z₁ + (OutPort.ṅ[2].*OutPort.z₂ + OutPort.ṅ[3].*OutPort.z₃) .+ collect(rᵥ[:, 2:end]*V[2:end]))...

        scalarize(rᵥ[:, 1] .~ sum(collect(rᵥ[:, 2:end]), dims = 2))...

        scalarize(rₐ[:, 1] .~ sum(collect(rₐ[:, 2:end]), dims = 2))...

        scalarize(Nᵢ .~ nᴸⱽ[1]*ControlVolumeState.z[:, 2] + nᴸⱽ[2]*ControlVolumeState.z[:, 3])...

        scalarize(sum(Nᵢ) ~ sum(nᴸⱽ))

        scalarize(sum(collect(ControlVolumeState.ϕ)) ~ 1.0)

        
        # Thermodynamic state equations

        scalarize(ControlVolumeState.z[:, 1] .~ Nᵢ ./ sum(collect(Nᵢ)))...

        ControlVolumeState.ϕ[1] ~ nᴸⱽ[1]/sum(collect(nᴸⱽ))

        
        # Control Volume properties

        U ~ (OutPort.h[1] - ControlVolumeState.p/ControlVolumeState.ρ[1])*sum(collect(Nᵢ)) 

        V[1] ~ sum(collect(V[2:end]))
        
        V[2]*ControlVolumeState.ρ[2] ~ nᴸⱽ[1]

        V[3]*ControlVolumeState.ρ[3] ~ nᴸⱽ[2]
        
        
        # Outlet port properties

        [OutPort.h[j] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[j], ControlVolumeState.T, collect(ControlVolumeState.z[:, j])) for j in 2:medium.Constants.nphases]...
        
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        OutPort.p ~ ControlVolumeState.p

        scalarize(OutPort.z₁ .~ ControlVolumeState.z[:, 1])...
        scalarize(OutPort.z₂ .~ ControlVolumeState.z[:, 2])...
        scalarize(OutPort.z₃ .~ ControlVolumeState.z[:, 3])...

        # Specifics for the two-phase system

        #= ControlVolumeState.p ~ p
        scalarize(ControlVolumeState.z[:, 3] ~ flash_mol_fractions_vapor(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1])))...
        scalarize(nᴸⱽ[1]/sum(nᴸⱽ) ~ flash_vaporized_fraction(medium.EoSModel, ControlVolumeState.p, ControlVolumeState.T, collect(ControlVolumeState.z[:, 1]))[1]) =#

        ]

        return ODESystem(eqs, t, collect(Iterators.flatten(vars)), pars; name, systems = [systems...])

end

@component function FixedBoundary_pTzn_(;medium, p, T, z, ṅ, name)

    systems = @named begin
        OutPort = PhZConnector_(medium = medium)
        ControlVolumeState = ρTz_ThermodynamicState_(medium = medium)
    end
    
    vars = []

    pars = []

    eqs = [
        # Port equations
        OutPort.p ~ p
        scalarize(OutPort.z₁ .~ z)...
        scalarize(OutPort.z₂ .~ ControlVolumeState.z[:, 2])...
        scalarize(OutPort.z₃ .~ ControlVolumeState.z[:, 3])...
        OutPort.ṅ[1] ~ ṅ
        OutPort.ṅ[2] ~ ControlVolumeState.ϕ[1]*ṅ
        OutPort.ṅ[3] ~ ControlVolumeState.ϕ[2]*ṅ
        [OutPort.h[i] ~ ρT_enthalpy(medium.EoSModel, ControlVolumeState.ρ[i], ControlVolumeState.T, ControlVolumeState.z[:, i]) for i in 2:medium.Constants.nphases]...
        OutPort.h[1] ~ dot(collect(OutPort.h[2:end]), collect(ControlVolumeState.ϕ))

        # State equations
        scalarize(ControlVolumeState.z[:, 1] .~ z)...
        scalarize(ControlVolumeState.z[:, 2] .~ flash_mol_fractions_liquid(medium.EoSModel, p, T, z))...
        scalarize(ControlVolumeState.z[:, 3] .~ flash_mol_fractions_vapor(medium.EoSModel, p, T, z))...
        ControlVolumeState.T ~ T
        OutPort.ṅ[2] + OutPort.ṅ[3] ~ OutPort.ṅ[1]
        ControlVolumeState.p ~ p
        scalarize(ControlVolumeState.ϕ[1] ~ flash_vaporized_fraction(medium.EoSModel, p, T, z)[1])
    ]

    return ODESystem(eqs, t, vars, collect(Iterators.flatten(pars)); name, systems = [systems...])


end




