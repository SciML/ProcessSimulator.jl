using UnicodePlots
using Printf

function Base.show(io::IO, valve::Valve)
    # Header with valve type and phase
    println(io, "Valve Model ($(valve.phase))")
    println(io, "=" ^ 60)
    
    # Flow coefficient and opening setpoint
    println(io, "Flow Coefficient (Cv): $(valve.Cv)")
    println(io, "Valve Opening Setpoint: $(valve.opening_setpoint * 100)% open")
    
    # Pressure drop function representation
    println(io, "\nPressure Drop Equation:")
    println(io, "  ΔP = f(Δp) where:")
    
    # Function representation with Unicode plot - with improved debugging
    if hasproperty(valve, :f)
        # Plot the function using UnicodePlots - only for positive pressure drops
        try 
            # Test if the function can be evaluated
            test_val = try 
                val = valve.f(100.0)
                println(io, "    • Sample evaluation: f(100.0) = $val")
                val
            catch e
                println(io, "    • Function evaluation failed: $(sprint(showerror, e))")
                nothing
            end
            
            if test_val !== nothing
                # Only positive pressure drops are physically meaningful for a valve
                x_range = range(1.0, 100_000.0, length=40)  # Reduced number of points for stability
                y_values = Float64[]
                x_valid = Float64[]
                
                # Evaluate points individually to identify issues
                for x in x_range
                    try
                        y = valve.f(x)
                        if !isnan(y) && !isinf(y)
                            push!(x_valid, x)
                            push!(y_values, y)
                        end
                    catch
                        # Skip problematic points
                    end
                end
                
                # Check if we have any valid points
                if !isempty(y_values)
                    # Create plot
                    plt = lineplot(
                        x_valid, 
                        y_values,
                        title="Valve Function Behavior",
                        xlabel="Pressure Drop (Pa)",
                        ylabel="Flow Factor",
                        width=30, 
                        height=10,
                        border=:ascii
                    )
                    
                    # Print the plot
                    println(io, "\n    • Function visualization (for positive pressure drops):")
                    show(io, plt)
                    println(io)
                else
                    println(io, "\n    • Function visualization unavailable (no valid points to plot)")
                end
            else
                println(io, "\n    • Function visualization unavailable (function evaluation failed)")
            end
        catch e
            println(io, "\n    • Function visualization unavailable: $(sprint(showerror, e))")
            println(io, "      $(typeof(e)): $(e)")
        end
    else
        println(io, "\n    • Function property not found on valve object")
    end
    
    println(io, "    • Flow equation: ṅ/ρ = opening·Cv·f(Pin - Pout)")
    
    # State guess section (thermodynamic state)
    println(io, "\nThermodynamic State Guess:")
    println(io, "    • Pressure: $(valve.state.p) Pa")
    println(io, "    • Temperature: $(valve.state.T) K")
    println(io, "    • Composition: $(valve.state.N ./ sum(valve.state.N))")
    
    # Density information from medium's guesses if available
    if hasproperty(valve.medium, :Guesses) && hasproperty(valve.medium.Guesses, :ρ)
        println(io, "    • Density:")
        for i in eachindex(valve.medium.Guesses.ρ)
            phase_name = i == 1 ? "overall" : (i == 2 ? "liquid" : "vapor")
            println(io, "        - $(phase_name): $(round(valve.medium.Guesses.ρ[i], digits=3)) mol/m³")
        end
    end
    
    # Molar flowrate guess
    println(io, "\nMolar Flowrate Guess: $(valve.molar_flowrate_guess) mol/s")
    
    # Medium information with detailed breakdown
    println(io, "\nMedium Details:")
    if isa(valve.medium, EoSBased)
        # EOS Model information
        eos_type = split(string(typeof(valve.medium.EoSModel)), "{")[1]
        println(io, "    • EoS Model: $eos_type")
        
        # Components information
        if hasproperty(valve.medium.Constants, :iupacName) && !isnothing(valve.medium.Constants.iupacName)
            println(io, "    • Components: $(valve.medium.Constants.iupacName)")
            if hasproperty(valve.medium.Constants, :molarMass) && !isnothing(valve.medium.Constants.molarMass)
                mw_str = join(["$(valve.medium.Constants.iupacName[i]): $(valve.medium.Constants.molarMass[i]) kg/mol" 
                               for i in 1:length(valve.medium.Constants.iupacName)], ", ")
                println(io, "    • Molar masses: $mw_str")
            end
        end
        
        # Transport model breakdown
        if hasproperty(valve.medium, :TransportModel) && !isnothing(valve.medium.TransportModel)
            println(io, "    • Transport Models:")
            
            # Mass transfer model
            if hasproperty(valve.medium.TransportModel, :MassTransferModel) && 
               !isnothing(valve.medium.TransportModel.MassTransferModel)
                mass_model = typeof(valve.medium.TransportModel.MassTransferModel)
                println(io, "        - Mass Transfer: $(split(string(mass_model), "{")[1])")
                if isa(valve.medium.TransportModel.MassTransferModel, ConstantMassTransferCoeff)
                    println(io, "          Coefficients: $(valve.medium.TransportModel.MassTransferModel.k)")
                end
            else
                println(io, "        - Mass Transfer: None")
            end
            
            # Heat transfer model
            if hasproperty(valve.medium.TransportModel, :HeatTransferModel) && 
               !isnothing(valve.medium.TransportModel.HeatTransferModel)
                heat_model = typeof(valve.medium.TransportModel.HeatTransferModel)
                println(io, "        - Heat Transfer: $(split(string(heat_model), "{")[1])")
                if isa(valve.medium.TransportModel.HeatTransferModel, ConstantHeatTransferCoeff)
                    println(io, "          Coefficient: $(valve.medium.TransportModel.HeatTransferModel.k)")
                end
            else
                println(io, "        - Heat Transfer: None")
            end
            
            # Viscosity model
            if hasproperty(valve.medium.TransportModel, :ViscosityModel) && 
               !isnothing(valve.medium.TransportModel.ViscosityModel)
                visc_model = typeof(valve.medium.TransportModel.ViscosityModel)
                println(io, "        - Viscosity: $(split(string(visc_model), "{")[1])")
            else
                println(io, "        - Viscosity: None")
            end
        end
    else
        # For non-EoSBased media
        medium_type = split(string(typeof(valve.medium)), "{")[1]
        println(io, "    • Type: $medium_type")
    end
    
    # End with a separator
    println(io, "=" ^ 60)
end




function Base.show(io::IO, model::SolidEoSModel)
    # Header with type name
    println(io, "SolidEoSModel")
    println(io, "=" ^ 50)
    
    # Density section
    println(io, "\nDensity Parameters:")
    println(io, "    • Base density (ρ0): $(round(model.ρ0, digits=3)) kg/m³")
    println(io, "    • Reference temperature (ρ_T0): $(round(model.ρ_T0, digits=3)) K")
    
    # Print density coefficients
    println(io, "    • Temperature coefficients:")
    if all(x -> x === nothing || x == 0, model.coeffs_ρ)
        println(io, "        - No temperature dependence")
    else
        for (i, coef) in enumerate(model.coeffs_ρ)
            if coef !== nothing
                println(io, "        - Order $(i): $(round(coef, digits=6))")
            end
        end
        
        # Show density equation
        println(io, "\n    • Density equation: ρ(T) = ρ0 + Σ coeff_i * ΔT^i")
        print(io, "      where ρ(T) = $(round(model.ρ0, digits=2))")
        
        for (i, coef) in enumerate(model.coeffs_ρ)
            if coef !== nothing && coef != 0
                sign = coef > 0 ? " + " : " - "
                print(io, "$(sign)$(round(abs(coef), digits=6)) * ΔT^$(i)")
            end
        end
        println(io)
        println(io, "      with ΔT = T - $(round(model.ρ_T0, digits=2))")
    end
    
    # Enthalpy section
    println(io, "\nEnthalpy Parameters:")
    println(io, "    • Base enthalpy (h0): $(round(model.h0, digits=3)) J/mol")
    println(io, "    • Reference temperature (h_T0): $(round(model.h_T0, digits=3)) K")
    
    # Print enthalpy coefficients
    println(io, "    • Temperature coefficients:")
    if all(x -> x === nothing || x == 0, model.coeffs_h)
        println(io, "        - No temperature dependence")
    else
        for (i, coef) in enumerate(model.coeffs_h)
            if coef !== nothing
                println(io, "        - Order $(i): $(round(coef, digits=6))")
            end
        end
        
        # Show enthalpy equation
        println(io, "\n    • Enthalpy equation: h(T) = h0 + Σ coeff_i * ΔT^i")
        print(io, "      where h(T) = $(round(model.h0, digits=2))")
        
        for (i, coef) in enumerate(model.coeffs_h)
            if coef !== nothing && coef != 0
                sign = coef > 0 ? " + " : " - "
                print(io, "$(sign)$(round(abs(coef), digits=6)) * ΔT^$(i)")
            end
        end
        println(io)
        println(io, "      with ΔT = T - $(round(model.h_T0, digits=2))")
    end
    
    # End with separator
    println(io, "=" ^ 50)
end


function Base.show(io::IO, fluid::EoSBased)
    # Header with type
    println(io, "EoSBased Fluid Medium")
    println(io, "=" ^ 50)
    
    # Constants section
    println(io, "\nFluid Constants:")
    println(io, "    • Number of components: $(fluid.Constants.Nc)")
    println(io, "    • Number of phases: $(fluid.Constants.nphases)")
    
    # Components information
    if !isnothing(fluid.Constants.iupacName)
        println(io, "    • Components: $(fluid.Constants.iupacName)")
        
        # Molar masses if available
        if !isnothing(fluid.Constants.molarMass)
            println(io, "    • Molar masses [kg/mol]:")
            for i in 1:fluid.Constants.Nc
                name = fluid.Constants.iupacName[i]
                mass = fluid.Constants.molarMass[i]
                println(io, "        - $(name): $(round(mass, digits=6))")
            end
        end
    end
    
    # EoS Model information
    println(io, "\nEquation of State:")
    eos_type = split(string(typeof(fluid.EoSModel)), "{")[1]
    println(io, "    • Model: $(eos_type)")
    
    # Transport model section
    println(io, "\nTransport Models:")
    
    # Mass transfer model
    if hasproperty(fluid.TransportModel, :MassTransferModel) && !isnothing(fluid.TransportModel.MassTransferModel)
        mass_model = typeof(fluid.TransportModel.MassTransferModel)
        println(io, "    • Mass Transfer: $(split(string(mass_model), "{")[1])")
        if isa(fluid.TransportModel.MassTransferModel, ConstantMassTransferCoeff)
            println(io, "        - Coefficients: $(fluid.TransportModel.MassTransferModel.k)")
        end
    end
    
    # Heat transfer model
    if hasproperty(fluid.TransportModel, :HeatTransferModel) && !isnothing(fluid.TransportModel.HeatTransferModel)
        heat_model = typeof(fluid.TransportModel.HeatTransferModel)
        println(io, "    • Heat Transfer: $(split(string(heat_model), "{")[1])")
        if isa(fluid.TransportModel.HeatTransferModel, ConstantHeatTransferCoeff)
            println(io, "        - Coefficient: $(fluid.TransportModel.HeatTransferModel.k)")
        end
    end
    
    # Viscosity model
    if hasproperty(fluid.TransportModel, :ViscosityModel) && !isnothing(fluid.TransportModel.ViscosityModel)
        visc_model = typeof(fluid.TransportModel.ViscosityModel)
        println(io, "    • Viscosity: $(split(string(visc_model), "{")[1])")
    end
    
    # Current state guesses
    println(io, "\nCurrent State Guesses:")
    println(io, "    • Pressure: $(round(fluid.Guesses.p, digits=2)) Pa")
    println(io, "    • Temperature: $(round(fluid.Guesses.T, digits=2)) K")
    
    # Vapor fraction
    if hasproperty(fluid.Guesses, :ϕ) && !isnothing(fluid.Guesses.ϕ)
        println(io, "    • Phase fractions: liquid=$(round(fluid.Guesses.ϕ[1], digits=4)), vapor=$(round(fluid.Guesses.ϕ[2], digits=4))")
    end
    
    # Densities
    println(io, "    • Densities [mol/m³]:")
    for i in eachindex(fluid.Guesses.ρ)
        phase_name = i == 1 ? "overall" : (i == 2 ? "liquid" : "vapor")
        println(io, "        - $(phase_name): $(round(fluid.Guesses.ρ[i], digits=3))")
    end
    
    # Enthalpies
    println(io, "    • Enthalpies [J/mol]:")
    for i in eachindex(fluid.Guesses.h)
        phase_name = i == 1 ? "overall" : (i == 2 ? "liquid" : "vapor")
        println(io, "        - $(phase_name): $(round(fluid.Guesses.h[i], digits=3))")
    end
    
    # Compositions
    println(io, "    • Compositions:")
    for i in 1:size(fluid.Guesses.x, 2)
        phase_name = i == 1 ? "overall" : (i == 2 ? "liquid" : "vapor")
        comp_str = join(["$(round(fluid.Guesses.x[j, i], digits=4))" for j in 1:fluid.Constants.Nc], ", ")
        println(io, "        - $(phase_name): [$(comp_str)]")
    end
    
    # End with separator
    println(io, "=" ^ 50)
end


# ============================================================================
# Pretty Printing for PowerLawReaction
# ============================================================================

"""
    Base.show(io::IO, rxn::PowerLawReaction)

Pretty print a PowerLawReaction showing the reaction equation, kinetics, and parameters.
"""
function Base.show(io::IO, rxn::PowerLawReaction)
    # Build reaction equation string
    reactants = String[]
    products = String[]

    for (i, sp) in enumerate(rxn.species)
        ν_val = rxn.ν[i]
        if ν_val < 0
            # Reactant
            coeff = abs(ν_val)
            if coeff ≈ 1.0
                push!(reactants, sp)
            else
                coeff_str = isinteger(coeff) ? string(Int(coeff)) : string(coeff)
                push!(reactants, "$(coeff_str)$(sp)")
            end
        elseif ν_val > 0
            # Product
            coeff = ν_val
            if coeff ≈ 1.0
                push!(products, sp)
            else
                coeff_str = isinteger(coeff) ? string(Int(coeff)) : string(coeff)
                push!(products, "$(coeff_str)$(sp)")
            end
        end
    end

    reactant_str = join(reactants, " + ")
    product_str = join(products, " + ")
    equation = isempty(products) ? reactant_str : "$(reactant_str) → $(product_str)"

    # Build rate law string
    rate_terms = String[]
    overall_order = 0.0

    for (i, sp) in enumerate(rxn.species)
        n_val = rxn.n[i]
        if n_val != 0
            overall_order += n_val
            if n_val ≈ 1.0
                push!(rate_terms, "[$(sp)]")
            else
                order_str = isinteger(n_val) ? string(Int(n_val)) : string(n_val)
                push!(rate_terms, "[$(sp)]^$(order_str)")
            end
        end
    end

    rate_law = isempty(rate_terms) ? "1" : join(rate_terms, " × ")

    # Format Arrhenius parameters
    A_str = @sprintf("%.3e", rxn.A)
    Eₐ_str = @sprintf("%.3e", rxn.Eₐ)

    # Print everything nicely
    println(io, "PowerLawReaction")
    println(io, "  Reaction: $(equation)")
    println(io, "  Rate Law: r = k(T) × $(rate_law)")
    println(io, "    where k(T) = A × exp(-Eₐ/RT)")
    println(io, "  ")
    println(io, "  Kinetic Parameters:")
    println(io, "    A  = $(A_str) (pre-exponential factor)")
    println(io, "    Eₐ = $(Eₐ_str) J/mol (activation energy)")
    println(io, "  ")
    println(io, "  Overall Order: $(overall_order)")
    println(io, "  ")
    println(io, "  Species Details:")

    # Print table header
    println(io, "    ┌────────────────┬──────────────┬──────────────┐")
    println(io, "    │ Species        │ Stoich. (ν)  │ Order (n)    │")
    println(io, "    ├────────────────┼──────────────┼──────────────┤")

    # Print each species
    for (i, sp) in enumerate(rxn.species)
        sp_padded = rpad(sp, 14)
        ν_str = @sprintf("%.2f", rxn.ν[i])
        n_str = @sprintf("%.2f", rxn.n[i])
        ν_padded = lpad(ν_str, 12)
        n_padded = lpad(n_str, 12)
        println(io, "    │ $(sp_padded) │ $(ν_padded) │ $(n_padded) │")
    end

    print(io, "    └────────────────┴──────────────┴──────────────┘")
end

"""
    Base.show(io::IO, ::MIME"text/plain", rxn::PowerLawReaction)

Compact display for PowerLawReaction in REPL.
"""
function Base.show(io::IO, ::MIME"text/plain", rxn::PowerLawReaction)
    show(io, rxn)
end


"""
    Base.show(io::IO, network::PowerLawReactionSet)

Pretty print a reaction network showing all reactions and the stoichiometric/order matrices.
"""
function Base.show(io::IO, network::PowerLawReactionSet)
    Nc = length(network.species)
    Nr = length(network.reactions)

    println(io, "PowerLawReactionSet")
    println(io, "  Number of species: $(Nc)")
    println(io, "  Number of reactions: $(Nr)")
    println(io, "  ")
    println(io, "  Reactions:")

    for j in 1:Nr
        # Build equation string
        reactants = String[]
        products = String[]

        for i in 1:Nc
            ν_val = network.ν[i, j]
            sp = network.species[i]

            if ν_val < 0
                coeff = abs(ν_val)
                if coeff ≈ 1.0
                    push!(reactants, sp)
                else
                    coeff_str = isinteger(coeff) ? string(Int(coeff)) : @sprintf("%.1f", coeff)
                    push!(reactants, "$(coeff_str)$(sp)")
                end
            elseif ν_val > 0
                coeff = ν_val
                if coeff ≈ 1.0
                    push!(products, sp)
                else
                    coeff_str = isinteger(coeff) ? string(Int(coeff)) : @sprintf("%.1f", coeff)
                    push!(products, "$(coeff_str)$(sp)")
                end
            end
        end

        reactant_str = join(reactants, " + ")
        product_str = join(products, " + ")
        equation = "$(reactant_str) → $(product_str)"

        A_str = @sprintf("%.2e", network.A[j])
        Eₐ_str = @sprintf("%.2e", network.Eₐ[j])

        println(io, "    R$(j): $(equation)")
        println(io, "        A = $(A_str), Eₐ = $(Eₐ_str) J/mol")
    end

    println(io, "  ")
    println(io, "  Stoichiometric Matrix ν (Nc × Nr):")
    print_matrix(io, network.ν, network.species, ["R$j" for j in 1:Nr], "    ")

    println(io, "  ")
    println(io, "  Reaction Order Matrix n (Nc × Nr):")
    print_matrix(io, network.n, network.species, ["R$j" for j in 1:Nr], "    ")
end

function Base.show(io::IO, ::MIME"text/plain", network::PowerLawReactionSet)
    show(io, network)
end

"""
    print_matrix(io::IO, mat::Matrix, row_labels::Vector{String}, col_labels::Vector{String}, indent::String)

Helper function to print a labeled matrix in a nice table format.
"""
function print_matrix(io::IO, mat::Matrix, row_labels::Vector{String}, col_labels::Vector{String}, indent::String)
    Nrows, Ncols = size(mat)

    # Determine column widths
    label_width = maximum(length.(row_labels))
    col_width = 10

    # Print header
    print(io, indent)
    print(io, rpad("", label_width + 2))
    for label in col_labels
        print(io, lpad(label, col_width))
    end
    println(io)

    # Print separator
    println(io, indent, "─"^(label_width + 2 + col_width * Ncols))

    # Print rows
    for i in 1:Nrows
        print(io, indent)
        print(io, rpad(row_labels[i], label_width + 2))
        for j in 1:Ncols
            val_str = @sprintf("%.2f", mat[i, j])
            print(io, lpad(val_str, col_width))
        end
        println(io)
    end
end