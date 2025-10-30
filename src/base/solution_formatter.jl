using Printf
using PrettyTables

"""
    extract_cstr_solution(sol, sys, cstr_name::Symbol; time_index=:end, components::Vector{String})

Extract CSTR solution data from solved system and format as tables.
"""
function extract_cstr_solution(sol, sys, cstr_name::Symbol; time_index=:end, components::Vector{String})
    # Handle both ODESolution (has time) and NonlinearSolution (steady-state, no time)
    t_idx = if hasproperty(sol, :t)
        time_index == :end ? length(sol.t) : time_index
    else
        1  # NonlinearSolution - just one solution point
    end

    Nc = length(components)

    result = Dict{Symbol, Any}()

    # STREAM TABLE
    stream_data = Matrix{Any}(undef, 7, 3)
    stream_data[:, 1] = ["Temperature", "Pressure", "Total Molar Flow",
                         "Liquid Molar Flow", "Vapor Molar Flow",
                         "Molar Enthalpy", "Liquid Fraction"]

    # Inlet - Note: Ports don't have T directly, use ControlVolumeState for temperature
    try
        reactor = getproperty(sys, cstr_name)
        stream_data[1, 2] = "N/A"  # Temperature not available at port
        stream_data[2, 2] = sol[reactor.InPort.p][t_idx] / 1e5
        stream_data[3, 2] = abs(sol[reactor.InPort.ṅ[1]][t_idx])
        stream_data[4, 2] = abs(sol[reactor.InPort.ṅ[2]][t_idx])
        stream_data[5, 2] = abs(sol[reactor.InPort.ṅ[3]][t_idx])
        stream_data[6, 2] = sol[reactor.InPort.h[1]][t_idx]
        stream_data[7, 2] = "-"
    catch e
        @warn "Could not extract inlet stream data for $cstr_name" exception=(e, catch_backtrace())
        stream_data[:, 2] .= "N/A"
    end

    # Outlet
    try
        reactor = getproperty(sys, cstr_name)
        stream_data[1, 3] = sol[reactor.ControlVolumeState.T][t_idx]  # Use reactor temperature
        stream_data[2, 3] = sol[reactor.OutPort.p][t_idx] / 1e5
        stream_data[3, 3] = abs(sol[reactor.OutPort.ṅ[1]][t_idx])
        stream_data[4, 3] = abs(sol[reactor.OutPort.ṅ[2]][t_idx])
        stream_data[5, 3] = abs(sol[reactor.OutPort.ṅ[3]][t_idx])
        stream_data[6, 3] = sol[reactor.OutPort.h[1]][t_idx]
        stream_data[7, 3] = "-"
    catch e
        @warn "Could not extract outlet stream data for $cstr_name" exception=(e, catch_backtrace())
        stream_data[:, 3] .= "N/A"
    end

    result[:streams] = stream_data

    # Units
    result[:stream_units] = ["K", "bar", "kmol/hr", "kmol/hr", "kmol/hr", "kJ/kmol", ""]

    # COMPOSITION TABLE
    comp_data = Matrix{Any}(undef, Nc, 3)
    reactor = getproperty(sys, cstr_name)
    for (i, comp) in enumerate(components)
        comp_data[i, 1] = comp
        try
            comp_data[i, 2] = sol[reactor.InPort.z[i,1]][t_idx]
        catch
            comp_data[i, 2] = "N/A"
        end
        try
            comp_data[i, 3] = sol[reactor.OutPort.z[i,1]][t_idx]
        catch
            comp_data[i, 3] = "N/A"
        end
    end
    result[:composition] = comp_data

    # EQUIPMENT TABLE
    equip_data = Matrix{Any}(undef, 9, 2)
    equip_data[:, 1] = ["Total Volume", "Liquid Volume", "Vapor Volume",
                        "Temperature", "Pressure", "Liquid Fraction",
                        "Internal Energy", "Heat Duty", "Shaft Work"]

    try
        reactor = getproperty(sys, cstr_name)
        equip_data[1, 2] = sol[reactor.V[1]][t_idx]
        equip_data[2, 2] = sol[reactor.V[2]][t_idx]
        equip_data[3, 2] = sol[reactor.V[3]][t_idx]
        equip_data[4, 2] = sol[reactor.ControlVolumeState.T][t_idx]
        equip_data[5, 2] = sol[reactor.ControlVolumeState.p][t_idx] / 1e5
        equip_data[6, 2] = sol[reactor.ControlVolumeState.ϕ[1]][t_idx]
        equip_data[7, 2] = sol[reactor.U][t_idx]
        equip_data[8, 2] = sol[reactor.Q][t_idx]
        equip_data[9, 2] = sol[reactor.Wₛ][t_idx]
    catch e
        @warn "Could not extract equipment data for $cstr_name" exception=(e, catch_backtrace())
        equip_data[:, 2] .= "N/A"
    end
    result[:equipment] = equip_data
    result[:equip_units] = ["m3", "m3", "m3", "K", "bar", "", "kJ", "kW", "kW"]

    # HOLDUP TABLE
    holdup_data = Matrix{Any}(undef, Nc, 2)
    reactor = getproperty(sys, cstr_name)
    for (i, comp) in enumerate(components)
        holdup_data[i, 1] = comp
        try
            holdup_data[i, 2] = sol[reactor.Nᵢ[i]][t_idx]
        catch
            holdup_data[i, 2] = "N/A"
        end
    end
    result[:holdup] = holdup_data

    # CONVERSION TABLE
    conversion_data = Matrix{Any}(undef, Nc, 2)
    reactor = getproperty(sys, cstr_name)
    for (i, comp) in enumerate(components)
        conversion_data[i, 1] = comp
        try
            X_val = sol[reactor.X[i]][t_idx]
            conversion_data[i, 2] = X_val * 100
        catch
            conversion_data[i, 2] = "N/A"
        end
    end
    result[:conversion] = conversion_data

    return result
end

"""
    extract_boundary_solution(sol, sys, boundary_name::Symbol; time_index=:end, components::Vector{String})

Extract fixed boundary solution data from solved system.
"""
function extract_boundary_solution(sol, sys, boundary_name::Symbol; time_index=:end, components::Vector{String})
    # Handle both ODESolution (has time) and NonlinearSolution (steady-state, no time)
    t_idx = if hasproperty(sol, :t)
        time_index == :end ? length(sol.t) : time_index
    else
        1  # NonlinearSolution - just one solution point
    end

    Nc = length(components)
    result = Dict{Symbol, Any}()

    # STREAM TABLE - Only outlet for boundary
    stream_data = Matrix{Any}(undef, 7, 2)
    stream_data[:, 1] = ["Temperature", "Pressure", "Total Molar Flow",
                         "Liquid Molar Flow", "Vapor Molar Flow",
                         "Molar Enthalpy", "Liquid Fraction"]

    try
        boundary = getproperty(sys, boundary_name)
        stream_data[1, 2] = sol[boundary.ControlVolumeState.T][t_idx]
        stream_data[2, 2] = sol[boundary.OutPort.p][t_idx] / 1e5
        stream_data[3, 2] = abs(sol[boundary.OutPort.ṅ[1]][t_idx])
        stream_data[4, 2] = abs(sol[boundary.OutPort.ṅ[2]][t_idx])
        stream_data[5, 2] = abs(sol[boundary.OutPort.ṅ[3]][t_idx])
        stream_data[6, 2] = sol[boundary.OutPort.h[1]][t_idx]
        stream_data[7, 2] = sol[boundary.ControlVolumeState.ϕ[1]][t_idx]
    catch e
        @warn "Could not extract boundary stream data for $boundary_name" exception=(e, catch_backtrace())
        stream_data[:, 2] .= "N/A"
    end

    result[:streams] = stream_data
    result[:stream_units] = ["K", "bar", "kmol/hr", "kmol/hr", "kmol/hr", "kJ/kmol", ""]

    # COMPOSITION TABLE
    comp_data = Matrix{Any}(undef, Nc, 2)
    boundary = getproperty(sys, boundary_name)
    for (i, comp) in enumerate(components)
        comp_data[i, 1] = comp
        try
            comp_data[i, 2] = sol[boundary.OutPort.z[i,1]][t_idx]
        catch
            comp_data[i, 2] = "N/A"
        end
    end
    result[:composition] = comp_data

    return result
end

"""
    print_boundary_report(sol, sys, boundary_name::Symbol, components::Vector{String}; time_index=:end)

Print boundary stream report in Aspen Plus style.
"""
function print_boundary_report(sol, sys, boundary_name::Symbol, components::Vector{String}; time_index=:end)
    data = extract_boundary_solution(sol, sys, boundary_name; time_index=time_index, components=components)
    io = stdout

    # Title Block
    println(io, "\n")
    println(io, " " * "="^78)
    println(io, " " * " "^78)
    println(io, center_string("BLOCK:  $(uppercase(String(boundary_name)))  MODEL: BOUNDARY", 78))
    println(io, " " * " "^78)
    println(io, " " * "="^78)
    println(io)

    # OUTLET STREAM SECTION
    println(io, "\n OUTLET MATERIAL STREAM")
    println(io, " " * "-"^78)
    pretty_table(io, data[:streams],
                 header = ["", "OUTLET"],
                 alignment = [:l, :r],
                 formatters = (ft_printf("%.4f", 2),),
                 tf = tf_simple)

    # COMPOSITION SECTION
    println(io, "\n COMPOSITION (MOLE FRACTION)")
    println(io, " " * "-"^78)
    pretty_table(io, data[:composition],
                 header = ["Component", "Outlet"],
                 alignment = [:l, :r],
                 formatters = (ft_printf("%.6f", 2),),
                 tf = tf_simple)

    println(io, "\n" * " " * "="^78 * "\n")
end

"""
    print_cstr_report(sol, sys, cstr_name::Symbol, components::Vector{String}; time_index=:end)

Print CSTR report in Aspen Plus style using PrettyTables.
"""
function print_cstr_report(sol, sys, cstr_name::Symbol, components::Vector{String}; time_index=:end)
    data = extract_cstr_solution(sol, sys, cstr_name; time_index=time_index, components=components)
    io = stdout

    # Title Block
    println(io, "\n")
    println(io, " " * "="^78)
    println(io, " " * " "^78)
    println(io, center_string("BLOCK:  $(uppercase(String(cstr_name)))  MODEL: RCSTR", 78))
    println(io, " " * " "^78)
    println(io, " " * "="^78)
    println(io)

    # INLET STREAMS SECTION
    println(io, "\n INLET MATERIAL STREAM")
    println(io, " " * "-"^78)
    pretty_table(io, data[:streams],
                 header = ["", "INLET", "OUTLET"],
                 alignment = [:l, :r, :r],
                 formatters = (ft_printf("%.4f", [2, 3]),),
                 tf = tf_simple)

    # COMPOSITION SECTION
    println(io, "\n COMPOSITION (MOLE FRACTION)")
    println(io, " " * "-"^78)
    pretty_table(io, data[:composition],
                 header = ["Component", "Inlet", "Outlet"],
                 alignment = [:l, :r, :r],
                 formatters = (ft_printf("%.6f", [2, 3]),),
                 tf = tf_simple)

    # EQUIPMENT SPECIFICATIONS
    println(io, "\n EQUIPMENT SPECIFICATIONS")
    println(io, " " * "-"^78)
    pretty_table(io, data[:equipment],
                 header = ["Property", "Value"],
                 alignment = [:l, :r],
                 formatters = (ft_printf("%.6e", 2),),
                 tf = tf_simple)

    # COMPONENT HOLDUP
    println(io, "\n COMPONENT HOLDUP (mol)")
    println(io, " " * "-"^78)
    pretty_table(io, data[:holdup],
                 header = ["Component", "Holdup"],
                 alignment = [:l, :r],
                 formatters = (ft_printf("%.4f", 2),),
                 tf = tf_simple)

    # CONVERSION
    println(io, "\n CONVERSION")
    println(io, " " * "-"^78)
    pretty_table(io, data[:conversion],
                 header = ["Component", "Conversion (%)"],
                 alignment = [:l, :r],
                 formatters = (ft_printf("%.2f", 2),),
                 tf = tf_simple)

    println(io, "\n" * " " * "="^78 * "\n")
end

"""
    print_flowsheet_summary(sol, sys, unit_ops::Dict{Symbol, Symbol}, components::Vector{String})

Print complete flowsheet summary report (Aspen Plus style).

# Arguments
- `sol`: Solution from solver
- `sys`: Simplified/compiled system
- `unit_ops`: Dict mapping unit names to model types (e.g., Dict(:R1 => :CSTR, :S1 => :Feed))
- `components`: Vector of component names
"""
function print_flowsheet_summary(sol, sys, unit_ops::Dict{Symbol, Symbol}, components::Vector{String})
    io = stdout

    # Main Header
    println(io, "\n\n")
    println(io, " " * "="^78)
    println(io, center_string("PROCESS SIMULATION RESULTS", 78))
    println(io, center_string("ProcessSimulator.jl - Equation-Based Flowsheet Simulator", 78))
    println(io, " " * "="^78)

    # Simulation Summary
    println(io, "\n\n SIMULATION SUMMARY")
    println(io, " " * "-"^78)

    summary_data = Matrix{Any}(undef, 4, 2)
    summary_data[1, :] = ["Run Status", string(sol.retcode)]
    summary_data[2, :] = ["Solution Time", "Steady State"]
    summary_data[3, :] = ["Number of Components", length(components)]
    summary_data[4, :] = ["Components", join(components, ", ")]

    pretty_table(io, summary_data,
                 header = ["Parameter", "Value"],
                 alignment = [:l, :l],
                 tf = tf_simple)

    # Unit Operations List
    println(io, "\n\n UNIT OPERATIONS")
    println(io, " " * "-"^78)

    units_data = Matrix{Any}(undef, length(unit_ops), 2)
    for (idx, (name, type)) in enumerate(unit_ops)
        units_data[idx, 1] = String(name)
        units_data[idx, 2] = String(type)
    end

    pretty_table(io, units_data,
                 header = ["Block ID", "Model"],
                 alignment = [:l, :l],
                 tf = tf_simple)

    # Print each unit operation report
    for (unit_name, unit_type) in unit_ops
        if unit_type == :CSTR
            print_cstr_report(sol, sys, unit_name, components)
        elseif unit_type == :Boundary || unit_type == :Feed
            print_boundary_report(sol, sys, unit_name, components)
        elseif unit_type == :FlashDrum
            println(io, "\n BLOCK:  $(uppercase(String(unit_name)))  MODEL: FLASH")
            println(io, " " * "-"^78)
            println(io, " [Flash drum report not yet implemented]")
            println(io)
        else
            println(io, "\n BLOCK:  $(uppercase(String(unit_name)))  MODEL: $(uppercase(String(unit_type)))")
            println(io, " " * "-"^78)
            println(io, " [Report formatter not yet implemented]")
            println(io)
        end
    end

    # Footer
    println(io, "\n" * " " * "="^78)
    println(io, center_string("END OF REPORT", 78))
    println(io, " " * "="^78 * "\n")
end

"""
    center_string(s::String, width::Int)

Center string within specified width.
"""
function center_string(s::String, width::Int)
    len = length(s)
    if len >= width
        return s[1:width]
    end
    left_pad = div(width - len, 2)
    right_pad = width - len - left_pad
    return " "^left_pad * s * " "^right_pad
end

"""
    print_flowsheet_summary(sol, sys, components::Vector{String}, units...)

Print complete flowsheet summary report (Aspen Plus style) with automatic model type detection.

Automatically detects model types from unit operation structs that have a `model_type` field.

# Arguments
- `sol`: Solution from solver
- `sys`: Simplified/compiled system
- `components`: Vector of component names
- `units...`: Variable number of unit operation objects (e.g., R1, S1, etc.)

# Example
```julia
@named R1 = FixedVolumeSteadyStateCSTR(...)
@named S1 = FixedBoundary_pTzn_(...)

# Automatically detects R1 is CSTR, S1 is boundary
print_flowsheet_summary(sol, simplified_sys, components, R1, S1)
```
"""
function print_flowsheet_summary(sol, sys, components::Vector{String}, units...)
    # Build unit_ops dictionary automatically from model_type field
    unit_ops = Dict{Symbol, Symbol}()

    for unit in units
        unit_name = Symbol(nameof(unit.odesystem))

        # Try to get model_type from struct
        if hasproperty(unit, :model_type)
            unit_ops[unit_name] = unit.model_type
        else
            # Fallback: try to infer from type name
            type_str = string(typeof(unit))
            if occursin("CSTR", type_str)
                unit_ops[unit_name] = :CSTR
            elseif occursin("Boundary", type_str) || occursin("Feed", type_str)
                unit_ops[unit_name] = :Feed
            elseif occursin("Flash", type_str)
                unit_ops[unit_name] = :FlashDrum
            else
                unit_ops[unit_name] = :Unknown
            end
        end
    end

    # Call the main formatter
    print_flowsheet_summary(sol, sys, unit_ops, components)
end

export extract_cstr_solution, print_cstr_report, extract_boundary_solution, print_boundary_report, print_flowsheet_summary
