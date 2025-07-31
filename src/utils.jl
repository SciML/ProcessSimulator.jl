"""
    load_component_properties(filepath::String)

Load component properties from a file.

# Arguments

  - `filepath::String`: Path to the file containing the component properties.

# Returns

  - `properties::Dict`: A dictionary containing the component properties.
"""
function load_component_properties(component_name::String)
    file_path = abspath(joinpath(@__DIR__, "database/$(component_name).json"))
    if isfile(file_path)
        return JSON.parsefile(file_path)
    else
        error("Component file does not exist: $file_path")
    end
end

function read_reidcp(substances::Vector{String})
    _size = length(substances)
    a = Vector{Float64}(undef, _size)
    b = Vector{Float64}(undef, _size)
    c = Vector{Float64}(undef, _size)
    d = Vector{Float64}(undef, _size)
    e = Vector{Float64}(undef, _size)

    data = Dict(subs => load_component_properties(subs) for subs in substances)

    for i in eachindex(substances)
        reidcp = data[substances[i]]["ReidCp"]
        a[i] = reidcp[1]
        b[i] = reidcp[2]
        c[i] = reidcp[3]
        d[i] = reidcp[4]
        e[i] = reidcp[5]
    end

    return (a = a, b = b, c = c, d = d, e = e)
end

export load_component_properties, read_reidcp
