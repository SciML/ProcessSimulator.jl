using ProcessSimulator
using Documenter

DocMeta.setdocmeta!(ProcessSimulator, :DocTestSetup, :(using ProcessSimulator); recursive = true)

makedocs(;
    modules = [ProcessSimulator],
    authors = "SciML",
    repo = "https://github.com/SciML/ProcessSimulator.jl/blob/{commit}{path}#{line}",
    sitename = "ProcessSimulator.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/ProcessSimulator/stable/",
        edit_link = "main",
        assets = String[],
        size_threshold = 512000,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Media & Thermodynamics" => "guide/media.md",
            "Components" => "guide/components.md",
            "Reactors" => "guide/reactors.md",
            "Separation Units" => "guide/separation.md",
        ],
        "Examples" => [
            "Flash Drum" => "examples/flash_drum.md",
            "CSTR" => "examples/cstr.md",
        ],
        "API Reference" => [
            "Thermodynamics" => "api/thermodynamics.md",
            "Base Components" => "api/base_components.md",
            "Reactors" => "api/reactors.md",
            "Separation" => "api/separation.md",
        ],
    ]
)

deploydocs(;
    repo = "github.com/SciML/ProcessSimulator.jl",
    devbranch = "main",
    push_preview = true
)
