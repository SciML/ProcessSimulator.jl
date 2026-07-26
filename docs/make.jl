using ProcessSimulator
using Documenter

DocMeta.setdocmeta!(
    ProcessSimulator,
    :DocTestSetup,
    quote
        using ProcessSimulator: CSTR, MaterialSource, MaterialStream, Port, Reaction,
            SimpleAdiabaticCompressor, SimpleIsobaricHeatExchanger
    end;
    recursive = true
)

makedocs(;
    modules = [ProcessSimulator],
    authors = "SciML",
    sitename = "ProcessSimulator.jl",
    format = Documenter.HTML(;
        canonical = "https://docs.sciml.ai/ProcessSimulator/stable/",
        edit_link = "main",
        assets = String[]
    ),
    doctest = true,
    checkdocs = :public,
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo = "github.com/SciML/ProcessSimulator.jl",
    devbranch = "main"
)
