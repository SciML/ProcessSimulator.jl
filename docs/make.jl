using ProcessSimulator
using Documenter

DocMeta.setdocmeta!(
    ProcessSimulator, :DocTestSetup, :(using ProcessSimulator); recursive = true)

makedocs(;
    modules = [ProcessSimulator],
    authors = "SciML",
    sitename = "ProcessSimulator.jl",
    format = Documenter.HTML(;
        canonical = "https://docs.sciml.ai/ProcessSimulator/stable/",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(;
    repo = "github.com/SciML/ProcessSimulator.jl",
    devbranch = "main"
)
