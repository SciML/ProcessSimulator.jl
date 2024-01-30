using ProcessSimulator
using Documenter

DocMeta.setdocmeta!(ProcessSimulator, :DocTestSetup, :(using ProcessSimulator); recursive=true)

makedocs(;
    modules=[ProcessSimulator],
    authors="Avinash Subramanian",
    sitename="ProcessSimulator.jl",
    format=Documenter.HTML(;
        canonical="https://avinashresearch1.github.io/ProcessSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/avinashresearch1/ProcessSimulator.jl",
    devbranch="main",
)
