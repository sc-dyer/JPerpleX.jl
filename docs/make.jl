using JPerplex
using Documenter

DocMeta.setdocmeta!(JPerplex, :DocTestSetup, :(using JPerplex); recursive=true)

makedocs(;
    modules=[JPerplex],
    authors="Sabastien Dyer <scdyer@uwaterloo.ca> and contributors",
    repo="https://github.com/sc-dyer/JPerplex.jl/blob/{commit}{path}#{line}",
    sitename="JPerplex.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sc-dyer.github.io/JPerplex.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sc-dyer/JPerplex.jl",
    devbranch="main",
)
