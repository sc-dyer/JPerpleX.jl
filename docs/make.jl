using JPerpleX
using Documenter

DocMeta.setdocmeta!(JPerpleX, :DocTestSetup, :(using JPerpleX); recursive=true)

makedocs(;
    modules=[JPerpleX],
    authors="Sabastien Dyer <scdyer@uwaterloo.ca> and contributors",
    repo="https://github.com/sc-dyer/JPerpleX.jl/blob/{commit}{path}#{line}",
    sitename="JPerpleX.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sc-dyer.github.io/JPerpleX.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sc-dyer/JPerpleX.jl",
    devbranch="main",
)
