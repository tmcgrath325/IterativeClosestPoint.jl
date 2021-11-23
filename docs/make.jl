using IterativeClosestPoint
using Documenter

DocMeta.setdocmeta!(IterativeClosestPoint, :DocTestSetup, :(using IterativeClosestPoint); recursive=true)

makedocs(;
    modules=[IterativeClosestPoint],
    authors="Tom McGrath <tmcgrath325@gmail.com> and contributors",
    repo="https://github.com/tmcgrath325/IterativeClosestPoint.jl/blob/{commit}{path}#{line}",
    sitename="IterativeClosestPoint.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmcgrath325.github.io/IterativeClosestPoint.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tmcgrath325/IterativeClosestPoint.jl",
    devbranch="main",
)
