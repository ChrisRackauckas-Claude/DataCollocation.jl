using Documenter
using SmoothedCollocation

makedocs(
    sitename = "SmoothedCollocation.jl",
    authors = "SciML Contributors", 
    modules = [SmoothedCollocation],
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Smoothed Collocation for Fast Two-Stage Training" => "examples/collocation.md",
        ],
        "API Reference" => [
            "Smoothed Collocation" => "utilities/Collocation.md",
        ],
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo = "github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl.git"
)