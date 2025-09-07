using Documenter
using SmoothedCollocation

makedocs(
    sitename = "SmoothedCollocation.jl",
    authors = "SciML Contributors", 
    modules = [SmoothedCollocation],
    pages = [
        "Home" => "index.md",
        "Smoothed Collocation for Fast Two-Stage Training" => "collocation.md",
        "API Reference" => "Collocation.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo = "github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl.git"
)