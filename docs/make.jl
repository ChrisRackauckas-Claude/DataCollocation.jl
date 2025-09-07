using Documenter
using DataCollocation

makedocs(
    sitename = "DataCollocation.jl",
    authors = "SciML Contributors", 
    modules = [DataCollocation],
    pages = [
        "Home" => "index.md",
        "Kernel Smoothing for Neural ODE Training" => "collocation.md",
        "DataInterpolations Methods for Clean Data" => "datainterpolations.md",
        "API Reference" => "Collocation.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo = "github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl.git"
)