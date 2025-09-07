using Documenter
using DataCollocation

makedocs(
    sitename = "DataCollocation.jl",
    authors = "SciML Contributors", 
    modules = [DataCollocation],
    pages = [
        "Home" => "index.md",
        "Neural ODE Training with Kernel Smoothing" => "neural_ode_training.md",
        "DataInterpolations Methods for Clean Data" => "datainterpolations.md", 
        "Performance Optimization Techniques" => "optimization_tutorials.md",
        "API Reference" => "Collocation.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo = "github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl.git"
)