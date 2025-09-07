using Documenter
using DataCollocations

makedocs(
    sitename = "DataCollocations.jl",
    authors = "SciML Contributors", 
    modules = [DataCollocations],
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
    repo = "github.com/SciML/DataCollocations.jl.git"
)