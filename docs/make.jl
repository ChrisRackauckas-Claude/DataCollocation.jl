using Documenter
using SmoothedCollocation

makedocs(
    sitename = "SmoothedCollocation.jl",
    authors = "SciML Contributors",
    modules = [SmoothedCollocation],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Kernel Selection Guide" => "tutorials/kernel_selection.md",
        ],
        "Examples" => [
            "Neural ODE Training" => "examples/neural_ode_training.md",
            "Noisy Data Smoothing" => "examples/noisy_data_smoothing.md",
            "DataInterpolations Integration" => "examples/data_interpolations.md",
        ],
        "API Reference" => [
            "Core Functions" => "api/core.md",
            "Kernels" => "api/kernels.md",
        ],
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

deploydocs(
    repo = "github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl.git"
)