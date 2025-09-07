# SmoothedCollocation.jl

**Non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives**

SmoothedCollocation.jl provides non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives. This package was extracted from DiffEqFlux.jl to provide a lightweight, standalone solution for data collocation tasks.

## Installation

Since this package is not yet registered, you can install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl")
```

Once registered in the General registry:

```julia
using Pkg
Pkg.add("SmoothedCollocation")
```

## Quick Example

```julia
using SmoothedCollocation
using OrdinaryDiffEq

# Generate some noisy data from an ODE
function trueODEfunc(du, u, p, t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u .^ 3)'true_A)'
end

u0 = [2.0; 0.0]
tspan = (0.0, 1.5)
prob = ODEProblem(trueODEfunc, u0, tspan)
tsteps = range(tspan[1], tspan[2]; length = 300)
data = Array(solve(prob, Tsit5(); saveat = tsteps)) .+ 0.1*randn(2, 300)

# Perform smoothed collocation
du, u = collocate_data(data, tsteps, EpanechnikovKernel())

# du contains estimated derivatives
# u contains smoothed data
```

## Citing SmoothedCollocation.jl

If you use SmoothedCollocation.jl in your research, please cite the collocation methodology paper:

```bibtex
@article{roesch2021collocation,
  title={Collocation based training of neural ordinary differential equations},
  author={Roesch, Elisabeth and Rackauckas, Christopher and Stumpf, Michael P. H.},
  journal={Statistical Applications in Genetics and Molecular Biology},
  volume={20},
  number={2},
  pages={37--49},
  year={2021},
  publisher={De Gruyter},
  doi={10.1515/sagmb-2020-0025},
  url={https://doi.org/10.1515/sagmb-2020-0025}
}
```