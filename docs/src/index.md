# DataCollocation.jl

**Non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives**

DataCollocation.jl provides non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives.

## Two Approaches for Data Collocation

DataCollocation.jl provides two different approaches for handling data collocation, each suited for different noise characteristics:

### 1. Kernel Smoothing Methods (For Noisy Data)
These methods use kernel-based smoothing to handle noisy measurements:
- Multiple kernel functions (Epanechnikov, Triangular, Gaussian, etc.)
- Robust to measurement noise and outliers
- Automatic bandwidth selection for optimal smoothing
- Best for: Experimental data with significant measurement noise

### 2. DataInterpolations.jl Approach (For Clean Data)  
Uses standard polynomial interpolations like cubic splines:
- Exact interpolation through data points (minimal noise assumption)
- Methods: CubicSpline, QuadraticInterpolation, BSpline, etc.
- Efficient and precise for clean data
- Best for: Simulation data or high-quality measurements with minimal noise

For datasets with significant noise, consider [NoiseRobustDifferentiation.jl](https://adrianhill.de/NoiseRobustDifferentiation.jl/dev/examples/) as an alternative specialized library for estimating derivatives from noisy data.

## Installation

Since this package is not yet registered, you can install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/ChrisRackauckas-Claude/DataCollocation.jl")
```

Once registered in the General registry:

```julia
using Pkg
Pkg.add("DataCollocation")
```

## Quick Example

```julia
using DataCollocation
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

# Method 1: Kernel smoothing (for noisy data)
du, u = collocate_data(data, tsteps, EpanechnikovKernel())

# Method 2: DataInterpolations approach (for clean data)
# using DataInterpolations
# du_interp, u_interp = collocate_data(data, tsteps, new_timepoints, CubicSpline)

# du contains estimated derivatives
# u contains smoothed data
```

## Citing DataCollocation.jl

If you use DataCollocation.jl in your research, please cite the collocation methodology paper:

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