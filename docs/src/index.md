# SmoothedCollocation.jl

**Non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives**

SmoothedCollocation.jl provides non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives. This package was extracted from DiffEqFlux.jl to provide a lightweight, standalone solution for data collocation tasks.

## What is Smoothed Collocation?

Smoothed collocation, also referred to as the two-stage method, allows for fitting differential equations to time series data without relying on a numerical differential equation solver by building a smoothed collocating polynomial and using this to estimate the true `(u',u)` pairs, at which point `u'-f(u,p,t)` can be directly estimated as a loss to determine the correct parameters `p`. 

This method can be extremely fast and robust to noise, though, because it does not accumulate through time, is not as exact as other methods.

## Key Features

- **Multiple kernel functions** for data smoothing with both bounded and unbounded support
- **Automatic bandwidth selection** for optimal smoothing
- **DataInterpolations.jl integration** via package extensions
- **Derivative estimation** from noisy data
- **Efficient implementation** with pre-allocated arrays
- **Fast neural ODE training** through two-stage optimization

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

## Getting Started

- [Getting Started Tutorial](tutorials/getting_started.md) - Learn the basics of smoothed collocation
- [Kernel Selection Guide](tutorials/kernel_selection.md) - Choose the right kernel for your data
- [Neural ODE Training Example](examples/neural_ode_training.md) - See how to use collocation for fast neural ODE training

## Citing SmoothedCollocation.jl

If you use SmoothedCollocation.jl in your research, please cite:

```bibtex
@article{rackauckas2020universal,
  title={Universal differential equations for scientific machine learning},
  author={Rackauckas, Christopher and Ma, Yingbo and Martensen, Julius and Warner, Collin and Zubov, Kirill and Supekar, Rohit and Skinner, Dominic and Ramadhan, Ali and Edelman, Alan},
  journal={arXiv preprint arXiv:2001.04385},
  year={2020}
}
```