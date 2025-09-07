# SmoothedCollocation.jl

[![CI](https://github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/ChrisRackauckas-Claude/SmoothedCollocation.jl/actions/workflows/documentation.yml/badge.svg)](https://chrismrackauckas-claude.github.io/SmoothedCollocation.jl/)

SmoothedCollocation.jl provides non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives. This package was extracted from DiffEqFlux.jl to provide a lightweight, standalone solution for data collocation tasks.

## Features

- Multiple kernel functions for data smoothing
- Automatic bandwidth selection
- Support for DataInterpolations.jl integration
- Derivative estimation from noisy data
- Efficient implementation with pre-allocated arrays

## Installation

```julia
using Pkg
Pkg.add("SmoothedCollocation")
```

## Quick Start

```julia
using SmoothedCollocation
using OrdinaryDiffEq

# Generate some sample data
f(u, p, t) = p .* u
prob = ODEProblem(f, [1.0], (0.0, 10.0), [-0.1])
t = collect(0.0:0.1:10.0)
data = Array(solve(prob, Tsit5(); saveat=t))

# Perform collocation to estimate derivatives and smooth data
u′, u = collocate_data(data, t, TriangularKernel(), 0.1)

# u′ contains the estimated derivatives
# u contains the smoothed data
```

## Available Kernels

SmoothedCollocation.jl supports multiple kernel functions:

**Bounded Support Kernels (support on [-1, 1]):**
- `EpanechnikovKernel()`
- `UniformKernel()`
- `TriangularKernel()` (default)
- `QuarticKernel()`
- `TriweightKernel()`
- `TricubeKernel()`
- `CosineKernel()`

**Unbounded Support Kernels:**
- `GaussianKernel()`
- `LogisticKernel()`
- `SigmoidKernel()`
- `SilvermanKernel()`

## DataInterpolations.jl Integration

With DataInterpolations.jl loaded, you can use interpolation methods:

```julia
using DataInterpolations

# Use interpolation to generate data at intermediate timepoints
tpoints_sample = 0.05:0.1:9.95
u′, u = collocate_data(data, t, tpoints_sample, LinearInterpolation)
```

## API Reference

### `collocate_data`

```julia
u′, u = collocate_data(data, tpoints, kernel=TriangularKernel(), bandwidth=nothing)
u′, u = collocate_data(data, tpoints, tpoints_sample, interp, args...)
```

**Arguments:**
- `data`: Matrix where each column is a snapshot of the timeseries
- `tpoints`: Time points corresponding to data columns
- `kernel`: Kernel function for smoothing (default: `TriangularKernel()`)
- `bandwidth`: Smoothing bandwidth (auto-selected if `nothing`)
- `tpoints_sample`: Sample points for interpolation method
- `interp`: Interpolation method from DataInterpolations.jl

**Returns:**
- `u′`: Estimated derivatives
- `u`: Smoothed data

## Citation

If you use SmoothedCollocation.jl in your research, please cite:

```bibtex
@article{rackauckas2020universal,
  title={Universal differential equations for scientific machine learning},
  author={Rackauckas, Christopher and Ma, Yingbo and Martensen, Julius and Warner, Collin and Zubov, Kirill and Supekar, Rohit and Skinner, Dominic and Ramadhan, Ali and Edelman, Alan},
  journal={arXiv preprint arXiv:2001.04385},
  year={2020}
}
```

## Contributing

Contributions are welcome! Please see the [contributing guidelines](CONTRIBUTING.md) for more information.

## Related Packages

- [DiffEqFlux.jl](https://github.com/SciML/DiffEqFlux.jl) - Neural differential equations
- [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) - Interpolation methods
- [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) - ODE solvers

## Reference

The collocation methodology is based on:

Roesch, E., Rackauckas, C., & Stumpf, M. P. H. (2021). Collocation based training of neural ordinary differential equations. *Statistical Applications in Genetics and Molecular Biology*, 20(2), 37-49. https://doi.org/10.1515/sagmb-2020-0025