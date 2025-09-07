# DataCollocations.jl

[![CI](https://github.com/SciML/DataCollocations.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SciML/DataCollocations.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/SciML/DataCollocations.jl/actions/workflows/documentation.yml/badge.svg)](https://sciml.github.io/DataCollocations.jl/)

DataCollocations.jl provides non-parametric data collocation functionality for smoothing timeseries data and estimating derivatives.

## Two Approaches for Data Collocation

DataCollocations.jl offers two distinct methodologies for data collocation, each optimized for different data characteristics:

### 1. Kernel Smoothing Methods (For Noisy Data)
**Robust regression-based approach for handling noisy measurements:**
- **Multiple kernel functions**: Epanechnikov, Triangular, Gaussian, Quartic, Triweight, Tricube, Cosine, Logistic, Sigmoid, Silverman
- **Automatic bandwidth selection**: Optimally balances bias-variance tradeoff
- **Noise robustness**: Designed to handle measurement noise and outliers
- **Regression splines**: Smoothed fit that doesn't necessarily pass through data points
- **Best for**: Experimental data with significant measurement noise

### 2. DataInterpolations.jl Integration (For Clean Data)
**Exact interpolation approach for high-quality data:**
- **Standard interpolation methods**: CubicSpline, QuadraticInterpolation, BSpline, Akima, etc.
- **Exact fitting**: Interpolation curves pass exactly through data points
- **Minimal noise assumption**: Assumes data points are accurate measurements
- **High efficiency**: Fast computation for clean, well-sampled data
- **Best for**: Simulation data or high-precision measurements with minimal noise

### When to Use Each Approach

| Data Characteristics | Recommended Method | Reason |
|---------------------|-------------------|---------|
| Experimental measurements with noise | Kernel smoothing | Robust to noise, provides smoothed estimates |
| Simulation results | DataInterpolations | Exact, efficient, preserves accuracy |
| Sparse, clean data | DataInterpolations (CubicSpline) | Exact interpolation between points |
| Dense, noisy data | Kernel smoothing (Epanechnikov) | Optimal noise handling |
| Very noisy data | [NoiseRobustDifferentiation.jl](https://adrianhill.de/NoiseRobustDifferentiation.jl/dev/examples/) | Specialized for heavy noise |

## Features

- Multiple kernel functions for data smoothing
- Automatic bandwidth selection
- Support for DataInterpolations.jl integration
- Derivative estimation from noisy data
- Efficient implementation with pre-allocated arrays

## Installation

```julia
using Pkg
Pkg.add("DataCollocations")
```

## Quick Start

```julia
using DataCollocations
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

DataCollocations.jl supports multiple kernel functions for noisy data:

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
- [NoiseRobustDifferentiation.jl](https://adrianhill.de/NoiseRobustDifferentiation.jl/dev/examples/) - Specialized library for estimating derivatives from very noisy data

## Cite Us

If you use DataCollocations.jl in your research, please cite the collocation methodology paper:

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