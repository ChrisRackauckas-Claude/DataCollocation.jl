# Smoothed Collocation

Smoothed collocation, also referred to as the two-stage method, allows
for fitting differential equations to time series data without relying
on a numerical differential equation solver by building a smoothed
collocating polynomial and using this to estimate the true `(u',u)`
pairs, at which point `u'-f(u,p,t)` can be directly estimated as a
loss to determine the correct parameters `p`. This method can be
extremely fast and robust to noise, though, because it does not
accumulate through time, is not as exact as other methods.

!!! note
    
    This is one of many methods for calculating the collocation coefficients
    for the training process. For a more comprehensive set of collocation
    methods, see [JuliaSimModelOptimizer](https://help.juliahub.com/jsmo/stable/manual/collocation/).

```@docs
collocate_data
```

## Kernel Choice

The smoothed kernels are regression splines that provide robust estimates for noisy data. Each kernel has different mathematical properties:

- **`EpanechnikovKernel()`**: Optimal in terms of mean squared error, recommended default for most noisy data applications
- **`TriangularKernel()`**: Simple linear weighting, fast computation, good general purpose choice
- **`UniformKernel()`**: Constant weighting within bandwidth, minimal smoothing assumptions
- **`GaussianKernel()`**: Infinitely smooth results, best for very smooth underlying functions
- **`QuarticKernel()`, `TriweightKernel()`, `TricubeKernel()`**: Higher-order smoothness, good for applications requiring smooth derivatives

**For noisy data**: Use `EpanechnikovKernel()` or `GaussianKernel()` with appropriate bandwidth.
**For fast computation**: Use `TriangularKernel()` or `UniformKernel()`.
**For very smooth results**: Use `GaussianKernel()` or higher-order kernels.

