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

Note that the kernel choices of DataInterpolations.jl, such as `CubicSpline()`,
are exact, i.e. go through the data points, while the smoothed kernels are
regression splines. Thus `CubicSpline()` is preferred if the data is not too
noisy or is relatively sparse. If data is sparse and very noisy, a `BSpline()`
can be the best regression spline, otherwise one of the other kernels such as as
`EpanechnikovKernel`.

