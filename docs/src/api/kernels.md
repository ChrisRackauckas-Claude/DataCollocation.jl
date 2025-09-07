# Kernel Functions

SmoothedCollocation.jl provides 11 different kernel functions for data smoothing. Each kernel has different mathematical properties and is suitable for different applications.

## Kernel Types

All kernels are subtypes of `CollocationKernel`:

```julia
abstract type CollocationKernel end
```

## Bounded Support Kernels

These kernels have support on the interval [-1, 1], meaning they are zero outside this range.

### `EpanechnikovKernel`

```@docs
EpanechnikovKernel
```

The Epanechnikov kernel is optimal in terms of minimizing mean squared error and is the recommended default choice.

**Mathematical Form:**
```
K(t) = 3/4 * (1 - t²)  for |t| ≤ 1
K(t) = 0               for |t| > 1
```

**Properties:**
- Optimal statistical properties
- Smooth at the center
- Bounded support
- Good balance of bias and variance

**Example:**
```julia
using SmoothedCollocation, Plots

# Visualize kernel shape
t = -1.5:0.01:1.5
kernel = EpanechnikovKernel()
weights = [SmoothedCollocation.calckernel(kernel, ti) for ti in t]

plot(t, weights, label="Epanechnikov Kernel", lw=3)
plot!(xlim=(-1.5, 1.5), xlabel="Distance", ylabel="Weight")
```

### `TriangularKernel`

```@docs
TriangularKernel
```

Simple triangular weighting function.

**Mathematical Form:**
```
K(t) = 1 - |t|  for |t| ≤ 1
K(t) = 0        for |t| > 1
```

**Properties:**
- Simple and intuitive
- Fast computation
- Linear weight decay
- Good general-purpose choice

### `UniformKernel`

```@docs
UniformKernel
```

Constant weighting within the bandwidth.

**Mathematical Form:**
```
K(t) = 1/2  for |t| ≤ 1
K(t) = 0    for |t| > 1
```

**Properties:**
- Simplest kernel
- Fastest computation
- Can produce less smooth results
- Minimal smoothing assumptions

### `QuarticKernel`

```@docs
QuarticKernel
```

Fourth-order polynomial kernel.

**Mathematical Form:**
```
K(t) = 15/16 * (1 - t²)²  for |t| ≤ 1
K(t) = 0                  for |t| > 1
```

**Properties:**
- Higher-order smoothness
- Good for applications requiring smoothness
- Moderate computational cost

### `TriweightKernel`

```@docs
TriweightKernel
```

Sixth-order polynomial kernel.

**Mathematical Form:**
```
K(t) = 35/32 * (1 - t²)³  for |t| ≤ 1
K(t) = 0                  for |t| > 1
```

**Properties:**
- Very smooth results
- Higher-order derivatives exist
- More expensive computation

### `TricubeKernel`

```@docs
TricubeKernel
```

Tricube kernel with cubic terms.

**Mathematical Form:**
```
K(t) = 70/81 * (1 - |t|³)³  for |t| ≤ 1
K(t) = 0                    for |t| > 1
```

**Properties:**
- Smooth transitions
- Robust to outliers
- Popular in robust regression

### `CosineKernel`

```@docs
CosineKernel
```

Cosine-based kernel function.

**Mathematical Form:**
```
K(t) = π/4 * cos(πt/2)  for |t| ≤ 1
K(t) = 0                for |t| > 1
```

**Properties:**
- Smooth, oscillation-free
- Good for periodic-like data
- Moderate computational cost

## Unbounded Support Kernels

These kernels have infinite support, using information from all data points with decreasing weights.

### `GaussianKernel`

```@docs
GaussianKernel
```

Standard normal kernel with infinite support.

**Mathematical Form:**
```
K(t) = 1/√(2π) * exp(-t²/2)  for all t
```

**Properties:**
- Infinitely smooth (all derivatives exist)
- Global influence
- Excellent smoothing properties
- More expensive computation

**Example:**
```julia
# Compare Gaussian with bounded kernel
t = -3:0.01:3
gaussian_weights = [SmoothedCollocation.calckernel(GaussianKernel(), ti) for ti in t]
epan_weights = [SmoothedCollocation.calckernel(EpanechnikovKernel(), ti) for ti in t]

plot(t, gaussian_weights, label="Gaussian", lw=3)
plot!(t, epan_weights, label="Epanechnikov", lw=3)
plot!(xlabel="Distance", ylabel="Weight", title="Bounded vs Unbounded Kernels")
```

### `LogisticKernel`

```@docs
LogisticKernel
```

Logistic function kernel.

**Mathematical Form:**
```
K(t) = 1/(exp(t) + 2 + exp(-t))  for all t
```

**Properties:**
- S-shaped weight function
- Unbounded support
- Moderate smoothness

### `SigmoidKernel`

```@docs
SigmoidKernel
```

Sigmoid-based kernel function.

**Mathematical Form:**
```
K(t) = 2/(π(exp(t) + exp(-t)))  for all t
```

**Properties:**
- Sigmoid-shaped weights
- Unbounded support
- Used in specialized applications

### `SilvermanKernel`

```@docs
SilvermanKernel
```

Silverman kernel with exponential decay.

**Mathematical Form:**
```
K(t) = 1/2 * exp(-|t|/√2) * sin(|t|/2 + π/4)  for all t
```

**Properties:**
- Oscillatory weights
- Specialized statistical properties
- Used in density estimation

## Kernel Comparison

### Performance Characteristics

| Kernel | Support | Computation | Smoothness | Best Use |
|--------|---------|-------------|------------|----------|
| Uniform | [-1,1] | Fastest | Low | Quick estimates |
| Triangular | [-1,1] | Very Fast | Moderate | General purpose |
| Epanechnikov | [-1,1] | Fast | High | Optimal MSE |
| Quartic | [-1,1] | Moderate | Very High | Smooth results |
| Gaussian | Unbounded | Slower | Infinite | Maximum smoothness |

### Visual Comparison

```julia
using SmoothedCollocation, Plots

# Compare kernel shapes
x = -2:0.01:2
kernels = [
    ("Uniform", UniformKernel()),
    ("Triangular", TriangularKernel()),
    ("Epanechnikov", EpanechnikovKernel()),
    ("Quartic", QuarticKernel()),
    ("Gaussian", GaussianKernel())
]

plot()
for (name, kernel) in kernels
    weights = [SmoothedCollocation.calckernel(kernel, xi) for xi in x]
    plot!(x, weights, label=name, lw=3)
end
plot!(xlabel="Distance from Center", ylabel="Weight", 
      title="Kernel Function Comparison")
```

### Bandwidth Sensitivity

Different kernels may require different bandwidth adjustments:

```julia
# Test kernel sensitivity to bandwidth
using Random
Random.seed!(42)

t_test = 0:0.1:2π
true_func = sin.(t_test)
noisy_data = true_func + 0.2*randn(length(t_test))

kernels_test = [TriangularKernel(), EpanechnikovKernel(), GaussianKernel()]
kernel_names = ["Triangular", "Epanechnikov", "Gaussian"]
bandwidths = [0.1, 0.2, 0.4]

# Plot results for different bandwidth/kernel combinations
plots = []
for (i, (kernel, name)) in enumerate(zip(kernels_test, kernel_names))
    p = plot(title=name)
    scatter!(p, t_test[1:3:end], noisy_data[1:3:end], alpha=0.4, label="Data")
    plot!(p, t_test, true_func, label="True", lw=2, ls=:dash, color=:black)
    
    for bw in bandwidths
        _, u_smooth = collocate_data(noisy_data, t_test, kernel, bw)
        plot!(p, t_test, u_smooth, label="bw=$bw", lw=2)
    end
    push!(plots, p)
end

plot(plots..., layout=(3,1), size=(800,900))
```

## Custom Kernel Implementation

Advanced users can implement custom kernels by subtyping `CollocationKernel`:

```julia
# Example: Custom exponential kernel
struct ExponentialKernel <: SmoothedCollocation.CollocationKernel end

function SmoothedCollocation.calckernel(::ExponentialKernel, t::T) where T
    return exp(-abs(t))  # Unbounded support
end

# For bounded support kernels, implement the bounded version:
function SmoothedCollocation.calckernel(::ExponentialKernel, t::T, abst::T) where T
    return exp(-abst)  # Only called when |t| ≤ 1
end

# Usage
custom_kernel = ExponentialKernel()
du, u = collocate_data(data, time_points, custom_kernel, 0.2)
```

## Kernel Selection Guidelines

1. **Default choice**: Use `EpanechnikovKernel()` for most applications
2. **Speed priority**: Use `TriangularKernel()` or `UniformKernel()` 
3. **Maximum smoothness**: Use `GaussianKernel()`
4. **Noisy data**: Use `EpanechnikovKernel()` or `GaussianKernel()`
5. **Sparse data**: Use `GaussianKernel()` (unbounded support helps)
6. **Outlier robust**: Consider `TricubeKernel()`

The choice of kernel significantly affects the smoothing behavior, so experiment with different options to find what works best for your specific application.