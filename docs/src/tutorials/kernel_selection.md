# Kernel Selection Guide

Choosing the right kernel is crucial for optimal smoothing performance. This guide will help you understand the different kernels available and when to use each one.

## Available Kernels

SmoothedCollocation.jl provides 11 different kernel functions, divided into two categories:

### Bounded Support Kernels (Support on [-1, 1])

These kernels are zero outside the interval [-1, 1], which means they only use local data for smoothing.

#### `EpanechnikovKernel()` - **Recommended Default**
- **Mathematical optimality**: Optimal kernel in terms of mean squared error
- **Smooth results**: Good balance between smoothing and detail preservation
- **Computational efficiency**: Relatively fast to compute
- **Best for**: Most general applications, especially when you want optimal statistical properties

```julia
using SmoothedCollocation, Plots, Random
Random.seed!(123)

# Generate test data
t = 0:0.05:5
true_func = t -> sin(t) + 0.2*sin(5*t)
data = true_func.(t) + 0.1*randn(length(t))

# Apply Epanechnikov kernel
du, u = collocate_data(data, t, EpanechnikovKernel(), 0.2)
plot(t, u, label="Epanechnikov", lw=3)
```

#### `TriangularKernel()` 
- **Simple and intuitive**: Linear weight decay from center
- **Fast computation**: Very efficient
- **Moderate smoothing**: Good balance for most applications
- **Best for**: When you need fast computation and don't require optimal statistical properties

#### `UniformKernel()`
- **Simple**: Constant weight within bandwidth
- **Less smooth**: Can produce blocky results
- **Fast**: Very efficient to compute
- **Best for**: When you want minimal smoothing assumptions

#### `QuarticKernel()`, `TriweightKernel()`, `TricubeKernel()`
- **Higher-order smoothness**: More derivatives exist at boundaries
- **Smoother results**: Better for applications requiring smoothness
- **More computation**: Slightly more expensive
- **Best for**: When smoothness is critical

#### `CosineKernel()`
- **Smooth transitions**: Cosine-based weight function
- **Moderate computation**: Balanced performance
- **Best for**: When you want smooth, oscillation-free results

### Unbounded Support Kernels

These kernels have infinite support, meaning they use information from all data points (though with decreasing weight).

#### `GaussianKernel()`
- **Infinitely smooth**: All derivatives exist everywhere
- **Global influence**: Uses all data points
- **Excellent smoothing**: Very smooth results
- **Best for**: When you need the smoothest possible results and have dense data

```julia
# Compare Gaussian with bounded kernel
du_gauss, u_gauss = collocate_data(data, t, GaussianKernel(), 0.2)
du_tri, u_tri = collocate_data(data, t, TriangularKernel(), 0.2)

plot(t, u_gauss, label="Gaussian", lw=3)
plot!(t, u_tri, label="Triangular", lw=3)
scatter!(t, data, label="Data", alpha=0.4)
```

#### `LogisticKernel()`, `SigmoidKernel()`, `SilvermanKernel()`
- **Specialized applications**: Each has unique mathematical properties
- **Research use**: Often used in specific statistical contexts
- **Best for**: Specialized applications where these specific kernel shapes are desired

## Kernel Comparison

Let's compare different kernels on the same dataset:

```julia
using SmoothedCollocation, Plots, Random
Random.seed!(42)

# Create test data with different noise levels
t = 0:0.02:3
clean_data = sin.(2*t) .* exp.(-0.3*t)
noisy_data = clean_data + 0.15*randn(length(t))

# Test different kernels
kernels = [
    ("Epanechnikov", EpanechnikovKernel()),
    ("Triangular", TriangularKernel()), 
    ("Gaussian", GaussianKernel()),
    ("Uniform", UniformKernel()),
    ("Quartic", QuarticKernel())
]

bandwidth = 0.15
plots = []

for (name, kernel) in kernels
    du, u = collocate_data(noisy_data, t, kernel, bandwidth)
    p = plot(t, u, label="$name Kernel", lw=3)
    plot!(p, t, clean_data, label="True Function", lw=2, ls=:dash)
    scatter!(p, t[1:5:end], noisy_data[1:5:end], label="Noisy Data", alpha=0.6, ms=2)
    push!(plots, p)
end

plot(plots..., layout=(3,2), size=(800,900))
```

## Kernel Properties Summary

| Kernel | Support | Smoothness | Speed | Best Use Case |
|--------|---------|------------|-------|---------------|
| Epanechnikov | [-1,1] | High | Fast | General purpose, optimal MSE |
| Triangular | [-1,1] | Moderate | Very Fast | Quick computations |
| Gaussian | Unbounded | Infinite | Moderate | Maximum smoothness |
| Uniform | [-1,1] | Low | Very Fast | Minimal assumptions |
| Quartic/Triweight | [-1,1] | Very High | Moderate | High smoothness needs |

## Choosing the Right Kernel

### For most applications: `EpanechnikovKernel()`
This is the statistically optimal choice and provides the best trade-off between bias and variance.

### For maximum speed: `TriangularKernel()` or `UniformKernel()`
When computational efficiency is critical.

### For maximum smoothness: `GaussianKernel()`
When you need the smoothest possible results and have sufficient data density.

### For noisy data: `EpanechnikovKernel()` or `GaussianKernel()`
These provide better noise reduction.

### For sparse data: `GaussianKernel()`
The unbounded support helps when data points are far apart.

## Bandwidth Considerations by Kernel

Different kernels may require different bandwidth adjustments:

```julia
# Same bandwidth, different kernels
bandwidth = 0.2
kernels = [EpanechnikovKernel(), TriangularKernel(), GaussianKernel()]
names = ["Epanechnikov", "Triangular", "Gaussian"]

plot()
for (i, (kernel, name)) in enumerate(zip(kernels, names))
    du, u = collocate_data(noisy_data, t, kernel, bandwidth)
    plot!(t, u, label=name, lw=3)
end
scatter!(t[1:5:end], noisy_data[1:5:end], label="Data", alpha=0.6, ms=2)
plot!(title="Same Bandwidth, Different Kernels")
```

### Bandwidth Guidelines:
- **Gaussian kernels** typically need smaller bandwidths due to their unbounded support
- **Uniform kernels** may need larger bandwidths for adequate smoothing
- **Epanechnikov/Triangular** kernels work well with moderate bandwidths

## Advanced: Kernel Shape Visualization

Understanding kernel shapes helps with selection:

```julia
# Visualize kernel shapes
x = -2:0.01:2
kernels = [
    ("Epanechnikov", EpanechnikovKernel()),
    ("Triangular", TriangularKernel()),
    ("Gaussian", GaussianKernel()),
    ("Uniform", UniformKernel())
]

plot()
for (name, kernel) in kernels
    weights = [SmoothedCollocation.calckernel(kernel, xi) for xi in x]
    plot!(x, weights, label=name, lw=3)
end
plot!(title="Kernel Shapes", xlabel="Distance from Center", ylabel="Weight")
```

This visualization shows how each kernel weights nearby points, which directly affects the smoothing behavior.