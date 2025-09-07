# DataInterpolations.jl Integration

SmoothedCollocation.jl seamlessly integrates with [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) through package extensions, providing access to a wide variety of interpolation methods for generating intermediate data points.

## Basic DataInterpolations Usage

When DataInterpolations.jl is loaded, `collocate_data` gains additional functionality:

```julia
using SmoothedCollocation
using DataInterpolations
using Plots

# Generate some sparse data
t_sparse = [0.0, 0.5, 1.2, 2.1, 3.0, 4.5, 5.0]
data_sparse = sin.(t_sparse) + 0.1*randn(length(t_sparse))

# Define dense sampling points
t_dense = 0:0.1:5.0

# Use DataInterpolations for intermediate points
du_interp, u_interp = collocate_data(
    data_sparse, t_sparse, collect(t_dense), LinearInterpolation
)

# Plot results
scatter(t_sparse, data_sparse, label="Original Data", ms=6)
plot!(t_dense, u_interp, label="Linear Interpolation", lw=2)
plot!(xlabel="Time", ylabel="Value", title="DataInterpolations Integration")
```

## Available Interpolation Methods

DataInterpolations.jl provides many interpolation methods. Here are the most commonly used:

### Linear Interpolation
Simple and fast, good for data with linear trends:

```julia
# Linear interpolation - connects points with straight lines
du, u = collocate_data(data_sparse, t_sparse, collect(t_dense), LinearInterpolation)
```

### Cubic Splines
Smooth interpolation, good for smooth underlying functions:

```julia
# Cubic spline - smooth curves through data points
du_cubic, u_cubic = collocate_data(data_sparse, t_sparse, collect(t_dense), CubicSpline)
```

### Quadratic Interpolation
Balance between linear and cubic:

```julia
# Quadratic interpolation
du_quad, u_quad = collocate_data(data_sparse, t_sparse, collect(t_dense), QuadraticInterpolation)
```

### B-Splines
Flexible splines with controlled smoothness:

```julia
# B-spline with degree 3 (cubic B-spline)
du_bspline, u_bspline = collocate_data(data_sparse, t_sparse, collect(t_dense), 
                                       BSplineInterpolation(3, 2, :uniform, :uniform))
```

## Comprehensive Comparison

Let's compare different interpolation methods on the same dataset:

```julia
using SmoothedCollocation, DataInterpolations, Plots
using Random
Random.seed!(42)

# Create test data with some curvature
t_original = 0:0.3:4π
true_func = t -> sin(t) * exp(-0.1*t) + 0.2*cos(3*t)
data_original = true_func.(t_original) + 0.05*randn(length(t_original))

# Dense evaluation points
t_eval = 0:0.05:4π
true_values = true_func.(t_eval)

# Test different interpolation methods
methods = [
    ("Linear", LinearInterpolation),
    ("Quadratic", QuadraticInterpolation),
    ("Cubic Spline", CubicSpline),
    ("Akima", AkimaInterpolation),
]

plots = []
for (name, method) in methods
    try
        du, u = collocate_data(data_original, t_original, collect(t_eval), method)
        
        p = plot(t_eval, true_values, label="True Function", lw=2, color=:blue, ls=:dash)
        scatter!(p, t_original, data_original, label="Data Points", ms=4, alpha=0.7)
        plot!(p, t_eval, u, label="$name", lw=2, color=:red)
        plot!(p, title="$name Interpolation")
        
        # Calculate and display error
        rmse = sqrt(mean((u - true_values).^2))
        plot!(p, title="$name (RMSE: $(round(rmse, digits=4)))")
        
        push!(plots, p)
    catch e
        println("Method $name failed: $e")
    end
end

plot(plots..., layout=(2,2), size=(800,600))
```

## Working with Multidimensional Data

DataInterpolations works seamlessly with multidimensional data:

```julia
# 2D system example - damped oscillator
t_2d_sparse = 0:0.5:3π
true_x = exp.(-0.2*t_2d_sparse) .* cos.(t_2d_sparse)
true_y = exp.(-0.2*t_2d_sparse) .* sin.(t_2d_sparse)
data_2d = [true_x + 0.05*randn(length(t_2d_sparse)), 
           true_y + 0.05*randn(length(t_2d_sparse))]
data_matrix_2d = hcat(data_2d...)'

# Dense evaluation
t_2d_dense = 0:0.1:3π

# Apply interpolation
du_2d, u_2d = collocate_data(data_matrix_2d, t_2d_sparse, collect(t_2d_dense), CubicSpline)

# Plot phase portrait
scatter(data_2d[1], data_2d[2], label="Sparse Data", ms=6, alpha=0.7)
plot!(u_2d[1,:], u_2d[2,:], label="Interpolated Trajectory", lw=2)
plot!(xlabel="x", ylabel="y", title="2D Interpolation", aspect_ratio=:equal)

# Time series plots
p1 = scatter(t_2d_sparse, data_2d[1], label="Sparse x", ms=4)
plot!(p1, t_2d_dense, u_2d[1,:], label="Interpolated x", lw=2)

p2 = scatter(t_2d_sparse, data_2d[2], label="Sparse y", ms=4)
plot!(p2, t_2d_dense, u_2d[2,:], label="Interpolated y", lw=2)

plot(p1, p2, layout=(2,1), size=(800,500))
```

## Advanced: Custom Interpolation Parameters

Many interpolation methods accept additional parameters:

```julia
# B-spline with custom parameters
# BSplineInterpolation(degree, order, boundary_conditions)
methods_bspline = [
    ("B-spline (deg=1)", BSplineInterpolation(1, 2, :uniform, :uniform)),
    ("B-spline (deg=2)", BSplineInterpolation(2, 2, :uniform, :uniform)),
    ("B-spline (deg=3)", BSplineInterpolation(3, 2, :uniform, :uniform)),
]

# Test with same data
t_test = 0:0.4:2π
data_test = sin.(t_test) + 0.1*randn(length(t_test))
t_eval_test = 0:0.1:2π

plots_bspline = []
for (name, method) in methods_bspline
    du, u = collocate_data(data_test, t_test, collect(t_eval_test), method)
    
    p = plot(t_eval_test, sin.(t_eval_test), label="True", lw=2, ls=:dash)
    scatter!(p, t_test, data_test, label="Data", ms=4)
    plot!(p, t_eval_test, u, label=name, lw=2)
    plot!(p, title=name)
    
    push!(plots_bspline, p)
end

plot(plots_bspline..., layout=(3,1), size=(800,700))
```

## Practical Applications

### 1. Upsampling Time Series Data

Convert low-frequency measurements to high-frequency estimates:

```julia
# Simulate hourly temperature measurements
hours_measured = 0:3:24  # Every 3 hours
temp_measured = 20 .+ 5*sin.(2π*hours_measured/24) + randn(length(hours_measured))

# Upsample to every 15 minutes
minutes_dense = 0:0.25:24
du_temp, temp_dense = collocate_data(temp_measured, hours_measured, 
                                    minutes_dense, CubicSpline)

scatter(hours_measured, temp_measured, label="3-hour measurements", ms=6)
plot!(minutes_dense, temp_dense, label="15-min interpolation", lw=2)
plot!(xlabel="Hour", ylabel="Temperature (°C)", title="Temperature Upsampling")
```

### 2. Irregular to Regular Grid Conversion

Convert irregularly spaced measurements to regular grid:

```julia
# Irregular measurement times (simulating sensor data with missed readings)
t_irregular = sort([0; 0.5; 1.2; 1.8; 2.9; 3.1; 4.2; 4.8; 5.0])
data_irregular = exp.(-0.3*t_irregular) .* sin.(2*t_irregular) + 0.1*randn(length(t_irregular))

# Regular grid
t_regular = 0:0.2:5.0
du_reg, u_reg = collocate_data(data_irregular, t_irregular, collect(t_regular), 
                               AkimaInterpolation())

scatter(t_irregular, data_irregular, label="Irregular Data", ms=6)
plot!(t_regular, u_reg, label="Regular Grid", lw=2, marker=:circle, ms=2)
plot!(xlabel="Time", ylabel="Value", title="Irregular to Regular Grid")
```

### 3. Derivative Estimation on Interpolated Data

Use interpolation to get derivatives at arbitrary points:

```julia
# Sparse derivative estimation
t_sparse_deriv = [0, 1, 2.5, 4, 6]
data_sparse_deriv = t_sparse_deriv.^2 .* exp.(-0.2*t_sparse_deriv) + 0.1*randn(length(t_sparse_deriv))

# Evaluate derivatives at many points
t_eval_deriv = 0:0.1:6
du_eval, u_eval = collocate_data(data_sparse_deriv, t_sparse_deriv, 
                                 collect(t_eval_deriv), CubicSpline)

# True derivative for comparison
true_deriv = (2*t_eval_deriv - 0.2*t_eval_deriv.^2) .* exp.(-0.2*t_eval_deriv)

scatter(t_sparse_deriv, data_sparse_deriv, label="Data Points", ms=6)
plot!(t_eval_deriv, u_eval, label="Interpolated Function", lw=2)
plot!(t_eval_deriv, du_eval, label="Estimated Derivative", lw=2)
plot!(t_eval_deriv, true_deriv, label="True Derivative", lw=2, ls=:dash)
plot!(xlabel="Time", ylabel="Value", title="Derivative Estimation via Interpolation")
```

## Performance Considerations

### Speed Comparison

```julia
using BenchmarkTools

# Setup test data
n_points = 100
t_bench = sort(rand(n_points) * 10)
data_bench = sin.(t_bench) + 0.1*randn(n_points)
t_eval_bench = 0:0.1:10

# Benchmark different methods
methods_bench = [
    LinearInterpolation,
    QuadraticInterpolation,
    CubicSpline,
    AkimaInterpolation
]

println("Interpolation Method Benchmarks:")
for method in methods_bench
    time = @belapsed collocate_data($data_bench, $t_bench, $t_eval_bench, $method)
    println("$(string(method)): $(round(time*1000, digits=2)) ms")
end
```

### Memory Efficiency

For large datasets, consider memory usage:

```julia
# For very large datasets, process in chunks
function chunk_interpolate(data, t_data, t_eval, method; chunk_size=1000)
    n_chunks = div(length(t_eval), chunk_size) + 1
    results = []
    
    for i in 1:n_chunks
        start_idx = (i-1) * chunk_size + 1
        end_idx = min(i * chunk_size, length(t_eval))
        
        if start_idx <= length(t_eval)
            chunk_t = t_eval[start_idx:end_idx]
            _, u_chunk = collocate_data(data, t_data, chunk_t, method)
            push!(results, u_chunk)
        end
    end
    
    return vcat(results...)
end
```

## Error Handling and Edge Cases

```julia
# Handle edge cases gracefully
function safe_interpolate(data, t_data, t_eval, method)
    try
        return collocate_data(data, t_data, t_eval, method)
    catch e
        if isa(e, ArgumentError)
            println("Interpolation failed: $e")
            println("Falling back to linear interpolation")
            return collocate_data(data, t_data, t_eval, LinearInterpolation)
        else
            rethrow(e)
        end
    end
end

# Example usage
t_problematic = [0, 0, 1, 2]  # Duplicate points
data_problematic = [1, 1.1, 2, 3]
t_eval_safe = 0:0.1:2

du_safe, u_safe = safe_interpolate(data_problematic, t_problematic, t_eval_safe, CubicSpline)
```

The integration with DataInterpolations.jl makes SmoothedCollocation.jl extremely versatile for handling various data interpolation and smoothing tasks across different domains and applications.