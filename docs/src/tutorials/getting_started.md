# Getting Started with SmoothedCollocation.jl

This tutorial will walk you through the basics of using SmoothedCollocation.jl for data smoothing and derivative estimation.

## What is Smoothed Collocation?

Smoothed collocation is a technique used to:
1. **Smooth noisy time series data** using kernel-based methods
2. **Estimate derivatives** of the underlying smooth function
3. **Enable fast neural ODE training** by avoiding expensive ODE solvers during optimization

The key insight is that instead of solving differential equations numerically, we can estimate the derivative directly from smoothed data and then fit our model to match these derivative estimates.

## Basic Usage

Let's start with a simple example using synthetic data:

```julia
using SmoothedCollocation
using OrdinaryDiffEq, Plots
using Random
Random.seed!(123)

# Create some synthetic ODE data with noise
function lotka_volterra!(du, u, p, t)
    α, β, δ, γ = p
    du[1] = α*u[1] - β*u[1]*u[2]  # prey
    du[2] = -δ*u[2] + γ*u[1]*u[2]  # predator
end

u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [1.5, 1.0, 3.0, 1.0]
prob = ODEProblem(lotka_volterra!, u0, tspan, p)

# Generate data with noise
tsteps = range(0.0, 10.0, length=200)
sol = solve(prob, Tsit5(), saveat=tsteps)
data = Array(sol) + 0.1 * randn(2, length(tsteps))

# Plot the noisy data
scatter(tsteps, data[1, :], label="Noisy Prey", alpha=0.6)
scatter!(tsteps, data[2, :], label="Noisy Predator", alpha=0.6)
plot!(sol, vars=(0,1), label="True Prey", lw=2)
plot!(sol, vars=(0,2), label="True Predator", lw=2)
```

Now let's apply smoothed collocation:

```julia
# Apply smoothed collocation
du_est, u_smooth = collocate_data(data, tsteps, TriangularKernel(), 0.2)

# Plot the results
plot(tsteps, u_smooth[1, :], label="Smoothed Prey", lw=3)
plot!(tsteps, u_smooth[2, :], label="Smoothed Predator", lw=3)
scatter!(tsteps, data[1, :], label="Noisy Prey", alpha=0.4)
scatter!(tsteps, data[2, :], label="Noisy Predator", alpha=0.4)
plot!(sol, vars=(0,1), label="True Prey", lw=2, ls=:dash)
plot!(sol, vars=(0,2), label="True Predator", lw=2, ls=:dash)
```

And visualize the estimated derivatives:

```julia
# Plot estimated vs true derivatives
du_true = similar(data)
for i in 1:length(tsteps)
    lotka_volterra!(view(du_true, :, i), view(Array(sol), :, i), p, tsteps[i])
end

plot(tsteps, du_est[1, :], label="Estimated dPrey/dt", lw=3)
plot!(tsteps, du_est[2, :], label="Estimated dPredator/dt", lw=3)
plot!(tsteps, du_true[1, :], label="True dPrey/dt", lw=2, ls=:dash)
plot!(tsteps, du_true[2, :], label="True dPredator/dt", lw=2, ls=:dash)
```

## Understanding the Parameters

### Kernels
The choice of kernel affects how the smoothing is performed:
- `TriangularKernel()` - Good general-purpose choice, simple implementation
- `EpanechnikovKernel()` - Optimal in terms of mean squared error
- `GaussianKernel()` - Smooth results, unbounded support

### Bandwidth
The bandwidth parameter controls the amount of smoothing:
- **Smaller bandwidth**: Less smoothing, preserves more detail, may be noisier
- **Larger bandwidth**: More smoothing, removes noise but may lose important features

```julia
# Compare different bandwidths
bandwidths = [0.05, 0.2, 0.5]
plot()
for (i, bw) in enumerate(bandwidths)
    du_est, u_smooth = collocate_data(data, tsteps, TriangularKernel(), bw)
    plot!(tsteps, u_smooth[1, :], label="Bandwidth = $bw", lw=2)
end
scatter!(tsteps, data[1, :], label="Noisy Data", alpha=0.4, color=:black)
plot!(title="Effect of Bandwidth on Smoothing")
```

## Automatic Bandwidth Selection

If you don't specify a bandwidth, SmoothedCollocation.jl will automatically select one based on the data size:

```julia
# Let the algorithm choose the bandwidth automatically
du_est, u_smooth = collocate_data(data, tsteps, TriangularKernel())

plot(tsteps, u_smooth', lw=3, label=["Auto-smoothed Prey" "Auto-smoothed Predator"])
scatter!(tsteps, data', alpha=0.4, label=["Noisy Prey" "Noisy Predator"])
```

The automatic bandwidth selection uses the formula:
```
bandwidth = n^(-1/5) * n^(-3/35) * (log(n))^(-1/16)
```
where `n` is the number of data points.

## Working with Vector Data

For single-variable time series, you can work directly with vectors:

```julia
# Single variable example
t_single = 0:0.1:5
data_single = sin.(t_single) .+ 0.1*randn(length(t_single))

# collocate_data automatically handles vector input
du_single, u_single = collocate_data(data_single, t_single, EpanechnikovKernel(), 0.3)

plot(t_single, u_single, label="Smoothed", lw=3)
scatter!(t_single, data_single, label="Noisy Data", alpha=0.6)
plot!(t_single, sin.(t_single), label="True Function", lw=2, ls=:dash)
```

## Error Handling

If the bandwidth is too small, the algorithm may fail:

```julia
try
    du_est, u_smooth = collocate_data(data, tsteps, TriangularKernel(), 0.001)
catch e
    println("Error: ", e)
    println("Try using a larger bandwidth!")
end
```

## Next Steps

- Learn about [Kernel Selection](kernel_selection.md) to choose the best kernel for your data
- See [Neural ODE Training](../examples/neural_ode_training.md) for advanced applications
- Explore [DataInterpolations Integration](../examples/data_interpolations.md) for additional interpolation methods