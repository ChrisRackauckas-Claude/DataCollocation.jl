# Noisy Data Smoothing Examples

This page demonstrates various applications of SmoothedCollocation.jl for cleaning and analyzing noisy time series data.

## Basic Smoothing Example

Let's start with a simple sinusoidal signal corrupted by noise:

```julia
using SmoothedCollocation, Plots, Random
Random.seed!(42)

# Create a noisy signal
t = 0:0.02:4π
clean_signal = sin.(t) + 0.5*sin.(3*t)
noisy_signal = clean_signal + 0.3*randn(length(t))

# Apply smoothed collocation
du, u_smooth = collocate_data(noisy_signal, t, EpanechnikovKernel(), 0.15)

# Plot results
plot(t, clean_signal, label="True Signal", lw=3, color=:blue)
scatter!(t[1:5:end], noisy_signal[1:5:end], label="Noisy Data", alpha=0.6, color=:red, ms=2)
plot!(t, u_smooth, label="Smoothed Signal", lw=2, color=:green)
plot!(title="Basic Signal Smoothing", xlabel="Time", ylabel="Amplitude")
```

## Comparing Different Kernels for Noisy Data

Different kernels perform differently depending on the noise characteristics:

```julia
# Create test signal with different noise levels
t = 0:0.01:2π
true_signal = exp.(-0.3*t) .* sin.(2*t)

# Light noise
light_noise = true_signal + 0.05*randn(length(t))

# Heavy noise  
heavy_noise = true_signal + 0.2*randn(length(t))

# Test different kernels
kernels = [
    ("Triangular", TriangularKernel()),
    ("Epanechnikov", EpanechnikovKernel()),
    ("Gaussian", GaussianKernel()),
    ("Quartic", QuarticKernel())
]

bandwidth = 0.1
plots_light = []
plots_heavy = []

# Light noise comparison
for (name, kernel) in kernels
    _, u_smooth = collocate_data(light_noise, t, kernel, bandwidth)
    p = plot(t, u_smooth, label="$name", lw=2)
    plot!(p, t, true_signal, label="True", lw=1, ls=:dash, color=:black)
    scatter!(p, t[1:10:end], light_noise[1:10:end], alpha=0.4, ms=1, label="Noisy")
    plot!(p, title="Light Noise - $name Kernel")
    push!(plots_light, p)
end

# Heavy noise comparison
for (name, kernel) in kernels
    _, u_smooth = collocate_data(heavy_noise, t, kernel, bandwidth)
    p = plot(t, u_smooth, label="$name", lw=2)
    plot!(p, t, true_signal, label="True", lw=1, ls=:dash, color=:black)
    scatter!(p, t[1:10:end], heavy_noise[1:10:end], alpha=0.4, ms=1, label="Noisy")
    plot!(p, title="Heavy Noise - $name Kernel")
    push!(plots_heavy, p)
end

plot(plots_light..., layout=(2,2), size=(800,600))
plot(plots_heavy..., layout=(2,2), size=(800,600))
```

## Bandwidth Selection for Different Noise Levels

The optimal bandwidth depends on the noise level:

```julia
# Create signals with varying noise levels
noise_levels = [0.05, 0.1, 0.2, 0.4]
bandwidths = [0.05, 0.1, 0.2, 0.3]

# Performance metric: Mean Squared Error
function mse(predicted, true_values)
    return mean((predicted - true_values).^2)
end

# Test matrix
results = zeros(length(noise_levels), length(bandwidths))

for (i, noise_level) in enumerate(noise_levels)
    noisy_data = true_signal + noise_level * randn(length(t))
    
    for (j, bw) in enumerate(bandwidths)
        _, u_smooth = collocate_data(noisy_data, t, EpanechnikovKernel(), bw)
        results[i, j] = mse(u_smooth, true_signal)
    end
end

# Plot results
heatmap(bandwidths, noise_levels, results, 
        xlabel="Bandwidth", ylabel="Noise Level", 
        title="MSE vs Bandwidth and Noise Level",
        color=:viridis)
```

## Derivative Estimation from Noisy Data

Estimating derivatives from noisy data is challenging, but smoothed collocation handles it well:

```julia
# Create a function with known derivative
t_deriv = 0:0.01:3
true_func = t_deriv.^2 .* exp.(-0.5*t_deriv)
true_derivative = (2*t_deriv - 0.5*t_deriv.^2) .* exp.(-0.5*t_deriv)

# Add noise
noisy_func = true_func + 0.1*randn(length(t_deriv))

# Estimate derivative using smoothed collocation
du_est, u_smooth = collocate_data(noisy_func, t_deriv, EpanechnikovKernel(), 0.08)

# Compare different bandwidth choices for derivative estimation
bandwidths_deriv = [0.03, 0.08, 0.15]
plot_derivs = []

for bw in bandwidths_deriv
    du_est, u_smooth = collocate_data(noisy_func, t_deriv, EpanechnikovKernel(), bw)
    
    p = plot(t_deriv, true_derivative, label="True Derivative", lw=3, color=:blue)
    plot!(p, t_deriv, du_est, label="Estimated (bw=$bw)", lw=2)
    scatter!(p, t_deriv[1:10:end], noisy_func[1:10:end], label="Noisy Data", 
             alpha=0.4, ms=2, color=:red)
    plot!(p, title="Derivative Estimation (bandwidth = $bw)")
    push!(plot_derivs, p)
end

plot(plot_derivs..., layout=(3,1), size=(800,900))
```

## Real-World Example: Temperature Data

Let's simulate processing temperature sensor data:

```julia
# Simulate temperature sensor data with measurement noise
hours = 0:0.25:24  # Every 15 minutes for 24 hours
base_temp = 20     # Base temperature

# Simulate daily temperature cycle with some random variation
true_temperature = base_temp .+ 8*sin.(2π*hours/24 - π/2) .+ 
                   2*sin.(4π*hours/24) .+ 
                   1*sin.(6π*hours/24 - π/3)

# Add sensor noise and occasional outliers
sensor_noise = 0.5*randn(length(hours))
outliers = zeros(length(hours))
outliers[rand(1:length(hours), 5)] .= 5*randn(5)  # 5 random outliers

measured_temp = true_temperature + sensor_noise + outliers

# Apply smoothed collocation
du_temp, smoothed_temp = collocate_data(measured_temp, hours, 
                                       EpanechnikovKernel(), 0.8)

# Plot results
scatter(hours, measured_temp, label="Sensor Readings", alpha=0.6, ms=2)
plot!(hours, true_temperature, label="True Temperature", lw=3, color=:blue)
plot!(hours, smoothed_temp, label="Smoothed Temperature", lw=2, color=:red)
plot!(xlabel="Hour of Day", ylabel="Temperature (°C)", 
      title="Temperature Sensor Data Smoothing")

# Plot temperature rate of change
plot(hours, du_temp, label="Temperature Rate of Change", lw=2, color=:green)
plot!(xlabel="Hour of Day", ylabel="°C per hour", 
      title="Estimated Temperature Rate of Change")
```

## Multi-dimensional Smoothing: 2D Oscillator

For systems with multiple variables:

```julia
# Create a 2D damped oscillator with noise
t_2d = 0:0.02:4π
ω = 1.5  # frequency
γ = 0.1  # damping

# True solution
true_x = exp.(-γ*t_2d) .* cos.(ω*t_2d)
true_y = exp.(-γ*t_2d) .* sin.(ω*t_2d)

# Add correlated noise
noise_x = 0.1*randn(length(t_2d))
noise_y = 0.1*randn(length(t_2d)) + 0.3*noise_x  # Correlated noise

noisy_data = [true_x + noise_x, true_y + noise_y]
data_matrix = hcat(noisy_data...)'

# Apply smoothed collocation
du_2d, u_smooth_2d = collocate_data(data_matrix, t_2d, EpanechnikovKernel(), 0.15)

# Plot phase portrait
scatter(noisy_data[1], noisy_data[2], label="Noisy Data", alpha=0.4, ms=1)
plot!(true_x, true_y, label="True Trajectory", lw=3, color=:blue)
plot!(u_smooth_2d[1,:], u_smooth_2d[2,:], label="Smoothed Trajectory", lw=2, color=:red)
plot!(xlabel="x", ylabel="y", title="2D Oscillator Phase Portrait", aspect_ratio=:equal)

# Time series comparison
p1 = plot(t_2d, true_x, label="True x", lw=2, color=:blue)
scatter!(p1, t_2d[1:5:end], noisy_data[1][1:5:end], label="Noisy x", alpha=0.6, ms=2)
plot!(p1, t_2d, u_smooth_2d[1,:], label="Smoothed x", lw=2, color=:red)

p2 = plot(t_2d, true_y, label="True y", lw=2, color=:blue)
scatter!(p2, t_2d[1:5:end], noisy_data[2][1:5:end], label="Noisy y", alpha=0.6, ms=2)
plot!(p2, t_2d, u_smooth_2d[2,:], label="Smoothed y", lw=2, color=:red)

plot(p1, p2, layout=(2,1), size=(800,600))
```

## Performance Tips for Noisy Data

### 1. Bandwidth Selection Guidelines

```julia
# Rule of thumb for bandwidth selection based on noise level
function suggest_bandwidth(data, noise_estimate)
    n = length(data)
    data_range = maximum(data) - minimum(data)
    noise_ratio = noise_estimate / data_range
    
    # Base bandwidth from automatic selection
    base_bw = n^(-1/5) * n^(-3/35) * (log(n))^(-1/16)
    
    # Adjust based on noise level
    if noise_ratio < 0.05
        return base_bw * 0.5      # Less smoothing for clean data
    elseif noise_ratio > 0.2
        return base_bw * 2.0      # More smoothing for very noisy data
    else
        return base_bw            # Standard smoothing
    end
end

# Example usage
estimated_noise = std(diff(noisy_signal))  # Rough noise estimate
suggested_bw = suggest_bandwidth(noisy_signal, estimated_noise)
println("Suggested bandwidth: $suggested_bw")
```

### 2. Outlier Handling

```julia
# Smoothed collocation is robust to outliers, but extreme outliers can be pre-filtered
function filter_outliers(data, threshold=3.0)
    μ = median(data)
    σ = 1.4826 * median(abs.(data .- μ))  # Robust std estimate
    
    mask = abs.(data .- μ) .< threshold * σ
    return data[mask], mask
end

# Apply outlier filtering before smoothing for extreme cases
filtered_data, mask = filter_outliers(measured_temp, 2.5)
filtered_hours = hours[mask]

du_filtered, u_filtered = collocate_data(filtered_data, filtered_hours, 
                                        EpanechnikovKernel(), 0.5)
```

### 3. Cross-Validation for Bandwidth Selection

```julia
function cross_validate_bandwidth(data, t, kernel, bandwidths; k_folds=5)
    n = length(data)
    fold_size = div(n, k_folds)
    cv_errors = zeros(length(bandwidths))
    
    for (i, bw) in enumerate(bandwidths)
        fold_errors = []
        
        for fold in 1:k_folds
            # Create train/test split
            start_idx = (fold - 1) * fold_size + 1
            end_idx = min(fold * fold_size, n)
            
            test_indices = start_idx:end_idx
            train_indices = setdiff(1:n, test_indices)
            
            # Train on subset
            _, u_smooth = collocate_data(data[train_indices], t[train_indices], kernel, bw)
            
            # Simple interpolation for test points (basic validation)
            # In practice, you might use more sophisticated validation
            if length(u_smooth) > 0
                push!(fold_errors, 0.0)  # Placeholder
            end
        end
        
        cv_errors[i] = mean(fold_errors)
    end
    
    return cv_errors
end
```

These examples demonstrate the versatility and robustness of SmoothedCollocation.jl for handling various types of noisy data across different domains.