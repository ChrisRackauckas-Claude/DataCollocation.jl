# Neural ODE Training with Smoothed Collocation

This example demonstrates how to use smoothed collocation for fast two-stage neural ODE training, based on the methodology from DiffEqFlux.jl.

## Overview

Traditional neural ODE training requires solving the ODE at each optimization step, which can be computationally expensive. Smoothed collocation provides a faster alternative by:

1. **Stage 1**: Fit the neural network to match estimated derivatives from smoothed data
2. **Stage 2**: Fine-tune using traditional neural ODE loss (optional)

This two-stage approach provides better initial conditions and faster convergence.

## Complete Example

```julia
using SmoothedCollocation, OrdinaryDiffEq, Optimization, OptimizationOptimisers
using ComponentArrays, Lux, Random, Plots
using Zygote  # for automatic differentiation

# Set random seed for reproducibility
rng = Xoshiro(0)

# Define the true ODE system
u0 = Float32[2.0; 0.0]
datasize = 300
tspan = (0.0f0, 1.5f0)
tsteps = range(tspan[1], tspan[2]; length = datasize)

function trueODEfunc(du, u, p, t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u .^ 3)'true_A)'
end

# Generate noisy training data
prob_trueode = ODEProblem(trueODEfunc, u0, tspan)
data = Array(solve(prob_trueode, Tsit5(); saveat = tsteps)) .+ 0.1*randn(rng, 2, 300)

# Visualize the noisy data
scatter(tsteps, data[1,:], label="Noisy u₁", alpha=0.6)
scatter!(tsteps, data[2,:], label="Noisy u₂", alpha=0.6)
plot!(title="Noisy Training Data")
```

## Stage 1: Collocation-Based Training

First, we use smoothed collocation to estimate derivatives and train the neural network:

```julia
# Apply smoothed collocation
du, u = collocate_data(data, tsteps, EpanechnikovKernel())

# Visualize the smoothed results
scatter(tsteps, data[1,:], label="Noisy u₁", alpha=0.4)
scatter!(tsteps, data[2,:], label="Noisy u₂", alpha=0.4)
plot!(tsteps, u[1,:], label="Smoothed u₁", lw=3)
plot!(tsteps, u[2,:], label="Smoothed u₂", lw=3)
plot!(title="Original vs Smoothed Data")
```

```julia
# Plot the estimated derivatives
plot(tsteps, du[1,:], label="du₁/dt", lw=3)
plot!(tsteps, du[2,:], label="du₂/dt", lw=3)
plot!(title="Estimated Derivatives")
```

Now define and train the neural network:

```julia
# Define the neural network architecture
dudt_net = Chain(
    x -> x .^ 3,           # Preprocessing (known from problem structure)
    Dense(2, 50, tanh),    # Hidden layer
    Dense(50, 2)           # Output layer
)

# Initialize network parameters
pinit, st = Lux.setup(rng, dudt_net)

# Define the collocation loss function
function collocation_loss(p)
    cost = zero(first(p))
    for i in 1:size(du, 2)
        # Forward pass through network
        _du, _ = dudt_net(u[:, i], p, st)
        
        # Compare with estimated derivative
        dui = du[:, i]
        cost += sum(abs2, dui .- _du)
    end
    return sqrt(cost)
end

# Set up optimization
callback = function (p, l)
    println("Current loss: ", l)
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> collocation_loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentArray(pinit))

# Train using collocation loss
println("Stage 1: Training with collocation loss...")
result_collocation = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(0.05); 
    callback, maxiters = 5000
)
```

Test the collocation-trained network:

```julia
# Create a neural ODE with the trained parameters
prob_neuralode = NeuralODE(dudt_net, tspan, Tsit5(); saveat = tsteps)
nn_sol_stage1, st = prob_neuralode(u0, result_collocation.u, st)

# Plot Stage 1 results
scatter(tsteps, data[1,:], label="Noisy Data u₁", alpha=0.4)
scatter!(tsteps, data[2,:], label="Noisy Data u₂", alpha=0.4)
plot!(nn_sol_stage1, vars=(0,1), label="Stage 1 u₁", lw=3)
plot!(nn_sol_stage1, vars=(0,2), label="Stage 1 u₂", lw=3)
plot!(title="Stage 1: Collocation Training Results")
```

## Stage 2: Traditional Neural ODE Fine-tuning (Optional)

For even better accuracy, we can fine-tune using traditional neural ODE loss:

```julia
# Define traditional neural ODE prediction and loss
function predict_neuralode(p)
    Array(prob_neuralode(u0, p, st)[1])
end

function neuralode_loss(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, data .- pred)
    return loss
end

# Set up Stage 2 optimization
optf_stage2 = Optimization.OptimizationFunction((x, p) -> neuralode_loss(x), adtype)
optprob_stage2 = Optimization.OptimizationProblem(optf_stage2, result_collocation.u)

# Stage 2 training (fewer iterations needed due to good initialization)
println("Stage 2: Fine-tuning with neural ODE loss...")
result_final = Optimization.solve(
    optprob_stage2, OptimizationOptimisers.Adam(0.01); 
    callback, maxiters = 300
)

# Plot final results
nn_sol_final, st = prob_neuralode(u0, result_final.u, st)

scatter(tsteps, data[1,:], label="Noisy Data u₁", alpha=0.4)
scatter!(tsteps, data[2,:], label="Noisy Data u₂", alpha=0.4)
plot!(nn_sol_stage1, vars=(0,1), label="Stage 1 u₁", lw=2, ls=:dash)
plot!(nn_sol_stage1, vars=(0,2), label="Stage 1 u₂", lw=2, ls=:dash)
plot!(nn_sol_final, vars=(0,1), label="Final u₁", lw=3)
plot!(nn_sol_final, vars=(0,2), label="Final u₂", lw=3)
plot!(title="Two-Stage Training Results")
```

## Performance Comparison

Let's compare the two-stage approach with training from scratch:

```julia
# Train from scratch for comparison (warning: this takes longer!)
function train_from_scratch()
    # Fresh parameters
    p_scratch, st_scratch = Lux.setup(rng, dudt_net)
    
    # Traditional loss only
    optf_scratch = Optimization.OptimizationFunction((x, p) -> neuralode_loss(x), adtype)
    optprob_scratch = Optimization.OptimizationProblem(optf_scratch, ComponentArray(p_scratch))
    
    result_scratch = Optimization.solve(
        optprob_scratch, OptimizationOptimisers.Adam(0.05); 
        callback, maxiters = 5000  # Same total iterations as two-stage
    )
    
    return result_scratch
end

# Uncomment to run comparison (takes time!)
# result_scratch = train_from_scratch()
# nn_sol_scratch, _ = prob_neuralode(u0, result_scratch.u, st)
```

## Advanced: Non-Allocating Loss Function

For better performance, especially with larger systems, you can use a non-allocating loss function:

```julia
using PreallocationTools

# Create pre-allocated arrays
du_temp = PreallocationTools.dualcache(similar(u0))
preview_est_sol = [@view u[:, i] for i in 1:size(u, 2)]
preview_est_deriv = [@view du[:, i] for i in 1:size(du, 2)]

function efficient_collocation_loss(p)
    _du = PreallocationTools.get_tmp(du_temp, p)
    vecdu = vec(_du)
    cost = zero(first(p))
    
    for i in 1:length(preview_est_sol)
        est_sol = preview_est_sol[i]
        est_deriv = preview_est_deriv[i]
        
        # Forward pass
        _du, _ = dudt_net(est_sol, p, st)
        
        # Compute residual
        vecdu .= vec(est_deriv) .- vec(_du)
        cost += sum(abs2, vecdu)
    end
    
    return sqrt(cost)
end
```

## Key Advantages

1. **Faster convergence**: Good initial parameters from Stage 1
2. **Robust to noise**: Smoothed collocation handles noisy data well
3. **Less prone to local minima**: Better initialization landscape
4. **Flexible**: Can use Stage 1 results directly or fine-tune with Stage 2

## When to Use This Approach

- **Noisy data**: Smoothed collocation is robust to measurement noise
- **Fast prototyping**: Quick initial results from Stage 1
- **Complex dynamics**: When traditional neural ODE training struggles
- **Limited computational budget**: Stage 1 alone often gives good results

## Tips for Success

1. **Kernel selection**: `EpanechnikovKernel()` works well for most cases
2. **Bandwidth tuning**: Automatic selection usually works, but manual tuning can help
3. **Network architecture**: Include known problem structure (like the `x .^ 3` preprocessing)
4. **Stage 2 iterations**: Usually need fewer iterations than Stage 1
5. **Learning rates**: Often need to reduce learning rate for Stage 2