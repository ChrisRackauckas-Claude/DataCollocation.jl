# Performance Optimization Techniques

This page covers advanced techniques for optimizing performance when working with large datasets or when automatic differentiation performance is critical.

## Non-Allocating Forward-Mode L2 Collocation Loss

For large systems or when using forward-mode automatic differentiation, it's important to minimize allocations. The following example shows how to construct a non-allocating loss function over the collocation that is compatible with forward-mode automatic differentiation:

```julia
using DataCollocation, PreallocationTools
using OrdinaryDiffEq, Optimization, OptimizationOptimisers
using ComponentArrays, Lux, Random

# Setup example system
rng = Xoshiro(42)
u0 = [1.0, 0.5]
tspan = (0.0, 2.0)
tsteps = range(tspan[1], tspan[2], length=100)

function example_ode!(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = -p[3] * u[2] + p[4] * u[1] * u[2]
end

# Generate data and perform collocation
true_p = [1.5, 1.0, 3.0, 1.0]
prob = ODEProblem(example_ode!, u0, tspan, true_p)
data = Array(solve(prob, Tsit5(), saveat=tsteps)) + 0.05*randn(rng, 2, length(tsteps))
estimated_derivative, estimated_solution = collocate_data(data, tsteps, EpanechnikovKernel())

# Pre-allocate arrays for non-allocating loss
du = PreallocationTools.dualcache(similar(u0))
preview_est_sol = [@view estimated_solution[:, i] for i in 1:size(estimated_solution, 2)]
preview_est_deriv = [@view estimated_derivative[:, i] for i in 1:size(estimated_solution, 2)]

function construct_iip_cost_function(f, du, preview_est_sol, preview_est_deriv, tpoints)
    function (p)
        _du = PreallocationTools.get_tmp(du, p)
        vecdu = vec(_du)
        cost = zero(first(p))
        for i in 1:length(preview_est_sol)
            est_sol = preview_est_sol[i]
            f(_du, est_sol, p, tpoints[i])
            vecdu .= vec(preview_est_deriv[i]) .- vec(_du)
            cost += sum(abs2, vecdu)
        end
        sqrt(cost)
    end
end

# Create the non-allocating cost function
cost_function = construct_iip_cost_function(
    example_ode!, du, preview_est_sol, preview_est_deriv, tsteps)

# Test the cost function
p_test = [1.0, 1.0, 1.0, 1.0]
println("Cost with test parameters: ", cost_function(p_test))
```

## Benefits of Non-Allocating Loss

This approach provides several advantages:

1. **Memory efficiency**: No intermediate arrays allocated during optimization
2. **Forward-mode AD compatibility**: Works efficiently with ForwardDiff.jl
3. **Scalability**: Performance doesn't degrade with system size
4. **Cache efficiency**: Better CPU cache utilization

## Performance Comparison

Let's compare allocating vs non-allocating versions:

```julia
using BenchmarkTools

# Allocating version (simple but inefficient)
function allocating_loss(p, estimated_derivative, estimated_solution, tsteps, ode_func)
    cost = 0.0
    for i in 1:size(estimated_derivative, 2)
        u = estimated_solution[:, i]
        du_pred = zeros(eltype(p), length(u))
        ode_func(du_pred, u, p, tsteps[i])
        du_actual = estimated_derivative[:, i]
        cost += sum((du_actual - du_pred).^2)
    end
    return sqrt(cost)
end

# Benchmark comparison
p_bench = true_p
println("Allocating version:")
@btime allocating_loss($p_bench, $estimated_derivative, $estimated_solution, $tsteps, $example_ode!)

println("Non-allocating version:")
@btime cost_function($p_bench)
```

## Memory Usage Analysis

Monitor memory allocations during optimization:

```julia
# Track allocations
function track_allocations(cost_func, p)
    # Warm up
    cost_func(p)
    
    # Measure allocations
    allocated_before = Base.gc_bytes()
    result = cost_func(p)
    allocated_after = Base.gc_bytes()
    
    println("Allocated: $(allocated_after - allocated_before) bytes")
    return result
end

println("Testing allocation patterns:")
track_allocations(cost_function, p_test)
```

## Integration with Optimization

Using the non-allocating loss with optimization frameworks:

```julia
# Define neural network for parameter fitting
param_net = Dense(2, 4)  # Maps state to parameters
net_params, st = Lux.setup(rng, param_net)

function neural_ode_loss(net_p)
    # Predict parameters using neural network
    predicted_params, _ = param_net(u0, net_p, st)
    
    # Use non-allocating cost function
    return cost_function(predicted_params)
end

# Setup optimization
adtype = Optimization.AutoForwardDiff()  # Use forward-mode AD
optf = Optimization.OptimizationFunction((x, p) -> neural_ode_loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentArray(net_params))

# Optimize with reduced allocations
result = Optimization.solve(optprob, OptimizationOptimisers.Adam(0.01), maxiters=1000)
```

## Advanced: Multi-System Batch Processing

For processing multiple related systems efficiently:

```julia
function batch_collocation_loss(systems_data, systems_tsteps, ode_functions)
    total_cost = 0.0
    
    for (data, tsteps, ode_func) in zip(systems_data, systems_tsteps, ode_functions)
        # Perform collocation
        est_deriv, est_sol = collocate_data(data, tsteps, EpanechnikovKernel())
        
        # Create non-allocating loss for this system
        du_cache = PreallocationTools.dualcache(similar(data[:, 1]))
        preview_sol = [@view est_sol[:, i] for i in 1:size(est_sol, 2)]
        preview_deriv = [@view est_deriv[:, i] for i in 1:size(est_deriv, 2)]
        
        cost_func = construct_iip_cost_function(ode_func, du_cache, preview_sol, preview_deriv, tsteps)
        
        # Add to total cost
        total_cost += cost_func(p)
    end
    
    return total_cost
end
```

## Thread-Safe Implementation

For parallel processing of multiple datasets:

```julia
using Base.Threads

function threaded_collocation_loss(datasets, parameters)
    costs = Vector{Float64}(undef, length(datasets))
    
    Threads.@threads for i in 1:length(datasets)
        data, tsteps, ode_func = datasets[i]
        
        # Each thread gets its own pre-allocated arrays
        est_deriv, est_sol = collocate_data(data, tsteps, EpanechnikovKernel())
        du_local = PreallocationTools.dualcache(similar(data[:, 1]))
        preview_sol = [@view est_sol[:, j] for j in 1:size(est_sol, 2)]
        preview_deriv = [@view est_deriv[:, j] for j in 1:size(est_deriv, 2)]
        
        cost_func = construct_iip_cost_function(ode_func, du_local, preview_sol, preview_deriv, tsteps)
        costs[i] = cost_func(parameters)
    end
    
    return sum(costs)
end
```

## Best Practices for Performance

1. **Pre-allocate all working arrays** using PreallocationTools.jl
2. **Use views instead of copies** when accessing data slices
3. **Choose appropriate AD backend** (ForwardDiff for small systems, Zygote for large)
4. **Profile your code** to identify remaining bottlenecks
5. **Consider threading** for independent parallel computations

These techniques can significantly improve performance for large-scale parameter estimation problems using collocation methods.