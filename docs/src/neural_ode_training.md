# Smoothed Collocation for Fast Two-Stage Training

!!! note
    
    This is one of many methods for calculating the collocation coefficients
    for the training process. For a more comprehensive set of collocation
    methods, see [JuliaSimModelOptimizer](https://help.juliahub.com/jsmo/stable/manual/collocation/).

One can avoid a lot of the computational cost of the ODE solver by
pretraining the neural network against a smoothed collocation of the
data. First the example and then an explanation.

## Copy-Pasteable Code

Before getting to the explanation, here's some code to start with. We will follow a full explanation of the definition and training process:

```@example collocation_cp
using ComponentArrays, Lux, DataCollocation, DiffEqFlux, OrdinaryDiffEq, SciMLSensitivity, Optimization,
      OptimizationOptimisers, Plots
using DataCollocation: EpanechnikovKernel

using Random
rng = Xoshiro(0)

u0 = Float32[2.0; 0.0]
datasize = 300
tspan = (0.0f0, 1.5f0)
tsteps = range(tspan[1], tspan[2]; length = datasize)

function trueODEfunc(du, u, p, t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u .^ 3)'true_A)'
end

prob_trueode = ODEProblem(trueODEfunc, u0, tspan)
data = Array(solve(prob_trueode, Tsit5(); saveat = tsteps)) .+ 0.1randn(2, 300)

du, u = collocate_data(data, tsteps, EpanechnikovKernel())

scatter(tsteps, data')
plot!(tsteps, u'; lw = 5)
savefig("colloc.png")
plot(tsteps, du')
savefig("colloc_du.png")

dudt2 = Chain(x -> x .^ 3, Dense(2, 50, tanh), Dense(50, 2))

function loss(p)
    cost = zero(first(p))
    for i in 1:size(du, 2)
        _du, _ = dudt2(@view(u[:, i]), p, st)
        dui = @view du[:, i]
        cost += sum(abs2, dui .- _du)
    end
    sqrt(cost)
end

pinit, st = Lux.setup(rng, dudt2)

callback = function (p, l)
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentArray(pinit))

result_neuralode = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(0.05); callback, maxiters = 10000)

prob_neuralode = NeuralODE(dudt2, tspan, Tsit5(); saveat = tsteps)
nn_sol, st = prob_neuralode(u0, result_neuralode.u, st)
scatter(tsteps, data')
plot!(nn_sol)
savefig("colloc_trained.png")

function predict_neuralode(p)
    Array(prob_neuralode(u0, p, st)[1])
end

function loss_neuralode(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, data .- pred)
    return loss
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss_neuralode(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentArray(pinit))

numerical_neuralode = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(0.05); callback, maxiters = 300)

nn_sol, st = prob_neuralode(u0, numerical_neuralode.u, st)
scatter(tsteps, data')
plot!(nn_sol; lw = 5)
```

## Generating the Collocation

The smoothed collocation is a spline fit of the data points which allows
us to get an estimate of the approximate noiseless dynamics:

```@example collocation_cp
scatter(tsteps, data')
plot!(tsteps, u'; lw = 5)
```

We can then differentiate the smoothed function to get estimates of the
derivative at each data point:

```@example collocation_cp
plot(tsteps, du')
```

Because we have `(u',u)` pairs, we can write a loss function that
calculates the squared difference between `f(u,p,t)` and `u'` at each
point, and find the parameters which minimize this difference:

```@example collocation_cp
dudt2 = Chain(x -> x .^ 3, Dense(2, 50, tanh), Dense(50, 2))

function loss(p)
    cost = zero(first(p))
    for i in 1:size(du, 2)
        _du, _ = dudt2(@view(u[:, i]), p, st)
        dui = @view du[:, i]
        cost += sum(abs2, dui .- _du)
    end
    sqrt(cost)
end

pinit, st = Lux.setup(rng, dudt2)

callback = function (p, l)
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentArray(pinit))

result_neuralode = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(0.05); callback, maxiters = 10000)

prob_neuralode = NeuralODE(dudt2, tspan, Tsit5(); saveat = tsteps)
nn_sol, st = prob_neuralode(u0, result_neuralode.u, st)
scatter(tsteps, data')
plot!(nn_sol)
```

While this doesn't look great, it has the characteristics of the
full solution all throughout the timeseries, but it does have a drift.
We can continue to optimize like this, or we can use this as the
initial condition to the next phase of our fitting:

```@example collocation_cp
function predict_neuralode(p)
    Array(prob_neuralode(u0, p, st)[1])
end

function loss_neuralode(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, data .- pred)
    return loss
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss_neuralode(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentArray(pinit))

numerical_neuralode = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(0.05); callback, maxiters = 300)

nn_sol, st = prob_neuralode(u0, numerical_neuralode.u, st)
scatter(tsteps, data')
plot!(nn_sol; lw = 5)
```

This method then has a good global starting position, making it less
prone to local minima, and this method is thus a great method to mix in with other
fitting methods for neural ODEs.
