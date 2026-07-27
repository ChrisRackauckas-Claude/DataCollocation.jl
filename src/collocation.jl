"""
    CollocationKernel

Abstract interface for smoothing kernels used by [`collocate_data`](@ref).

To define a custom kernel, subtype `CollocationKernel` and implement
[`calckernel`](@ref) for the custom type and a scalar evaluation point. The
method must return the kernel weight without mutating either argument. Its
result must be compatible with the element type of the data and time points
passed to `collocate_data`.
"""
abstract type CollocationKernel end

"""
    EpanechnikovKernel()

Epanechnikov smoothing kernel for `collocate_data`.
"""
struct EpanechnikovKernel <: CollocationKernel end

"""
    UniformKernel()

Uniform smoothing kernel for `collocate_data`.
"""
struct UniformKernel <: CollocationKernel end

"""
    TriangularKernel()

Triangular smoothing kernel for `collocate_data`.
"""
struct TriangularKernel <: CollocationKernel end

"""
    QuarticKernel()

Quartic smoothing kernel for `collocate_data`.
"""
struct QuarticKernel <: CollocationKernel end

"""
    TriweightKernel()

Triweight smoothing kernel for `collocate_data`.
"""
struct TriweightKernel <: CollocationKernel end

"""
    TricubeKernel()

Tricube smoothing kernel for `collocate_data`.
"""
struct TricubeKernel <: CollocationKernel end

"""
    GaussianKernel()

Gaussian smoothing kernel for `collocate_data`.
"""
struct GaussianKernel <: CollocationKernel end

"""
    CosineKernel()

Cosine smoothing kernel for `collocate_data`.
"""
struct CosineKernel <: CollocationKernel end

"""
    LogisticKernel()

Logistic smoothing kernel for `collocate_data`.
"""
struct LogisticKernel <: CollocationKernel end

"""
    SigmoidKernel()

Sigmoid smoothing kernel for `collocate_data`.
"""
struct SigmoidKernel <: CollocationKernel end

"""
    SilvermanKernel()

Silverman smoothing kernel for `collocate_data`.
"""
struct SilvermanKernel <: CollocationKernel end

"""
    calckernel(kernel, t)

Evaluate a [`CollocationKernel`](@ref) at scalar point `t`.

Custom `CollocationKernel` subtypes must implement this generic function. The
method must return a scalar kernel weight, preserve compatibility with `t`'s
numeric type, and not mutate its arguments.
"""
function calckernel(kernel::CollocationKernel, t::T) where {T}
    abst = abs(t)
    return ifelse(abst > 1, T(0), calckernel(kernel, t, abst))
end
calckernel(::EpanechnikovKernel, t::T, abst::T) where {T} = T(0.75) * (T(1) - t^2)
calckernel(::UniformKernel, t::T, abst::T) where {T} = T(0.5)
calckernel(::TriangularKernel, t::T, abst::T) where {T} = T(1) - abst
calckernel(::QuarticKernel, t::T, abst::T) where {T} = T(15) * (T(1) - t^2)^2 / T(16)
calckernel(::TriweightKernel, t::T, abst::T) where {T} = T(35) * (T(1) - t^2)^3 / T(32)
calckernel(::TricubeKernel, t::T, abst::T) where {T} = T(70) * (T(1) - abst^3)^3 / T(81)
calckernel(::CosineKernel, t::T, abst::T) where {T} = T(π) * cospi(t / T(2)) / T(4)

calckernel(::GaussianKernel, t::T) where {T} = exp(-t^2 / T(2)) / sqrt(T(2) * π)
calckernel(::LogisticKernel, t::T) where {T} = T(1) / (exp(t) + T(2) + exp(-t))
calckernel(::SigmoidKernel, t::T) where {T} = T(2) / (π * (exp(t) + exp(-t)))
function calckernel(::SilvermanKernel, t::T) where {T}
    return sin(abs(t) / T(2) + π / T(4)) * T(0.5) * exp(-abs(t) / sqrt(T(2)))
end

function construct_t1(t::Number, tpoints::AbstractVector)
    return hcat(ones(eltype(tpoints), length(tpoints)), tpoints .- t)
end

function construct_t2(t::Number, tpoints::AbstractVector)
    return hcat(ones(eltype(tpoints), length(tpoints)), tpoints .- t, (tpoints .- t) .^ 2)
end

function construct_w(t::Number, tpoints::AbstractVector, h::Number, kernel::CollocationKernel)
    W = @. calckernel((kernel,), ((tpoints - t) / (tpoints[end] - tpoints[begin])) / h) / h
    return Diagonal(W)
end

"""
    du, u = collocate_data(data, tpoints, kernel = TriangularKernel(), bandwidth = nothing)
    du, u = collocate_data(data, tpoints, tpoints_sample, interp, args...)

Estimate smoothed derivative and state samples from timeseries data.

For kernel smoothing, each column of matrix-valued `data` is a snapshot at
`tpoints[i]`. Vector-valued `data` is treated as a one-dimensional timeseries. The
returned `du` and `u` have the same dimensionality convention as `data`.

## Arguments

- `data`: observations to smooth, either as a vector or as a matrix with observations
  stored by column.
- `tpoints`: timestamps corresponding to the observations in `data`.
- `kernel`: smoothing kernel instance, such as `TriangularKernel()` or
  `GaussianKernel()`.
- `bandwidth`: kernel bandwidth. When `nothing`, an adaptive default based on
  `length(tpoints)` is used.
- `tpoints_sample`: timestamps where the interpolation-based method evaluates `u`
  and `du`.
- `interp`: `DataInterpolations.jl` interpolation constructor used by the extension
  method.
- `args...`: additional arguments forwarded to `interp`.

Available kernels are `EpanechnikovKernel`, `UniformKernel`, `TriangularKernel`,
`QuarticKernel`, `TriweightKernel`, `TricubeKernel`, `GaussianKernel`,
`CosineKernel`, `LogisticKernel`, `SigmoidKernel`, and `SilvermanKernel`.

Additionally, we can use interpolation methods from
[DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) to generate
data from intermediate timesteps. In this case, pass any of the methods like
`QuadraticInterpolation` as `interp`, and the timestamps to sample from as `tpoints_sample`.

See also: [Kernel smoothing and collocation methods](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631937/).
"""
function collocate_data(
        data::AbstractMatrix, tpoints::AbstractVector,
        kernel::CollocationKernel = TriangularKernel(), bandwidth::Union{Number, Nothing} = nothing
    )
    if !fast_scalar_indexing(data) || !fast_scalar_indexing(tpoints)
        throw(
            ArgumentError(
                "collocate_data requires arrays that support fast scalar indexing. " *
                    "GPU arrays are not supported. Please use CPU arrays instead."
            )
        )
    end
    _one = oneunit(first(data))
    _zero = zero(first(data))
    e1 = [_one; _zero]
    e2 = [_zero; _one; _zero]
    n = length(tpoints)
    bandwidth = bandwidth === nothing ?
        (n^(-1 / 5)) * (n^(-3 / 35)) * ((log(n))^(-1 / 16)) : bandwidth

    Wd = similar(data, n, size(data, 1))
    WT1 = similar(data, n, 2)
    WT2 = similar(data, n, 3)
    T2WT2 = similar(data, 3, 3)
    T1WT1 = similar(data, 2, 2)
    x = map(tpoints) do _t
        T1 = construct_t1(_t, tpoints)
        T2 = construct_t2(_t, tpoints)
        W = construct_w(_t, tpoints, bandwidth, kernel)
        mul!(Wd, W, data')
        mul!(WT1, W, T1)
        mul!(WT2, W, T2)
        mul!(T2WT2, T2', WT2)
        mul!(T1WT1, T1', WT1)
        (det(T2WT2) ≈ 0.0 || det(T1WT1) ≈ 0.0) &&
            error("Collocation failed with bandwidth $bandwidth. Please choose a higher bandwidth")
        (e2' * ((T2' * WT2) \ T2')) * Wd, (e1' * ((T1' * WT1) \ T1')) * Wd
    end
    estimated_derivative = mapreduce(xᵢ -> transpose(first(xᵢ)), hcat, x)
    estimated_solution = mapreduce(xᵢ -> transpose(last(xᵢ)), hcat, x)
    return estimated_derivative, estimated_solution
end

# Convenience method for 1D vector data (kernel-based)
function collocate_data(
        data::AbstractVector, tpoints::AbstractVector,
        kernel::CollocationKernel = TriangularKernel(), bandwidth::Union{Number, Nothing} = nothing
    )
    du, u = collocate_data(reshape(data, 1, :), tpoints, kernel, bandwidth)
    return vec(du), vec(u)
end

# Convenience method for 1D vector data (interpolation-based)
@views function collocate_data(
        data::AbstractVector, tpoints::AbstractVector,
        tpoints_sample::AbstractVector, interp, args...
    )
    du, u = collocate_data(reshape(data, 1, :), tpoints, tpoints_sample, interp, args...)
    return du[1, :], u[1, :]
end
