# Core Functions

## Main Functions

```@docs
collocate_data
```

## Function Details

### `collocate_data`

The main function for performing smoothed collocation. It has two primary methods:

#### Method 1: Kernel-based Smoothing

```julia
collocate_data(data, tpoints, kernel=TriangularKernel(), bandwidth=nothing)
```

**Arguments:**
- `data::AbstractMatrix`: Input data where each column represents a snapshot at time `tpoints[i]`
- `data::AbstractVector`: For single-variable data, automatically reshaped to matrix form
- `tpoints::AbstractVector`: Time points corresponding to data columns
- `kernel::CollocationKernel`: Kernel function for smoothing (default: `TriangularKernel()`)
- `bandwidth::Union{Number,Nothing}`: Smoothing bandwidth. If `nothing`, automatically selected

**Returns:**
- `estimated_derivative`: Matrix/Vector of estimated derivatives `u'`
- `estimated_solution`: Matrix/Vector of smoothed data `u`

**Example:**
```julia
using SmoothedCollocation
using OrdinaryDiffEq

# Generate sample data
f(u, p, t) = -0.1 * u
prob = ODEProblem(f, [1.0], (0.0, 5.0))
t = 0:0.1:5
data = Array(solve(prob, Tsit5(), saveat=t)) + 0.1*randn(1, length(t))

# Apply smoothed collocation
du, u = collocate_data(data, t, EpanechnikovKernel(), 0.2)
```

#### Method 2: DataInterpolations Integration

```julia
collocate_data(data, tpoints, tpoints_sample, interp, args...)
```

**Arguments:**
- `data::AbstractMatrix`: Input data matrix
- `data::AbstractVector`: Single-variable input data  
- `tpoints::AbstractVector`: Original time points
- `tpoints_sample::AbstractVector`: New time points for interpolation
- `interp`: Interpolation method from DataInterpolations.jl
- `args...`: Additional arguments passed to the interpolation constructor

**Returns:**
- `estimated_derivative`: Derivatives estimated at `tpoints_sample`
- `estimated_solution`: Interpolated values at `tpoints_sample`

**Example:**
```julia
using SmoothedCollocation, DataInterpolations

# Sparse data
t_sparse = [0, 1, 2, 3, 4, 5]
data_sparse = sin.(t_sparse)

# Dense evaluation points
t_dense = 0:0.1:5
du, u = collocate_data(data_sparse, t_sparse, t_dense, CubicSpline)
```

## Algorithm Details

### Automatic Bandwidth Selection

When `bandwidth=nothing`, the algorithm uses the following formula:

```julia
bandwidth = n^(-1/5) * n^(-3/35) * (log(n))^(-1/16)
```

where `n` is the number of data points. This provides a good balance between bias and variance for most applications.

### Mathematical Foundation

The smoothed collocation method works by:

1. **Kernel Smoothing**: For each time point `t`, construct a local polynomial approximation using weighted data points
2. **Weight Calculation**: Weights are determined by the kernel function: `w_i = K((t_i - t)/h)/h`
3. **Local Regression**: Solve the weighted least squares problem to get coefficients
4. **Derivative Estimation**: Extract derivative from the polynomial coefficients

The method constructs matrices:
- `T1 = [1, (t_i - t)]` for function estimation
- `T2 = [1, (t_i - t), (t_i - t)²]` for derivative estimation

And solves:
- Function value: `u(t) = e₁ᵀ (T₁ᵀWT₁)⁻¹T₁ᵀW data`
- Derivative: `u'(t) = e₂ᵀ (T₂ᵀWT₂)⁻¹T₂ᵀW data`

where `W` is the diagonal weight matrix and `e₁ = [1, 0]`, `e₂ = [0, 1, 0]`.

### Error Conditions

The function will throw an error if:
- The bandwidth is too small, leading to singular matrices
- Data and time points have mismatched dimensions
- Time points are not properly ordered for interpolation methods

**Error Message Example:**
```
ERROR: Collocation failed with bandwidth 0.001. Please choose a higher bandwidth
```

### Performance Considerations

- **Memory Allocation**: The algorithm pre-allocates working arrays to minimize allocations
- **Computational Complexity**: O(n²) for each evaluation point, where n is the number of data points
- **Kernel Choice**: Bounded kernels (like Triangular) are faster than unbounded kernels (like Gaussian)

### Vector vs Matrix Input

```julia
# Vector input (single variable)
data_vec = [1.0, 2.0, 1.5, 0.8]
t = [0.0, 1.0, 2.0, 3.0]
du, u = collocate_data(data_vec, t)  # Returns vectors

# Matrix input (multiple variables)
data_mat = [1.0 2.0 1.5 0.8; 0.5 1.0 0.8 0.4]  # 2 variables, 4 time points
du, u = collocate_data(data_mat, t)  # Returns matrices
```

### Thread Safety

The core collocation algorithm is thread-safe and can be used in multithreaded environments. However, be aware that:
- Random number generation (if used for data generation) should use thread-local RNGs
- DataInterpolations.jl methods may have their own thread safety considerations