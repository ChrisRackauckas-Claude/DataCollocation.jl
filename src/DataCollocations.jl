module DataCollocations

using LinearAlgebra: Diagonal, det, mul!
using ArrayInterface: fast_scalar_indexing

export collocate_data
export CollocationKernel, calckernel
export EpanechnikovKernel, UniformKernel, TriangularKernel, QuarticKernel
export TriweightKernel, TricubeKernel, GaussianKernel, CosineKernel
export LogisticKernel, SigmoidKernel, SilvermanKernel

include("collocation.jl")

using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    @compile_workload begin
        kernel = TriangularKernel()
        calckernel(kernel, 0.25)
        calckernel(GaussianKernel(), 0.25)
        tpoints = collect(range(0.0, 1.0; length = 12))
        data = reshape(sin.(tpoints), 1, length(tpoints))
        collocate_data(data, tpoints, kernel, 0.2)
    end
end

end
