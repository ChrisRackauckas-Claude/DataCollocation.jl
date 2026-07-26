using DataCollocations
import DataCollocations: calckernel
using JLArrays
using Test

struct ConstantKernel <: CollocationKernel end

calckernel(::ConstantKernel, t) = one(t)

@testset "Interface Compatibility" begin
    @testset "Custom kernel extension" begin
        tpoints = collect(range(0.0, stop = 1.0, length = 30))
        data = sin.(tpoints)

        @test calckernel(ConstantKernel(), 2.0) == 1.0

        du, u = collocate_data(data, tpoints, ConstantKernel(), 1.0)
        @test length(du) == length(tpoints)
        @test length(u) == length(tpoints)
    end

    @testset "BigFloat support" begin
        # Test that BigFloat inputs are supported and eltype is preserved
        tpoints = BigFloat.(collect(range(0.0, stop = 10.0, length = 30)))
        data = BigFloat.(reshape(sin.(Float64.(tpoints)), 1, :))
        bandwidth = BigFloat(0.5)

        kernels = [
            EpanechnikovKernel(), UniformKernel(), TriangularKernel(),
            QuarticKernel(), TriweightKernel(), TricubeKernel(),
            GaussianKernel(), CosineKernel(), LogisticKernel(),
            SigmoidKernel(), SilvermanKernel(),
        ]

        @testset "$kernel" for kernel in kernels
            u_prime, u = collocate_data(data, tpoints, kernel, bandwidth)
            @test eltype(u_prime) == BigFloat
            @test eltype(u) == BigFloat
        end
    end

    @testset "GPU array error message" begin
        # Test that GPU-like arrays (JLArrays) give a clear error message
        # instead of cryptic scalar indexing errors
        tpoints = JLArray([1.0, 2.0, 3.0, 4.0, 5.0])
        data = JLArray([1.0 2.0 3.0 4.0 5.0])

        @test_throws ArgumentError collocate_data(data, tpoints)

        # Also test that only GPU tpoints triggers error
        tpoints_cpu = collect(1.0:5.0)
        @test_throws ArgumentError collocate_data(data, tpoints_cpu)

        # And only GPU data triggers error
        data_cpu = reshape(collect(1.0:5.0), 1, :)
        @test_throws ArgumentError collocate_data(data_cpu, tpoints)
    end
end
