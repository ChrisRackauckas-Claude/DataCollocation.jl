using DataCollocations, DataInterpolations, SciMLTesting, Test

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension modules actually exist rather than trusting a green run_qa.
@testset "Extensions loaded" begin
    @test Base.get_extension(DataCollocations, :DataCollocationsDataInterpolationsExt) !== nothing
end

run_qa(DataCollocations)
