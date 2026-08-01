# DataInterpolations is loaded so that DataCollocationsDataInterpolationsExt exists:
# ExplicitImports only scans an extension module once its trigger package is loaded.
using DataCollocations, DataInterpolations, SciMLTesting, Test

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension modules actually exist rather than trusting a green run_qa.
@testset "Extensions loaded" begin
    @test Base.get_extension(DataCollocations, :DataCollocationsDataInterpolationsExt) !== nothing
end

run_qa(
    DataCollocations;
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            # DataInterpolations: `derivative` is the documented way to differentiate an
            # interpolation object, but DataInterpolations neither exports it nor declares
            # it `public`, so there is no public spelling to use instead.
            ignore = (:derivative,),
        ),
    ),
)
