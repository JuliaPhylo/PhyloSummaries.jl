using PhyloSummaries
using Test
using Aqua

@testset "PhyloSummaries.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PhyloSummaries)
    end
    # Write your tests here.
end
