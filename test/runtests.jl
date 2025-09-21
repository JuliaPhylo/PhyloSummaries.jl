using PhyloSummaries
using Test
using Aqua

using PhyloNetworks

const PN = PhyloNetworks

@testset "PhyloSummaries Code quality (Aqua.jl)" begin
    Aqua.test_all(PhyloSummaries)
    #= Test.detect_ambiguities(PhyloTraits)
    Aqua.test_all(
        PhyloTraits;
        ambiguities = (broken=false),
        persistent_tasks = false,
    )
    =#
end

@testset "PhyloSummaries.jl" begin
    include("test_consensusTrees.jl")
end
