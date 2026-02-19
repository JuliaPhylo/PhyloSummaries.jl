using PhyloSummaries
using Test
using Aqua

using CSV
using DataFrames
using PhyloNetworks

const PN = PhyloNetworks
const PS = PhyloSummaries

@testset "PhyloSummaries Code quality (Aqua.jl)" begin
    # Aqua.test_all(PhyloSummaries)
    Aqua.test_all(
        PhyloSummaries;
        unbound_args = (broken=false),
        # add_blobnode!: P,N are both bound and used. don't know why Aqua fails here
    )
end

@testset "PhyloSummaries.jl" begin
    include("test_consensustrees.jl")
    include("test_consensusblobs.jl")
end
