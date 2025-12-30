@testset "consensus blobs" begin

nwk_strings = [
    "(A,(((B,(C)#H1:::0.7),(D,#H1:::0.3))#H0,#H0),E);",
    "(A,((C,(B)#H1:::0.65),(E,#H1:::0.35)),D);",
]
net = readnewick.(nwk_strings)
# sanity-check canonical order helper
@testset "canonicalorders" begin
    idxmap = [5, 1, 3, 2, 4]
    start_idx = 2 # position of partition[1] in idxmap
    hybrid_idx = 2
    cw, ccw = PhyloSummaries.canonicalorders(idxmap, start_idx, hybrid_idx)
    @test cw == (1, 3, 2, 4, 5)
    @test ccw == (1, 5, 4, 2, 3)
end

@testset "consensusblobs basics" begin
    blobs, bps = consensus_treeofblobs(net)
    @test !isempty(blobs)
    @test reduce(+, b.freq for b in blobs) == length(net)
    @test any(b.freq == length(net) for b in blobs)
end
# fixit: get consensus ToB etc. & test
#= to look at these networks locally:
using RCall, PhyloNetworks, PhyloPlots
R"layout"([1 2 3 4 5; 6 7 8 9 10]);
R"par"(mar=[0,0,0,0])
for i in 1:10
    plot(net[i], shownodenumber=true, tipoffset=0.2, nodecex=0.2)
    R"mtext"("net $i", side=1, line=-1)
end
=#
end
