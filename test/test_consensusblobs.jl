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

@testset "blobs & bipartitions with chains of 2-blobs" begin
chainnwk = "(((((((((e)#H3,(#H3,d)),c,b))#H2),#H2))#H1,#H1),a);"
net = readnewick.([chainnwk, chainnwk, chainnwk])
suppressroot!(net[1]); net[1].rooti = 14 # to root at leaf "a"
removedegree2nodes!(net[2]); # suppress root & another node
PN.deletehybridedge!(net[3], net[3].edge[3]) # shrinks the 3-cycle
taxa = sort(tiplabels(net[1]))
blobs, bps = PhyloSummaries.count_blobpartitions(net, taxa, 4)
@test isempty(blobs)
@test length(bps) == 1
@test bps[1].split == (true, true, true, false, false)
@test PhyloSummaries.freq(bps[1]) == 3
blobs, bps = PhyloSummaries.count_blobpartitions(net, taxa, 3)
@test length(blobs) == 1
# fixit: check correct 3-way partition, hybrid and circular order
@test PhyloSummaries.freq(blobs[1]) == 2
@test length(bps) == 1
@test bps[1].split == (true, true, true, false, false)
@test PhyloSummaries.freq(bps[1]) == 1
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
