@testset "consensus blobs" begin

nwk = [ # some rooted, some unrooted
  "(A,((((B,(C)#H1:::0.7),(#H1:::0.3,D)))#H0,#H0),E);",
  "(A,(((B,(C)#H1:::0.65),(#H1:::0.35,D))#H0,#H0),E);", # 2 bicomp = 1 blob
  "(A,((C,(B)#H1),(#H1,E)),D);",
  "(D,((C,(B)#H1),(#H1,(A,E))));",
  "(C,((B)#H1,((#H1,(A,E)),D)));", # same as previous but rooted at C
]
# sanity-check canonical order helper
@testset "canonicalorders" begin
    idxmap = [5, 1, 3, 2, 4]
    start_idx = 2 # position of partition[1] in idxmap
    hybrid_idx = 2
    cw, ccw = PS.canonicalorders(idxmap, start_idx, hybrid_idx)
    @test_broken cw == (1, 3, 2, 4, 5)
    @test_broken ccw == (1, 5, 4, 2, 3)
end

@testset "count blobs" begin
taxa = ["D","C","B","A","E"] # not what consensusblob would use
net = readnewick.(nwk)
blobs, bps = @test_logs (:warn, r"^non-binary articulation") PS.count_blobpartitions(net, taxa, 4)
@test [b.partition for b in blobs] == [
 ((0,0,0,1,1), (0,0,1,0,0), (0,1,0,0,0), (1,0,0,0,0)),
 ((1,0,0,1,0), (0,1,0,0,0), (0,0,1,0,0), (0,0,0,0,1)) ]
@test [PS.freq(b) for b in blobs] == [4,1]
@test [b.circorder for b in blobs] == [Dict((1,2,3,4)=>4), Dict((1,2,3,4)=>1)]
@test [b.hybrid for b in blobs] == [Dict(2=>2, 3=>2), Dict(3=>1)]

# blob of 3 biconnected components
net1 = readnewick("((a,(b,(c)#H0)),(#H0,(d)#H1,#H1,#H2,(e)#H2));")
taxa = ["a","b","c","d","e"]
bb, bp = @test_logs (:warn,r"^non-binary") (:warn,r"^n") (:warn,r"^n") PS.count_blobpartitions([net1], taxa, 2)
@test isempty(bp)
@test length(bb) == 1
@test all(sum(p) == 1 for p in bb[1].partition)
@test all(n.intn1==1 for n in net1.node if !n.leaf) # all nodes in same blob 1
@test unique(e.inte1 for e in net1.edge if !isexternal(e)) == [1,2,3] # 3 bicomps

# level-4 (or level-3) blob with 1 exit hybrid only
l4nwk = [
  # net1 has no circular order. net3 = net1 rooted inside the a-clade.
  # net2 = net1 rooted "more inside" blob and ≠ a-clade.
  "((c2,#H1),((((h)#H0,b1))#H2,(((c1,#H0))#H1,#H3)),(((#H2,b2))#H3,(a1,(a2,a3))));",
  "(((c1,#H0))#H1,#H3,((((h)#H0,b1))#H2,((c2,#H1),(((#H2,b2))#H3,((a1,a2),a3)))));",
  "(a2,a3,(a1,(((#H2,b2))#H3,((c2,#H1),((((h)#H0,b1))#H2,(((c1,#H0))#H1,#H3))))));",
  # after deleting edge 8: outer-labeled planar, so has a circular order
  "((c2,#H1),(((b2,((h)#H0,b1)))#H3,(a1,(a2,a3))),(((c1,#H0))#H1,#H3));",
  # after deleting edge 16: 2 exit hybrids
  "((c2,#H1),((((c1,#H0))#H1,#H3),((h)#H0,b1)),((b2)#H3,(a1,(a2,a3))));",
]
net = readnewick.(l4nwk)
taxa = ["a1","a2","a3", "c1","c2", "h", "b1","b2"]
bb, bp = @test_logs (:warn, r"single exit hybrid") PS.count_blobpartitions(net, taxa, 4)
@test [x.split for x in bp] == [(0,1,1,0,0,0,0,0), (1,1,0,0,0,0,0,0)]
@test [PS.freq(x) for x in bp] == [4,1]
@test [x.partition for x in bb] == [
  ((0,0,0,0,1,0,0,0),(0,0,0,1,0,0,0,0),(0,0,0,0,0,1,0,0),
   (0,0,0,0,0,0,1,0),(0,0,0,0,0,0,0,1),(1,1,1,0,0,0,0,0))]
@test_broken isempty(bb[1].circorder)
# fixit: we should turn off calculation of the circular order
# for blobs with a biconnected component of level >1
@test bb[1].hybrid == Dict(
  3 => 5,  # h hybrid in all 5 nets,
  5 => 1) # b2 hybrid in net[5] only
end

@testset "blobs & bipartitions with chains of 2-blobs" begin
chainnwk = "((((((((((e)#H3,(#H3,d))),c,b))#H2),#H2))#H1,#H1),a);"
net = readnewick.([chainnwk, chainnwk, chainnwk])
suppressroot!(net[1]); net[1].rooti = 14 # to root at leaf "a"
removedegree2nodes!(net[2]); # suppress root & another node
PN.deletehybridedge!(net[3], net[3].edge[3]) # shrinks the 3-cycle
taxa = ["a","b","c","d","e"]
blobs, bps = PS.count_blobpartitions(net, taxa, 4)
@test isempty(blobs)
@test length(bps) == 1
@test bps[1].split == (true, true, true, false, false)
@test PS.freq(bps[1]) == 3
blobs, bps = PS.count_blobpartitions(net, taxa, 3)
@test length(blobs) == 1
@test blobs[1].partition == ((1,1,1,0,0),(0,0,0,0,1),(0,0,0,1,0))
@test PS.freq(blobs[1]) == 2
@test blobs[1].hybrid == Dict(2 => 2)
@test blobs[1].circorder == Dict((1,2,3) => 2)
@test length(bps) == 1
@test bps[1].split == (true, true, true, false, false)
@test PS.freq(bps[1]) == 1
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
