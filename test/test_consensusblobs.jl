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
    cw, ccw = PS.canonicalorders([5, 1, 3, 2, 4], 2) # start at index 2
    @test cw == (1, 3, 2, 4, 5)
    @test ccw == (1, 5, 4, 2, 3)
end

@testset "blobcompatible" begin
    A = ((true,true,false,false,false),   # {A,B}
         (false,false,true,false,false),  # {C}
         (false,false,false,true,true))   # {D,E}
    B =  (true,true,true,false,false)     # ABC, not DE
    @test PS.splitblobcompatible(B, A)
    A2 = ((true,true,true,false,false,false),  # {A,B,C}
          (false,false,false,true,true,true))  # {D,E,F}
    B2 = ((true,false,false,true,false,false), # {A,D}
          (false,true,false,false,true,false), # {B,E}
          (false,false,true,false,false,true)) # {C,F}
    @test !PS.blobcompatible(A2, B2)
    A2 = (B2[1], B2[2], # AD|BE|C|F
          (false,false,true,false,false,false), (false,false,false,false,false,true))
    @test !PS.blobcompatible(A2, B2) # tree-compatible, but not blob-compatible
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
  # after deleting edge 16: 2 exit hybrids in same blob
  "((c2,#H1),((((c1,#H0))#H1,#H3),((h)#H0,b1)),((b2)#H3,(a1,(a2,a3))));",
]
net = readnewick.(l4nwk)
taxa = ["a1","a2","a3", "c1","c2", "h", "b1","b2"]
bb, bp = PS.count_blobpartitions(net, taxa, 4)
@test [x.split for x in bp] == [(0,1,1,0,0,0,0,0), (1,1,0,0,0,0,0,0)]
@test [PS.freq(x) for x in bp] == [4,1]
@test [x.partition for x in bb] == [
  ((0,0,0,0,1,0,0,0),(0,0,0,1,0,0,0,0),(0,0,0,0,0,1,0,0),
   (0,0,0,0,0,0,1,0),(0,0,0,0,0,0,0,1),(1,1,1,0,0,0,0,0))]
@test isempty(bb[1].circorder)
@test bb[1].hybrid == Dict( # 6 total in only 5 networks because one has 2 hybrids
  3 => 5,  # h hybrid in all 5 nets,
  5 => 1) # b2 hybrid in net[5] only
@test occursin(
  r"""BlobFreq on 8 taxa, partitioned into 6 blocks
  taxon blocks: 5|4|6|7|8|1,2,3
  frequency: 5
  hybrid block => frequency:""", repr("text/plain", bb[1]))
@test repr(bp[1]) == "SplitFreq on 8 taxa, taxa in split cluster: 2,3, frequency: 4"

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

tob = consensus_treeofblobs(net, minimumblobdegree=3)
@test writenewick(tob) == "(d,e,(c,b,a));"
@test [n.fvalue for n in tob.node] ≈ [-1,-1,-1,-1,-1, 2/3, -1]
@test [e.y for e in tob.edge] ≈ [-1,-1,-1,-1,-1, 1/3]

tob = consensus_treeofblobs(net)
@test writenewick(tob) == "(d,e,(c,b,a));"
@test all(n.fvalue == -1 for n in tob.node)
@test [e.y for e in tob.edge] ≈ [-1,-1,-1,-1,-1, 1]
end

@testset "consensus ToB" begin
# nets 1,2,4,5: same blob AE|B|C|D, same hybrid C, same circular order
# net 3; blob AD|C|B(hybrid)|E
net = readnewick.(nwk)
tob = @test_logs (:warn, r"^non-binary") consensus_treeofblobs(net)
@test writenewick(tob) == "(A,E,(D,C,B));"
@test tob.node[7].fvalue == 0.8 # blob in 4/5 input nets

nfile = joinpath(@__DIR__,"..","test","bootstrapnets_h1.nwk")
# nfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","bootstrapnets_h1.nwk")
net = readmultinewick(nfile)
# 4 blob partitions, each with a single circular order, freqs: 6,1,1,1,1
# 0 non-redundant bipartitions, 4 hybrid clades: t6 ×5, t3 ×2, t4 ×2, t8 ×1
#= to look at networks locally:
using RCall, PhyloPlots
R"layout"([1 2 3 4 5; 6 7 8 9 10]);
R"par"(mar=[0,0,0,0]);
for (i,n) in enumerate(net)
    plot(n, shownodenumber=true, tipoffset=0.2, nodecex=0.2)
    R"mtext"("net $i", side=1, line=-1)
end
=#
# consensus tree-of-blobs: star, blob support 60%:
tob = consensus_treeofblobs(net)
@test writenewick(tob) == "(t1,t2,t3,t4,t5,t6);"
@test getroot(tob).fvalue == 6/10
tob = consensus_treeofblobs(net[[2,3,4,9]])
@test writenewick(tob) == "(t5,t6,(t4,t3,t2,t1));" # blob from nets 4,9
@test tob.node[8].fvalue == 0.5
tob = consensus_treeofblobs(net[[4,9,3,2]])
@test writenewick(tob) == "(t3,t4,t5,t6,(t2,t1));" # blob from nets 2,3
@test getroot(tob).fvalue == 0.5
# fixit: make consensus_treeofblobs return a data table, and test it
conl1, _ = consensus_level1network(net)
@test writenewick(conl1) == "(t2,(t4,(t3,(t5,#H1))),(t1,(t6)#H1));"
# or different ismajor: "(t2,(t4,(t3,(t5,(t6)#H0))),(#H0,t1));"
@test [e.inte1 for e in conl1.edge] == repeat([-1,1], inner=6)
@test all(e -> e.y==-1, conl1.edge)
@test [n.fvalue for n in conl1.node if n.hybrid] == [.5] # hybrid support
@test getroot(conl1).fvalue == 0.6 # blob support
@test all(i -> conl1.node[i].fvalue == 0.6, [8,9,10,12]) # circular order support
# fixit: text returned data frames, here and below

nfile = joinpath(@__DIR__,"..","test","level1_7taxa_abc.nwk")
# nfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","level1_7taxa_abc.nwk")
net = readmultinewick(nfile)
con,_ = consensus_level1network(net, minimumblobdegree=3) # 2 blobs, 0 bps
@test writenewick(con) == "(a3,(a4,#H2),(a2,(a1,(((c2,(c1,(b1)#H1)),#H1))#H2)));"
@test [n.fvalue for n in con.node if n.hybrid] == [.6, .2]
@test [n.fvalue for n in con.node if !n.leaf && !n.hybrid] == [.4,.6,.6,.6,.6,.4,.4]
@test [con.edge[i].inte1 for i in 8:17] == [-1, 2,2,2,2,2, 1,1,1,1]
@test [n.intn1 for n in con.node if !n.leaf] == [1, 2,2,2,2,2, 1,1,1]
con,_ = consensus_level1network(net, outgroup="a2") # 1 blob, 1 bps, hyb edge reversed
@test writenewick(con) == "(a2,((a1,((b1,(c1,c2)))#H1),(a3,(a4,#H1))));"
@test [n.fvalue for n in con.node if !n.leaf] == [-1,.6,-1,-1,.6,.6,.6,.6]
@test [e.inte1 for e in con.edge if !isexternal(e)] == [-1,-1,-1, 1,1,1,1,1]
@test [n.intn1 for n in con.node if !n.leaf] == [-1, 1, -1,-1, 1,1,1,1]
con,_ = consensus_level1network(net[1:4], minimumblobdegree=3) # hedge of weight 0
@test writenewick(con) == "(a3,(a4,#H2),(a2,(a1,((((c1,c2),(b1)#H1),#H1))#H2)));"
@test [n.fvalue for n in con.node if n.hybrid] == [.5, 0]
con,_ = consensus_level1network(net[[1,3,4]], minimumblobdegree=3)
@test writenewick(con) == "(c1,c2,((b1,(((a2,a1),(a4,a3)))#H1),#H1));"
@test [n.fvalue for n in con.node if n.intn1 ≠ -1] ≈ [1,1,2]./3
# fixit: test returned data frames
end

@testset "blob compatibility filtering, show" begin
  b(x...) = NTuple{length(x),Bool}(Bool.(x))
  net1 = readnewick("((((A,B),(C)#H1),(#H1,D)),E);")  # blob: E|AB|C|D
  net2 = readnewick("((((A,B),(C)#H1),(#H1,D)),E);")  # same as net1
  net3 = readnewick("((((A,C),(B)#H1),(#H1,D)),E);")  # blob: E|AC|B|D
  taxa = ["A","B","C","D","E"]
  blobs, bps = PS.count_blobpartitions([net1, net2, net3], taxa, 4)
  s = IOBuffer(); show(s, MIME"text/plain"(), blobs)
  @test String(take!(s)) == """
2-element Vector{PhyloSummaries.BlobFreq{5}}:
 BlobFreq on 5 taxa, 4 blocks 5|1,2|3|4, frequency 2, 1 circular orders, 1 hybrid blocks
 BlobFreq on 5 taxa, 4 blocks 5|1,3|2|4, frequency 1, 1 circular orders, 1 hybrid blocks"""
  s = sprint(show, MIME"text/plain"(), blobs[1])
  @test occursin("1 entry:\n  (1, 2, 3, 4) => 2", s) # 1 circular order, freq 2
  @test occursin("1 entry:\n  3 => 2", s) # 1 hybrid block (C), freq 2
  PS.filter_sort_compatible_partitions!(blobs, bps, 3, 0)
  @test occursin("[BlobFreq on 5 taxa, 4 blocks 5|1,2|3|4, frequency 2, 1 circular orders, 1 hybrid blocks]",
      repr(blobs))
  net4 = readnewick("((((A,B),(C)#H1),(#H1,(D,E))));") # no interesting blob, bipartitions from AB|C|DE
  blobs4, bps4 = PS.count_blobpartitions([net1, net2, net4], taxa, 4)
  s = IOBuffer(); show(s, MIME"text/plain"(), bps4)
  @test String(take!(s)) == """
2-element Vector{PhyloSummaries.SplitFreq{5}}:
 SplitFreq on 5 taxa, taxa in split cluster: 1,2,3, frequency: 1
 SplitFreq on 5 taxa, taxa in split cluster: 1,2, frequency: 1"""
  # test filter_sort_compatible_partitions! with bipartitions
  PS.filter_sort_compatible_partitions!(blobs4, bps4, 3, 0)
  @test length(blobs4) == 1  
  @test length(bps4) == 1    
  @test bps4[1].split == (true, true, false, false, false)  
end

end
