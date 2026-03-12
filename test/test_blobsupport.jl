@testset "blobpartitions_support" begin

# level-1 networks in level1_7taxa_abc.nwk jldoctest-ed in docstring

@testset "level2 nets, 7 taxa" begin
nfile = joinpath(@__DIR__,"..","test","level2_7taxa_abc.nwk")
# nfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","level2_7taxa_abc.nwk")
net = readmultinewick(nfile)
refnet = readnewick("(((a3,(a4,#H1)),a2),(((c2,(#H2,c1)),(b1)#H2))#H1,a1);") # 5th in level1_7taxa_abc
res = blobpartitions_support(net, refnet)
@test keys(res) == (:blob_table, :circorder_table, :hybrid_table, :bipartition_table, :taxa)
@test res[:taxa] == ["a1","a2","a3","a4","b1","c1","c2"]
@test res[:blob_table] == (blob=[2,1], degree=[4,5], node=[-7,-2], support_partition=[0,.6],
    partition_num=["1,2,3,4|7|6|5", "1|2|3|4|5,6,7"],
    partition = ["a1,a2,a3,a4|c2|c1|b1", "a1|a2|a3|a4|b1,c1,c2"])
@test res[:circorder_table].support_circorder == [0,0] # level-2 input: circular order not calculated
@test res[:hybrid_table] == (blob=[2,1,1], node_from=[6,3,-2], node_to=[8,-7,9],
    edge=[13,15,17], support_hybrid=[0,.6,.4], cluster_num=["5","5,6,7","1"],
    cluster=["b1","b1,c1,c2","a1"])
for k in keys(res[:bipartition_table]) @test isempty(res[:bipartition_table][k]); end
end

@testset "level1 nets, 2-cycles, some blob made of 2 bicomps" begin
nfile = joinpath(@__DIR__,"..","test","level1_5taxa_2cycles_1blob2bicomp.nwk")
# nfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","level1_5taxa_2cycles_1blob2bicomp.nwk")
net = readmultinewick(nfile) # 5 networks
refnet = deepcopy(net[2]) # has 2 bicomponents (one of them being a 2-cycle) in 1 blob
res = (@test_logs (:warn, r"non-binary articulation node") (:warn,) blobpartitions_support(net, refnet))
@test res[:blob_table] == (blob=[1], degree=[4], node=[-3], support_partition=[.8],
    partition_num=["1,5|2|3|4"], partition=["A,E|B|C|D"])
@test res[:circorder_table].support_circorder == [0.8]
@test res[:hybrid_table] == (blob=[1,1], node_from=[-5,4], node_to=[2,3],
    edge=[2,3], support_hybrid=[.6,.4], cluster_num=["2","3"], cluster=["B", "C"])
for k in(:node1, :node2, :edge, :support_nonredundant, :cluster_num, :cluster)
    @test isempty(res[:bipartition_table][k])
end
end

@testset "misc: many 2&3-cycles, use net weights, expected errors" begin
net = readnewick.(
["(((((((((e)#H3,(#H3,d))),c,b))#H2),#H2))#H1,#H1,a);",
 "(((#H2,((c,b,((e)#H3,(#H3,d))))#H2))#H1,#H1,a);",
 "(((((((((e,d)),(c,b)))#H2),#H2))#H1,#H1),a);"])
refnet = readnewick("(a,((c,b),((d)#H3,(#H3,e))));")
res = blobpartitions_support(net, refnet; netweight=[1,1,3])
@test res[:blob_table] == (blob=[1], degree=[3], node=[-5], support_partition=[.4],
    partition_num=["1,2,3|4|5"], partition=["a,b,c|d|e"])
@test res[:circorder_table].support_circorder ≈ [.4] # (1+1)/(1+1+3) from weights
@test res[:hybrid_table] == (blob=[1,1], node_from=[5,-7], node_to=[4,6],
    edge=[5,8], support_hybrid=[0.,.4], cluster_num=["4","5"], cluster=["d","e"])
@test res[:bipartition_table] == (node1=[-3], node2=[-4], edge=[4],
    support_nonredundant=[.6], cluster_num=["2,3"], cluster=["b,c"])
res = blobpartitions_support(net, refnet; netweight=[1,1,3], minimumblobdegree=4)
for k in keys(res[:blob_table]) @test isempty(res[:blob_table][k]); end
@test res[:bipartition_table] == (node1=[-3,-3], node2=[-5,-4], edge=[10,4],
    support_nonredundant=[1,.6], cluster_num=["1,2,3","2,3"], cluster=["a,b,c","b,c"])

@test_throws ArgumentError blobpartitions_support(net, refnet; minimumblobdegree=2)
@test_throws ArgumentError blobpartitions_support(HybridNetwork[], refnet)
@test_throws ErrorException blobpartitions_support(net, refnet; netweight=[1.0,1])
@test_throws ErrorException blobpartitions_support(net, refnet; netweight=[1,0.,-1])
end

end
