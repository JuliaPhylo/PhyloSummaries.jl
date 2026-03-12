@testset "blobpartitions_support" begin

nfile = joinpath(@__DIR__,"..","test","level1_7taxa_abc.nwk")
# nfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","level1_7taxa_abc.nwk")
nets = readmultinewick(nfile) # 5 networks
refnet = deepcopy(nets[5])

res = blobpartitions_support(nets, refnet)
@test res[:taxa] == ["a1", "a2", "a3", "a4", "b1", "c1", "c2"]

bt = res.blob_table
@test length(bt.blob) == 2
@test bt.degree == [4, 5]
@test bt.partition == ["a1,a2,a3,a4|c2|c1|b1", "a1|a2|a3|a4|b1,c1,c2"]
@test bt.support_partition ≈ [0.4, 0.6]
@test bt.node == [-7, -2]
@test bt.partition_num == ["1,2,3,4|7|6|5", "1|2|3|4|5,6,7"]

ct = res.circorder_table
@test length(ct.blob) == 2
@test ct.support_circorder ≈ [0.4, 0.6]
@test ct.order == [(1, 2, 3, 4), (1, 2, 3, 4, 5)]
@test ct.partition_num == ["1,2,3,4|7|6|5", "1|2|3|4|5,6,7"]
@test ct.partition == ["a1,a2,a3,a4|c2|c1|b1", "a1|a2|a3|a4|b1,c1,c2"]

ht = res.hybrid_table
@test length(ht.blob) == 3
@test ht.support_hybrid ≈ [0.2, 0.2, 0.6]
@test ht.cluster == ["b1", "a1,a2,a3,a4", "b1,c1,c2"]
@test ht.node_from == [6, -7, 3]
@test ht.node_to == [8, 3, -7]
@test ht.edge == [13, 15, 15]
@test ht.cluster_num == ["5", "1,2,3,4", "5,6,7"]

sdt = res.bipartition_table
@test isempty(sdt.edge)

@testset "errors" begin
    @test_throws ArgumentError blobpartitions_support(nets, refnet; minimumblobdegree=2)
    @test_throws ArgumentError blobpartitions_support(HybridNetwork[], refnet)
    @test_throws ErrorException blobpartitions_support(nets, refnet; netweight=[1.0, 1.0])
    @test_throws ErrorException blobpartitions_support(nets, refnet; netweight=fill(-1.0, 5))
end

@testset "Tree case (no blobs)" begin
    tre = readnewick("(a1,(a2,(a3,(a4,(b1,(c1,c2))))));")
    res_tre = blobpartitions_support(nets, tre)
    @test isempty(res_tre.blob_table.blob)
    @test isempty(res_tre.circorder_table.blob)
    @test isempty(res_tre.hybrid_table.blob)
    @test length(res_tre.bipartition_table.edge) > 0
end

end
