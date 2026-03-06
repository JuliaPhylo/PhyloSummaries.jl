# test blobpartitions_support using level1_7taxa_abc.nwk
@testset "blobpartitions_support" begin

nets = readmultinewick(joinpath(@__DIR__, "level1_7taxa_abc.nwk"))
@test length(nets) == 5

# network 5 has 2 blobs: interesting reference
# currently throws BoundsError: intn1 != blob index in refblobs (fixit)
refnet = nets[5]
# @test_throws BoundsError blobpartitions_support(nets, refnet)

 
res = blobpartitions_support(nets, refnet)
@test res.taxa == ["a1", "a2", "a3", "a4", "b1", "c1", "c2"]

bt = res.blob_table
@test length(bt.blob) == 2
@test bt.degree == [4, 5]
@test bt.partition == ["a1,a2,a3,a4|c2|c1|b1", "a1|a2|a3|a4|b1,c1,c2"]
@test bt.support_partition ≈ [0.4, 0.6]

ct = res.circorder_table
@test length(ct.blob) == 2
@test ct.support_circorder ≈ [0.4, 0.6]

ht = res.hybrid_table
@test length(ht.blob) == 3
@test ht.support_hybrid ≈ [0.2, 0.2, 0.6]
@test ht.cluster == ["b1", "a1,a2,a3,a4", "b1,c1,c2"]

sdt = res.bipartition_table
@test isempty(sdt.edge)

# 2-blob network: .intn1 != blob index in refblobs
# blob 2 has intn1=8 (bicomp index) but is refblobs[2].
# fails with BoundsError if using b_i (=intn1) directly to index refblobs.
refnet2 = readnewick("(((A,(B)#H1),(#H1,C)),((D,(E)#H2),(#H2,F)));")
taxa2 = sort(tiplabels(refnet2))
refblobs = PhyloSummaries.BlobFreq{6}[]
refbps = PhyloSummaries.SplitFreq{6}[]
hwmatrix, edgemap, blobdegree = PhyloSummaries.count_blobpartitions!(refblobs, refbps, refnet2, taxa2, 3, false, 0.0)
@info refblobs
@info refnet2.partition
# test isolated _blobnode_blobedges function for the BoundsError fix
bbn, bbei_nums = PhyloSummaries._blobnode_blobedges(refnet2, hwmatrix, edgemap, refblobs, taxa2, blobdegree, 3)
@test length(bbn) == 2
@test length(bbei_nums) == 2

end
