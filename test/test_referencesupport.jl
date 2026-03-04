# test blobpartitions_support using level1_7taxa_abc.nwk
@testset "blobpartitions_support" begin

nets = readmultinewick(joinpath(@__DIR__, "level1_7taxa_abc.nwk"))
@test length(nets) == 5

# network 5 has 2 blobs: interesting reference
refnet = nets[5]
res = blobpartitions_support(nets, refnet)

@test res.taxa == ["a1", "a2", "a3", "a4", "b1", "c1", "c2"]

# --- blob_table ---
bt = res.blob_table
@test length(bt.blob) == 2
@test bt.degree == [4, 5]
@test bt.partition == ["a1,a2,a3,a4|c2|c1|b1", "a1|a2|a3|a4|b1,c1,c2"]
@test bt.support_partition ≈ [0.4, 0.6]

# --- circorder_table ---
ct = res.circorder_table
@test length(ct.blob) == 2
@test ct.support_circorder ≈ [0.4, 0.6]

# --- hybrid_table ---
ht = res.hybrid_table
@test length(ht.blob) == 3
@test ht.support_hybrid ≈ [0.2, 0.2, 0.6]
@test ht.cluster == ["b1", "a1,a2,a3,a4", "b1,c1,c2"]

# --- bipartition_table ---
sdt = res.bipartition_table
@test isempty(sdt.edge)

end
