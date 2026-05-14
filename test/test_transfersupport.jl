@testset "blobsupport_transfer, 7 taxa level-1" begin

nfile = joinpath(@__DIR__, "..", "test", "level1_7taxa_abc.nwk")
# nfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","level1_7taxa_abc.nwk")
nets = readmultinewick(nfile)
refnet = deepcopy(nets[5]) # 2 blobs: a1|a2|a3|a4|b1,c1,c2 and a1,a2,a3,a4|c2|c1|b1

res = PS.blobsupport_transfer(nets, refnet)
@test length(res.refblobs) == 2
# refblob 1: a1|a2|a3|a4|b1,c1,c2 — matches nets 1,2 (dist 0), nets 3,4 differ
# refblob 2: a1,a2,a3,a4|c2|c1|b1 — matches nets 3 (close), nets 1,2 differ
@test res.transfer_index == [0 0 10 8 0; 10 10 0 2 0]

end
