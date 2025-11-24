@testset "consensus blobs" begin

nfile = joinpath(@__DIR__,"..","test","bootstrapnets_h1.nwk")
# tfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","bootstrapnets_h1.nwk")
net = readmultinewick(nfile)
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
