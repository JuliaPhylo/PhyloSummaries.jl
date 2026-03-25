using PhyloNetworks
using PhyloSummaries

taxa = ["a1", "a2", "a3", "a4", "b1", "c1", "c2"]

nwk = String[]
for _ in 1:9 push!(nwk, "((((a1,a2,a3,a4))#H1, c1), (c2, (b1, #H1)));") end
for _ in 1:6 push!(nwk, "(((b1)#H1, c1), (c2, (((a1,a2,a3,a4)), #H1)));") end

for _ in 1:8 push!(nwk, "(((((b1,c1,c2))#H2, a1), a2), (a3, (a4, #H2)));") end
for _ in 1:2 push!(nwk, "((((a2)#H2, a1), ((b1,c1,c2))), (a3, (a4, #H2)));") end

net = readnewick.(nwk)
res = consensus_level1network(net, minimumblobdegree=3, proportion=0.01, suppressinfo=true)

println("Hybrid clusters found in consensus network: ", res[:blob_table].hybrid_cluster)
