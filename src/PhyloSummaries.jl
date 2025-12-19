module PhyloSummaries

using Dictionaries
using PhyloNetworks

const PN = PhyloNetworks

export
consensustree,
consensus_treeofblobs

include("consensustrees.jl")
include("consensusblobs.jl")
end
