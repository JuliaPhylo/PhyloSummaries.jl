module PhyloSummaries

using PhyloNetworks

# import: to extend methods from another packages
import Base: show

const PN = PhyloNetworks

export
consensustree,
consensus_treeofblobs,
consensus_level1network

include("utils.jl")
include("consensustrees.jl")
include("consensusblobs.jl")
end
