#= treeofblobs
exitnodes_preindex ( maybe loop through edges for each exit node)
process_biconnectedcinoibebts\
descendentweight / hardwiredclusters for getting descendants of outgoing edges
warn if one exit node has 2 edge

ig exit node does not have two try to match the cycle ordering
add a note saying that in the future we will add functionality to maintain circular order if the net is level 1
=#

#TODO: parse in circular order and store info that way. check funcs after canonicalizepartition

"""
    consensusblobs(networks; proportion=0)

Construct a consensus network from a collection of networks by identifying common blob structures.
Returns a single network representing the consensus of the input networks.
"""
function consensusblobs(
    networks::AbstractVector{PN.HybridNetwork};
    proportion::Number=0,
)
    isempty(networks) &&
        throw(ArgumentError("consensusblobs requires at least one network"))
    
    if length(networks) == 1
        net = deepcopy(networks[1])
        return net
    end
    
    taxa = sort!(tiplabels(networks[1]))
    blobcounts = Dictionary{Tuple{Vararg{BitSet}},Int}()
    
    for net in networks
        length(net.leaf) == length(taxa) ||
            throw(ArgumentError("input networks do not share the same taxon set"))
        count_blobs!(blobcounts, net, taxa)
    end
    
    nnets = length(networks)
    consensus_blobs!(blobcounts, proportion, nnets)
    return network_from_blobs(taxa, blobcounts, nnets)
end


struct BlobFreq{N,P}
    partition::NTuple{P,NTuple{N,Int}}
    count::Int
    circorder::Dict{NTuple{P,Int},Int}
    hybrid::Vector{Int}
end

# array of blob frequencies
blobs = Vector{BlobFreq{N,P} where P}

function parseblob(
    blob,
    blobarray,
    hw_matrix
)

end
function find_matching_blob(blobs, splits)
    for (i, bf_any) in pairs(blobs)
        partition = bf_any.partition

        # same number of blocks?
        length(partition) == length(splits) || continue

        # same tuples, ignoring order
        if Set(partition) == Set(splits)
            return 
        end
    end
    return nothing
end
"""
    traverse_blob_nodes(
        net, blob, blob_index, edge_map, hw_matrix, taxa_cols, rooted
    ) 
"""
function traverse_blob_splits(
    net::PN.HybridNetwork,
    blob,
    blob_index::Int,
    edge_map::Dictionary{Int,Int},  
    hw_matrix,                      
    taxa_cols,                      
)

    for v in net.node
        v.booln5 = false
    end

    entry_idx  = PN.entrynode_preindex(blob)
    entry_node = net.vec_node[entry_idx]

    splits = Tuple[]
    hybrids = Tuple[]
    # Traverse the blob starting at entry_node
    visit!(entry_node,blob_index,edge_map,hw_matrix,taxa_cols,splits, hybrids)
    return splits, hybrids
end
function visit!(
    node::PN.Node,
    blob_index::Int,
    edge_map,
    hw_matrix,
    taxa_cols,
    splits, 
    hybrids,
)
    
    if node.hybrid
        if node.booln5
            return
        end
        node.booln5 = true

    end
    for e in node.edge
        PN.getparent(e) === node || continue

        if e.inte1 != blob_index
            edgenum     = Int(e.number)
            row_idx     = edge_map[edgenum]              # row in hw_matrix

            split = Tuple(view(hw_matrix, row_idx, taxa_cols))
            push!(splits, split)
            if node.hybrid
                push!(hybrids, length(splits))
            end

        end
    end


    for e in node.edge
        PN.getparent(e) === node || continue


        e.inte1 == blob_index || continue

        child = PN.getchild(e)
        visit!(child, blob_index, edge_map, hw_matrix, taxa_cols, splits, hybrids)
    end

    return nothing
end
