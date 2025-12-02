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
    hybrid::Dict{Int,Int}
end

# fixit: problem below, we do not want this global variable in PhyloSummaries' namespace
# array of blob frequencies
blobs = Vector{BlobFreq{N,P} where P}

function parseblob(
    blob,
    blobarray,
    hw_matrix,
    edge_map,
)

    splits, hybrids = traverse_blob_splits(
        net,
        blob,
        blob_index,
        edge_map,
        hw_matrix,
        taxa_cols,
    )
    # canonicalize splits to get partition

    # check if this partition already exists in blobarray
    #currently assuming only level 1
    blob_idx, idxmap = find_matching_blob(blobarray, splits)
    if blob_idx == -1
        # new blob
        circorder = Dict{NTuple{length(partition),Int},Int}()
        for (i, split) in pairs(splits)
            circorder[split] = i
        end
        hybridmap = Dict{Int,Int}()
        hybridmap[hybrids[1]] = 1
        newblob = BlobFreq{length(partition),length(splits)}(
            partition,
            1,
            circorder,
            hybridmap,
        )
        push!(blobarray, newblob)
    else
        # existing blob, increment count
        blobarray[blob_idx].count += 1

        # update circorder and hybrid
    end


end
function find_matching_blob(blobs, splits)
    for (i, blob) in pairs(blobs)
        partition = blob.partition
        length(partition) == length(splits) || continue

        idxmap = Vector{Int}(undef, length(splits))
        used = falses(length(partition))              # tracks which partition slots are taken
        equalblob = true

        for (k, s) in pairs(splits)
            pos = findfirst(isequal(s), partition)
            if pos === nothing || used[pos]
                equalblob = false
                break
            end
            idxmap[k] = pos
            used[pos] = true
        end

        equalblob && return i, idxmap
    end
    return -1, Int[]
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
