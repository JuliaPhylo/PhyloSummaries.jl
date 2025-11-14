#= treeofblobs
exitnodes_preindex ( maybe loop through edges for each exit node)
process_biconnectedcinoibebts\
descendentweight / hardwiredclusters for getting descendants of outgoing edges
warn if one exit node has 2 edge

ig exit node does not have two try to match the cycle ordering    add a note saying that in the future we will add functionality to maintain circular order if the net is level 1
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
    partition::NTuple{P,NTuple{N,Bool}}
    count::Int
    circorder::Dict{NTuple{P,Int},Int}
    hybrid::Vector{Int}
end

"""
    processblobs!(::Val{N}, blobs, net, taxa)

Extract blob information from a network and add it to the blobs collection.
Processes biconnected components and records partition patterns with their frequencies.
"""
function processblobs!(
    ::Val{N},
    blobs::Vector{BlobFreq{N,P} where P},
    net::PN.HybridNetwork,
    taxa::Vector{String},
) where {N}
    PN.directedges!(net)
    M = PN.descendenceweight(net)
    PN.process_biconnectedcomponents!(net)

    leafrow = Vector{Int}(undef, N)
    for i in 1:N
        name = taxa[i]
        idx = findfirst(j -> (net.vec_node[j].leaf && net.vec_node[j].name == name), eachindex(net.vec_node))
        @assert idx !== nothing "taxon $(name) not found as a leaf node"
        leafrow[i] = Int(idx)
    end

    for (ii, blob) in pairs(net.partition)
        PN.istrivial(blob) && continue

        exit_edges = blobexitedges(net, blob, ii)  # already circular; first is the hybrid clade
        @assert !isempty(exit_edges) "blob $ii has no exit edge"

        masks_vec = NTuple{N,Bool}[]
        for e in exit_edges
            en = Int(e.number)
            push!(masks_vec, ntuple(i -> M[leafrow[i], en] > 0, N))
        end

        part_tuple, old2new = canonicalizepartition(masks_vec)
        P = length(part_tuple)

        ord_tuple = NTuple{P,Int}(old2new[i] for i in 1:P)
        h_new     = old2new[1]

        insertblob!(blobs, part_tuple, ord_tuple, h_new)
    end
    return blobs
end

"""
    blobexitedges(net, blob, blob_index)

Get the ordered list of edges exiting a blob, starting with the hybrid clade edge.
Returns edges in circular order around the blob.
"""
function blobexitedges(net::PN.HybridNetwork, blob, blob_index::Int)
    out_nodes = blobexitnodes(net, blob, blob_index)  
    res = PN.Edge[]

    for n in out_nodes
        # edge that goes out of the blob from this exit node
        candidates = [e for e in n.edge if e.inte1 != blob_index && PN.getparent(e) === n]
        if length(candidates) != 1
            @warn "exit node $(n.number) has $(length(candidates)) external outgoing edges for blob $blob_index"
        end
        if !isempty(candidates)
            push!(res, candidates[1])
        end
    end
    return res
end


"""
    blobexitnodes(net, blob, blob_index)

Get the ordered list of nodes that are exit points from a blob.
Traverses the cycle starting from the hybrid clade exit node.
"""
function blobexitnodes(net::PN.HybridNetwork, blob, blob_index::Int)
    outs = Set(net.vec_node[i] for i in PN.exitnodes_preindex(blob))  # exit nodes (as Nodes)
    res  = PN.Node[]  # ordered exit nodes


    hn = nothing
    for e in blob.edges
        if e.hybrid
            hn = PN.getchild(e)
            break
        end
    end
    @assert hn !== nothing "no hybrid edge found in blob"
    @assert hn in outs "assumption violated: hybrid clade is not an exit node"


    start_parent = nothing
    for e in hn.edge
        if PN.getchild(e) === hn && e.inte1 == blob_index && e.ismajor
            start_parent = PN.getparent(e)
            break
        end
    end
    @assert start_parent !== nothing "no major parent on the blob's cycle"


    push!(res, hn)               # start at the hybrid-clade exit (your assumption)
    prev_node = hn
    cur_node  = start_parent


    while cur_node !== hn
        if cur_node in outs
            if isempty(res) || res[end] !== cur_node
                push!(res, cur_node)
            end
        end

        next_node = nothing
        for e in cur_node.edge
            if e.inte1 == blob_index
                nb = (PN.getchild(e) === cur_node) ? PN.getparent(e) : PN.getchild(e)
                if nb !== prev_node
                    next_node = nb
                    break
                end
            end
        end
        @assert next_node !== nothing "broken cycle traversal in blob $blob_index"
        prev_node, cur_node = cur_node, next_node

    end
    return res
end

"""
    canonicalizepartition(masks)

Canonicalize a partition by sorting masks and returning the sorted partition 
along with the mapping from old to new positions.
"""
function canonicalizepartition(masks::Vector{NTuple{N,Bool}}) where {N}
    @assert all(any, masks) "exit mask with no true taxa"
    idxs = collect(1:length(masks))
    sort!(idxs; by = i -> (findfirst(true, masks[i])::Int, masks[i]))
    newpos = zeros(Int, length(masks))
    for (newi, oldi) in pairs(idxs)
        newpos[oldi] = newi
    end
    return (ntuple(i -> masks[idxs[i]], length(masks)), newpos)
end


"""
    insertblob!(blobs, partition, ord, hybrid_idx)

Insert or update a blob frequency record in the blobs collection.
Tracks partition patterns, circular orderings, and hybrid clade positions.
"""
function  insertblob!(
    blobs::Vector{BlobFreq{N,P} where P},
    partition::NTuple{P,NTuple{N,Bool}},
    ord::NTuple{P,Int},
    hybrid_idx::Union{Nothing,Int},
) where {N,P}
    for i in eachindex(blobs)
        if blobs[i] isa BlobFreq{N,P} && (blobs[i]::BlobFreq{N,P}).partition == partition
            bf = blobs[i]::BlobFreq{N,P}
            cc = bf.count + 1
            oh = bf.circorder
            oh[ord] = get(oh, ord, 0) + 1
            hv = bf.hybrid
            if hybrid_idx !== nothing
                idx = Int(hybrid_idx)
                if idx >= 1 && idx <= length(hv)
                    hv[idx] += 1
                end
            end
            blobs[i] = BlobFreq{N,P}(partition, cc, oh, hv)
            return
        end
    end
    oh = Dict{NTuple{P,Int},Int}()
    oh[ord] = 1
    hv = fill(0, P)
    if hybrid_idx !== nothing
        idx = Int(hybrid_idx)
        if idx >= 1 && idx <= P
            hv[idx] = 1
        end
    end
    push!(blobs, BlobFreq{N,P}(partition, 1, oh, hv))
    return
end

