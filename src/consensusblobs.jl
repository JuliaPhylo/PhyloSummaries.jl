#= treeofblobs
exitnodes_preindex ( maybe loop through edges for each exit node)
process_biconnectedcinoibebts\
descendentweight / hardwiredclusters for getting descendants of outgoing edges
warn if one exit node has 2 edge

ig exit node does not have two try to match the cycle ordering    add a note saying that in the future we will add functionality to maintain circular order if the net is level 1
=#

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

function ingest_blobs!(
    ::Val{N},
    blobs::Vector{BlobFreq{N,P} where P},
    net::PN.HybridNetwork,
    taxa::Vector{String};
    order_fn::Union{Nothing,Function}=nothing,
    hybrid_fn::Union{Nothing,Function}=nothing,
    reflect_invariant::Bool=true,
) where {N}
    PN.directedges!(net)
    M = PN.descendenceweight(net)              # may reset edge numbers; use after this point
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

        exit_edges = blobexitedges(net, blob, ii)
        @assert !isempty(exit_edges) "blob $ii has no exit edge"

        masks_vec = NTuple{N,Bool}[]
        for e in exit_edges
            en = Int(e.number)
            push!(masks_vec, ntuple(i -> M[leafrow[i], en] > 0, N))
        end

        part_tuple, old2new = canonicalizepartition(masks_vec)
        P = length(part_tuple)

        ord_tuple = if order_fn === nothing
            ntuple(i->i, P)
        else
            _canon_cycle(_map_order(order_fn(net, blob, exit_edges, taxa), old2new), reflect_invariant)
        end
        hidx = hybrid_fn === nothing ? nothing :
               _map_hybrid(hybrid_fn(net, blob, exit_edges, taxa), old2new)

        _upsert_blob!(blobs, part_tuple, ord_tuple, hidx)
    end
    return blobs
end

function blobexitedges(net::PN.HybridNetwork, blob, blob_index::Int)
    outs = PN.exitnodes_preindex(blob)
    eds = PN.Edge[]
    for ni in outs
        v = net.vec_node[ni]
        candidates = PN.Edge[]
        for e in v.edge
            if e.inte1 != blob_index && PN.getparent(e) === v
                push!(candidates, e)
            end
        end
        if length(candidates) != 1
            @warn "exit node $(v.number) has $(length(candidates)) outgoing edges outside blob $blob_index"
        end
        !isempty(candidates) && push!(eds, candidates[1])
    end
    return eds
end

_edge_mask(::Val{N}, net::PN.HybridNetwork, e, taxa::Vector{String}) where {N} =
    ntuple(i->Bool(PN.hardwiredcluster(e, taxa)[i]), N)

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


_map_order(order_vec::AbstractVector{<:Integer}, old2new::Vector{Int}) =
    NTuple{length(order_vec),Int}(old2new[i] for i in order_vec)

_map_hybrid(hidx::Integer, old2new::Vector{Int}) = old2new[Int(hidx)]

function _canon_cycle(ord::NTuple{P,Int}, reflect::Bool) where {P}
    best = _min_rotation(ord)
    if reflect
        rev = ntuple(i->ord[P - i + 1], P)
        best = min(best, _min_rotation(rev))
    end
    return best
end

function _min_rotation(t::NTuple{P,Int}) where {P}
    best = t
    for s in 1:P-1
        rot = ntuple(i->t[(i + s - 1) % P + 1], P)
        if rot < best
            best = rot
        end
    end
    return best
end

function _upsert_blob!(
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

function find_hybrid_idx(net::PN.HybridNetwork, blob, exit_edges, taxa)::Int
    hn = nothing
    for e in blob.edges
        if e.hybrid
            hn = PN.getchild(e)
            break
        end
    end
    @assert hn !== nothing "no hybrid edge found in blob"
    clade_edge = nothing
    for e in hn.edge
        if PN.getparent(e) === hn
            clade_edge = e
            break
        end
    end
    @assert clade_edge !== nothing "no child edge found for hybrid node"
    hc = PN.hardwiredcluster(clade_edge, taxa)
    for i in 1:length(exit_edges)
        if PN.hardwiredcluster(exit_edges[i], taxa) == hc
            return i
        end
    end
    error("could not match hybrid clade to any exit")
end

function order_from_major(net::PN.HybridNetwork, blob, exit_edges, taxa)
    P = length(exit_edges)
    h = find_hybrid_idx(net, blob, exit_edges, taxa)
    return vcat(h:P, 1:h-1)
end

function consensusblobs(networks::AbstractVector{PN.HybridNetwork})
    taxa = sort!(PN.tiplabels(networks[1]))
    N = length(taxa)
    blobs = Vector{BlobFreq{N,P} where P}()
    for net in networks
        ingest_blobs!(Val(N), blobs, net, taxa; order_fn=order_from_major, hybrid_fn=find_hybrid_idx, reflect_invariant=false)
    end
    return blobs
end
