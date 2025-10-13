
"""
    consensus_tree(networks::AbstractVector{PN.HybridNetwork}, rooted::Bool)

Construct a consensus network from the supplied phylogenetic networks. Throws an
`ArgumentError` when no networks are provided or there exists a network with different tiplables. When a 
single network is given,returns a deep copy of it.
# Assumptions:
- All networks share the same taxa set.
- All inputs are trees.
# Arguments:
- `networks::AbstractVector{PN.HybridNetwork}`: A vector of phylogenetic networks.
- `rooted::Bool`: Whether to treat the networks as rooted or unrooted when
  extracting bipartitions.
"""
function consensus_tree(networks::AbstractVector{PN.HybridNetwork},  rooted::Bool )
    isempty(networks) && throw(ArgumentError("consensus_tree requires at least one network"))
    if length(networks) == 1
        return deepcopy(networks[1])
    end


    taxa = sort(PN.tiplabels(networks[1]))

    counts = Dict{Tuple{Vararg{Int}},Int}()

    for net in networks
        if length(PN.tiplabels(net)) != length(taxa)
            throw(ArgumentError("networks do not share the same taxa set"))
        end
        extract_bipartitions!(counts, net, taxa, rooted)
    end

    return consensus_bipartition() #TODO
end

"""
    extract_bipartitions!(counts, net, taxa, rooted)

Accumulate bipartitions from `net` into `counts`. Hardwired clusters are
generated from a defensive copy of the network so the original structure stays
untouched. Only valid splits contribute to the count totals. If tiplabels do not
match the taxa set, PN.hardwiredclusters will throw an error.

"""
function extract_bipartitions!(
    counts::Dict{Tuple{Vararg{Int}},Int},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    rooted::Bool
    )
    
    net_copy = deepcopy(net)
    PN.directedges!(net_copy)

    hw_matrix = PN.hardwiredclusters(net_copy, taxa)
    taxa_len = length(taxa)
    taxa_cols = 2:(taxa_len + 1)

    for row_idx in axes(hw_matrix, 1)
        edge_number = hw_matrix[row_idx, 1]
        edge_number > 0 || continue
        edge_kind = hw_matrix[row_idx, end]
        edge_kind == 10 || continue # only consider tree edges

        split = bipartition_from_cluster(view(hw_matrix, row_idx, taxa_cols), rooted, taxa)
        split === nothing && continue

        counts[split] = get(counts, split, 0) + 1

    end

    return nothing
end

"""
    bipartition_from_cluster(cluster, rooted, taxa)

Transform a hardwired cluster indicator vector into a tuple key usable inside
`counts`. Returns `nothing` when the cluster does not represent a proper split.
When `rooted` is false, if the root is included in the cluster, the
partition is inverted to exclude it.
"""
function bipartition_from_cluster(cluster::AbstractVector, rooted::Bool, taxa)
    partition = copy(cluster)
    if !rooted && partition[length(taxa)] == 1
        partition .= 1 .- partition
    end
    
    ones = count(==(1), partition)
    zeros = length(partition) - ones
    (ones <= 1 || zeros <= 1) && return nothing
    return Tuple(partition)
end

"""
    consensus_bipartition()

Placeholder for assembling the consensus network from accumulated bipartition
frequencies. 
"""
function consensus_bipartition()
    return nothing
end
