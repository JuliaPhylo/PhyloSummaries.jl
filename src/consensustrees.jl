
"""
    consensustree(networks::AbstractVector{PN.HybridNetwork};
                  rooted=false, proportion=0.5)

Consensus tree summarizing the bipartitions (or clades) shared by more than
the required `proportion` of input `networks`.
An `ArgumentError` is thrown if one input network is not a tree, or the list of
input trees is empty, or if the input trees do not all have the same tip labels.
When a single tree is given, the result is a deep copy of it.

By default, input trees are considered unrooted, and bipartitions are considered.
Use `rooted=true` to consider all input trees as rooted, in which case clades
(rather than bipartitions) are used to build the output rooted consensus tree.

By default, the majority-rule consensus is calculated: the tree is built from
the bipartitions (or clades) present in more than 50% of the input trees.
The greedy consensus tree can be obtained by using `proportion=0`.

# example

fixit / todo: add a jldoctest block, to serve both as an example and as a unit test
(during documentation build).
"""
function consensustree(
    networks::AbstractVector{PN.HybridNetwork};
    rooted::Bool=false,
    proportion::Number=0.5,
)
    isempty(networks) &&
        throw(ArgumentError("consensustree requires at least one network"))
    all(n.numhybrids==0 for n in networks) ||
        throw(ArgumentError("consensustree requires input trees (without reticulations)"))
    if length(networks) == 1
        return deepcopy(networks[1])
    end
    taxa = sort!(tiplabels(networks[1]))

    splitcounts = Dict{Tuple{Vararg{Int}},Int}()
    for net in networks
        length(net.leaf) == length(taxa) ||
            throw(ArgumentError("input trees do not share the same taxon set"))
        count_bipartitions!(splitcounts, net, taxa, rooted)
    end

    return consensus_bipartition(splitcounts, proportion) #TODO
end

"""
    count_bipartitions!(counts, net, taxa, rooted)

Count bipartitions in `net`, and add them to `counts`.
If a new bipartition is found, a new key is added to `counts` with value 1.
If `net` has a bipartition already present as a key in `counts`, then the
corresponding value is incremented by 1.
`net` is not modified.

todo / fixit: explain `rooted`.
Only valid splits contribute to the count totals. (fixit: what does 'valid' mean here?)
If the tip labels in `net` do not match those in `taxa` (as a set), then an
error will be thrown indirectly (via `PhyloNetworks.hardwiredclusters`).
"""
function count_bipartitions!(
    counts::Dict{Tuple{Vararg{Int}},Int},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    rooted::Bool
)
    net_copy = deepcopy(net) # can we avoid this? this will cause a slow down and extra memory usage
    directedges!(net_copy) # fixit: add an argument to have the option to avoid doing this

    hw_matrix = hardwiredclusters(net_copy, taxa)
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
