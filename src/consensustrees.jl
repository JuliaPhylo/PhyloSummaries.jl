
"""
    consensustree(trees::AbstractVector{PN.HybridNetwork};
                  rooted=false, proportion=0.5)

Consensus tree summarizing the bipartitions (or clades) shared by more than
the required `proportion` of input `trees`.
An `ArgumentError` is thrown if one input network is not a tree, or the list of
input trees is empty, or if the input trees do not all have the same tip labels.
When a single tree is given, the result is a deep copy of it.

By default, input trees are considered unrooted, and bipartitions are considered.
Use `rooted=true` to consider all input trees as rooted, in which case clades
(rather than bipartitions) are used to build the output rooted consensus tree.

By default, the majority-rule consensus is calculated: the tree is built from
the bipartitions (or clades) present in more than 50% of the input trees.
The greedy consensus tree can be obtained by using `proportion=0`.

assumptions and **warnings**:
- Input trees are assumed to have their edges correctly directed.
  If unsure, run `directedges!.(trees)` prior.
- Input trees should not have degree-2 nodes other than the root
  (nodes with 1 only parent and 1 child).
  If unsure, run `removedegree2nodes!.(trees, true))` to keep their root even
  of degree 2 or `removedegree2nodes!.(trees, false))` unroot them also.

# example

fixit / todo: add a jldoctest block, to serve both as an example and as a unit test
(during documentation build).
"""
function consensustree(
    trees::AbstractVector{PN.HybridNetwork};
    rooted::Bool=false,
    proportion::Number=0.5,
)
    isempty(trees) &&
        throw(ArgumentError("consensustree requires at least one network"))
    all(n.numhybrids==0 for n in trees) ||
        throw(ArgumentError("consensustree requires input trees (without reticulations)"))
    if length(trees) == 1
        return deepcopy(trees[1])
    end
    taxa = sort!(tiplabels(trees[1]))

    splitcounts = Dict{NTuple{length(taxa),Int},Int}()
    for net in trees
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
    counts::Dict{NTuple{N,Int},Int},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    rooted::Bool,
) where N
    if !rooted
        rootdegree = length(getroot(net).edge)
        if rootdegree == 2
            net = deepcopy(net) # re-binds the variable 'net'
            PN.fuseedgesat!(net.rooti, net)
        elseif rootdegree == 1 # should almost never happen
            net = deepcopy(net)
            deleteleaf!(net, net.rooti; index=true)
            PN.fuseedgesat!(net.rooti, net)
        end
    end
    hw_matrix = hardwiredclusters(net, taxa)
    taxa_cols = 2:(length(taxa) + 1)

    for row_idx in axes(hw_matrix, 1)
        edge_number = hw_matrix[row_idx, 1]
        edge_number > 0 || continue # fixit: why? edge numbers are arbitrary. could be negative. they just need to be unique
        edge_kind = hw_matrix[row_idx, end]
        edge_kind == 10 || continue # only consider tree edges. fixit: consensustree already checked that 'net' is a tree

        split = tuple_from_clustervector(view(hw_matrix, row_idx, taxa_cols), rooted)
        isnothing(split) && continue
        counts[split] = get(counts, split, 0) + 1
    end
    return counts
end

"""
    tuple_from_clustervector(cluster01vector, rooted)

Tuple of booleans `t` with `t[i]` true / false if `taxa[i]` does / does not
belong in the hardwired cluster vector `cluster01vector` of 0/1 integers;
or `nothing` if the cluster is trivial (all 0s or all 1s).

If `rooted` is false, then clusters are considered as bipartitions and the last
taxon is used as outgroup with a `false` entry.
For example, clusters `0011` and `1100` represent the same bipartition, and
both would return tuple `(true,true,false,false)`.
"""
function tuple_from_clustervector(cluster01vector::AbstractVector, rooted::Bool)
    if all(isequal(1), cluster01vector) || all(isequal(0), cluster01vector)
        return nothing
    end
    if !rooted && cluster01vector[end] == 1
        return Tuple(x==0 for x in cluster01vector)
    end
    return Tuple(x==1 for x in cluster01vector)
end

"""
    consensus_bipartition()

Placeholder for assembling the consensus network from accumulated bipartition
frequencies. 
"""
function consensus_bipartition()
    return nothing
end
