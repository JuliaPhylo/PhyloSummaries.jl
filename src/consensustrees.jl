
"""
    consensustree(trees::AbstractVector{PN.HybridNetwork};
                  rooted=false, proportion=0.5)

Consensus tree summarizing the bipartitions (or clades) shared by more than
the required `proportion` of input `trees`.
An `ArgumentError` is thrown if one input network is not a tree, or the list of
input trees is empty, or if the input trees do not all have the same tip labels.
Input trees are not modified.

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
        net = deepcopy(trees[1])
        # fixit: suppressroot!(net) when ready, from PN
        return 
    end
    taxa = sort!(tiplabels(trees[1]))

    splitcounts = Dict{BitVector,Int}()
    for net in trees
        length(net.leaf) == length(taxa) ||
            throw(ArgumentError("input trees do not share the same taxon set"))
        count_bipartitions!(splitcounts, net, taxa, rooted)
    end

    return consensus_bipartition(splitcounts, proportion, length(trees), taxa) #TODO
end

"""
    count_bipartitions!(counts, net, taxa, rooted)

Count bipartitions in `net`, and add them to `counts`.
If a new bipartition is found, a new key is added to `counts` with value 1.
If `net` has a bipartition already present as a key in `counts`, then the
corresponding value is incremented by 1.
`net` is not modified.

By default, input trees are considered unrooted, and bipartitions are considered.
Use `rooted=true` to consider all input trees as rooted, in which case clades
(rather than bipartitions) are used to build the output rooted consensus tree.

Only non-trivial (not all 0s or all 1s) splits contribute to the count totals. 

If the tip labels in `net` do not match those in `taxa` (as a set), then an
error will be thrown indirectly (via `PhyloNetworks.hardwiredclusters`).
"""
function count_bipartitions!(
    counts::Dict{BitVector,Int},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    rooted::Bool,
) 
    if !rooted
        rootdegree = length(getroot(net).edge)
        # fixit: replace code below by suppressroot!(net), after deepcopy if needed
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
        return (x==0 for x in cluster01vector)
    end
    return (x==1 for x in cluster01vector)
end

"""
    consensus_bipartition(splitcounts, proportion, numtrees, taxa)

Select all bipartitions that meet the frequency threshold and
assemble the final consensus tree.

Each bipartition in `splitcounts` is included if it occurs in at least
`ceil(proportion * numtrees)` of the input trees. The selected bipartitions
are passed to [`create_tree_from_bipartition_set`](@ref) to construct
the resulting consensus topology.

By default `proportion = 0.5`, this yields a majority-rule consensus tree.
Setting `proportion < 0.5` produces a greedy consensus tree that includes
all compatible bipartitions with frequency higher than.

# Arguments
- `splitcounts::Dict{BitVector,Int}`: Mapping of bipartition vectors to their counts.
- `proportion::Number`: Inclusion threshold (e.g., `0.5` for majority rule).
- `numtrees::Number`: Total number of input trees.
- `taxa::Vector{String}`: Ordered list of taxon labels corresponding to bipartition indices.

# Returns
A `PhyloNetworks.HybridNetwork` object representing the consensus tree.

# Example
```julia
splitcounts = Dict(BitVector([1,1,0,0]) => 3, BitVector([0,0,1,1]) => 3)
taxa = ["A","B","C","D"]
consensus_bipartition(splitcounts, 0.5, 4, taxa)
# → HybridNetwork for ((A,B),(C,D));
```
"""
function consensus_bipartition( splitcounts::Dict{BitVector,Int},
    proportion::Number, 
    numtrees::Number,
    taxa::Vector{String},
) 

    p = max(0.5, proportion)
    threshold = ceil(p * numtrees)

    final_bipartitions = Vector{BitVector}()


    for (bipartition, frequency) in splitcounts
        if frequency >= threshold
            push!(final_bipartitions, bipartition)
        end
    end
    for (bipartition, frequency) in splitcounts
        if frequency < threshold & frequency > ceil(proportion * numtrees)
            for final_bipartition in final_bipartitions
                if is_compatible(bipartition, final_bipartition)
                    push!(final_bipartitions, bipartition)
                end
            end
            
        end
    end



    return create_tree_from_bipartition_set(taxa, final_bipartitions)

end

"""
    is_compatible(a, b)
    Check if two bipartitions `a` and `b` (as `BitVector`s) are compatible.

"""
function is_compatible(a::BitVector, b::BitVector)::Bool
    @assert length(a) == length(b)
    inter = a .& b
    if !any(inter)                     
        return true
    end

    if inter == a || inter == b       
        return true
    end
    return false
end
"""
    create_tree_from_bipartition_set(taxa, bipartitions)

Construct a consensus tree topology from a compatible set of bipartitions.

Each bipartition (represented as a `BitVector`) is converted into a cluster
of taxa and combined hierarchically to reconstruct the tree structure.
Clusters are processed in ascending order of size, ensuring that smaller
subtrees are nested before larger clusters.  
Trivial clusters (size 1 or equal to the number of taxa) are ignored.

This function assumes that the bipartitions provided are mutually compatible
and collectively represent a valid tree. 

# Arguments
- `taxa::Vector{String}`: Ordered list of taxon labels corresponding to the bit positions.
- `bipartitions::Vector{BitVector}`: Compatible bipartitions representing clusters of taxa.

# Returns
A `PhyloNetworks.HybridNetwork` object representing the reconstructed consensus tree.

# Example
```julia
taxa = ["A","B","C","D"]
bipartitions = [BitVector([1,1,0,0]), BitVector([0,0,1,1])]
tree = create_tree_from_bipartition_set(taxa, bipartitions)
# → HybridNetwork for ((A,B),(C,D));
```
"""
function create_tree_from_bipartition_set(taxa::Vector{String}, bipartitions::Vector{BitVector})
    n = length(taxa)


    net = PN.HybridNetwork()
    root = PN.Node(-2,false) # root has number -2
    PN.pushNode!(net, root)
    net.rooti = 1

    # leaves have numbers 1:n
    leaf_nodes = Dict{String, PN.Node}()
    for (i,t) in enumerate(taxa)
        edge = PN.Edge(i) # ischild1 is true by default
        leaf = PN.Node(i,true,false, # true: leaf
            -1.,edge,false,false,false,false,false,false,-1,nothing,-1,-1,t)
        PN.pushNode!(net, leaf)
        leaf_nodes[t] = leaf
        edge.node = [leaf, root] # to match ischild1 is true
        PN.pushEdge!(net, edge)
        push!(root.edge, edge)
    end

    node_counter = n+1 # alternatively: start at -3 and go down -= 1
    for bv in bipartitions
        internal = PN.Node(node_counter,false)
        node_counter += 1
        PN.pushNode!(net, internal)


        for i in eachindex(bv)
            if bv[i] == 1
                leaf = leaf_nodes[taxa[i]]
                ## fix dont remove parent edge. change at LCA
                parent_edge = getparentedge(leaf)
                if parent_edge !== nothing
                    parent = parent_edge.node[1]
                    PN.deleteEdge!(net, parent_edge)
                    edge_counter -= 1
                end

                e = PN.Edge(edge_counter)
                e.node = [leaf, internal]
                e.ischild1 = true
                PN.pushEdge!(net, e)
                push!(internal.edge, e)
                push!(leaf.edge, e)
            end
        end


    end

    PN.directEdges!(net)
    return net
end


