"""
    startree(taxa)

Star tree `t` such that `t.node` lists all leaves named after the N `taxa`
in that order: for `i` in 1 through N: node `n = t.node[i]` has `n.number = i`
and `n.name = taxa[i]`.
The root node is last, `t.node[N+1]`, and has number `-2`.
"""
function startree(taxa::Vector{String})
    net = PN.HybridNetwork()
    root = PN.Node(-2,false) # root has number -2
    # push leaves first, numbered 1:n. we will have taxa[i] = label of net.node[i]
    for (i,t) in enumerate(taxa)
        edge = PN.Edge(i,-1.0) # ischild1 is true by default. length=-1 for NA
        leaf = PN.Node(i,true,false, # true: leaf
            -1.,[edge],false,false,false,false,false,false,-1,nothing,-1,-1,t)
        PN.pushNode!(net, leaf)
        edge.node = [leaf, root] # to match ischild1 is true
        PN.pushEdge!(net, edge)
        push!(root.edge, edge)
    end
    PN.pushNode!(net, root)
    net.rooti = length(taxa) + 1 # n leaves listed first, root listed next in net.node
    return net
end

"""
    istrivialsplit(v)

true/false if `v` does / does not represent a trivial split. `v` should contain
booleans or 0/1 values. A split is trivial if it has:
all 0s, or all 1s, or a single 0, or a single 1.
"""
function istrivialsplit(v)
    N = length(v)
    n1 = sum(v)
    n0 = N - n1
    return n1 <= 1 || n0 <= 1
end

"""
    isredundantsplit(v, b)

true/false if `v` does / does not represent one taxon block of blob
(or multi-partition) `b`, with `v` considered as unrooted.
"""
function isredundantsplit(v, b)
    return (v ∈ b || .!v ∈ b) # any(isequal(v), b) || ...
end


"""
    treecompatible(a::NTuple{N,Bool}, b::NTuple{N,Bool})

true / false if two clusters `a` and `b` are / are not tree-compatible.

If A is the cluster of descendants of `a` (with `true` entries in `a`) and
if B is the cluster of descendant of `b`, these 2 clusters are tree-compatible
if there exists some tree that has both clusters.
This can be checked by the condition: A∩B is empty, or A⊆B, or B⊆A.
"""
function treecompatible(a::NTuple{N,Bool}, b::NTuple{N,Bool})::Bool where N
    inter = ntuple(i -> a[i] & b[i], length(a))
    if !any(inter)
        return true
    end
    if inter == a || inter == b
        return true
    end
    return false
end

"""
    blobcompatible(A, B)

Check if two partitions A and B are blob-compatible (unrooted).

Two partitions are compatible if there exists a pair (i, j) such that all parts
of A except Aᵢ are subsets of Bⱼ. This reduces to tree-compatibility when both
A and B are bipartitions.


"""
function blobcompatible(
    A::NTuple{Ka, NTuple{N,Bool}},
    B::NTuple{Kb, NTuple{N,Bool}}
) where {N, Ka, Kb}
    for i in 1:Ka
        for j in 1:Kb
            all_included = true
            for q in 1:Ka
                q == i && continue  # skip A_i
                # Check if A_q ⊆ B_j (every taxon in A_q is also in B_j)
                if !issubsetsplit(A[q], B[j])
                    all_included = false
                    break
                end
            end
            if all_included
                return true
            end
        end
    end
    return false
end

# Edge case: empty partitions are trivially compatible
blobcompatible(::Tuple{}, ::Tuple{}) = true
blobcompatible(::Tuple{}, ::NTuple{Kb, NTuple{N,Bool}}) where {N, Kb} = true
blobcompatible(::NTuple{Ka, NTuple{N,Bool}}, ::Tuple{}) where {N, Ka} = true

"""
    issubsetsplit(a, b)

Check if split `a` is a subset of split `b` (all taxa in `a` are also in `b`).
"""
function issubsetsplit(a::NTuple{N,Bool}, b::NTuple{N,Bool}) where N
    inter = ntuple(i -> a[i] & b[i], N)
    # return !any(inter) || inter == a # should we check disjoint as  unrooted?
    return inter == a 
end

function blobcompatible(
    split::NTuple{N,Bool},
    B::NTuple{Kb, NTuple{N,Bool}}
) where {N, Kb}
    splitcomplement = ntuple(i -> !split[i], N)
    A = (split, splitcomplement)
    return blobcompatible(A, B)
end


function blobcompatible(
    A::NTuple{Ka, NTuple{N,Bool}},
    split::NTuple{N,Bool}
) where {N, Ka}
    splitcomplement = ntuple(i -> !split[i], N)
    B = (split, splitcomplement)
    return blobcompatible(A, B)
end

function blobcompatible(
    split1::NTuple{N,Bool},
    split2::NTuple{N,Bool}
) where N
    return treecompatible(split1, split2)
end
