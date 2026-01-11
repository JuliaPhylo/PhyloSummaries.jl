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
    blobcompatible(a::NTuple{P,NTuple{N,Bool}}, b::NTuple{N,Bool})
    blobcompatible(a::NTuple{P,NTuple{N,Bool}}, b::NTuple{P,NTuple{N,Bool}})

fixit / todo: implement
fixit: explain. unrooted.
"""
blobcompatible(p1,p2) = true

