const SplitTuple = NTuple{N,Bool} where N # tuple used to represent a bipartition
splitstring(obj::SplitTuple) = join(findall(obj),",")
splitstring_names(v, taxa::Vector) =
    join((taxa[j] for (j,b) in enumerate(v) if b), ",")

"""
    check_nonnumericnames(taxa)

Check that the vector of `taxa` contains names that do not start with "_"
followed by a positive integer, such as "_12".
The exception is that the ith taxon `taxa[i]` is allowed to be "_i".

Output: `nothing` if the conditions are met, otherwise the first index
`i` such that `taxa[i]` is `_n` with `n` ≠ `i`.
"""
function check_nonnumericnames(taxa::AbstractVector{<:AbstractString})
    r = r"^_([1-9]\d*)$"
    i_bad = nothing
    for (i,label) in enumerate(taxa)
        m = match(r, label)
        (isnothing(m) || tryparse(Int, m.captures[1]) == i) && continue
        i_bad = i
        break
    end
    return i_bad
end

"""
    startree(taxa)

Star tree `t` such that `t.node` lists all leaves named after the N `taxa`
in that order: for `i` in 1 through N: node `n = t.node[i]` has `n.number = i`
and `n.name = taxa[i]`.
It is incident to edge `e = t.edge[i]` which has `e.number = i`.
The root node is last, `t.node[N+1]`, has number `N+1` and name "_`N+1`".
"""
function startree(taxa::Vector{String})
    i = check_nonnumericnames(taxa)
    isnothing(i) ||
        error("""taxon names need be non-numeric (or ith taxon should be "i").
        $(i)th taxon is $(taxa[i]).
        Numeric names will be used for internal node labels here.""")
    net = PN.HybridNetwork()
    ntax = length(taxa)
    root = PN.Node(ntax+1, false) # root has number N+1
    root.name = "_$(ntax+1)"
    # push leaves first, numbered 1:N. we will have taxa[i] = label of net.node[i]
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
    net.rooti = ntax + 1 # N leaves listed first, root listed next in net.node
    return net
end

"""
    resetnodenumbers_fromnames!(net)

Reset node numbers so that:
- leaves are numbered 1 through the number of taxa, in alphabetical order of
  taxon names
- any internal node named `_n*` is numbered `n`, which should be an integer
  that is not starting with 0, and `*` is what comes after
- any hybrid node named `Hn` is numbered `n`
- any node with an empty name is given some other number, so that node numbers
  are unique.

An error occurs if an internal node has a name that cannot be parsed as
described above, or results in a number that is ≤ the number of taxa,
or results in duplicated numbers.

See also: `PhyloNetworks.resetedgenumbers!` and
`PhyloNetworks.resetnodenumbers!`.

# example

```jldoctest
julia> net = readnewick("((t2,(t1,(t4,#H11)_10))_8,(t3)#H11)_7_blob1;");

julia> printnodes(net)
node leaf  hybrid name     i_cycle edges'numbers
1    true  false  t2       -1      1   
2    true  false  t1       -1      2   
3    true  false  t4       -1      3   
5    false false  _10      -1      3    4    5   
-4   false false           -1      2    5    6   
6    false false  _8       -1      1    6    7   
7    true  false  t3       -1      8   
4    false true   H11      -1      8    4    9   
-2   false false  _7_blob1 -1      7    9   

julia> resetnodenumbers_fromnames!(net);

julia> printnodes(net)
node leaf  hybrid name     i_cycle edges'numbers
2    true  false  t2       -1      1   
1    true  false  t1       -1      2   
4    true  false  t4       -1      3   
10   false false  _10      -1      3    4    5   
5    false false           -1      2    5    6   
8    false false  _8       -1      1    6    7   
3    true  false  t3       -1      8   
11   false true   H11      -1      8    4    9   
7    false false  _7_blob1 -1      7    9   
```
"""
function resetnodenumbers_fromnames!(net::PN.HybridNetwork)
    taxa = sort!(tiplabels(net))
    leafname2num = Dict(t => i for (i,t) in enumerate(taxa))
    ntax = length(taxa)
    nnod = length(net.node)
    rx = r"^_([1-9]\d*)"
    rh = r"^H([1-9]\d*)"
    newnum = zeros(Int, nnod)
    for (i,n) in enumerate(net.node)
        if n.leaf
            newnum[i] = leafname2num[n.name]
        elseif n.name != ""
            m = match((n.hybrid ? rh : rx), n.name)
            isnothing(m) &&
                error("node $(n.name): cannot find desired number")
            nn = parse(Int, m.captures[1])
            nn > ntax || error("node $(n.name): $nn is ≤ number of taxa $ntax")
            any(newnum[j]==nn for j in 1:(i-1)) &&
                error("multiple nodes with desired number $nn")
            newnum[i] = nn
        end
    end
    # if some internal nodes had no names: assign them unused numbers
    i = findfirst(isequal(0), newnum)
    if !isnothing(i)
        unused = setdiff(1:nnod, newnum)
        nextj = (isempty(unused) ? nnod+1 : popfirst!(unused))
        while !isnothing(i)
            newnum[i] = nextj
            i = findnext(isequal(0), newnum, i+1)
            nextj = (isempty(unused) ? nextj+1 : popfirst!(unused))
        end
    end
    for i in 1:nnod
        net.node[i].number = newnum[i]
    end
    return net
end

"""
    edgenumbers_fromnodenumbers(n1_nums::Vector, n2_nums::Vector, net::PN.HybridNetwork)

Vector of edge numbers, of edges in `net` between pairs of nodes of specified numbers.
An error is thrown if, in one pair, the two nodes are not adjacent.
This is to avoid unknowingly using wrong edge numbers for the other pairs of nodes.

Uses `PhyloNetworks.getconnectingedge`.
"""
function edgenumbers_fromnodenumbers(n1num::AbstractVector, n2num::AbstractVector, net::PN.HybridNetwork)
    ne = length(n1num)
    length(n2num) == ne ||
        error("node number lists of different lengths: $ne and $(length(n2num))")
    ne == 0 && return Int[]
    enum = map(zip(n1num, n2num)) do (n1,n2)
        e = PN.getconnectingedge(n1, n2, net)
        isnothing(e) && error("no edge between nodes $n1 and $n2")
        return e.number
    end
    return enum
end
"""
    edgenumbers_fromnodenumbers(table, net::PN.HybridNetwork)

Vector of edge numbers in `net`, with element `i` corresponding to the edge connecting
(incident to) the two nodes listed in row `i` of the input table. This `table` may be
a data frame or a named tuple, such as produced by [`consensus_treeofblobs`](@ref) and
[`consensus_level1network`](@ref). It should have 2 columns with node numbers, either named
`:node_from` and `:node_to`, or `:node1` and `:node2`.
"""
function edgenumbers_fromnodenumbers(df, net::PN.HybridNetwork)
    colnames = propertynames(df) # works for NamedTuples and DataFrames
    if (:node_from in colnames) && (:node_to in colnames)
        edge = edgenumbers_fromnodenumbers(df.node_from, df.node_to, net)
    elseif (:node1 in colnames) && (:node2 in colnames)
        edge = edgenumbers_fromnodenumbers(df.node1, df.node2, net)
    else
        error("could not find the 2 columns with the node numbers")
    end
    return edge
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
    issubsetsplit(a, b)

true if `a ⊆ b` (all taxa in `a` are also in `b`), false otherwise.
`a` and `b` are splits (or clusters) on N taxa, encoded as tuples of N booleans,
with `a[i]` true to mean that taxon i ∈ `a`.
"""
issubsetsplit(a::NTuple{N,Bool}, b::NTuple{N,Bool}) where {N} =
    all(x -> !x[1] | x[2], zip(a,b)) # 1 iteration, no storage of intersection
"""
    isnotdisjointsplit(a, b)

true if `a ∩ b` is non-empty (`a` and `b` share 1 or more taxa), false otherwise.
`a` and `b` are splits (or clusters) on N taxa, encoded as tuples of N booleans.
"""
isnotdisjointsplit(a::NTuple{N,Bool}, b::NTuple{N,Bool}) where {N} =
    any(x -> x[1] & x[2], zip(a,b))

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

true if two partitions A and B of N taxa are blob-compatible; false otherwise.

A and B are compatible if there exists i and j such that all parts
of B except Bⱼ are subsets of Aᵢ. Then, since A and B are partitions of
the same set (not checked!), all parts of A except Aᵢ are subsets of Bⱼ.
This reduces to tree-compatibility when both A and B are bipartitions.

If A or B is given as a single cluster, it is considered as unrooted, as one
part of the bipartition (A,Aᶜ) or (B, Bᶜ) where Aᶜ denotes the complement of A.
"""
function blobcompatible(
    A::AbstractVector{NTuple{N,Bool}},
    B::AbstractVector{NTuple{N,Bool}}
) where N
    Ka = length(A)
    Kb = length(B)
    if Ka == 0 # a partition must have at least 1 part, assuming N > 0
        ArgumentError("$Ka parts in partition A")
    end
    for i in 1:Ka
        j0 = 0 # j0 such that B[j0] ⊇ all but one A[i]
        found_i0 = true # i0 "good" if A[i0] ⊇ B[j] for j ≠ j0
        Ai = A[i]
        for j in 1:Kb
            if !issubsetsplit(B[j], Ai)
                if j0 > 0
                    found_i0 = false
                    break
                else
                    j0 = j
                end
            end
        end
        if found_i0
            return true
        end
    end
    return false # assumes A non-empty
end

function splitblobcompatible(
    split::NTuple{N,Bool},
    B::AbstractVector{NTuple{N,Bool}}
) where N
    splitcomplement = ntuple(i -> !split[i], N)
    A = [split, splitcomplement]
    return blobcompatible(A, B)
end
