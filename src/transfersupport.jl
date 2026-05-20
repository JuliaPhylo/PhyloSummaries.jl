"""
    transferdistance!(C::AbstractMatrix{Int},
        β1::AbstractVector{NTuple{N,Bool}},
        β2::AbstractVector{NTuple{N,Bool}})

Twice the transfer distance `2 * d_transfer(β1,β2)` between partitions `β1`
and `β2` of the same set of `N` taxa.
The transfer distance `d_transfer(β1,β2)` is the minimum number of taxa to
transfer from one block to another for the two partitions to be equal
(one partition possibly padded with extra empty taxon blocks).
Alternatively, this is the minimum number of taxa to exclude, for the two
partitions to be equal (possibly resulting in some empty taxon blocks).

The output here is twice the transfer distance: counting each transfer
with a cost of 2 for removing a taxon from one block, then adding it back
to another block.

`C[i,j]` is used (and mutated) to store the cost of associating
taxon block `i` from `β1` to taxon block `j` from `β2`.
This cost is |β1[i] Δ β2[j]| where Δ is the symmetric difference:
the number of taxa to be removed (or transferred to another block).

The Hungarian algorithm is then used with this cost matrix, from
[Hungarian.jl](https://github.com/Gnimuc/Hungarian.jl)

Assumptions:
- `N` ≥ 1, and β1, β2 are both correct partition of the `N` taxa into
  `k1` and `k2` taxon blocks respectively, hence with `k1 ≥ 1` and k2 ≥ 1`.
- `C` has size `m1 × m2` with `m1,m2 ≥ k` where `k = max{k1,k2}`.
  The first `k` rows and columns of `C` are used to store the costs.

# example

```jldoctest
julia> const PS = PhyloSummaries; # to use internals with less typing

julia> # below: partitions β1 = 123|5|4 and β2 = 12|345

julia> β1 = [Bool.(v) for v in [(1,1,1,0,0), (0,0,0,0,1), (0,0,0,1,0)]]
3-element Vector{NTuple{5, Bool}}:
 (1, 1, 1, 0, 0)
 (0, 0, 0, 0, 1)
 (0, 0, 0, 1, 0)

julia> β2 = [(true,true,false,false,false), (false,false,true,true,true)]
2-element Vector{NTuple{5, Bool}}:
 (1, 1, 0, 0, 0)
 (0, 0, 1, 1, 1)

julia> C = Matrix{Int}(undef, 3, 3); # max 3 taxon blocks

julia> PS.transferdistance!(C, β1, β2) # transfer 3 and 5 (or 3 and 4)
4

julia> PS.transferdistance!(C, β2[1], β1) # same, from 1 side of bipartition β2
4
```
"""
function transferdistance!(
    C::AbstractMatrix{Int},
    blob1::AbstractVector{NTuple{N,Bool}},
    blob2::AbstractVector{NTuple{N,Bool}},
) where N
    k1 = length(blob1)
    k2 = length(blob2)
    # (k1 == 0 || k2 == 0) && return 0
    if k2 < k1
        blob1,blob2 = blob2,blob1
        k1,k2 = k2,k1
    end # now k1 <= k2
    Cv = view(C, 1:k2, 1:k2)
    fill!(Cv, 0)
    for j in 1:k2
        for i in 1:k1
            Cv[i,j] = count(t -> xor(blob1[i][t], blob2[j][t]), 1:N)
        end
        Cv[(k1+1):k2,j] .= count(blob2[j])
    end
    _, mincost = hungarian(Cv)
    return mincost
end
"""
    transferdistance!(C::AbstractMatrix{Int},
        s1::NTuple{N,Bool},
        β2::AbstractVector{NTuple{N,Bool}})

Twice the transfer distance `2 * d_transfer(β1,β2)` between partitions `β1`
and `β2` of the same set of `N` taxa, where β1 is the bipartition `s1|s̄1`.
"""
function transferdistance!(
    C::AbstractMatrix{Int},
    split::NTuple{N,Bool},
    blob2::AbstractVector{NTuple{N,Bool}},
) where N
    k2 = length(blob2)
    @assert k2 > 1
    Cv = view(C, 1:k2, 1:k2)
    fill!(Cv, 0)
    for j in 1:k2
        Cv[1,j] = count(t -> xor(split[t], blob2[j][t]), 1:N)
        Cv[2,j] = N - Cv[1,j] # split complement ↔ blob2[j]
        Cv[3:k2,j] .= count(blob2[j])
    end
    _, mincost = hungarian(Cv)
    return mincost
end

"""
    transferindex!(C, dist_cache, β_idx, β,
                   blb_idx, blobs, bip_idx, biparts)

Twice the minimum transfer distance between a reference blob β,
at row index `β_idx`, and all partitions β* in a sample network:
blobs at indices `blb_idx` in the `blobs` vector, and
non-redundant bipartitions at indices `bip_idx` in the `biparts` vector.
This is `2 * TI` with
`TI = min{ d_transfer(β,β*); β* ∈ blobs[blb_idx] ∪ biparts[bip_idx]}`.

Distances are lazily computed via [`transferdistance!`](@ref) and
stored in a dictionary `dist_cache[(β_idx, j, blob/bipart)]`
for reuse across networks.

For each calculation of the transfer distance between β and one of the
network's partition β*, `C` is used to store the costs to map the taxon blocks
in β to those in β*.
"""
function transferindex!(
    C::AbstractMatrix{Int},
    dist_cache::Dict{Tuple{Int,Int,Bool},Int},
    ref_idx::Int,
    refblob::AbstractVector{NTuple{N,Bool}},
    blb_idx::AbstractVector{Int},
    blobs::Vector{BlobFreq{N}},
    bip_idx::AbstractVector{Int},
    biparts::Vector{SplitFreq{N}},
) where N
    mindist = N # max: transfer all N taxa. No need for typemax(Int)
    for (idx, isblob) in zip((blb_idx,bip_idx), (true,false)), j in idx
        key = (ref_idx, j, isblob)
        βnet = (isblob ? blobs[j].partition : biparts[j].split)
        d = get!(dist_cache, key) do # enters block only if key not in dict yet
            transferdistance!(C, βnet, refblob)
        end
        if d < mindist  mindist = d; end
        mindist == 0 && break
    end
    return mindist
end

"""
    transferindex(refblobs, blobs, biparts, netweights, nnets, bbidx, bpidx)

Vector listing the transfer index of each blob in `refblobs` with respect to
the sample of networks, whose blobs are stored in `blobs` and whose
non-redundant bipartitions are stored in `biparts`. Network `i` has weight
`netweights[i]`, blobs at indices `bbidx` in `blobs`, and bipartitions at
indices `bpidx` in `biparts`.

The transfer index of a blob β is the average minimum transfer distance d(β,N),
with the average (or weighted average) taken over all networks N in the sample.
`nnets` is a normalization factor: the total weight of all networks.

See [`transferindex!`](@ref), which calculates twice the transfer index of
one reference partition with respect to one network.
"""
function transferindex(
    refblobs::Vector{BlobFreq{N}},
    blobs::Vector{BlobFreq{N}},
    biparts::Vector{SplitFreq{N}},
    netweights::Union{Nothing,AbstractVector},
    nnets::Number,
    bbidx::Vector{Vector{Int}},
    bpidx::Vector{Vector{Int}},
) where N
    transfercost = Matrix{Int}(undef, N, N)
    transferdist = Dict{Tuple{Int,Int,Bool},Int}() # d(ref β i, β* j, blob/bipart)
    ti = zeros(length(refblobs)) # sum of TI(refblob i, N) over N
    usenw = !isnothing(netweights)
    for n_i in eachindex(bbidx) # for each network
        for (rb_i,rb) in enumerate(refblobs)
            d = transferindex!(transfercost, transferdist, rb_i, rb.partition,
                bbidx[n_i], blobs, bpidx[n_i], biparts)
            ti[rb_i] += (usenw ? netweights[n_i] * d : d)
        end
    end
    ti ./= (2 * nnets) # 2: bc transferindex! gets twice the transfer index
    return ti
end

"""
Table listing each blob β in the reference network and its transfer support
across a sample of `networks`. This support is the blob's average transfer index
with respect to a network, averaged across all networks in the sample.
For each network N, the transfer index TI(β,N) of a blob β with respect to N
is the transfer distance between β (as a partition into taxon blocks) and the
best-matching blob in N:
`TI(β,N) = min{ d_transfer(β,β*), β* blob partition in N}`

*Warning*: this is **experimental**.
The partitions being considered here are restricted to
(1) partitions into 3+ taxon blocks from non-trivial 'interesting' blobs, and
(2) non-redundant bipartitions.
Not considered are partitions into 3+ taxon blocks from trivial blobs
(which are not 'interesting').

When taking the average, sample networks are weighted equally by default,
unless a vector of `netweights` is provided.

Note that high index means low support, because the transfer index is the average
number of taxa to remove to match the given partition β with a partition in
a sample network. β is found in N if there are 0 taxa to remove (or transfer),
so a transfer index of 0 means perfect support.
"""
function blobdata_onToB_transferindex(
    blobparts::Vector{BlobFreq{N}},
    blobnode::Vector{PN.Node},
    ti::AbstractVector,
    taxa::Vector{<:AbstractString},
) where N
    nB = length(blobparts)
    @assert nB == length(blobnode)
    bitr = ((i,blobparts[i]) for i in nB:-1:1) # from most to least frequent blob
    blob_data = (blob = [i for (i,b) in bitr],
        node = [blobnode[i].number for (i,b) in bitr],
        transferindex = [ti[i] for (i,b) in bitr],
        partition_num = [partitionstring(b) for (i,b) in bitr],
        partition = [partitionstring_names(b.partition, taxa) for (i,b) in bitr])
    return blob_data
end
