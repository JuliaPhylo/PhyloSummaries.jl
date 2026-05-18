"""
    blobpartitions_support(networks, referencenet;
        minimumblobdegree=3, netweights=nothing)

Support for blob partitions and related features (circular orders,
hybrid clades, bipartitions non-redundant with a blob) for the blobs present
in a reference network `referencenet`,
based on their frequency in a sample of `networks`.

Sample networks are weighted equally by default, unless a vector of
`netweights` is provided (of same length as the vector of sample networks).

Output: a `NamedTuple` of tables, named as follows
- `:blob_table`: support for each blob partition in the reference network,
  including support for the circular order of the taxon blocks if the sample
  networks are of level 1.
- `:transfer_table`: transfer index for each blob partition β in the reference
  network: average TI(β,N) over networks N in the sample. This is the average
  number of taxa to remove (from β and from N) to make β match some taxon block
  in N (either from a blob or from a non-redundant bipartition).
- `:hybrid_table`: support for each taxon block to be below a blob's
  lowest hybrid node in a sample network.
- `:bipartition_table`: support for each non-redundant bipartition (cut-edge)
  in the reference network.
- `:taxa`: sorted taxon labels shared by all networks. `taxa[i]` is indicated
  as taxon `i` in the other tables.

See also:
[`PhyloNetworks.treeedges_support`](@extref),
[`PhyloNetworks.hybridclades_support`](@extref),
[`consensus_treeofblobs`](@ref),
[`PhyloSummaries.count_blobpartitions`](@ref).

# example

```jldoctest
julia> netfile = joinpath(dirname(pathof(PhyloSummaries)), "..", "test","level1_7taxa_abc.nwk");

julia> bootnet = readmultinewick(netfile); # could be bootstrap networks

julia> nwk = "(((a3,(a4,#H1)),a2),(((c2,(#H2,c1)),(b1)#H2))#H1,a1);";

julia> refnet = readnewick(nwk); # same as 5th bootnet

julia> res = blobpartitions_support(bootnet, refnet);

julia> keys(res)
(:blob_table, :transfer_table, :circorder_table, :hybrid_table, :bipartition_table, :taxa)

julia> using DataFrames; DataFrame(res[:blob_table])
2×6 DataFrame
 Row │ blob   degree  node   support_partition  partition_num  partition            
     │ Int64  Int64   Int64  Float64            String         String               
─────┼──────────────────────────────────────────────────────────────────────────────
   1 │     2       4     -7                0.4  1,2,3,4|7|6|5  a1,a2,a3,a4|c2|c1|b1
   2 │     1       5     -2                0.6  1|2|3|4|5,6,7  a1|a2|a3|a4|b1,c1,c2

julia> DataFrame(res[:transfer_table])
2×5 DataFrame
 Row │ blob   node   transferindex  partition_num  partition            
     │ Int64  Int64  Float64        String         String               
─────┼──────────────────────────────────────────────────────────────────
   1 │     2     -7            1.0  1,2,3,4|7|6|5  a1,a2,a3,a4|c2|c1|b1
   2 │     1     -2            1.2  1|2|3|4|5,6,7  a1|a2|a3|a4|b1,c1,c2

julia> DataFrame(res[:circorder_table])
2×5 DataFrame
 Row │ blob   order            support_circorder  partition_num  partition            
     │ Int64  Tuple…           Float64            String         String               
─────┼────────────────────────────────────────────────────────────────────────────────
   1 │     2  (1, 2, 3, 4)                   0.4  1,2,3,4|7|6|5  a1,a2,a3,a4|c2|c1|b1
   2 │     1  (1, 2, 3, 4, 5)                0.6  1|2|3|4|5,6,7  a1|a2|a3|a4|b1,c1,c2

julia> DataFrame(res[:hybrid_table])
3×7 DataFrame
 Row │ blob   node_from  node_to  edge   support_hybrid  cluster_num  cluster     
     │ Int64  Int64      Int64    Int64  Float64         String       String      
─────┼────────────────────────────────────────────────────────────────────────────
   1 │     2          6        8     13             0.2  5            b1
   2 │     2         -7        3     15             0.4  1,2,3,4      a1,a2,a3,a4
   3 │     1          3       -7     15             0.6  5,6,7        b1,c1,c2

julia> DataFrame(res[:bipartition_table]) # refnet has 0 non-redundant bipartitions
0×6 DataFrame
 Row │ node1  node2  edge   support_nonredundant  cluster_num  cluster 
     │ Int64  Int64  Int64  Float64               String       String  
─────┴─────────────────────────────────────────────────────────────────

```
To plot these support values onto the reference network,
see examples in the package manual.
"""
function blobpartitions_support(
    networks::AbstractVector{PN.HybridNetwork},
    referencenet::PN.HybridNetwork;
    minimumblobdegree::Int=3,
    netweights::Union{Nothing,AbstractVector}=nothing,
)
    isempty(networks) &&
        throw(ArgumentError("No input networks: cannot compute reference support"))
    minimumblobdegree ≥ 3 ||
        throw(ArgumentError("minimumblobdegree should be 3 or higher, not $minimumblobdegree"))
    taxa = sort(tiplabels(referencenet))
    ntaxa = length(taxa)
    nnets = length(networks)
    if !isnothing(netweights)
      length(netweights) == nnets ||
          error("there should be $nnets network weights, got $(length(netweights))")
      nnets = sum(netweights)
      nnets>0 && all(netweights .>= 0) || error("network weights should be ≥ 0")
    end
    # collect info on sample networks. do *not* filter them
    bbidx = [Int[] for _ in eachindex(networks)]
    bpidx = [Int[] for _ in eachindex(networks)]
    blobvec, bpvec = count_blobpartitions(networks, taxa, minimumblobdegree,
        false, netweights, bbidx, bpidx)
    hybdict = count_hybridclusters(blobvec)
    # collect info on reference blobs
    refblobs = BlobFreq{ntaxa}[]
    refbps = SplitFreq{ntaxa}[]
    hwmatrix, edgemap, blobdegree = count_blobpartitions!(refblobs, refbps, referencenet,
        taxa, minimumblobdegree, false, 0.0) # frequencies initialized at 0
    # compare reference and sample blobs
    transfer = transferindex(refblobs, blobvec, bpvec, netweights, nnets, bbidx, bpidx)
    update_blobcircorderfrequency!(refblobs, blobvec)
    update_hybridclusterfrequency!(refblobs, hybdict)
    update_bipartitionfrequency!(refbps, bpvec)
    bbn, bbe_nums = _blobnode_blobedges(referencenet,
        hwmatrix, edgemap, refblobs, taxa, blobdegree, minimumblobdegree)
    edgenum2idx = Dict(e.number => i for (i, e) in pairs(referencenet.edge))
    # specify type in case of empty iterator
    bbei = Vector{Int}[[edgenum2idx[n] for n in nums] for nums in bbe_nums]
    bpei, bpe_nums = _bipartition_edgeindices(refbps, hwmatrix, edgenum2idx)
    # build output tables
    bdat, odat = blobdata_onToB(refblobs, bbn, nnets, taxa)
    tidat = blobdata_onToB_transferindex(refblobs, bbn, transfer, taxa)
    hdat = hybriddata_onToB(refblobs, bbn, bbei, nnets, referencenet.edge, taxa, bbe_nums)
    sdat = bipartdata_onToB(refbps, bpei, nnets, referencenet.edge, taxa, bpe_nums)
    return (blob_table=bdat, transfer_table=tidat, circorder_table=odat,
        hybrid_table=hdat, bipartition_table=sdat, taxa=taxa)
end

"""
    _bipartition_edgeindices(refbps, hwmatrix, edgenum2idx)

For each non-redundant bipartition in `refbps`, find the corresponding edge
index and edge number by matching the split against the `hwmatrix`.
"""
function _bipartition_edgeindices(
    refbps::Vector{SplitFreq{N}},
    hwmatrix::AbstractMatrix,
    edgenum2idx::Dict{Int,Int},
) where N
    bpei = Int[]
    bpe_nums = Int[]
    for rb in refbps
        matched = false
        csplit = .!rb.split
        for row in axes(hwmatrix, 1)
            rowsplit = cluster_fromHmatrix(hwmatrix, row, N)
            if rowsplit == rb.split || rowsplit == csplit
                edge_num = hwmatrix[row, 1]
                push!(bpe_nums, edge_num)
                push!(bpei, edgenum2idx[edge_num])
                matched = true
                break
            end
        end
        matched || error("bipartition $(rb.split) not found in hwmatrix")
    end
    return bpei, bpe_nums
end

"""
    update_blobcircorderfrequency!(refblobs, sampleblobs)

Update the frequency and circular order dictionary of each "reference" blob
partition blob (item in `refblobs`), to match that from `sampleblobs`.
A reference blob partition not found in `sampleblobs` is not modified --
which assumes that its frequencies had been initialized to 0.

Note the blobs hybrid dictionaries are *not* updated. For this, see instead
[`update_hybridclusterfrequency!`](@ref) to use the frequencies of hybrid
clades aggregated over all sample blobs
(a hybrid clade can originate below different blobs, of different partitions).
"""
function update_blobcircorderfrequency!(
    refblobs::Vector{BlobFreq{N}},
    sampleblobs::Vector{BlobFreq{N}},
) where N
    for rb in refblobs
        matchidx, idxmap = findmatchingblob(sampleblobs, rb.partition)
        isnothing(matchidx) && continue # all frequencies were initialized to 0
        bf = sampleblobs[matchidx]
        freq!(rb, freq(bf))
        # copy circular orders from bf into rb, remapping block indices.
        # rb.partition[k] = bf.partition[idxmap[k]], so the inverse gives us:
        # bf block index j → rb block index invmap[j]
        P = nblocks(rb) # also = length(idxmap)
        invmap = Vector{Int}(undef, P)
        for k in 1:P
            invmap[idxmap[k]] = k
        end
        for (co_key, co_freq) in bf.circorder
            remapped = [invmap[j] for j in co_key]
            add_canonical_circularorder!(rb.circorder, remapped, co_freq)
        end
    end
    return nothing
end

"""
    update_bipartitionfrequency!(refbps, samplebps)

Update the frequency of each bipartition in `refbps` to match that from
`samplebps`. A bipartition not found in `samplebps` is not modified --
which assumes that its frequencies had been initialized to 0.

Bipartitions are considered unrooted: so a bipartition in `samplebps` is a
match if either its clade or its complement matches that from `refbps`.
"""
function update_bipartitionfrequency!(
    refbps::Vector{SplitFreq{N}},
    samplebps::Vector{SplitFreq{N}},
) where N
    for rb in refbps
        i = findfirst(bp -> bp.split == rb.split || bp.split == .!rb.split, samplebps)
        isnothing(i) && continue
        freq!(rb, freq(samplebps[i]))
    end
    return nothing
end

"""
    _blobnode_blobedges

Output: `(blobnode, blobedge_numbers)` where
- `blobnode` is the vector of entry nodes, one per blob and
- `blobedge_numbers` has one vector per blob, each being a vector with
  one edge number per taxon block.

**Warning**: uses `net.partition` and internal fields `.inte1` and `.intn1`.
These fields should store the edge's bicomponent number, and the node's blob
number (if non-trivial, 0 if the node is a trivial blob).
They do if `net` was already traversed by [`count_blobpartitions`](@ref),
which calls `process_biconnectedcomponents!` to build `net.partition`.
The blob number should also be its index in the input `refblobs`.

assumptions about input arguments, met after calling `count_blobpartitions`
with the same `minBdegree`:
- `blobdegree`: vector of length `net.partition`, with 1 degree per bicomponent:
  blob degree if the bicomponent is at the top (entry) of a blob, 0 if it is
  trivial (1 edge) or below another bicomponent in the same blob.
- `blobnode` has 1 node for each blob of degree ≥ `minBdegree`, in the same
  order as they appear in `blobdegree` (as in a `net.partition`)
"""
function _blobnode_blobedges(
    net::PN.HybridNetwork,
    hwmatrix::AbstractMatrix,
    edgemap::Dict{<:Integer,<:Integer},
    refblobs::Vector{BlobFreq{N}},
    taxa::AbstractVector{<:String},
    blobdegree::Vector{Int},
    minBdegree::Int
) where N
    nblobs = length(refblobs)
    blobedge_nums = [zeros(Int, nblocks(bb)) for bb in refblobs]
    # 1. build blobnode, and map entry node intn1 → blob index in refblobs
    blobnode = PN.Node[]
    intn1_to_blobidx = Dict{Int,Int}()
    idx = 1
    for (bidx, bbd) in pairs(blobdegree)
        bbd < minBdegree && continue
        nn = net.vec_node[PN.entrynode_preindex(net.partition[bidx])]
        nn.intn1 == bidx ||
            error("blob starting at bicomponent $bidx ≠ entry node intn1 $(nn.intn1)")
        push!(blobnode, nn)
        intn1_to_blobidx[nn.intn1] = idx
        idx += 1
    end
    length(blobnode) == nblobs ||
        error("found $(length(blobnode)) interesting blobs instead of $nblobs")

    # 2. Assign cut-edges to their interesting blobs
    for p in net.partition
        if PN.istrivial(p) # add cut-edge to blobedge_nums?
            ee = p.edges[1]
            b1_i = ee.node[1].intn1
            b2_i = ee.node[2].intn1
            (b1_i == 0 && b2_i == 0) && continue # non-redundant bipart, could be external
            enum = ee.number
            row = get(edgemap, enum, nothing)
            if isnothing(row)
                cn = getchild(ee)
                cn.leaf || error("cut-edge without a row in hwmatrix, yet child isn't a leaf")
                i0 = findfirst(isequal(cn.name), taxa)
                split = ntuple(isequal(i0), N)
            else
                split = cluster_fromHmatrix(hwmatrix, row, N)
            end
            vsplitmatch(v) = v == split || all(v .!== split)
            for b_i in (b1_i, b2_i)
                b_i == 0 && continue # adjacent to trivial blob
                bi = get(intn1_to_blobidx, b_i, nothing)
                if isnothing(bi)
                    blobdegree[b_i] < minBdegree ||
                        error("interesting blob at node with intn1=$b_i yet index not found")
                    continue
                end
                block_j = findfirst(vsplitmatch, refblobs[bi].partition)
                isnothing(block_j) && error("cut-edge $enum, split $split: no matching taxon block in $(refblobs[bi])")
                blobedge_nums[bi][block_j] = enum
            end
        end
    end
    any(any(enums .== 0) for enums in blobedge_nums) &&
        error("some blob taxon blocks have no matching edge")
    return blobnode, blobedge_nums
end
