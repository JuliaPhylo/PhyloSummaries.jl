using Base: BitSet


function consensus_tree(networks::AbstractVector{PN.HybridNetwork};)
    isempty(networks) && throw(ArgumentError("consensus_tree requires at least one network"))
    if length(networks) == 1
        return deepcopy(networks[1])
    end

    

    taxa = sort(PN.tiplabels(networks[1]))
    validate_taxa_match!(networks, taxa)

    counts = Dict{Tuple{Vararg{Int}},Int}()

    for net in networks
        root_idx = if net.isrooted
            root_name = net.node[net.rooti].name
            idx = root_name === "" ? nothing : findfirst(==(root_name), taxa)
            something(idx, length(taxa))
        else
            length(taxa)
        end
        extract_bipartitions!(counts, net, taxa, root_idx)
    end

    return consensus_bipartition() #TODO
end

function validate_taxa_match!(networks::AbstractVector{PN.HybridNetwork}, reference_taxa::Vector{String})
    refsorted = reference_taxa
    for (idx, net) in pairs(networks)
        taxa = sort(PN.tiplabels(net))
        taxa == refsorted || throw(ArgumentError("network $idx does not share the same taxa set"))
    end
    return nothing
end

function extract_bipartitions!(
    counts::Dict{Tuple{Vararg{Int}},Int},
    net::PN.HybridNetwork,
    taxa::Vector{String},
    root_idx::Int
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
        edge_kind == 10 || continue  
        split = bipartition_from_cluster(view(hw_matrix, row_idx, taxa_cols), root_idx)
        split === nothing && continue

        counts[split] = get(counts, split, 0) + 1

    end

    return nothing
end

function bipartition_from_cluster(cluster::AbstractVector, root_idx::Int)
    partition = copy(cluster)
    if partition[root_idx] == 1
        partition .= .!partition
    end
    
    ones = count(==(1), partition)
    zeros = length(partition) - ones
    (ones <= 1 || zeros <= 1) && return nothing
    return Tuple(partition)
end


function consensus_bipartition()
    return nothing
end


#     counts::Dict{Tuple{Vararg{Int}},Int},
#     subset_sets::Dict{Tuple{Vararg{Int}},BitSet}
# )
    
#     sorted_splits = sort!(
#         collect(keys(counts));
#         by = split -> (-counts[split], length(split), split)
#     )

#     selected = Tuple{Vararg{Int}}[]
#     selected_sets = BitSet[]

#     for split in sorted_splits
#         bs = subset_sets[split]
#         compatible = true
#         for existing in selected_sets
#             compatible = issubset(bs, existing) ||
#                           issubset(existing, bs) ||
#                           isempty(intersect(bs, existing))
#             compatible || break
#         end
#         if compatible
#             push!(selected, split)
#             push!(selected_sets, bs)
#         end
#     end

#     return selected
# end

# change to use the matrix
# change to tuple of boolean and then use that as a key.
# once we hace n-3 stop