using Test
using PhyloNetworks
using PhyloSummaries
using Base: BitSet

const PN = PhyloNetworks

function direct_tree(newick::AbstractString)
    net = PN.readnewick(newick, false)
    PN.directedges!(net)
    return net
end

@testset "extract_edge_bipartition " begin
    tree = direct_tree("((A,B),(C,D));")
    taxa = sort(PN.tiplabels(tree))

    splits = Tuple{Vararg{Int}}[]
    for edge in tree.edge
        split = PhyloSummaries.extract_edge_bipartition(edge, taxa)
        split === nothing && continue
        push!(splits, split)
    end

    @test !isempty(splits)
    @test length(unique(splits)) == 1
    @test unique(splits)[1] == (1, 1, 0, 0)

end

@testset "consensus_bipartition " begin
    counts = Dict{Tuple{Vararg{Int}},Int}()
    subsets = Dict{Tuple{Vararg{Int}},Vector{Int}}()
    subset_sets = Dict{Tuple{Vararg{Int}},BitSet}()

    tree1 = direct_tree("((A,B),(C,D));")
    taxa = sort(PN.tiplabels(tree1))
    PhyloSummaries.extract_bipartitions!(counts, subsets, subset_sets, tree1, taxa)

    tree2 = direct_tree("((A,C),(B,D));")
    PhyloSummaries.extract_bipartitions!(counts, subsets, subset_sets, tree2, taxa)

    selected = PhyloSummaries.consensus_bipartition(counts, subset_sets)

    @test selected == [(1, 1, 0, 0)]
    @test counts[(1, 1, 0, 0)] == 2

end
