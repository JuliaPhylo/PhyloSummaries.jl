
@testset "extract_bipartitions! canonicalises unrooted splits" begin
    first_tree = readnewick("((A,B),(C,D));")
    directedges!(first_tree)
    taxa = sort(tiplabels(first_tree))

    counts = Dict{Tuple{Vararg{Int}}, Int}()
    PhyloSummaries.extract_bipartitions!(counts, first_tree, taxa, false)

    second_tree = readnewick("((A,C),(B,D));")
    directedges!(second_tree)
    PhyloSummaries.extract_bipartitions!(counts, second_tree, taxa, false)

    expected = Dict(
        (1, 1, 0, 0) => 2,  
        (1, 0, 1, 0) => 2   
    )

    @test counts == expected
end


@testset "extract_bipartitions! counts expected splits" begin

    first_tree = readnewick("((A,B),(C,D));")
    directedges!(first_tree)
    taxa = sort(tiplabels(first_tree))

    counts = Dict{Tuple{Vararg{Int}}, Int}()
    PhyloSummaries.extract_bipartitions!(counts, first_tree, taxa, true)

    second_tree = readnewick("((A,C),(B,D));")
    directedges!(second_tree)
    PhyloSummaries.extract_bipartitions!(counts, second_tree, taxa, true)

    expected = Dict(
        (1, 1, 0, 0) => 1,  
        (0, 0, 1, 1) => 1, 
        (1, 0, 1, 0) => 1,  
        (0, 1, 0, 1) => 1  
    )

    @test counts == expected
end
