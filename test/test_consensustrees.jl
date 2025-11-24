@testset "consensus trees" begin

tree1 = readnewick("((A,B),(C,D));")
tree2 = readnewick("((C,D),(B,A));") # ≈ tree1
tree3 = readnewick("((A,C),(B,D));")
tree4 = readnewick("(D,(C,(A,B)));") # ≈ trees 1,2 if unrooted
tree5 = readnewick("(D, C,(A,B)) ;") # ≈ tree4 but root suppressed
taxa = ["A","B","C","D"]

@testset "count_bipartitions!" begin
# unrooted
counts = Dictionary{BitVector,Int}()
PhyloSummaries.count_bipartitions!(counts, tree1, taxa, false)
PhyloSummaries.count_bipartitions!(counts, tree3, taxa, false)
expected = Dict(
    BitVector([1, 1, 0, 0]) => 1,
    BitVector([1, 0, 1, 0]) => 1,
)
@test length(counts) == length(expected)
for bp in keys(expected)
    @test counts[bp] == expected[bp]
end
# rooted
empty!(counts)
PhyloSummaries.count_bipartitions!(counts, tree1, taxa, true)
PhyloSummaries.count_bipartitions!(counts, tree3, taxa, true)
expected = Dict(
    BitVector([1, 1, 0, 0]) => 1,
    BitVector([0, 0, 1, 1]) => 1,
    BitVector([1, 0, 1, 0]) => 1,
    BitVector([0, 1, 0, 1]) => 1,
)
@test length(counts) == length(expected)
for bp in keys(expected)
    @test counts[bp] == expected[bp]
end
end

@testset "consensustree" begin

@test_throws ArgumentError consensustree(PN.HybridNetwork[])
@test_throws ArgumentError consensustree([readnewick("((A,#H1),(B)#H1);")])
@test_throws "not share the same taxon set" consensustree([tree1, readnewick("((A,B),C);")])
@test_throws "not in taxon list" consensustree([tree1, readnewick("((A,B),C,E);")])

# single-tree input
@test writenewick(consensustree([tree1])) == "(C,D,(A,B));"
@test writenewick(consensustree([tree1]; rooted=true)) == "((A,B),(C,D));"

# 4 taxa, 5 trees, missing edge lengths
trees = [tree2,tree3,tree1,tree4,tree5]
con = consensustree(trees; rooted=true, proportion=0.8)
writenewick(con) == "(A,B,C,D);"
con = consensustree(trees; proportion=0.7)
writenewick(con,round=true,support=true) == "(C,D,(B,A)::0.8);"
con = consensustree(trees; rooted=true, supportaslength=true) # greedy
@test writenewick(con,round=true) == "((D,C):0.4,(B,A):0.8);"
@test [n.fvalue for n in con.node if !n.leaf] == [-1,.4,.8]
@test [e.y for e in con.edge if !isexternal(e)] == [.4,.8]

tfile = joinpath(@__DIR__,"..","test","raxmltrees.tre")
# tfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","raxmltrees.tre")
trees = readmultinewick(tfile)
@test writenewick(consensustree(trees, proportion=1)) == "(A,B,E,O,(D,C));"
@test writenewick(consensustree(trees),round=true,support=true) == "(E,O,((A,B)::0.833,(C,D)::1.0)::0.533);"
con = consensustree(trees; rooted=true, supportaslength=true)
@test writenewick(con, round=true) == "((O,E):0.033,((A,B):0.767,(C,D):1.0):0.5);"
#=
checked correctness with ape::consensus
```r
library(ape)
mytrees <- read.tree("test/raxmltrees.tre")
write.tree(consensus(mytrees, p=0.5)) #  "((E,O)0.533,(A,B)0.833,(C,D)1)1;"
write.tree(consensus(mytrees, p=0.5, rooted=T)) # "((A,B)0.7667,(C,D)1,E,O)1;"
also by plotting them: tree #2 has (O,E). 15 trees have A-D: #1,4-7,11,16-17,20-21,25-28,30
```
=#

end # of sub-testset


end
