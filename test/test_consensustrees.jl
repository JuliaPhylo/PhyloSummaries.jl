@testset "consensus trees" begin

tree1 = readnewick("((A,B),(C,D));")
tree2 = readnewick("((C,D),(B,A));") # ≈ tree1
tree3 = readnewick("((A,C),(B,D));")
tree4 = readnewick("(D,(C,(A,B)));") # ≈ trees 1,2 if unrooted
tree5 = readnewick("(D, C,(A,B)) ;") # ≈ tree4 but root suppressed
taxa = ["A","B","C","D"]

@testset "count_bipartitions!" begin
# unrooted
counts = Dictionary{NTuple{4,Bool},Int}()
PS.count_bipartitions!(counts, tree1, taxa, false)
PS.count_bipartitions!(counts, tree3, taxa, false)
expected = Dict(
    (true, true, false, false) => 1,
    (true, false, true, false) => 1,
)
@test length(counts) == length(expected)
for bp in keys(expected)
    @test counts[bp] == expected[bp]
end
# rooted
empty!(counts)
PS.count_bipartitions!(counts, tree1, taxa, true)
PS.count_bipartitions!(counts, tree3, taxa, true)
expected = Dict(
    (true, true, false, false) => 1,
    (false, false, true, true) => 1,
    (true, false, true, false) => 1,
    (false, true, false, true) => 1,
)
@test length(counts) == length(expected)
for bp in keys(expected)
    @test counts[bp] == expected[bp]
end
end

@testset "utilities" begin
t = PS.startree(["t1","t2","t3"])
ni = Ref(30); ei = Ref(4)
ne = @test_logs (:warn, r"^will skip trivial clade") PS.add_clusteredge!(t,ni,ei,(true,true,true),10,false)
@test isnothing(ne)
ne = @test_logs (:warn, r"^clade already in tree") PS.add_clusteredge!(t,ni,ei,(true,false,false),10,false)
@test isnothing(ne)
ne = PS.add_clusteredge!(t,ni,ei,(true,true,false),10,false)
@test getchild(ne).number == 30
@test t.edge[4].y == 10
## only warning from an incompatible cluster:
# PS.add_clusteredge!(t,ni,ei,(false,true,true),0,false)

@test PS.isredundantsplit((true,true,true,false),
  ((false,false,false,true), (true,true,false,false), (false,false,true,false)))
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
@test writenewick(con, internallabel=false) == "(A,B,C,D);"
con = consensustree(trees; proportion=0.7)
@test writenewick(con,round=true,support=true,internallabel=false) == "(C,D,(B,A)::0.8);"
con = consensustree(trees; rooted=true, supportaslength=true) # greedy
@test writenewick(con,round=true) == "((D,C):0.4,(B,A):0.8);"
@test all(n.fvalue == -1 for n in con.node)
@test [e.y for e in con.edge if !isexternal(e)] == [.4,.8]

tfile = joinpath(@__DIR__,"..","test","raxmltrees.tre")
# tfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","raxmltrees.tre")
trees = readmultinewick(tfile)
@test writenewick(consensustree(trees, proportion=1)) == "(A,B,E,O,(D,C));"
@test writenewick(consensustree(trees),round=true,support=true,
    internallabel=false) == "(E,O,((A,B)::0.833,(C,D)::1.0)::0.533);"
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

tfile = joinpath(@__DIR__,"..","test","tobs20_15taxa.tre")
# tfile = joinpath(dirname(pathof(PhyloSummaries)), "..","test","tobs20_15taxa.tre")
trees = readmultinewick(tfile)
con = (@test_logs consensustree(trees)) # fixit: broken
# @test writenewick(con, support=:y) == fixit

end # of sub-testset


end
