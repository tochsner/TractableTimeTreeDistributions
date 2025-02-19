@testset "cladify tree with one taxa" begin
    newick = "(A:1);"
    tree = load_trees_from_newick(newick)

    cladified_tree = cladify_tree(tree[1])

    @test cladified_tree.clades == Set([Leaf(1, 1, "A")])
    @test cladified_tree.splits == Set()
end

@testset "cladify tree with two taxa" begin
    newick = "(A:1,B:1);"
    tree = load_trees_from_newick(newick)

    cladified_tree = cladify_tree(tree[1])

    @test length(cladified_tree.clades) == 3
    @test Leaf(1, 2, "A") in cladified_tree.clades
    @test Leaf(2, 2, "B") in cladified_tree.clades
    @test Clade([1, 2], 2) in cladified_tree.clades
    @test Clade([1, 2], 2) in cladified_tree.clades

    @test length(cladified_tree.splits) == 1
    @test CladeSplit(Clade(1, 2), Clade(2, 2)) in cladified_tree.splits
end

@testset "cladify caterpillar tree with four taxa" begin
    newick = "(A:1,(B:1,(C:1,D:1):1));"
    tree = load_trees_from_newick(newick)

    cladified_tree = cladify_tree(tree[1])

    @test length(cladified_tree.clades) == 7

    @test Leaf(1, 4, "A") in cladified_tree.clades
    @test Leaf(2, 4, "B") in cladified_tree.clades
    @test Leaf(3, 4, "C") in cladified_tree.clades
    @test Leaf(4, 4, "D") in cladified_tree.clades
    
    @test Clade([3, 4], 4) in cladified_tree.clades
    @test Clade([2, 3, 4], 4) in cladified_tree.clades
    @test Clade([1, 2, 3, 4], 4) in cladified_tree.clades

    @test length(cladified_tree.splits) == 3
    @test CladeSplit(Clade(3, 4), Clade(4, 4)) in cladified_tree.splits
    @test CladeSplit(Clade([3, 4], 4), Clade(2, 4)) in cladified_tree.splits
    @test CladeSplit(Clade([2, 3, 4], 4), Clade(1, 4)) in cladified_tree.splits
end

@testset "cladify balanced tree with four taxa" begin
    newick = "((A:1,B:1):1,(C:1,D:1):1):1;"
    tree = load_trees_from_newick(newick)

    cladified_tree = cladify_tree(tree[1])

    @test length(cladified_tree.clades) == 7

    @test Leaf(1, 4, "A") in cladified_tree.clades
    @test Leaf(2, 4, "B") in cladified_tree.clades
    @test Leaf(3, 4, "C") in cladified_tree.clades
    @test Leaf(4, 4, "D") in cladified_tree.clades
    
    @test Clade([1, 2], 4) in cladified_tree.clades
    @test Clade([3, 4], 4) in cladified_tree.clades
    @test Clade([1, 2, 3, 4], 4) in cladified_tree.clades

    @test length(cladified_tree.splits) == 3
    @test CladeSplit(Clade(1, 4), Clade(2, 4)) in cladified_tree.splits
    @test CladeSplit(Clade(3, 4), Clade(4, 4)) in cladified_tree.splits
    @test CladeSplit(Clade([1, 2], 4), Clade([3, 4], 4)) in cladified_tree.splits
end