@testset "parameterize tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):1);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    trees_with_short_branches = transform_short_branches.(cladified_trees)
    @test trees_with_short_branches[1].parameters[Clade(1:4, 4)] ≈ 1.0
    @test trees_with_short_branches[1].parameters[Clade(2:4, 4)] ≈ 1.0
    @test trees_with_short_branches[1].parameters[Clade(3:4, 4)] ≈ 1.0
end

@testset "parameterize another tree with four taxa" begin
    newick = "((A:1,B:1):2,(C:2,D:2):1):0;" 
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    trees_with_short_branches = transform_short_branches.(cladified_trees)
    @test trees_with_short_branches[1].parameters[Clade(1:4, 4)] ≈ 1
    @test trees_with_short_branches[1].parameters[Clade(1:2, 4)] ≈ 1
    @test trees_with_short_branches[1].parameters[Clade(3:4, 4)] ≈ 2
end

@testset "parameterize a tree with six taxa" begin
    newick = "((((A:1,B:1):2,(C:2,D:2):1):4,(E:5, F:5):2):0);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    trees_with_short_branches = transform_short_branches.(cladified_trees)

    @test trees_with_short_branches[1].parameters[Clade(1:6, 6)] ≈ 2
    @test trees_with_short_branches[1].parameters[Clade(1:4, 6)] ≈ 1
    @test trees_with_short_branches[1].parameters[Clade(1:2, 6)] ≈ 1
    @test trees_with_short_branches[1].parameters[Clade(3:4, 6)] ≈ 2
    @test trees_with_short_branches[1].parameters[Clade(5:6, 6)] ≈ 5
end

@testset "set heights for tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):1);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    trees_with_short_branches = transform_short_branches.(cladified_trees)
    recovered_trees = invert_short_branches.(trees_with_short_branches)

    @test recovered_trees[1].parameters[recovered_trees[1].root] ≈ 3
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:4, 4)].clade2] ≈ 2
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(2:4, 4)].clade2] ≈ 1
end

@testset "set heights for another tree with four taxa" begin
    newick = "((A:1,B:1):2,(C:2,D:2):1):0;" 
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    trees_with_short_branches = transform_short_branches.(cladified_trees)
    recovered_trees = invert_short_branches.(trees_with_short_branches)

    @test recovered_trees[1].parameters[recovered_trees[1].root] ≈ 3
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:4, 4)].clade1] ≈ 1
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:4, 4)].clade2] ≈ 2
end

@testset "set heights for a tree with six taxa" begin
    newick = "((((A:1,B:1):2,(C:2,D:2):1):4,(E:5, F:5):2):0);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    trees_with_short_branches = transform_short_branches.(cladified_trees)
    recovered_trees = invert_short_branches.(trees_with_short_branches)

    @test recovered_trees[1].parameters[recovered_trees[1].root] ≈ 7
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:6, 6)].clade1] ≈ 3
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:6, 6)].clade2] ≈ 5
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:4, 6)].clade1] ≈ 1
    @test recovered_trees[1].parameters[recovered_trees[1].splits[Clade(1:4, 6)].clade2] ≈ 2
end