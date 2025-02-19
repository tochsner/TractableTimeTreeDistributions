@testset "parameterize tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):1);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)

    @test parameterized_trees[1].parameters[Clade(2:4, 4)] ≈ 1 / 3
    @test parameterized_trees[1].parameters[Clade(3:4, 4)] ≈ 0.5
    @test parameterized_trees[1].parameters[Clade(1:4, 4)] ≈ 3
end

@testset "parameterize another tree with four taxa" begin
    newick = "((A:1,B:1):3,(C:2,D:2):2):0;" 
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)

    @test parameterized_trees[1].parameters[Clade(1:2, 4)] ≈ 0.75
    @test parameterized_trees[1].parameters[Clade(3:4, 4)] ≈ 0.5
    @test parameterized_trees[1].parameters[Clade(1:4, 4)] ≈ 4
end

@testset "set heights for tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):1);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)
    recovered_trees = set_heights(HeightRatioParameterization(), parameterized_trees)

    for clade in keys(recovered_trees[1].splits)
        @test recovered_trees[1].splits[clade].clade1.height == cladified_trees[1].splits[clade].clade1.height
        @test recovered_trees[1].splits[clade].clade2.height == cladified_trees[1].splits[clade].clade2.height
        @test recovered_trees[1].splits[clade].parent.height == cladified_trees[1].splits[clade].parent.height
    end
end