@testset "parameterize tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):0);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)

    @test parameterized_trees[1].parameters[Clade(2:4, 4)] ≈ 1 / 3
    @test parameterized_trees[1].parameters[Clade(3:4, 4)] ≈ 0.5
    @test parameterized_trees[1].parameters[Clade(1:4, 4)] ≈ 3
end

@testset "parameterize another tree with four taxa" begin
    newick = "((A:1,B:1):2,(C:2,D:2):1):0;" 
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)

    @test parameterized_trees[1].parameters[Clade(1:2, 4)] ≈ 2 / 3
    @test parameterized_trees[1].parameters[Clade(3:4, 4)] ≈ 1 / 3
    @test parameterized_trees[1].parameters[Clade(1:4, 4)] ≈ 3
end

@testset "parameterize a tree with six taxa" begin
    newick = "((((A:1,B:1):2,(C:2,D:2):1):4,(E:5, F:5):2):0);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)

    @test parameterized_trees[1].parameters[Clade(1:4, 6)] ≈ 4 / 7
    @test parameterized_trees[1].parameters[Clade(1:2, 6)] ≈ 2 / 3
    @test parameterized_trees[1].parameters[Clade(3:4, 6)] ≈ 1 / 3
    @test parameterized_trees[1].parameters[Clade(5:6, 6)] ≈ 2  /7
end

@testset "set heights for tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):0);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)
    recovered_trees = set_heights!(HeightRatioParameterization(), parameterized_trees)

    @test recovered_trees[1].root.height ≈ 3
    @test recovered_trees[1].splits[Clade(1:4, 4)].clade1.height ≈ 0
    @test recovered_trees[1].splits[Clade(1:4, 4)].clade2.height ≈ 2
    @test recovered_trees[1].splits[Clade(2:4, 4)].clade1.height ≈ 0
    @test recovered_trees[1].splits[Clade(2:4, 4)].clade2.height ≈ 1
    @test recovered_trees[1].splits[Clade(3:4, 4)].clade1.height ≈ 0 
    @test recovered_trees[1].splits[Clade(3:4, 4)].clade2.height ≈ 0 
end

@testset "set heights for another tree with four taxa" begin
    newick = "((A:1,B:1):2,(C:2,D:2):1):0;" 
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)
    recovered_trees = set_heights!(HeightRatioParameterization(), parameterized_trees)

    @test recovered_trees[1].root.height ≈ 3
    @test recovered_trees[1].splits[Clade(1:4, 4)].clade1.height ≈ 1
    @test recovered_trees[1].splits[Clade(1:4, 4)].clade2.height ≈ 2
    @test recovered_trees[1].splits[Clade(1:2, 4)].clade1.height ≈ 0
    @test recovered_trees[1].splits[Clade(1:2, 4)].clade2.height ≈ 0
    @test recovered_trees[1].splits[Clade(3:4, 4)].clade1.height ≈ 0 
    @test recovered_trees[1].splits[Clade(3:4, 4)].clade2.height ≈ 0 
end

@testset "set heights for a tree with six taxa" begin
    newick = "((((A:1,B:1):2,(C:2,D:2):1):4,(E:5, F:5):2):0);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterized_trees = parameterize(HeightRatioParameterization(), cladified_trees)
    recovered_trees = set_heights!(HeightRatioParameterization(), parameterized_trees)

    @test recovered_trees[1].root.height ≈ 7
    @test recovered_trees[1].splits[Clade(1:6, 6)].clade1.height ≈ 3
    @test recovered_trees[1].splits[Clade(1:6, 6)].clade2.height ≈ 5
    @test recovered_trees[1].splits[Clade(1:4, 6)].clade1.height ≈ 1
    @test recovered_trees[1].splits[Clade(1:4, 6)].clade2.height ≈ 2
    @test recovered_trees[1].splits[Clade(1:2, 6)].clade1.height ≈ 0
    @test recovered_trees[1].splits[Clade(1:2, 6)].clade2.height ≈ 0
    @test recovered_trees[1].splits[Clade(3:4, 6)].clade1.height ≈ 0
    @test recovered_trees[1].splits[Clade(3:4, 6)].clade2.height ≈ 0
    @test recovered_trees[1].splits[Clade(5:6, 6)].clade1.height ≈ 0
    @test recovered_trees[1].splits[Clade(5:6, 6)].clade1.height ≈ 0
end