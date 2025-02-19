@testset "embed tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):1);"
    trees = load_trees_from_newick(newick)
    cladified_trees = map(tree -> cladify_tree(tree).clades, trees)

    embeddings = get_embeddings(HeightRatioParameterization(), cladified_trees)

    @test embeddings[1][Clade(1, 4)] ≈ 1.0
    @test embeddings[1][Clade(2, 4)] ≈ 1.0
    @test embeddings[1][Clade(3, 4)] ≈ 1.0
    @test embeddings[1][Clade(4, 4)] ≈ 1.0
    @test embeddings[1][Clade(1:3, 4)] ≈ 2 / 3
    @test embeddings[1][Clade(2:3, 4)] ≈ 0.5
end