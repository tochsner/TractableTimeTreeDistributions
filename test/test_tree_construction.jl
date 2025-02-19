using Phylo

@testset "construct tree with one taxa" begin
    newick = "(A:1);"
    tree = load_trees_from_newick(newick)[1]
    cladified_tree = cladify_tree(tree)

    reconstructed_tree = construct_tree(cladified_tree.clades)
    @test nleaves(reconstructed_tree) == 1
    @test nbranches(reconstructed_tree) == 0
end

@testset "construct tree with four taxa" begin
    newick = "(A:3,(B:2,(C:1,D:1):1):1);"
    tree = load_trees_from_newick(newick)[1]
    cladified_tree = cladify_tree(tree)

    reconstructed_tree = construct_tree(cladified_tree.clades)
    @test nleaves(reconstructed_tree) == 4
    @test nbranches(reconstructed_tree) == 6
end