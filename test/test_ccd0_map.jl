@testset "test ccd map on bigger Yule-10 treeset" begin
    reference_trees = load_trees("ref_trees.trees") .|> cladify_tree
    true_map_tree = load_trees("ccd0_map_tree.trees") |> first |> cladify_tree

    ccd = CCD0(reference_trees)
    map_tree = point_estimate(ccd)

    @test map_tree.tip_names == true_map_tree.tip_names
    @test map_tree.root == true_map_tree.root
    @test map_tree.splits == true_map_tree.splits
end
