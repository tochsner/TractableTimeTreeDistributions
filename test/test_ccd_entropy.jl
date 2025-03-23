@testset "test ccd0 entropy on Yule-10 treeset" begin
    reference_trees = load_trees("ref_trees.trees") .|> cladify_tree
    true_map_tree = load_trees("ccd0_map_tree.trees") |> first |> cladify_tree

    ccd = CCD0(reference_trees)
    cdd_entropy = entropy(ccd)

    @test cdd_entropy ≈ 2.0346400956014645
end

@testset "test ccd1 entropy on Yule-10 treeset" begin
    reference_trees = load_trees("ref_trees.trees") .|> cladify_tree
    true_map_tree = load_trees("ccd0_map_tree.trees") |> first |> cladify_tree

    ccd = CCD1(reference_trees)
    cdd_entropy = entropy(ccd)

    @test cdd_entropy ≈ 2.0703082105171133
end
