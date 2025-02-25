@testset "fit on bigger Yule-10 treeset" begin
    trees = load_trees("ref_trees.trees")
    query_trees = load_trees("query_trees.trees")

    distribution = TractableTimeTreeDist{
        CCD1,
        ShorterBranchDist{IndependentDist{LogNormal}},
    }(trees)

    density = log_density(distribution, query_trees[1])
    @test density ≈ 34.97404368489742 # obtained using the java implementation

    density = log_density(distribution, query_trees[2])
    @test density ≈ 35.74017239860576 # obtained using the java implementation
end

@testset "fit on bigger Yule-10 treeset and sample" begin
    trees = load_trees("ref_trees.trees")

    distribution = TractableTimeTreeDist{
        CCD1,
        ShorterBranchDist{IndependentDist{LogNormal}}
    }(trees)

    sample_tree(distribution)
end
