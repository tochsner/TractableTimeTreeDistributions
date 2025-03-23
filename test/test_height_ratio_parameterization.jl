@testset "fit on trees with four taxa" begin
    newick = ["(A:3,(B:2,(C:1,D:1):1):1);", "(A:5,(B:2.5,(C:2,D:2):0.5):2.5);"]
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterization = HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}(cladified_trees)

    @test parameterization.height.distributions[Clade(1:4, 4)].μ ≈ (log(3) + log(5)) / 2
    @test parameterization.height.distributions[Clade(1:4, 4)].σ ≈
          sqrt(((log(3) - (log(3) + log(5)) / 2)^2 + (log(5) - (log(3) + log(5)) / 2)^2) / 1)
end

@testset "fit on bigger Yule-10 treeset" begin
    trees = load_trees("ref_trees.trees")
    query_trees = load_trees("query_trees.trees")

    distribution = TractableTimeTreeDist{CCD1,HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}}(trees)

    density = log_density(distribution, query_trees[1])
    @test density ≈ 34.464263095568896 # obtained using the java implementation

    density = log_density(distribution, query_trees[2])
    @test density ≈ 36.63401964278324 # obtained using the java implementation
end

@testset "fit on bigger Yule-10 treeset and sample" begin
    trees = load_trees("ref_trees.trees")

    distribution = TractableTimeTreeDist{CCD1,HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}}(trees)

    sample_tree(distribution)
end
