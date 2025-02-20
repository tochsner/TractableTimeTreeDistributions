@testset "fit on trees with four taxa" begin
    newick = [
        "(A:3,(B:2,(C:1,D:1):1):1);",
        "(A:5,(B:2.5,(C:2,D:2):0.5):2.5);"
    ]
    trees = load_trees_from_newick(newick)
    cladified_trees = map(cladify_tree, trees)

    parameterization = HeightRatioTransform(
        IndependentDist{LogNormal}(),
        IndependentDist{Beta}()
    )

    fit!(parameterization, cladified_trees)

    @test parameterization.height.distributions[Clade(1:4, 4)].μ ≈ (log(3) + log(5)) / 2
    @test parameterization.height.distributions[Clade(1:4, 4)].σ ≈ sqrt(((log(3) - (log(3) + log(5)) / 2)^2 + (log(5) - (log(3) + log(5)) / 2)^2) / 1)
end

@testset "fit on bigger Yule-10 treeset" begin
    trees = load_trees("ref_trees.trees")
    cladified_trees = map(cladify_tree, trees)

    ccd = CCD1(trees)

    parameterization = HeightRatioTransform(
        IndependentDist{LogNormal}(),
        IndependentDist{Beta}()
    )
    fit!(parameterization, cladified_trees)

    query_trees = load_trees("query_trees.trees")
    cladified_query_trees = map(cladify_tree, query_trees)

    weight_density = log_density(parameterization, cladified_query_trees[1])
    ccd_density = log_density(ccd, query_trees[1])

    @test weight_density + ccd_density ≈ 34.464263095568896 # obtained using the java implementation
    
    weight_density = log_density(parameterization, cladified_query_trees[2])
    ccd_density = log_density(ccd, query_trees[2])
    
    @test weight_density + ccd_density ≈ 36.63401964278324 # obtained using the java implementation
end