@testset "test dirichlet parameter estimation on Yule-10 treeset" begin
    reference_trees = load_trees("ref_trees.trees") .|> cladify_tree

    height_ratio_dist = HeightRatioDist{IndependentDist{LogNormal},TreeDirichletDist}(reference_trees)
    dirichlet = height_ratio_dist.ratios

    # we compare it to the Java implementation

    @test dirichlet.K ≈ 18.2405410265742

    @test isapprox(dirichlet.alphas[Clade([1, 2, 3, 5, 6, 7, 8, 9, 10], 10)], 17.010542676191466; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 3, 4, 5, 6, 8, 9, 10], 10)], 16.23784594347937; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 3, 5, 6, 8, 9, 10], 10)], 12.990276754783498; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 2, 3, 4, 7, 8], 10)], 17.228913270834177; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([4, 5, 6, 9, 10], 10)], 12.591788864694154; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 2, 3, 7, 8], 10)], 15.77264747191974; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 3, 4, 8], 10)], 12.990276754783498; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 3, 8], 10)], 3.723821963541383; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([6, 9, 10], 10)], 6.650088226252324; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([5, 9], 10)], 6.649316703366856; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([2, 7], 10)], 5.570124712797397; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 8], 10)], 0.9911430255665317; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([3, 8], 10)], 3.044780490026429; rtol=0.01)
    @test isapprox(dirichlet.alphas[Clade([1, 3], 10)], 2.9790575708331066; rtol=0.01)
end

@testset "test dirichlet log density on Yule-10 treeset" begin
    reference_trees = load_trees("ref_trees.trees")
    query_trees = load_trees("query_trees.trees")

    distribution = TractableTimeTreeDist{CCD1,HeightRatioDist{IndependentDist{LogNormal},TreeDirichletDist}}(reference_trees)

    # we compare it to the Java implementation

    @test isapprox(log_density(distribution, query_trees[1]), 33.869891395562924; rtol=0.01)
    @test isapprox(log_density(distribution, query_trees[2]), 34.80513327135303; rtol=0.01)
    @test isapprox(log_density(distribution, query_trees[3]), 35.201671510034664; rtol=0.01)
    @test isapprox(log_density(distribution, query_trees[4]), 32.8378870900746; rtol=0.01)
    @test isapprox(log_density(distribution, query_trees[4]), 32.94259841913006; rtol=0.01)
end
