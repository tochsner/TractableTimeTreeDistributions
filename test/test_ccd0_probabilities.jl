@testset "test probabilities on bigger Yule-10 treeset" begin
    reference_trees = load_trees("ref_trees.trees") .|> cladify_tree
    query_trees = load_trees("query_trees.trees")

    ccd = CCD0(reference_trees)

    # compare with values from the CCD java implementation 
    @test isapprox(log_density(ccd, query_trees[1]), -2.45875626064967)
    @test isapprox(log_density(ccd, query_trees[2]), -1.3102171524213377)
    @test isapprox(log_density(ccd, query_trees[3]), -1.3102171524213377)
    @test isapprox(log_density(ccd, query_trees[4]), -2.478810462065421)
    @test isapprox(log_density(ccd, query_trees[5]), -1.9936326853413378)
end
