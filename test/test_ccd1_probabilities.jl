@testset "test ccd probabilities with two trees and four taxa" begin
    newicks = [
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "(A:1,(B:1,(C:1,D:1):1));",
    ]
    trees = load_trees_from_newick(newicks)

    ccd = CCD1(trees)

    @test isapprox(get_log_probability(ccd, trees[1]), log(0.5))
    @test isapprox(get_log_probability(ccd, trees[2]), log(0.5))
end

@testset "test ccd probabilities with three trees and five taxa" begin
    newicks = [
        "(((A, B), C), (D, E));",
        "(((A, B), C), (D, E));",
        "(((A, B), C), (D, E));",
        "(((A, (B, C)), D), E);",
        "(((A, (B, C)), D), E);",
        "(((A, (B, C)), E), D);",
        "(((A, (B, C)), E), D);",
    ]
    trees = load_trees_from_newick(newicks)

    ccd = CCD1(trees)

    @test isapprox(get_log_probability(ccd, CladeSplit(Clade(1:2, 5), Clade(3, 5))), log(3 / 7))
    @test isapprox(get_log_probability(ccd, CladeSplit(Clade(1, 5), Clade(2:3, 5))), log(4 / 7))

    @test isapprox(get_log_probability(ccd, CladeSplit(Clade(1:3, 5), Clade(4:5, 5))), log(3 / 7))
    @test isapprox(get_log_probability(ccd, CladeSplit(Clade(1:4, 5), Clade(5, 5))), log(2 / 7))
    @test isapprox(get_log_probability(ccd, CladeSplit(Clade([1, 2, 3, 5], 5), Clade(4, 5))), log(2 / 7))
    
    @test isapprox(get_log_probability(ccd, trees[1]), log(9 / 49))
    @test isapprox(get_log_probability(ccd, trees[4]), log(8 / 49))
    @test isapprox(get_log_probability(ccd, trees[6]), log(8 / 49))

    other_newicks = [
        "((A, (B, C)), (D, E));",
        "((((A, B), C), D), E);",
        "((((A, B), C), E), D);",
    ]
    other_trees = load_trees_from_newick(other_newicks)

    @test isapprox(get_log_probability(ccd, other_trees[1]), log(12 / 49))
    @test isapprox(get_log_probability(ccd, other_trees[2]), log(6/ 49))
    @test isapprox(get_log_probability(ccd, other_trees[3]), log(6 / 49))
end