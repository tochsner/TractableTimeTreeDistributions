@testset "test ccd map with two trees and four taxa" begin
    newicks = [
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "(A:1,(B:1,(C:1,D:1):1));",
    ]
    trees = load_trees_from_newick(newicks)

    ccd = CCD1(trees)

    map_tree = get_most_likely_tree(ccd)

    @test CladeSplit(Clade(1, 4), Clade(2, 4)) in map_tree
    @test CladeSplit(Clade(3, 4), Clade(4, 4)) in map_tree
    @test CladeSplit(Clade(1:2, 4), Clade(3:4, 4)) in map_tree
end


@testset "test ccd map with two other trees and four taxa" begin
    newicks = [
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "(A:1,(B:1,(C:1,D:1):1));",
        "(A:1,(B:1,(C:1,D:1):1));",
    ]
    trees = load_trees_from_newick(newicks)

    ccd = CCD1(trees)

    map_tree = get_most_likely_tree(ccd)

    @test CladeSplit(Clade(3, 4), Clade(4, 4)) in map_tree
    @test CladeSplit(Clade(2, 4), Clade(3:4, 4)) in map_tree
    @test CladeSplit(Clade(1, 4), Clade(2:4, 4)) in map_tree
end


@testset "test ccd map with two other trees and four taxa" begin
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

    map_tree = get_most_likely_tree(ccd)

    @test CladeSplit(Clade(2, 5), Clade(3, 5)) in map_tree
    @test CladeSplit(Clade(4, 5), Clade(5, 5)) in map_tree
    @test CladeSplit(Clade(1, 5), Clade(2:3, 5)) in map_tree
    @test CladeSplit(Clade(1:3, 5), Clade(4:5, 5)) in map_tree
end
