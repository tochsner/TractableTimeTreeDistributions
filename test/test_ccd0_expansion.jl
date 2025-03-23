@testset "test expansion of ccd0 with two trees and five taxa" begin
    newicks = ["(((A, B), C), (D, E));", "((A, (B, (C, D))), E);"]
    trees = load_trees_from_newick(newicks) .|> cladify_tree

    ccd = CCD0(trees)

    @test ccd.num_taxa == 5
    @test ccd.root_clade == Clade(1:5, 5)

    @test Clade([1, 2], 5) in ccd.clades
    @test Clade([1, 2, 3], 5) in ccd.clades
    @test Clade([4, 5], 5) in ccd.clades
    @test Clade([3, 4], 5) in ccd.clades
    @test Clade([2, 3, 4], 5) in ccd.clades
    @test Clade([1, 2, 3, 4], 5) in ccd.clades

    # observed splits

    @test Split(Clade(1, 5), Clade(2, 5)) in ccd.splits
    @test Split(Clade(4, 5), Clade(5, 5)) in ccd.splits
    @test Split(Clade(3, 5), Clade(4, 5)) in ccd.splits
    @test Split(Clade([1, 2], 5), Clade(3, 5)) in ccd.splits
    @test Split(Clade(2, 5), Clade([3, 4], 5)) in ccd.splits
    @test Split(Clade([1, 2, 3], 5), Clade([4, 5], 5)) in ccd.splits
    @test Split(Clade(1, 5), Clade([2, 3, 4], 5)) in ccd.splits
    @test Split(Clade([1, 2, 3, 4], 5), Clade(5, 5)) in ccd.splits

    # expanded splits

    @test Split(Clade([1, 2], 5), Clade([3, 4], 5)) in ccd.splits
    @test Split(Clade([1, 2, 3], 5), Clade(4, 5)) in ccd.splits

    @test length(ccd.splits) == 10
end

@testset "test expansion of ccd0 with three trees and five taxa" begin
    newicks = ["(((A, B), C), (D, E));", "(((A, (B, C)), D), E);", "(((A, (B, C)), E), D);"]
    trees = load_trees_from_newick(newicks) .|> cladify_tree

    ccd = CCD0(trees)

    @test ccd.num_taxa == 5
    @test ccd.root_clade == Clade(1:5, 5)

    @test Clade([1, 2], 5) in ccd.clades
    @test Clade([4, 5], 5) in ccd.clades
    @test Clade([1, 2, 3], 5) in ccd.clades
    @test Clade([2, 3], 5) in ccd.clades
    @test Clade([1, 2, 3, 4], 5) in ccd.clades
    @test Clade([1, 2, 3, 5], 5) in ccd.clades

    # observed splits

    @test Split(Clade(1, 5), Clade(2, 5)) in ccd.splits
    @test Split(Clade(2, 5), Clade(3, 5)) in ccd.splits
    @test Split(Clade(4, 5), Clade(5, 5)) in ccd.splits
    @test Split(Clade(1, 5), Clade([2, 3], 5)) in ccd.splits
    @test Split(Clade([1, 2], 5), Clade(3, 5)) in ccd.splits
    @test Split(Clade([1, 2, 3], 5), Clade([4, 5], 5)) in ccd.splits
    @test Split(Clade([1, 2, 3], 5), Clade(4, 5)) in ccd.splits
    @test Split(Clade([1, 2, 3], 5), Clade(5, 5)) in ccd.splits
    @test Split(Clade([1, 2, 3, 4], 5), Clade(5, 5)) in ccd.splits
    @test Split(Clade([1, 2, 3, 5], 5), Clade(4, 5)) in ccd.splits

    # there are no expanded splits

    @test length(ccd.splits) == 10
end
