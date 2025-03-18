@testset "create ccd with one tree and two taxa" begin
    newick = "(A:1,B:1);"
    trees = load_trees_from_newick(newick)

    ccd = CCD1(trees)

    @test ccd.num_taxa == 2
    @test ccd.root_clade == Clade(1:2, 2)
    
    @test Clade(1:2, 2) in ccd.clades

    @test Split(Clade(1, 2), Clade(2, 2)) in ccd.splits

    @test ccd.num_clade_occurrences[Clade(1:2, 2)] == 1
    
    @test ccd.num_split_occurrences[Split(Clade(1, 2), Clade(2, 2))] == 1
end

@testset "create ccd with two trees and four taxa" begin
    newicks = [
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "(A:1,(B:1,(C:1,D:1):1));",
    ]
    trees = load_trees_from_newick(newicks)

    ccd = CCD1(trees)

    @test ccd.num_taxa == 4
    @test ccd.root_clade == Clade(1:4, 4)
    
    @test Clade([1, 2], 4) in ccd.clades
    @test Clade([3, 4], 4) in ccd.clades
    @test Clade([2, 3, 4], 4) in ccd.clades
    @test Clade([1, 2, 3, 4], 4) in ccd.clades

    @test Split(Clade(1, 4), Clade(2, 4)) in ccd.splits
    @test Split(Clade(3, 4), Clade(4, 4)) in ccd.splits
    @test Split(Clade([1, 2], 4), Clade([3, 4], 4)) in ccd.splits
    @test Split(Clade([2, 3, 4], 4), Clade(1, 4)) in ccd.splits

    @test ccd.num_clade_occurrences[Clade([1, 2], 4)] == 1
    @test ccd.num_clade_occurrences[Clade([3, 4], 4)] == 2
    @test ccd.num_clade_occurrences[Clade([2, 3, 4], 4)] == 1
    @test ccd.num_clade_occurrences[Clade([1, 2, 3, 4], 4)] == 2
end