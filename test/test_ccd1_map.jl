@testset "test ccd map with two trees and four taxa" begin
    newicks = [
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "((A:1,B:1):1,(C:1,D:1):1):1;",
        "(A:1,(B:1,(C:1,D:1):1));",
    ]
    trees = load_trees_from_newick(newicks)

    ccd = CCD1(trees)

    map_tree = get_most_likely_tree(ccd)

    # @test Clade(1:2, 4) in map_tree
    # @test Clade(3:4, 4) in map_tree
end

# @testset "test ccd map with two other trees and four taxa" begin
#     newicks = [
#         "((A:1,B:1):1,(C:1,D:1):1):1;",
#         "(A:1,(B:1,(C:1,D:1):1));",
#         "(A:1,(B:1,(C:1,D:1):1));",
#     ]
#     trees = load_trees_from_newick(newicks)

#     ccd = CCD1(trees)

#     map_tree = get_most_likely_tree(ccd)

#     @test Clade(3:4, 4) in map_tree
#     @test Clade(2:4, 4) in map_tree
# end

# @testset "test ccd map with two other trees and four taxa" begin
#     newicks = [
#         "(((A, B), C), (D, E));",
#         "(((A, B), C), (D, E));",
#         "(((A, B), C), (D, E));",
#         "(((A, (B, C)), D), E);",
#         "(((A, (B, C)), D), E);",
#         "(((A, (B, C)), E), D);",
#         "(((A, (B, C)), E), D);",
#     ]
#     trees = load_trees_from_newick(newicks)

#     ccd = CCD1(trees)

#     map_tree = get_most_likely_tree(ccd)

#     @test Clade(2:3, 5) in map_tree
#     @test Clade(4:5, 5) in map_tree
#     @test Clade(1:3, 5) in map_tree
# end

# @testset "test ccd map on bigger Yule-10 treeset" begin
#     reference_trees = load_trees("ref_trees.trees")
#     true_map_tree = load_trees("map_tree.trees") |> first |> cladify_tree

#     ccd = CCD1(reference_trees)
#     map_tree = get_most_likely_tree(ccd)

#     @test map_tree == true_map_tree.clades
# end