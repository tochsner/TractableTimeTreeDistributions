@testset "leaf clade only contains leaf" begin
    leaf_clade = Clade(2, 10)
    @test 2 in leaf_clade
    @test (3 in leaf_clade) == false
end

@testset "subset is contained in clade" begin
    clade_1 = Clade([1, 2], 10)
    clade_2 = Clade([1, 2, 3], 10)
    clade_3 = Clade([1, 2, 3, 4], 10)

    @test clade_1 in clade_1
    @test clade_1 in clade_2
    @test clade_1 in clade_3

    @test (clade_2 in clade_1) == false
    @test clade_2 in clade_2
    @test clade_2 in clade_3

    @test (clade_3 in clade_1) == false
    @test (clade_3 in clade_2) == false
    @test clade_3 in clade_3
end

@testset "clade union of two clades" begin
    leaf_clade_1 = Clade(1, 10)
    leaf_clade_2 = Clade(2, 10)
    leaf_clade_3 = Clade(3, 10)

    clade_1_2 = union(leaf_clade_1, leaf_clade_2)
    @test 1 in clade_1_2
    @test 2 in clade_1_2
    @test (3 in clade_1_2) == false

    clade_1_2_3 = union(clade_1_2, leaf_clade_3)
    @test 1 in clade_1_2_3
    @test 2 in clade_1_2_3
    @test 3 in clade_1_2_3
    @test (4 in clade_1_2_3) == false
end
