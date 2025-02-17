@testset "clade split contains clades" begin
    clade1 = Clade([1, 2], 10)
    clade2 = Clade([3, 4], 10)
    clade3 = Clade([3, 4, 5], 10)
    
    clade_split = CladeSplit(clade1, clade2)
    
    @test clade1 in clade_split
    @test clade2 in clade_split
    @test (clade3 in clade_split) == false
end
