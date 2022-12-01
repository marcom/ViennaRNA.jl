@testset "prob_of_basepairs" begin
    showtestset()
    sequence  = "GGGAACCUC"
    structure = "(((...)))"
    pt = Pairtable(structure)
    fc = FoldCompound(sequence)
    partfn(fc)
    bppm = bpp(fc)
    @test prob_of_basepairs(bppm, pt) ==
        prob_of_basepairs(bppm, structure) ==
        prob_of_basepairs(fc, pt) ==
        prob_of_basepairs(fc, structure) ==
        prob_of_basepairs(sequence, pt) ==
        prob_of_basepairs(sequence, structure)
end
