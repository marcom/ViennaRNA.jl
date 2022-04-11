using Test
using ViennaRNA
using Unitful

@testset "FoldCompound" begin
    seq = "GGGAAACCC"
    @test length(FoldCompound(seq)) == 9
    fc = FoldCompound(seq; uniq_ML=1)
    @test length(fc) == 9
    # TODO
    # @test fc.ptr.params.model_details.uniq_ML[] == 1
    fc = FoldCompound(seq; options=ViennaRNA.LibRNA.VRNA_OPTION_MFE)
    @test length(fc) == 9
end

@testset "Pairtable" begin
    pt = ViennaRNA.Pairtable("(((...)))")
    @test length(pt) == 9
    # getindex
    @test pt[1] == 9
    @test pt[2] == 8
    @test pt[3] == 7
    @test pt[4] == 0
    @test pt[5] == 0
    @test pt[6] == 0
    @test pt[7] == 3
    @test pt[8] == 2
    @test pt[9] == 1
    # setindex!
    pt[3] = 4
    @test pt[3] == 4
    @test pt[4] == 3
    @test pt[7] == 0
    # show
    buf = IOBuffer()
    show(buf, MIME("text/plain"), pt)
    @test length(String(take!(buf))) > 0
end

@testset "bp_distance" begin
    @test ViennaRNA.bp_distance("((..))", "(....)") == 1
    @test_throws ArgumentError ViennaRNA.bp_distance("()", ".")
end

@testset "tree_edit_dist" begin
    @test ViennaRNA.tree_edit_dist("...", "()") == 3.0f0
end

@testset "mean_bp_distance" begin
    seq = "GGGAAACCC"
    fc = FoldCompound(seq)
    @test_throws ArgumentError ViennaRNA.mean_bp_distance(fc)
    ViennaRNA.partfn(fc)
    @test ViennaRNA.mean_bp_distance(fc) isa AbstractFloat
    @test ViennaRNA.mean_bp_distance(seq) isa AbstractFloat
end

@testset "ensemble_defect" begin
    seq = "GGGAAACCCC"
    str = "(((....)))"
    fc = FoldCompound(seq)
    ViennaRNA.partfn(fc)
    @test ViennaRNA.ensemble_defect(fc, str) isa AbstractFloat
    @test ViennaRNA.ensemble_defect(fc, ViennaRNA.Pairtable(str)) isa AbstractFloat
    @test_throws ArgumentError ViennaRNA.ensemble_defect(fc, ".")
    @test_throws ArgumentError ViennaRNA.ensemble_defect(fc, ViennaRNA.Pairtable("."))
    @test ViennaRNA.ensemble_defect(seq, ViennaRNA.Pairtable(str)) isa AbstractFloat
    @test ViennaRNA.ensemble_defect(seq, str) isa AbstractFloat
end

@testset "prob_of_structure" begin
    seq = "GGGAAACCCC"
    str = "(((....)))"
    fc = FoldCompound(seq)
    ViennaRNA.partfn(fc)
    @test ViennaRNA.prob_of_structure(fc, str) isa AbstractFloat
    @test_throws ArgumentError ViennaRNA.prob_of_structure(fc, ".")
    @test ViennaRNA.prob_of_structure(seq, str) isa AbstractFloat
end

@testset "energy" begin
    seq = "GGGAAAACCCC"
    str = "(((....)))."
    fc = FoldCompound(seq)
    @test ViennaRNA.energy(fc, str) isa Unitful.Quantity
    @test ViennaRNA.energy(seq, str) isa Unitful.Quantity
    @test_throws ArgumentError ViennaRNA.energy(fc, ".")
    @test_throws ArgumentError ViennaRNA.energy("(...)", ".")
end

@testset "mfe" begin
    seq = "GGGAAACCCC"
    fc = FoldCompound(seq)
    @test ViennaRNA.mfe(fc) isa Tuple{String,Unitful.Quantity}
    @test ViennaRNA.mfe(seq) isa Tuple{String,Unitful.Quantity}
end

@testset "partfn" begin
    seq = "GGGAAACCCC"
    fc = FoldCompound(seq)
    @test ViennaRNA.partfn(fc) isa Tuple{String,Unitful.Quantity}
    @test ViennaRNA.partfn(seq) isa Tuple{String,Unitful.Quantity}
end

@testset "bpp" begin
    seq = "GGGAAACCCC"
    n = length(seq)

    fc = FoldCompound(seq)
    @test_throws ArgumentError ViennaRNA.bpp(fc)
    ViennaRNA.partfn(fc)
    p = ViennaRNA.bpp(fc)
    @test eltype(p) <: AbstractFloat
    @test size(p) == (n,n)

    p = ViennaRNA.bpp(seq)
    @test eltype(p) <: AbstractFloat
    @test size(p) == (n,n)
end

@testset "pbacktrack" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq; uniq_ML = 1)
    ViennaRNA.partfn(fc)
    @test length(ViennaRNA.pbacktrack(fc)) == 1
    @test length(ViennaRNA.pbacktrack(seq)) == 1
    s = ViennaRNA.pbacktrack(fc; num_samples=10)
    @test length(s) == 10
    s = ViennaRNA.pbacktrack(fc; num_samples=5,
                             options=ViennaRNA.LibRNA.VRNA_PBACKTRACK_NON_REDUNDANT)
    @test length(s) == 5
end

@testset "mea" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    ViennaRNA.partfn(fc)
    @test ViennaRNA.mea(fc) isa Tuple{String,AbstractFloat}
    @test ViennaRNA.mea(fc; gamma=1.0) isa Tuple{String,AbstractFloat}
    @test ViennaRNA.mea(fc; gamma=0.5) isa Tuple{String,AbstractFloat}
    @test ViennaRNA.mea(seq) isa Tuple{String,AbstractFloat}
    @test ViennaRNA.mea(seq; gamma=1.0) isa Tuple{String,AbstractFloat}
    @test ViennaRNA.mea(seq; gamma=0.5) isa Tuple{String,AbstractFloat}
    fc = FoldCompound(seq)
    @test_throws ArgumentError ViennaRNA.mea(fc)
end

@testset "centroid" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    ViennaRNA.partfn(fc)
    @test ViennaRNA.centroid(fc) isa Tuple{String,AbstractFloat}
    @test ViennaRNA.centroid(seq) isa Tuple{String,AbstractFloat}
    fc = FoldCompound(seq)
    @test_throws ArgumentError ViennaRNA.centroid(fc)
end

@testset "subopt" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq; uniq_ML = 1)
    s = ViennaRNA.subopt(fc; delta=5u"kcal/mol")
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = ViennaRNA.subopt(fc; delta=5u"kcal/mol", sorted=true)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = ViennaRNA.subopt(seq; delta=5u"kcal/mol")
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = ViennaRNA.subopt(seq; delta=5u"kcal/mol", sorted=true)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
end

@testset "subopt_zuker" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    s = ViennaRNA.subopt_zuker(fc)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = ViennaRNA.subopt_zuker(seq)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
end

@testset "inverse_fold" begin
    s = "((((....))))"
    @test ViennaRNA.inverse_fold("A"^length(s), s) isa Tuple{String,AbstractFloat}
    @test_throws ArgumentError ViennaRNA.inverse_fold("A", "()")
    @test ViennaRNA.inverse_pf_fold("A"^length(s), s) isa Tuple{String,Unitful.Quantity}
    @test_throws ArgumentError ViennaRNA.inverse_pf_fold("A", "()")
end

@testset "neighbors" begin
    fc = FoldCompound("GGGAAACCC")
    pt = ViennaRNA.Pairtable(".((...)).")
    @test ViennaRNA.neighbors(fc, pt) == [[(-2, -8)], [(-3, -7)], [(1, 9)]]
end

@testset "plot_coords" begin
    function test_plot_xy(s, x, y)
        @test eltype(x) <: AbstractFloat
        @test length(x) == length(s)
        @test eltype(y) <: AbstractFloat
        @test length(y) == length(s)
    end
    s = "(((...)))"
    pt = ViennaRNA.Pairtable(s)

    # test: plot_coords(::String), plot_coords(::ViennaRNA.Pairtable)
    x, y = ViennaRNA.plot_coords(s)
    test_plot_xy(s, x, y)
    x, y = ViennaRNA.plot_coords(pt)
    test_plot_xy(pt, x, y)
    for plot_type in (:simple, :naview, :circular, :turtle, :puzzler)
        x, y = ViennaRNA.plot_coords(s; plot_type)
        test_plot_xy(s, x, y)
        x, y = ViennaRNA.plot_coords(pt; plot_type)
        test_plot_xy(pt, x, y)
    end
    @test_throws ArgumentError x, y = ViennaRNA.plot_coords(s; plot_type = :unknown)
    @test_throws ArgumentError x, y = ViennaRNA.plot_coords(pt; plot_type = :unknown)
end
