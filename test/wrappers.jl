using Test
using ViennaRNA
using ViennaRNA: LibRNA
using Unitful

@testset "FoldCompound" begin
    for s in ["GGGAAACCC",
              ["GGG-A-CC-", "G-GGC-G-G"],
              ["GGG-A-CC-", "G-GGC-G-G", "CGC----CG"]]
        fc = FoldCompound(s)
        @test length(fc) == 9
        @test sum(size(fc)) == length(fc)

        # test FoldCompound properties
        @test fc.circular         == LibRNA.VRNA_MODEL_DEFAULT_CIRC
        @test fc.dangles          == LibRNA.VRNA_MODEL_DEFAULT_DANGLES
        @test fc.gquadruplex      == LibRNA.VRNA_MODEL_DEFAULT_GQUAD
        @test fc.has_exp_matrices == false
        @test fc.log_ML           == LibRNA.VRNA_MODEL_DEFAULT_LOG_ML
        @test fc.min_loop_size    == LibRNA.TURN
        @test fc.no_GU_basepairs  == LibRNA.VRNA_MODEL_DEFAULT_NO_GU
        @test fc.no_GU_closure    == LibRNA.VRNA_MODEL_DEFAULT_NO_GU_CLOSURE
        @test fc.no_lonely_pairs  == LibRNA.VRNA_MODEL_DEFAULT_NO_LP
        @test fc.nstrands         == 1
        @test fc.params_name      == "RNA - Turner 2004"
        @test fc.special_hairpins == LibRNA.VRNA_MODEL_DEFAULT_SPECIAL_HP
        @test fc.temperature      == 37u"°C"
        @test fc.uniq_ML          == LibRNA.VRNA_MODEL_DEFAULT_UNIQ_ML

        # options
        fc = FoldCompound(s; options=ViennaRNA.LibRNA.VRNA_OPTION_MFE)
        @test length(fc) == 9
        # use different energy parameter sets
        fc = FoldCompound(s; params=:RNA_Turner1999)
        @test length(fc) == 9
        @test_throws ArgumentError FoldCompound(s; params=:UNKNOWN_ENERGY_PARAMS)
        # use different temperatures to evaluate energies
        fc = FoldCompound(s; temperature=55u"°C")
        @test length(fc) == 9
        @test fc.temperature == 55u"°C"
        fc = FoldCompound(s; temperature=310u"K")
        @test length(fc) == 9
        @test fc.temperature == 310u"K"
        @test_throws Unitful.DimensionError FoldCompound(s; temperature=100u"m")

        # circular
        fc = FoldCompound(s; circular=true)
        @test length(fc) == 9
        @test fc.circular == true
        # dangles
        fc = FoldCompound(s; dangles=1)
        @test length(fc) == 9
        @test fc.dangles == 1
        @test_throws ArgumentError FoldCompound(s; dangles=-1)
        @test_throws ArgumentError FoldCompound(s; dangles=4)
        # gquadruplex
        fc = FoldCompound(s; gquadruplex=true)
        @test length(fc) == 9
        @test fc.gquadruplex == true
        # log_ML
        fc = FoldCompound(s; log_ML=true)
        @test length(fc) == 9
        @test fc.log_ML == true
        # min_loop_size
        fc = FoldCompound(s; min_loop_size=2)
        @test length(fc) == 9
        @test fc.dangles == 2
        @test_throws ArgumentError FoldCompound(s; min_loop_size=-1)
        # no_GU_basepairs
        fc = FoldCompound(s; no_GU_basepairs=true)
        @test length(fc) == 9
        @test fc.no_GU_basepairs == true
        # no_GU_closure
        fc = FoldCompound(s; no_GU_closure=true)
        @test length(fc) == 9
        @test fc.no_GU_closure == true
        # no_lonely_pairs
        fc = FoldCompound(s; no_lonely_pairs=true)
        @test length(fc) == 9
        @test fc.no_lonely_pairs == true
        # special_hairpins
        fc = FoldCompound(s; special_hairpins=false)
        @test length(fc) == 9
        @test fc.special_hairpins == false
        # uniq_ML
        fc = FoldCompound(s; uniq_ML=true)
        @test length(fc) == 9
        @test fc.uniq_ML == true
        # show
        buf = IOBuffer()
        show(buf, MIME("text/plain"), fc)
        @test length(String(take!(buf))) > 0
    end

    # DNA
    fc = FoldCompound("GATTACA"; params=:DNA_Mathews1999)
    @test length(fc) == 7

    # multiple strands
    fc = FoldCompound("GGG&AAA&CCC")
    @test length(fc) == 9
    @test size(fc) == (3, 3, 3)
    @test fc.nstrands == 3

    # multistrand alifold (not supported yet)
    @test_throws ArgumentError fc = FoldCompound(["GGG&AAA&CCC", "G-G&AA-&CCA"])
    # @test length(fc) == 9
    # @test size(fc) == (3, 3, 3)
    # @test fc.nstrands == 3

    # multiple strands don't work in comparative (alifold) mode
    @test_throws ArgumentError FoldCompound(["GGGA", "GCGG", "G"])
    @test_throws ArgumentError FoldCompound(["GG&C", "AA&U"])

    # test setting of min_loop_size with single and multiple strands
    # - this didn't work in ViennaRNA-2.5.0 with multiple strands
    #   (e.g. "GGG&CCC") and h=3 (min_loop_size was always set to 0 by
    #   ViennaRNA). Works as of ViennaRNA-2.5.1.
    for h = 0:5
        fc = FoldCompound("GGG"; min_loop_size=h)
        @test fc.min_loop_size == h
        fc = FoldCompound("GGG&CCC"; min_loop_size=h)
        @test fc.min_loop_size == h
        fc = FoldCompound("GGG&CCC&AAAA"; min_loop_size=h)
        @test fc.min_loop_size == h
        fc = FoldCompound(["GGG", "G-C"]; min_loop_size=h)
        @test fc.min_loop_size == h
        fc = FoldCompound(["GGG", "G-C", "-AG"]; min_loop_size=h)
        @test fc.min_loop_size == h
    end
end

@testset "Pairtable" begin
    pt = Pairtable("(((...)))")
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
    @test_throws BoundsError pt[0]
    @test_throws BoundsError pt[length(pt) + 1]
    # setindex!
    pt[3] = 4
    @test pt[3] == 4
    @test pt[4] == 3
    @test pt[7] == 0
    @test_throws ArgumentError pt[1] = 1
    @test_throws ArgumentError pt[0] = 0
    @test_throws ArgumentError pt[length(pt) + 1] = 0
    @test_throws ArgumentError pt[1] = -1
    @test_throws ArgumentError pt[1] = length(pt) + 1
    # show
    buf = IOBuffer()
    show(buf, MIME("text/plain"), pt)
    @test length(String(take!(buf))) > 0
    # basepairs
    @test Set(basepairs(Pairtable("(((...)))"))) == Set([(1,9), (2,8), (3,7)])
end

@testset "bp_distance" begin
    @test bp_distance("((..))", "(....)") == 1
    @test_throws ArgumentError bp_distance("()", ".")
end

@testset "tree_edit_dist" begin
    @test tree_edit_dist("...", "()") == 3.0f0
end

@testset "mean_bp_distance" begin
    seq = "GGGAAACCC"
    fc = FoldCompound(seq)
    @test_throws ArgumentError mean_bp_distance(fc)
    partfn(fc)
    @test mean_bp_distance(fc) isa AbstractFloat
    @test mean_bp_distance(seq) isa AbstractFloat
end

@testset "ensemble_defect" begin
    seq = "GGGAAACCCC"
    str = "(((....)))"
    fc = FoldCompound(seq)
    partfn(fc)
    @test ensemble_defect(fc, str) isa AbstractFloat
    @test ensemble_defect(fc, Pairtable(str)) isa AbstractFloat
    @test_throws ArgumentError ensemble_defect(fc, ".")
    @test_throws ArgumentError ensemble_defect(fc, Pairtable("."))
    @test ensemble_defect(seq, Pairtable(str)) isa AbstractFloat
    @test ensemble_defect(seq, str) isa AbstractFloat
end

@testset "prob_of_structure" begin
    seq = "GGGAAACCCC"
    str = "(((....)))"
    fc = FoldCompound(seq)
    partfn(fc)
    @test prob_of_structure(fc, str) isa AbstractFloat
    @test_throws ArgumentError prob_of_structure(fc, ".")
    @test prob_of_structure(seq, str) isa AbstractFloat
end

@testset "energy" begin
    seq = "GGGAAAACCCC"
    str = "(((....)))."
    fc = FoldCompound(seq)
    @test energy(fc, str) isa Unitful.Quantity
    @test energy(seq, str) isa Unitful.Quantity
    redirect_stdout(devnull) do
        # suppress stdout output from ViennaRNA for verbose option
        @test energy(fc, str; verbose=true) isa Unitful.Quantity
        @test energy(seq, str; verbose=true) isa Unitful.Quantity
    end
    @test_throws ArgumentError energy(fc, ".")
    @test_throws ArgumentError energy("(...)", ".")
    # test different energy parameter sets
    @test energy(FoldCompound(seq; params=:RNA_Turner1999), str) ≈ -2.90u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; params=:RNA_Andronescu2007), str) ≈ -2.15u"kcal/mol" atol=1e-3u"kcal/mol"
    # test different temperatures
    @test energy(FoldCompound(seq; temperature=30u"°C"), str) ≈ -3.52u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; temperature=37u"°C"), str) ≈ -2.90u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; temperature=55u"°C"), str) ≈ -1.24u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; temperature=323.15u"K"), str) ≈ -1.69u"kcal/mol" atol=1e-3u"kcal/mol" # 323.15u"K" is 50u"°C"
end

@testset "mfe" begin
    seq = "GGGAAACCCC"
    fc = FoldCompound(seq)
    @test mfe(fc) isa Tuple{String,Unitful.Quantity}
    @test mfe(seq) isa Tuple{String,Unitful.Quantity}
end

@testset "partfn" begin
    seq = "GGGAAACCCC"
    fc = FoldCompound(seq)
    @test partfn(fc) isa Tuple{String,Unitful.Quantity}
    @test partfn(seq) isa Tuple{String,Unitful.Quantity}
end

@testset "bpp" begin
    seq = "GGGAAACCCC"
    n = length(seq)

    fc = FoldCompound(seq)
    @test_throws ArgumentError bpp(fc)
    partfn(fc)
    p = bpp(fc)
    @test eltype(p) <: AbstractFloat
    @test size(p) == (n,n)

    p = bpp(seq)
    @test eltype(p) <: AbstractFloat
    @test size(p) == (n,n)
end

@testset "pbacktrack" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq; uniq_ML=true)
    partfn(fc)
    @test length(pbacktrack(fc)) == 1
    @test length(pbacktrack(seq)) == 1
    s = pbacktrack(fc; num_samples=10)
    @test length(s) == 10
    s = pbacktrack(fc; num_samples=5,
                   options=ViennaRNA.LibRNA.VRNA_PBACKTRACK_NON_REDUNDANT)
    @test length(s) == 5
end

@testset "mea" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    partfn(fc)
    @test mea(fc) isa Tuple{String,AbstractFloat}
    @test mea(fc; gamma=1.0) isa Tuple{String,AbstractFloat}
    @test mea(fc; gamma=0.5) isa Tuple{String,AbstractFloat}
    @test mea(seq) isa Tuple{String,AbstractFloat}
    @test mea(seq; gamma=1.0) isa Tuple{String,AbstractFloat}
    @test mea(seq; gamma=0.5) isa Tuple{String,AbstractFloat}
    fc = FoldCompound(seq)
    @test_throws ArgumentError mea(fc)
end

@testset "centroid" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    partfn(fc)
    @test centroid(fc) isa Tuple{String,AbstractFloat}
    @test centroid(seq) isa Tuple{String,AbstractFloat}
    fc = FoldCompound(seq)
    @test_throws ArgumentError centroid(fc)
end

@testset "subopt" begin
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq; uniq_ML=true)
    s = subopt(fc; delta=5u"kcal/mol")
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = subopt(fc; delta=5u"kcal/mol", sorted=true)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = subopt(seq; delta=5u"kcal/mol")
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = subopt(seq; delta=5u"kcal/mol", sorted=true)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
end

@testset "subopt_zuker" begin
    seq = "GGGGAAAACCCC"
    # TODO: the following sequences generate warnings
    # seq = "GGGUAAACCCAUUCAC"
    # seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    s = subopt_zuker(fc)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
    s = subopt_zuker(seq)
    @test s isa Vector{Tuple{String,Unitful.Quantity}}
end

@testset "inverse_fold" begin
    s = "((((....))))"
    @test inverse_fold("A"^length(s), s) isa Tuple{String,AbstractFloat}
    @test_throws ArgumentError inverse_fold("A", "()")
    @test inverse_pf_fold("A"^length(s), s) isa Tuple{String,Unitful.Quantity}
    @test_throws ArgumentError inverse_pf_fold("A", "()")
end

@testset "neighbors" begin
    fc = FoldCompound("GGGAAACCC")
    pt = Pairtable(".((...)).")
    @test neighbors(fc, pt) == [[(-2, -8)], [(-3, -7)], [(1, 9)]]
end

@testset "plot_coords" begin
    function test_plot_xy(s, x, y)
        @test x isa Vector{Float32}
        @test length(x) == length(s)
        @test y isa Vector{Float32}
        @test length(y) == length(s)
    end
    s = "(((...)))"
    pt = Pairtable(s)

    # test: plot_coords(::String), plot_coords(::Pairtable)
    x, y = plot_coords(s)
    test_plot_xy(s, x, y)
    x, y = plot_coords(pt)
    test_plot_xy(pt, x, y)
    for plot_type in (:simple, :naview, :circular, :turtle, :puzzler)
        x, y = plot_coords(s; plot_type)
        test_plot_xy(s, x, y)
        x, y = plot_coords(pt; plot_type)
        test_plot_xy(pt, x, y)
    end
    @test_throws ArgumentError x, y = plot_coords(s; plot_type = :unknown)
    @test_throws ArgumentError x, y = plot_coords(pt; plot_type = :unknown)
    # zero-sized inputs
    x, y = plot_coords("")
    test_plot_xy("", x, y)
    x, y = plot_coords(Pairtable(""))
    test_plot_xy(Pairtable(""), x, y)
end

@testset "heat_capacity" begin
    s = "GGGAAACCC"
    fc = FoldCompound(s)
    Tmin = 10u"°C"
    Tmax = 60u"°C"
    Tincrement = 1u"°C"
    mpoints = 2
    # TODO: this range call doesn't work without ustrip
    n = length(range(start=ustrip(Tmin), stop=ustrip(Tmax), step=ustrip(Tincrement)))

    hcs = heat_capacity(fc, Tmin, Tmax, Tincrement; mpoints)
    @test hcs isa Vector{Tuple{typeof(1.0f0u"°C"),typeof(1.0f0u"kcal/mol/K")}}
    @test length(hcs) == n
    hcs = heat_capacity(fc, Tmin, Tmax; mpoints)
    @test hcs isa Vector{Tuple{typeof(1.0f0u"°C"),typeof(1.0f0u"kcal/mol/K")}}
    @test length(hcs) == n
    hcs = heat_capacity(s, Tmin, Tmax, Tincrement; mpoints)
    @test hcs isa Vector{Tuple{typeof(1.0f0u"°C"),typeof(1.0f0u"kcal/mol/K")}}
    @test length(hcs) == n
    hcs = heat_capacity(s, Tmin, Tmax; mpoints)
    @test hcs isa Vector{Tuple{typeof(1.0f0u"°C"),typeof(1.0f0u"kcal/mol/K")}}
    @test length(hcs) == n
end

@testset "FoldCompound dp matrices" begin
    @testset "mfe matrices" begin
        # test MFE matrices fields from fc.uptr.matrices[]
        seq = "GGGAAACCC"
        n = length(seq)
        fc = FoldCompound(seq; uniq_ML=true, circular=true)
        mfe(fc)
        for sym in (:c, :fML, :fM1)
            mat = getproperty(fc, Symbol("matrices_" * String(sym)))
            @test mat isa Matrix
            @test size(mat) == (n,n)
        end
        for sym in (:fM2, :f5)
            vec = getproperty(fc, Symbol("matrices_" * String(sym)))
            @test vec isa Vector
            @test length(vec) == n
        end
        @test fc.matrices_Fc  isa Int
        @test fc.matrices_FcH isa Int
        @test fc.matrices_FcI isa Int
        @test fc.matrices_FcM isa Int
        # TODO: when does f3 get filled in?
        @test fc.matrices_f3 isa Union{Nothing,Vector}
    end

    @testset "partfn matrices" begin
        # test partition function fields from fc.uptr.exp_matrices[]
        seq = "GGGAAACCC"
        n = length(seq)
        fc = FoldCompound(seq; uniq_ML=true, circular=true)
        partfn(fc)
        for sym in (:q, :qb, :qm, :qm1, :probs)
            mat = getproperty(fc, Symbol("exp_matrices_" * String(sym)))
            @test mat isa Matrix
            @test size(mat) == (n,n)
        end
        for sym in (:qm2, :scale, :expMLbase)
            vec = getproperty(fc, Symbol("exp_matrices_" * String(sym)))
            @test vec isa Vector
            @test length(vec) == n
        end
    end
end
