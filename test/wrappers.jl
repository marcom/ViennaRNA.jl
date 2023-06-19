using Test
using ViennaRNA
using ViennaRNA: LibRNA
using Unitful

@testset "params_{load*,save}" begin
    showtestset()
    for (sym, loadfn, params_name) in [
        (:RNA_Andronescu2007, ViennaRNA.params_load_RNA_Andronescu2007, "RNA - Andronescu 2007"),
        (:RNA_Langdon2018, ViennaRNA.params_load_RNA_Langdon2018, "RNA - Langdon 2018"),
        (:RNA_Turner1999, ViennaRNA.params_load_RNA_Turner1999, "RNA - Turner 1999"),
        (:RNA_Turner2004, ViennaRNA.params_load_RNA_Turner2004, "RNA - Turner 2004"),
        (:RNA_misc_special_hairpins, ViennaRNA.params_load_RNA_misc_special_hairpins, "RNA - Misc. Special Hairpins"),
        (:DNA_Mathews1999, ViennaRNA.params_load_DNA_Mathews1999, "DNA - Mathews 1999"),
        (:DNA_Mathews2004, ViennaRNA.params_load_DNA_Mathews2004, "DNA - Mathews 2004"),
        (:defaults, ViennaRNA.params_load_defaults, "RNA - Turner 2004"),
        ]

        loadfn()
        fc = FoldCompound("GGGAAACCC")
        @test fc.params_name == params_name
        finalize(fc)

        if sym ∉ (:defaults, :RNA_misc_special_hairpins)
            # ViennaRNA.params_load(sym) for sym ∈ [:defaults, :RNA_misc_special_hairpins] is not supported
            ViennaRNA.params_load(sym)
            fc = FoldCompound("GGGAAACCC")
            @test fc.params_name == params_name
            finalize(fc)
        end
    end
    @test_throws ArgumentError ViennaRNA.params_load(:UNKNOWN_PARAMS_NAME)

    # params_load_from_string(str::AbstractString, name::AbstractString)
    paramstr =
        """
        ## RNAfold parameter file v2.0

        # stack
        /*  CG    GC    GU    UG    AU    UA    @  */
          -200  -300  -200   -50  -200  -200     0
          -300  -300  -250   -50  -220  -200  -150
          -200  -250   100   -50     0  -130   100
           -50   -50   -50   -50   -50   -50   -50
          -200  -200     0   -50  -110   -90   -60
          -200  -200  -100   -50   -90  -130   -90
             0  -150   100   -50   -60   -90   100
        """
    ViennaRNA.params_load_from_string(paramstr, "My test params")
    fc = FoldCompound("GGGAAACCC")
    @test fc.params_name == "My test params"
    finalize(fc)

    # params_load(filename::AbstractString)
    paramfile = joinpath(LibRNA.ViennaRNA_jll.artifact_dir, "share", "ViennaRNA", "rna_turner1999.par")
    ViennaRNA.params_load(paramfile)
    fc = FoldCompound("GGGAAACCC")
    @test fc.params_name == "rna_turner1999.par"
    finalize(fc)

    # params_save(filename::AbstractString)
    mktemp() do path, _
        ViennaRNA.params_save(path)
        @test filesize(path) > 0
    end

    # load the default parameters again
    ViennaRNA.params_load_defaults()
end

@testset "FoldCompound" begin
    showtestset()
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
        @test fc.max_bp_span      == length(fc)  # ≠ LibRNA.VRNA_MODEL_DEFAULT_MAX_BP_SPAN == -1
        @test fc.min_loop_size    == LibRNA.TURN
        @test fc.no_GU_basepairs  == LibRNA.VRNA_MODEL_DEFAULT_NO_GU
        @test fc.no_GU_closure    == LibRNA.VRNA_MODEL_DEFAULT_NO_GU_CLOSURE
        @test fc.no_lonely_pairs  == LibRNA.VRNA_MODEL_DEFAULT_NO_LP
        @test fc.nstrands         == 1
        @test fc.params_name      == "RNA - Turner 2004"
        @test fc.special_hairpins == LibRNA.VRNA_MODEL_DEFAULT_SPECIAL_HP
        @test fc.temperature      == 37u"°C"
        @test fc.type             == (s isa String ? :single : :comparative)
        @test fc.uniq_ML          == LibRNA.VRNA_MODEL_DEFAULT_UNIQ_ML
        @test fc.window_size      == length(fc)  # ≠ LibRNA.VRNA_MODEL_DEFAULT_WINDOW_SIZE == -1
        if fc.type === :single
            @test fc.sequence     == s
        elseif fc.type === :comparative
            @test isnothing(fc.sequence)
        end

        # options
        fc = FoldCompound(s; options=[:mfe])
        @test length(fc) == 9

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

@testset "FoldCompound dp matrices" begin
    showtestset()
    @testset "mfe matrices" begin
        showtestset()
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
        showtestset()
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

@testset "Pairtable" begin
    showtestset()
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

@testset "bpp" begin
    showtestset()
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

@testset "bp_distance" begin
    showtestset()
    @test bp_distance("((..))", "(....)") == 1
    @test_throws ArgumentError bp_distance("()", ".")
end

@testset "centroid" begin
    showtestset()
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    fc = FoldCompound(seq)
    partfn(fc)
    @test centroid(fc) isa Tuple{String,AbstractFloat}
    @test centroid(seq) isa Tuple{String,AbstractFloat}
    fc = FoldCompound(seq)
    @test_throws ArgumentError centroid(fc)
end

@testset "energy" begin
    showtestset()
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
    # RNA_Turner1999
    ViennaRNA.params_load_RNA_Turner1999()
    @test energy(FoldCompound(seq), str) ≈ -2.90u"kcal/mol" atol=1e-3u"kcal/mol"
    ViennaRNA.params_load_defaults()
    # RNA_Andronescu2007
    ViennaRNA.params_load_RNA_Andronescu2007()
    @test energy(FoldCompound(seq), str) ≈ -2.15u"kcal/mol" atol=1e-3u"kcal/mol"
    ViennaRNA.params_load_defaults()
    # test different temperatures
    @test energy(FoldCompound(seq; temperature=30u"°C"), str) ≈ -3.52u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; temperature=37u"°C"), str) ≈ -2.90u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; temperature=55u"°C"), str) ≈ -1.24u"kcal/mol" atol=1e-3u"kcal/mol"
    @test energy(FoldCompound(seq; temperature=323.15u"K"), str) ≈ -1.69u"kcal/mol" atol=1e-3u"kcal/mol" # 323.15u"K" is 50u"°C"
end

@testset "ensemble_defect" begin
    showtestset()
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

@testset "heat_capacity" begin
    showtestset()
    s = "GGGAAACCC"
    fc = FoldCompound(s)
    Tmin = 10u"°C"
    Tmax = 60u"°C"
    Tincrement = 1u"°C"
    mpoints = 2
    # TODO: this range call doesn't work without ustrip
    n = length(range(ustrip(Tmin), ustrip(Tmax); step=ustrip(Tincrement)))
    # TODO: use this once using julia-1.8 or above
    #n = length(range(start=ustrip(Tmin), stop=ustrip(Tmax), step=ustrip(Tincrement)))

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

@testset "init_rand_seed" begin
    @test ViennaRNA.init_rand_seed(42) === nothing
end

@testset "inverse_fold" begin
    showtestset()
    s = "((((....))))"
    @test inverse_fold("A"^length(s), s) isa Tuple{String,AbstractFloat}
    @test_throws ArgumentError inverse_fold("A", "()")
    @test inverse_pf_fold("A"^length(s), s) isa Tuple{String,Unitful.Quantity}
    @test_throws ArgumentError inverse_pf_fold("A", "()")
end

@testset "mea" begin
    showtestset()
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

@testset "mean_bp_distance" begin
    showtestset()
    seq = "GGGAAACCC"
    fc = FoldCompound(seq)
    @test_throws ArgumentError mean_bp_distance(fc)
    partfn(fc)
    @test mean_bp_distance(fc) isa AbstractFloat
    @test mean_bp_distance(seq) isa AbstractFloat
end

@testset "mfe" begin
    showtestset()
    seq = "GGGAAACCCC"
    fc = FoldCompound(seq)
    @test mfe(fc) isa Tuple{String,Unitful.Quantity}
    @test mfe(seq) isa Tuple{String,Unitful.Quantity}
    @test_throws ArgumentError mfe(FoldCompound(seq; options=[:window]))
end

@testset "mfe_window" begin
    showtestset()
    Tres = Vector{ViennaRNA.ResultWindowMFE}
    seq = "G"^30 * "A"^3 * "C"^30
    window_size = 30
    fc = FoldCompound(seq; window_size, options=[:mfe, :window])
    res1 = mfe_window(fc)
    @test res1 isa Tres
    @test length(res1) > 0
    res2 = mfe_window(seq; window_size)
    @test res2 isa Tres
    @test res1 == res2
end

@testset "mfe_window_channel" begin
    showtestset()
    Tres = Vector{ViennaRNA.ResultWindowMFE}
    seq = "G"^30 * "A"^3 * "C"^30
    window_size = 30
    fc = FoldCompound(seq; window_size, options=[:mfe, :window])
    chan = mfe_window_channel(fc)
    res1 = collect(chan)
    @test res1 isa Tres
    @test length(res1) > 0
    chan = mfe_window_channel(seq; window_size)
    res2 = collect(chan)
    @test res2 isa Tres
    @test res1 == res2
end

@testset "neighbors" begin
    showtestset()
    fc = FoldCompound("GGGAAACCC")
    pt = Pairtable(".((...)).")
    @test neighbors(fc, pt) == [[(-2, -8)], [(-3, -7)], [(1, 9)]]
end

@testset "partfn" begin
    showtestset()
    seq = "GGGAAACCCC"
    fc = FoldCompound(seq)
    @test partfn(fc) isa Tuple{String,Unitful.Quantity}
    @test partfn(seq) isa Tuple{String,Unitful.Quantity}
end

@testset "plot_coords" begin
    showtestset()
    function test_plot(s::AbstractString, plot_type::Union{Nothing,Symbol}=nothing)
        x, y = isnothing(plot_type) ? plot_coords(s) : plot_coords(s; plot_type)
        @test x isa Vector{Float32}
        @test length(x) == length(s)
        @test y isa Vector{Float32}
        @test length(y) == length(s)
    end
    function test_plot(pt::Pairtable, plot_type::Union{Nothing,Symbol}=nothing)
        x, y = isnothing(plot_type) ? plot_coords(pt) : plot_coords(pt; plot_type)
        @test x isa Vector{Float32}
        @test length(x) == length(pt)
        @test y isa Vector{Float32}
        @test length(y) == length(pt)
    end

    # test: plot_coords(::String), plot_coords(::Pairtable)
    for s in [".", "(((...)))", "(((..)).((..).))(...).."]
        pt = Pairtable(s)
        test_plot(s)
        test_plot(pt)
        for plot_type in keys(ViennaRNA.Private.PLOT_COORDS_PLOT_TYPE)
            test_plot(s, plot_type)
            test_plot(pt, plot_type)
        end
    end
    @test_throws ArgumentError x, y = plot_coords("(...)"; plot_type = :unknown_plot_type)
    @test_throws ArgumentError x, y = plot_coords(Pairtable("(...)"); plot_type = :unknown_plot_type)

    # zero-sized inputs
    let s = "", pt = Pairtable(s)
        for plot_type in keys(ViennaRNA.Private.PLOT_COORDS_PLOT_TYPE)
            test_plot(s, plot_type)
            test_plot(pt, plot_type)
        end
    end

    # test zero-length hairpins
    for s in ["()", "(())", "()()", "(()())"]
        pt = Pairtable(s)
        for plot_type in keys(ViennaRNA.Private.PLOT_COORDS_PLOT_TYPE)
            if plot_type ∉ (:default, :turtle, :puzzler)
                test_plot(s, plot_type)
                test_plot(pt, plot_type)
            else
                @test_throws ArgumentError x, y = plot_coords(s; plot_type)
                @test_throws ArgumentError x, y = plot_coords(pt; plot_type)
            end
        end
    end
end

@testset "prob_of_structure" begin
    showtestset()
    seq = "GGGAAACCCC"
    str = "(((....)))"
    fc = FoldCompound(seq)
    partfn(fc)
    @test prob_of_structure(fc, str) isa AbstractFloat
    @test_throws ArgumentError prob_of_structure(fc, ".")
    @test prob_of_structure(seq, str) isa AbstractFloat
end

@testset "sample_structures" begin
    showtestset()
    seq = "GGGGGAAAAACCCCCCCCAUUCA"
    n = length(seq)
    fc = FoldCompound(seq; uniq_ML=true)
    partfn(fc)

    for input in [seq, fc]
        s = sample_structures(input)
        @test length(s) == 10  # default num_samples
        @test all(x -> length(x) == n, s)
        for num_samples in [1, 15]
            for options in keys(ViennaRNA.Private.SAMPLE_STRUCTURES_OPTIONS)
                s = sample_structures(input; num_samples, options)
                @test length(s) == num_samples
                @test all(x -> length(x) == n, s)
            end
        end
        @test_throws ArgumentError sample_structures(input; num_samples=-1)
        @test_throws ArgumentError sample_structures(input; options=:unknown_option)
    end

end

@testset "subopt" begin
    showtestset()
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
    showtestset()
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

@testset "tree_edit_dist" begin
    showtestset()
    @test tree_edit_dist("...", "()") == 3.0f0
end
