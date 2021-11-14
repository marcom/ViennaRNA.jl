module ViennaRNA

module LibRNA
    import ViennaRNA_jll
    using CBinding

    let
        incdir = joinpath(ViennaRNA_jll.artifact_dir, "include")
        libdir = dirname(ViennaRNA_jll.libRNA_path)
        c`-std=c11 -fparse-all-comments -I$(incdir) -L$(libdir) -lRNA
          -DVRNA_DISABLE_C11_FEATURES`
        #c`-fparse-all-comments -I$(incdir) -I$(incdir)/ViennaRNA
        #  -DVRNA_DISABLE_C11_FEATURES -pthread -fno-lto -Wl,-fno-lto
        #  -L${libdir} -lRNA -fopenmp -lgsl -lgslcblas -lpthread -lmpfr -lgmp
        #  -lstdc++`
    end

    const c"FILE" = Cvoid
    const c"va_list" = Cvoid
    const c"size_t" = Csize_t
    c"""
        #include "ViennaRNA/model.h"
        #include "ViennaRNA/sequence.h"
        #include "ViennaRNA/eval.h"
        #include "ViennaRNA/mfe.h"
        #include "ViennaRNA/part_func.h"
        #include "ViennaRNA/subopt.h"
        #include "ViennaRNA/inverse.h"
        #include "ViennaRNA/boltzmann_sampling.h"
        #include "ViennaRNA/MEA.h"
        #include "ViennaRNA/treedist.h"
        #include "ViennaRNA/utils/basic.h"
        #include "ViennaRNA/utils/structures.h"
        #include "ViennaRNA/params/basic.h"
        #include "ViennaRNA/params/io.h"
        #include "ViennaRNA/plotting/layouts.h"
    """ij

    # TODO: This workaround is needed because otherwise we get an
    #       error ('ERROR: LoadError: UndefVarError: ViennaRNA not
    #       defined') from inside a CBinding function. Might be a
    #       parsing problem of the macro because it's not a simple
    #       macro but an expression.
    macro VRNA_BRACKETS_ANY()
        :(
            (LibRNA.@VRNA_BRACKETS_RND) |
            (LibRNA.@VRNA_BRACKETS_SQR) |
            (LibRNA.@VRNA_BRACKETS_CLY) |
            (LibRNA.@VRNA_BRACKETS_ANG) |
            (LibRNA.@VRNA_BRACKETS_ALPHA)
        )
    end
end # module LibRNA


using CBinding
using .LibRNA
import Base
using Unitful: @u_str, Quantity, uconvert
export FoldCompound

const en_unit = 1.0u"kcal/mol"
const en_int_unit = 0.01u"kcal/mol"

# data types

mutable struct FoldCompound
    seq :: String
    ptr :: Cptr{LibRNA.vrna_fc_s}
    function FoldCompound(seq::AbstractString,
                          model_details::Cptr{LibRNA.vrna_md_s},
                          options::Unsigned)
        ptr = LibRNA.vrna_fold_compound(seq, model_details, options)
        ptr != C_NULL || throw(ErrorException("pointer == C_NULL"))
        fc = new(seq, ptr)
        finalizer(fc) do x
            # TODO: do we have to call vrna_mx_mfe_free or
            #       vrna_mx_pf_free here ourselves?
            LibRNA.vrna_fold_compound_free(x.ptr)
        end
    end
end

function FoldCompound(seq::AbstractString;
                      uniq_ML::Integer=0,
                      options::Unsigned=LibRNA.@VRNA_OPTION_DEFAULT)
    # TODO: who frees md? maybe vrna_fold_compound_free()
    #       possible memory leak
    md = Cptr{LibRNA.vrna_md_s}(C_NULL)
    if uniq_ML != 0
        # TODO: we assume default vrna_md_s has uniq_ML = 0
        md = Libc.malloc(LibRNA.vrna_md_s)
        LibRNA.vrna_md_set_default(md)
        md.uniq_ML[] = uniq_ML
    end
    return FoldCompound(seq, md, options)
end

Base.length(fc::FoldCompound) = length(fc.seq)

has_exp_matrices(fc::FoldCompound) = fc.ptr.exp_matrices[] != C_NULL

mutable struct Pairtable
    # Note: ptr[1] contains the number of elements, ptr[i+1] is the
    #       i-th element this means that functions like getindex or
    #       setindex! must use ptr[i+1] to access the i-th element
    ptr :: Cptr{Cshort}
    function Pairtable(structure::AbstractString)
        pt = new(LibRNA.vrna_ptable_from_string(structure,
                                                LibRNA.@VRNA_BRACKETS_ANY))
        finalizer(pt) do x
            Libc.free(x.ptr)
        end
    end
end

Base.length(pt::Pairtable) = pt.ptr[1]

function Base.show(io::IO, mime::MIME"text/plain", pt::Pairtable)
    show(io, mime, [pt[i] for i = 1:length(pt)])
end

Base.getindex(pt::Pairtable, i) = pt.ptr[i+1]

function Base.setindex!(pt::Pairtable, val, i)
    if pt.ptr[i+1] != 0
        pt.ptr[pt.ptr[i+1] + 1] = 0
    end
    if pt.ptr[val+1] != 0
        pt.ptr[pt.ptr[val+1] + 1] = 0
    end
    pt.ptr[i+1] = val
    pt.ptr[val+1] = i
end


# secondary structure distance measures

function bp_distance(a::AbstractString, b::AbstractString)
    length(a) == length(b) ||
        throw(ArgumentError("input structures have different lengths"))
    return LibRNA.vrna_bp_distance(a, b)
end

function tree_edit_dist(structure1::AbstractString, structure2::AbstractString;
                        hit_type::Integer = @LibRNA.VRNA_STRUCTURE_TREE_HIT)
    # vrna_db_to_tree_string prints result to stdout, which we discard
    ts1, ts2 = redirect_stdout(devnull) do
        return LibRNA.vrna_db_to_tree_string(structure1, hit_type),
               LibRNA.vrna_db_to_tree_string(structure2, hit_type)
    end
    t1 = LibRNA.make_tree(ts1)
    t2 = LibRNA.make_tree(ts2)
    d = LibRNA.tree_edit_distance(t1, t2)
    LibRNA.free_tree(t1)
    LibRNA.free_tree(t2)
    Libc.free(ts1)
    Libc.free(ts2)
    return d
end


# equilibrium properties

# mean basepair distance of ensemble to itself

function mean_bp_distance(fc::FoldCompound)
    has_exp_matrices(fc) ||
        throw(ArgumentError("partfn(fc::FoldCompound) must be run first"))
    return LibRNA.vrna_mean_bp_distance(fc.ptr)
end

function mean_bp_distance(sequence::AbstractString)
    fc = FoldCompound(sequence)
    partfn(fc)
    return mean_bp_distance(fc)
end

# ensemble defect: mean distance of ensemble to a target structure
# TODO: check if this is exactly expected basepair distance to target struct

function ensemble_defect(fc::FoldCompound, pt::Pairtable)
    has_exp_matrices(fc) ||
        throw(ArgumentError("partfn(::FoldCompound) must be run first"))
    length(fc) == length(pt) ||
        throw(ArgumentError("FoldCompound and Pairtable must have equal length"))
    return LibRNA.vrna_ensemble_defect_pt(fc.ptr, pt.ptr)
end

function ensemble_defect(sequence::AbstractString, pt::Pairtable)
    fc = FoldCompound(sequence)
    partfn(fc)
    return ensemble_defect(fc, pt)
end

ensemble_defect(fc::FoldCompound, structure::AbstractString) =
    ensemble_defect(fc, Pairtable(structure))

function ensemble_defect(sequence::AbstractString, structure::AbstractString)
    fc = FoldCompound(sequence)
    partfn(fc)
    return ensemble_defect(fc, Pairtable(structure))
end

# probability of a structure in the ensemble

function prob_of_structure(fc::FoldCompound, structure::AbstractString)
    has_exp_matrices(fc) ||
        throw(ArgumentError("must first call ViennaRNA.partfn(::FoldCompound)"))
    length(fc) == length(structure) ||
        throw(ArgumentError("FoldCompound and structure must have equal length"))
    return LibRNA.vrna_pr_structure(fc.ptr, structure)
end

function prob_of_structure(sequence::AbstractString, structure::AbstractString)
    fc = FoldCompound(sequence)
    partfn(fc)
    return prob_of_structure(fc, structure)
end

# free energy change of folding

function energy(fc::FoldCompound,
                structure::AbstractString;
                verbose::Bool=false,
                verbosity_level::Integer=LibRNA.@VRNA_VERBOSITY_DEFAULT)
    length(fc) == length(structure) ||
        throw(ArgumentError("sequence and structure must have equal length"))
    if verbose
        en = LibRNA.vrna_eval_structure_v(fc.ptr, structure,
                                          verbosity_level, C_NULL)
    else
        en = LibRNA.vrna_eval_structure(fc.ptr, structure)
    end
    return en * en_unit
end

function energy(sequence::AbstractString,
                structure::AbstractString;
                verbose::Bool=false,
                verbosity_level::Integer=LibRNA.@VRNA_VERBOSITY_DEFAULT)
    fc = FoldCompound(sequence; options=LibRNA.@VRNA_OPTION_EVAL_ONLY)
    return energy(fc, structure; verbose, verbosity_level)
end

# minimum free energy (mfe) structure

function mfe(fc::FoldCompound)
    seqlen = length(fc)
    cstr_structure = Ptr{Cchar}(LibRNA.vrna_alloc(seqlen + 1))
    en_mfe = LibRNA.vrna_mfe(fc.ptr, cstr_structure)
    structure = unsafe_string(cstr_structure)
    Libc.free(cstr_structure)
    return structure, en_mfe * en_unit
end

mfe(sequence::AbstractString) =
    mfe(FoldCompound(sequence; options=LibRNA.@VRNA_OPTION_MFE))

# partition function

function partfn(fc::FoldCompound)
    seqlen = length(fc)
    cstr_structure = Ptr{Cchar}(LibRNA.vrna_alloc(seqlen + 1))
    RTlogZ = LibRNA.vrna_pf(fc.ptr, cstr_structure)
    structure = unsafe_string(cstr_structure)
    Libc.free(cstr_structure)
    return structure, RTlogZ * en_unit
end

partfn(sequence::AbstractString) = partfn(FoldCompound(sequence))

# basepair probabilities

function bpp(fc::FoldCompound)
    # TODO: return upper triangular matrix type
    has_exp_matrices(fc) ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    n = length(fc)
    p = zeros(n,n)
    index = fc.ptr.iindx[]
    probs = fc.ptr.exp_matrices.probs[]
    for i = 1:n-1
        for j = i+1:n
            p[i,j] = probs[index[i + 1] - j + 1]
        end
    end
    return p
end

function bpp(sequence::AbstractString)
    fc = FoldCompound(sequence)
    partfn(fc)
    return bpp(fc)
end

# stochastic backtrack

function pbacktrack(fc::FoldCompound;
                    num_samples::Integer=1,
                    options::Integer=LibRNA.@VRNA_PBACKTRACK_DEFAULT)
    has_exp_matrices(fc) ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    fc.ptr.params.model_details.uniq_ML[] == 1 ||
        throw(ArgumentError("must have fc.params.model_details.uniq_ML == 1"))
    s = LibRNA.vrna_pbacktrack_num(fc.ptr, num_samples, options)
    samples = String[]
    i = 1
    while s[i] != C_NULL
        push!(samples, unsafe_string(s[i]))
        i += 1
    end
    return samples
end

function pbacktrack(sequence::AbstractString; num_samples::Integer=1, options::Integer=LibRNA.@VRNA_PBACKTRACK_DEFAULT)
    fc = FoldCompound(sequence; uniq_ML=1)
    partfn(fc)
    return pbacktrack(fc; num_samples, options)
end

# maximum expected accuracy (MEA) structure

function mea(fc::FoldCompound; gamma=1.0)
    has_exp_matrices(fc) ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    ptr_mea   = Libc.malloc(Cfloat)
    ptr_str   = LibRNA.vrna_MEA(fc.ptr, gamma, ptr_mea)
    mea_val   = ptr_mea[]
    mea_struc = unsafe_string(ptr_str)
    Libc.free(ptr_mea)
    Libc.free(ptr_str)
    return mea_struc, mea_val
end

function mea(sequence::AbstractString; gamma=1.0)
    fc = FoldCompound(sequence)
    partfn(fc)
    return mea(fc; gamma)
end

# centroid structure of ensemble

function centroid(fc::FoldCompound)
    has_exp_matrices(fc) ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    ptr_dist = Libc.malloc(Cdouble)
    ptr_str = LibRNA.vrna_centroid(fc.ptr, ptr_dist)
    cen_dist = ptr_dist[]
    cen_struc = unsafe_string(ptr_str)
    Libc.free(ptr_dist)
    Libc.free(ptr_str)
    return cen_struc, cen_dist
end

function centroid(sequence::AbstractString)
    fc = FoldCompound(sequence)
    partfn(fc)
    return centroid(fc)
end


# suboptimal structures

function subopt(fc::FoldCompound; delta::Quantity, sorted::Bool=true)
    fc.ptr.params.model_details.uniq_ML[] == 1 ||
        throw(ArgumentError("must have fc.params.model_details.uniq_ML == 1"))
    delta_nounit = uconvert(u"kcal/mol", delta) / en_int_unit
    delta_int = round(Int, delta_nounit)
    # TODO: warn if this float -> int conversion leads to loss of
    # precision, e.g. in the case that delta == 0.0001 kcal/mol where
    # we would have delta_int == 0
    subs = Tuple{String,Quantity}[]
    s = LibRNA.vrna_subopt(fc.ptr, delta_int, sorted, C_NULL)
    i = 1
    while s[i].structure != C_NULL
        push!(subs, (unsafe_string(s[i].structure), s[i].energy * en_unit))
        i += 1
    end
    return subs
end

function subopt(sequence::AbstractString; delta::Quantity, sorted::Bool=true)
    fc = FoldCompound(sequence; uniq_ML=1, options=LibRNA.@VRNA_OPTION_MFE)
    return subopt(fc; delta, sorted)
end

function subopt_zuker(fc::FoldCompound)
    subs = Tuple{String,Quantity}[]
    s = LibRNA.vrna_subopt_zuker(fc.ptr)
    i = 1
    while s[i].structure != C_NULL
        push!(subs, (unsafe_string(s[i].structure), s[i].energy * en_unit))
        i += 1
    end
    return subs
end

function subopt_zuker(sequence::AbstractString)
    fc = FoldCompound(sequence; options=LibRNA.@VRNA_OPTION_MFE)
    return subopt_zuker(fc)
end

# inverse folding / sequence design

function inverse_fold(start::AbstractString, target::AbstractString)
    length(start) == length(target) ||
        throw(ArgumentError("start and target must have same size"))
    # TODO: there is no copy(::String), so this is a substitute
    #       needed because inverse_fold mutates its first argument
    design = join(collect(start))
    dist = LibRNA.inverse_fold(design, target)
    return design, dist
end

function inverse_pf_fold(start::AbstractString, target::AbstractString)
    length(start) == length(target) ||
        throw(ArgumentError("start and target must have same size"))
    # TODO: there is no copy(::String), so this is a substitute
    #       needed because inverse_pf_fold mutates its first argument
    design = join(collect(start))
    RT_log_p = LibRNA.inverse_pf_fold(design, target)
    return design, RT_log_p * en_unit
end


# plotting secondary structures

function plot_coords(structure::AbstractString;
                     plot_type::Symbol=:naview)
    if plot_type === :simple
        type = LibRNA.@VRNA_PLOT_TYPE_SIMPLE
    elseif plot_type === :naview
        type = LibRNA.@VRNA_PLOT_TYPE_NAVIEW
    elseif plot_type === :circular
        type = LibRNA.@VRNA_PLOT_TYPE_CIRCULAR
    elseif plot_type === :turtle
        type = LibRNA.@VRNA_PLOT_TYPE_TURTLE
    elseif plot_type === :puzzler
        type = LibRNA.@VRNA_PLOT_TYPE_PUZZLER
    else
        throw(ArgumentError("unknown plot_type $plot_type, options are: :simple, :naview, :circular, :turtle, :puzzler"))
    end
    cx = Libc.malloc(Cptr{Cfloat})
    cy = Libc.malloc(Cptr{Cfloat})
    n = LibRNA.vrna_plot_coords(structure, cx, cy, type)
    x = [cx[][i] for i = 1:n]
    y = [cy[][i] for i = 1:n]
    Libc.free(cx[])
    Libc.free(cy[])
    Libc.free(cx)
    Libc.free(cy)
    return x, y
end

end # module ViennaRNA
