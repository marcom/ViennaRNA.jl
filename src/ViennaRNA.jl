# TODO
# - replace Libc.malloc with LibRNA.vrna_alloc
# - has_exp_matrices checks
# - uniq_ML option to FoldCompound()
# - bpp()
# - pbacktrack needs uniq_ML != 0

# - change energy param set

# - change model details
#   vrna_md_set_default(*md) -> void
#   vrna_md_copy([*to_md], [*from_md]) -> *md
#   vrna_md_defaults_reset([*md]) -> void
#
#   md_new = vrna_md_copy(NULL, NULL)
#   

module ViennaRNA

import Base
using Unitful: @u_str, Quantity, uconvert
export FoldCompound, Pairtable, bp_distance, energy, mfe, partfn, subopt,
    subopt_zuker, inverse_fold, inverse_pf_fold, neighbors, plot_coords

const en_unit = 1.0u"kcal/mol"
const en_int_unit = 0.01u"kcal/mol"


module LibRNA

import ViennaRNA_jll
const LIB = ViennaRNA_jll.libRNA

# data types
const Cfile      = Cvoid
const vrna_fc_s  = Cvoid
const vrna_md_s  = Cvoid
const Tree       = Cvoid
const FLT_OR_DBL = Cfloat
struct vrna_subopt_sol_s
    energy::Cfloat
    structure::Cstring
end
Base.getindex(p::Ptr{vrna_subopt_sol_s}, i) = unsafe_load(p, i)
struct vrna_move_s
    pos_5::Cint
    pos_3::Cint
    next::Ptr{Cvoid}    # TODO: actually Ptr{vrna_move_s}
end
Base.getindex(p::Ptr{vrna_move_s}, i) = unsafe_load(p, i)

# constants
const VRNA_OPTION_DEFAULT     = 0x0
const VRNA_OPTION_MFE         = 0x1
const VRNA_OPTION_PF          = 0x2
const VRNA_OPTION_HYBRID      = 0x4
const VRNA_OPTION_EVAL_ONLY   = 0x8
const VRNA_OPTION_WINDOW      = 0x10
const VRNA_BRACKETS_ALPHA     = 0x4
const VRNA_BRACKETS_RND       = 0x8
const VRNA_BRACKETS_CLY       = 0x10
const VRNA_BRACKETS_ANG       = 0x20
const VRNA_BRACKETS_SQR       = 0x40
const VRNA_BRACKETS_DEFAULT   = VRNA_BRACKETS_RND | VRNA_BRACKETS_CLY | VRNA_BRACKETS_ANG | VRNA_BRACKETS_SQR
const VRNA_BRACKETS_ANY       = VRNA_BRACKETS_RND | VRNA_BRACKETS_CLY | VRNA_BRACKETS_ANG | VRNA_BRACKETS_SQR | VRNA_BRACKETS_ALPHA
const VRNA_VERBOSITY_QUIET    = -1
const VRNA_VERBOSITY_DEFAULT  = 1
const VRNA_PLOT_TYPE_SIMPLE   = 0
const VRNA_PLOT_TYPE_NAVIEW   = 1
const VRNA_PLOT_TYPE_CIRCULAR = 2
const VRNA_PLOT_TYPE_TURTLE   = 3
const VRNA_PLOT_TYPE_PUZZLER  = 4
const VRNA_MOVESET_INSERTION  = 4
const VRNA_MOVESET_DELETION   = 8
const VRNA_MOVESET_SHIFT      = 16
const VRNA_MOVESET_NO_LP      = 32
const VRNA_MOVESET_DEFAULT    = VRNA_MOVESET_INSERTION | VRNA_MOVESET_DELETION
const VRNA_STRUCTURE_TREE_HIT            = 0x1
const VRNA_STRUCTURE_TREE_SHAPIRO_SHORT  = 0x2
const VRNA_STRUCTURE_TREE_SHAPIRO        = 0x3
const VRNA_STRUCTURE_TREE_SHAPIRO_EXT    = 0x4
const VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT = 0x5
const VRNA_STRUCTURE_TREE_EXPANDED       = 0x6
const VRNA_PBACKTRACK_DEFAULT            = 0
const VRNA_PBACKTRACK_NON_REDUNDANT      = 1

# functions

vrna_alloc(s) = @ccall LIB.vrna_alloc(s::Csize_t)::Ptr{Cvoid}

vrna_fold_compound(sequence, md_p, options) =
    @ccall LIB.vrna_fold_compound(sequence::Cstring,
                                  md_p::Ptr{vrna_md_s},
                                  options::Cuint)::Ptr{vrna_fc_s}
vrna_fold_compound_free(fc) =
    @ccall LIB.vrna_fold_compound_free(fc::Ptr{vrna_fc_s})::Cvoid

vrna_ptable_from_string(string, options) =
    @ccall LIB.vrna_ptable_from_string(string::Cstring, options::Cuint)::Ptr{Cshort}

vrna_bp_distance(str1, str2) =
    @ccall LIB.vrna_bp_distance(str1::Cstring, str2::Cstring)::Cint

vrna_db_to_tree_string(structure, type) =
    @ccall LIB.vrna_db_to_tree_string(structure::Cstring, type::Cuint)::Cstring

make_tree(struc) =
    @ccall LIB.make_tree(struc::Cstring)::Ptr{Tree}

tree_edit_distance(t1, t2) =
    @ccall LIB.tree_edit_distance(t1::Ptr{Tree}, t2::Ptr{Tree})::Cfloat

free_tree(t) =
    @ccall LIB.free_tree(t::Ptr{Tree})::Cvoid

vrna_mean_bp_distance(fc) =
    @ccall LIB.vrna_mean_bp_distance(fc::Ptr{vrna_fc_s})::Cdouble

vrna_ensemble_defect_pt(fc, pt) =
    @ccall LIB.vrna_ensemble_defect_pt(fc::Ptr{vrna_fc_s}, pt::Ptr{Cshort})::Cdouble

vrna_pr_structure(fc, structure) =
    @ccall LIB.vrna_pr_structure(fc::Ptr{vrna_fc_s}, structure::Cstring)::Cdouble

vrna_eval_structure(fc, structure) =
    @ccall LIB.vrna_eval_structure(fc::Ptr{vrna_fc_s},
                                      structure::Cstring)::Cfloat

vrna_eval_structure_v(fc, structure, verbosity_level, file) =
    @ccall LIB.vrna_eval_structure_v(fc::Ptr{vrna_fc_s},
                                     structure::Cstring,
                                     verbosity_level::Cint,
                                     file::Ptr{Cfile})::Cfloat

vrna_mfe(fc, structure) =
    @ccall LIB.vrna_mfe(fc::Ptr{vrna_fc_s}, structure::Cstring)::Cfloat

vrna_pf(fc, structure) =
    @ccall LIB.vrna_pf(fc::Ptr{vrna_fc_s}, structure::Cstring)::FLT_OR_DBL

vrna_pbacktrack_num(fc, num_samples, options) =
    @ccall LIB.vrna_pbacktrack_num(fc::Ptr{vrna_fc_s}, num_samples::Cuint, options::Cuint)::Ptr{Cstring}

vrna_MEA(fc, gamma, mea) =
    @ccall LIB.vrna_MEA(fc::Ptr{vrna_fc_s}, gamma::Cdouble, mea::Ptr{Cfloat})::Cstring

vrna_centroid(fc, dist) =
    @ccall LIB.vrna_centroid(fc::Ptr{vrna_fc_s}, dist::Ptr{Cdouble})::Cstring

vrna_subopt(fc, delta, sorted, fp) =
    @ccall LIB.vrna_subopt(fc::Ptr{vrna_fc_s}, delta::Cint, sorted::Cint, fp::Ptr{Cfile})::Ptr{vrna_subopt_sol_s}

vrna_subopt_zuker(fc) =
    @ccall LIB.vrna_subopt_zuker(fc::Ptr{vrna_fc_s})::Ptr{vrna_subopt_sol_s}

inverse_fold(start, target) =
    @ccall LIB.inverse_fold(start::Cstring, target::Cstring)::Cfloat

inverse_pf_fold(start, target) =
    @ccall LIB.inverse_pf_fold(start::Cstring, target::Cstring)::Cfloat

vrna_neighbors(fc, pt, options) =
    @ccall LIB.vrna_neighbors(fc::Ptr{vrna_fc_s}, pt::Ptr{Cshort}, options::Cuint)::Ptr{vrna_move_s}

vrna_move_list_free(moves) =
    @ccall LIB.vrna_move_list_free(moves::Ptr{vrna_move_s})::Cvoid

vrna_plot_coords(structure, x, y, plot_type) =
    @ccall LIB.vrna_plot_coords(structure::Cstring, x::Ptr{Ptr{Cfloat}}, y::Ptr{Ptr{Cfloat}}, plot_type::Cint)::Cint

vrna_plot_coords_pt(pt, x, y, plot_type) =
    @ccall LIB.vrna_plot_coords_pt(pt::Ptr{Cshort}, x::Ptr{Ptr{Cfloat}}, y::Ptr{Ptr{Cfloat}}, plot_type::Cint)::Cint

end # module LibRNA


mutable struct FoldCompound
    seq :: String
    ptr :: Ptr{LibRNA.vrna_fc_s}
    function FoldCompound(seq::AbstractString,
                          model_details::Ptr{LibRNA.vrna_md_s},
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
                      options::Unsigned=LibRNA.VRNA_OPTION_DEFAULT)
    # TODO: who frees md? maybe vrna_fold_compound_free()
    #       possible memory leak
    #md = Ptr{vrna_md_s}(C_NULL)
    md = C_NULL
    # if uniq_ML != 0
    #     # TODO: we assume default vrna_md_s has uniq_ML = 0
    #     error("not implemented")
    #     #md = Libc.malloc(LibRNA.vrna_md_s)
    #     #LibRNA.vrna_md_set_default(md)
    #     #md.uniq_ML[] = uniq_ML
    # end
    return FoldCompound(seq, md, options)
end

Base.length(fc::FoldCompound) = length(fc.seq)

#has_exp_matrices(fc::FoldCompound) = fc.ptr.exp_matrices[] != C_NULL
# TODO
has_exp_matrices(fc::FoldCompound) = true

mutable struct Pairtable
    # Note: ptr[1] contains the number of elements, ptr[i+1] is the
    #       i-th element this means that functions like getindex or
    #       setindex! must use ptr[i+1] to access the i-th element
    ptr :: Ptr{Cshort}
    function Pairtable(structure::AbstractString)
        ptr = LibRNA.vrna_ptable_from_string(structure, LibRNA.VRNA_BRACKETS_ANY)
        pt = new(ptr)
        finalizer(pt) do x
            Libc.free(x.ptr)
        end
    end
end

Base.length(pt::Pairtable) = Int(unsafe_load(pt.ptr, 1))

function Base.show(io::IO, mime::MIME"text/plain", pt::Pairtable)
    show(io, mime, [pt[i] for i = 1:length(pt)])
end

function Base.getindex(pt::Pairtable, i)
    @boundscheck if i < 1 || i > length(pt)
        throw(BoundsError(pt, i))
    end
    unsafe_load(pt.ptr, i + 1)
end

function Base.setindex!(pt::Pairtable, val, i)
    @boundscheck begin
        len = length(pt)
        if i < 1 || val < 0 || i > len || val > len
            throw(ArgumentError("illegal basepair set: pt[$i] = $val"))
        end
    end
    if pt[i] != 0
        unsafe_store!(pt.ptr, 0, pt[i] + 1)
    end
    if val != 0
        if pt[val] != 0
            unsafe_store!(pt.ptr, 0, pt[val] + 1)
        end
        unsafe_store!(pt.ptr, i, val + 1)
    end
    unsafe_store!(pt.ptr, val, i + 1)
end


# secondary structure distance measures

function bp_distance(a::AbstractString, b::AbstractString)
    length(a) == length(b) ||
        throw(ArgumentError("input structures have different lengths"))
    return LibRNA.vrna_bp_distance(a, b)
end

function tree_edit_dist(structure1::AbstractString, structure2::AbstractString;
                        hit_type::Integer = LibRNA.VRNA_STRUCTURE_TREE_HIT)
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
                verbosity_level::Integer=LibRNA.VRNA_VERBOSITY_DEFAULT)
    length(fc) == length(structure) ||
        throw(ArgumentError("sequence and structure must have equal length"))
    if verbose
        en = LibRNA.vrna_eval_structure_v(fc.ptr,
                                          structure,
                                          verbosity_level,
                                          C_NULL)
    else
        en = LibRNA.vrna_eval_structure(fc.ptr, structure)
    end
    return en * en_unit
end

function energy(sequence::AbstractString,
                structure::AbstractString;
                verbose::Bool=false,
                verbosity_level::Integer=LibRNA.VRNA_VERBOSITY_DEFAULT)
    fc = FoldCompound(sequence; options=LibRNA.VRNA_OPTION_EVAL_ONLY)
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
    mfe(FoldCompound(sequence; options=LibRNA.VRNA_OPTION_MFE))

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

# function bpp(fc::FoldCompound)
#     # TODO: return upper triangular matrix type
#     has_exp_matrices(fc) ||
#         throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
#     n = length(fc)
#     p = zeros(n,n)
#     index = fc.ptr.iindx[]
#     probs = fc.ptr.exp_matrices.probs[]
#     for i = 1:n-1
#         for j = i+1:n
#             p[i,j] = probs[index[i + 1] - j + 1]
#         end
#     end
#     return p
# end

# function bpp(sequence::AbstractString)
#     fc = FoldCompound(sequence)
#     partfn(fc)
#     return bpp(fc)
# end

# stochastic backtrack

# function pbacktrack(fc::FoldCompound;
#                     num_samples::Integer=1,
#                     options::Integer=LibRNA.VRNA_PBACKTRACK_DEFAULT)
#     has_exp_matrices(fc) ||
#         throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
#     # fc.ptr.params.model_details.uniq_ML[] == 1 ||
#     #     throw(ArgumentError("must have fc.params.model_details.uniq_ML == 1"))
#     s = LibRNA.vrna_pbacktrack_num(fc.ptr, num_samples, options)
#     samples = String[]
#     s == C_NULL && return samples
#     i = 1
#     while unsafe_load(s, i) != C_NULL
#         push!(samples, unsafe_string(unsafe_load(s, i)))
#         i += 1
#     end
#     return samples
# end

# function pbacktrack(sequence::AbstractString; num_samples::Integer=1,
#                     options::Integer=LibRNA.VRNA_PBACKTRACK_DEFAULT)
#     fc = FoldCompound(sequence; uniq_ML=1)
#     partfn(fc)
#     return pbacktrack(fc; num_samples, options)
# end

# maximum expected accuracy (MEA) structure

function mea(fc::FoldCompound; gamma=1.0)
    has_exp_matrices(fc) ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    ptr_mea   = Ptr{Cfloat}(Libc.malloc(sizeof(Cfloat)))
    ptr_str   = LibRNA.vrna_MEA(fc.ptr, gamma, ptr_mea)
    mea_val   = unsafe_load(ptr_mea)
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
    ptr_dist = Ptr{Cdouble}(Libc.malloc(sizeof(Cdouble)))
    ptr_str = LibRNA.vrna_centroid(fc.ptr, ptr_dist)
    cen_dist = unsafe_load(ptr_dist)
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
#    fc.ptr.params.model_details.uniq_ML[] == 1 ||
#        throw(ArgumentError("must have fc.params.model_details.uniq_ML == 1"))
    delta_nounit = uconvert(u"kcal/mol", delta) / en_int_unit
    delta_int = round(Int, delta_nounit)
    # TODO: warn if this float -> int conversion leads to loss of
    # precision, e.g. in the case that delta == 0.0001 kcal/mol where
    # we would have delta_int == 0
    subs = Tuple{String,Quantity}[]
    s = LibRNA.vrna_subopt(fc.ptr, delta_int, sorted, C_NULL)
    s == C_NULL && return subs
    i = 1
    while s[i].structure != C_NULL
        push!(subs, (unsafe_string(s[i].structure), s[i].energy * en_unit))
        i += 1
    end
    # TODO: free s recursively
    return subs
end

function subopt(sequence::AbstractString; delta::Quantity, sorted::Bool=true)
    fc = FoldCompound(sequence; uniq_ML=1, options=LibRNA.VRNA_OPTION_MFE)
    return subopt(fc; delta, sorted)
end

function subopt_zuker(fc::FoldCompound)
    subs = Tuple{String,Quantity}[]
    s = LibRNA.vrna_subopt_zuker(fc.ptr)
    s == C_NULL && return subs
    i = 1
    while s[i].structure != C_NULL
        push!(subs, (unsafe_string(s[i].structure), s[i].energy * en_unit))
        i += 1
    end
    # TODO: free s recursively
    return subs
end

function subopt_zuker(sequence::AbstractString)
    fc = FoldCompound(sequence; options=LibRNA.VRNA_OPTION_MFE)
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

# kinetics and move sets

function neighbors(fc::FoldCompound, pt::Pairtable;
                   options::Int=LibRNA.VRNA_MOVESET_DEFAULT)
    mptr = LibRNA.vrna_neighbors(fc.ptr, pt.ptr, options)
    moves = Vector{Tuple{Int,Int}}[]
    mptr == C_NULL && return moves
    i = 1
    while mptr[i].pos_5 != 0 && mptr[i].pos_3 != 0
        mlist = [(mptr[i].pos_5, mptr[i].pos_3)]
        if mptr[i].next != C_NULL
            # this move is not an atomic move, so follow the next
            # pointer and collect all moves
            next = mptr[i].next
            j = 1
            while next[j].pos_5 != 0 && next[j].pos_3 != 0
                if next[j].next != C_NULL
                    @warn "unexpected extra branching in neighbors()"
                end
                push!(mlist, (next[j].pos_5, next[j].pos_3))
                j += 1
            end
        end
        push!(moves, mlist)
        i += 1
    end
    LibRNA.vrna_move_list_free(mptr)
    return moves
end

# plotting secondary structures

function plot_coords(structure; plot_type::Symbol=:naview)
    _vrna_plot_coords(structure::AbstractString, cx, cy, type) =
        LibRNA.vrna_plot_coords(structure, cx, cy, type)
    _vrna_plot_coords(structure::Pairtable, cx, cy, type) =
        LibRNA.vrna_plot_coords_pt(structure.ptr, cx, cy, type)

    if plot_type === :simple
        type = LibRNA.VRNA_PLOT_TYPE_SIMPLE
    elseif plot_type === :naview
        type = LibRNA.VRNA_PLOT_TYPE_NAVIEW
    elseif plot_type === :circular
        type = LibRNA.VRNA_PLOT_TYPE_CIRCULAR
    elseif plot_type === :turtle
        type = LibRNA.VRNA_PLOT_TYPE_TURTLE
    elseif plot_type === :puzzler
        type = LibRNA.VRNA_PLOT_TYPE_PUZZLER
    else
        throw(ArgumentError("unknown plot_type $plot_type, options are:" *
                            " :simple, :naview, :circular, :turtle, :puzzler"))
    end
    # TODO: use LibRNA.vrna_alloc
    ptr_cx = Ptr{Ptr{Cfloat}}(Libc.malloc(sizeof(Ptr{Cfloat})))
    ptr_cy = Ptr{Ptr{Cfloat}}(Libc.malloc(sizeof(Ptr{Cfloat})))
    n = _vrna_plot_coords(structure, ptr_cx, ptr_cy, type)
    cx = unsafe_load(ptr_cx)
    cy = unsafe_load(ptr_cy)
    x = [unsafe_load(cx, i) for i = 1:n]
    y = [unsafe_load(cy, i) for i = 1:n]
    Libc.free(cx)
    Libc.free(cy)
    Libc.free(ptr_cx)
    Libc.free(ptr_cy)
    return x, y
end

end # module ViennaRNA
