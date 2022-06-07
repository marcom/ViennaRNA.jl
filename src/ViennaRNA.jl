module ViennaRNA

import Base
using Unitful: @u_str, Quantity, uconvert, ustrip
using UnsafePointers: UnsafePtr

# types
export FoldCompound, Pairtable
# wrapper functions
export bp_distance, bpp, centroid, energy, ensemble_defect,
    heat_capacity, inverse_fold, inverse_pf_fold, mea,
    mean_bp_distance, mfe, neighbors, partfn, pbacktrack, plot_coords,
    prob_of_structure, subopt, subopt_zuker, tree_edit_dist
# additional utility functions
export basepairs, plot_structure, prob_of_basepairs

include("../lib/LibRNA.jl")
import .LibRNA

module LibRNA_Helper
function free_subopt_solutions(ptr::Ptr)
    ptr == C_NULL && return
    i = 1
    while true
        sol = unsafe_load(ptr, i)
        sol.structure == C_NULL && break
        Libc.free(sol.structure)
        i += 1
    end
    Libc.free(ptr)
end

function free_structure_list(ptr::Ptr, num::Integer)
    for i = 1:num
        Libc.free(unsafe_load(ptr, i))
    end
    Libc.free(ptr)
end
end
import .LibRNA_Helper

# TODO: rename en_unit to unit_energy, en_int_unit to unit_energy_int
const en_unit = 1.0u"kcal/mol"
const en_int_unit = 0.01u"kcal/mol"
const unit_temperature = u"°C"

"""
    FoldCompound(seq::AbstractString; [params, temperature, circular, dangles, min_loop_size, uniq_ML])
    FoldCompound(msa::Vector{<:AbstractString}; [params, temperature, circular, dangles, min_loop_size, uniq_ML])

A `FoldCompound` encapsulates nucleotide sequences, energy
parameters, and model details.

Input arguments:
- `seq`: nucleotide sequence, multiple strands are separated by an
  '&' character
- `msa`: multiple sequence alignment, for comparative folding
   (alifold). A vector of sequences which may contain multiple
   strands, denoted by '&', and gap '-' characters
- `params`: energy parameter set, legal values are `:RNA_Turner1999`,
  `:RNA_Turner2004`, `:RNA_Andronescu2007`, `:RNA_Langdon2018`.
  Default is `:RNA_Turner2004`.
- `temperature`: the temperature at which calculations are performed.
  Default is `37u"°C"`.

Model details:
- `circular`: determines if the RNA strand is circular, i.e. the
  5'-end and 3'-end are covalently bonded. Default is `false`.
- `dangles`: how to treat dangling base pairs in multiloops and the
  exterior loop. Can be 0, 1, 2, or 3. See ViennaRNA docs for
  details. Default is `2`.
- `min_loop_length`: the minimum size of a loop (without the closing
   base pair). Default is `3`.
- `uniq_ML`: use unique decomposition for multiloops, needed for
  `pbacktrack` and `subopt`
"""
mutable struct FoldCompound
    ptr         :: Ptr{LibRNA.vrna_fc_s}
    uptr        :: UnsafePtr{LibRNA.vrna_fc_s}
    msa         :: Vector{String}
    msa_strands :: Vector{Vector{SubString{String}}}

    FoldCompound(seq::AbstractString; kwargs...) = FoldCompound([seq]; kwargs...)
    function FoldCompound(msa::Vector{<:AbstractString};
                          model_details::Ptr{LibRNA.vrna_md_s}=Ptr{LibRNA.vrna_md_s}(C_NULL),
                          options::Unsigned=LibRNA.VRNA_OPTION_DEFAULT,
                          params::Symbol=:RNA_Turner2004,
                          temperature::Quantity=37.0u"°C",
                          # model_details options
                          circular::Bool=false,
                          dangles::Int=2,
                          min_loop_size::Int=3,
                          uniq_ML::Bool=false)
        if dangles < 0 || dangles > 3
            throw(ArgumentError("dangles must be 0, 1, 2, or 3"))
        end
        if min_loop_size < 0
            throw(ArgumentError("min_loop_size must be >= 0"))
        end
        if length(msa) > 1
            firstlen = length(first(msa))
            if any(s -> length(s) != firstlen, msa)
                throw(ArgumentError("all sequences in msa must have same length"))
            end
            if any(s -> contains(s, '&'), msa)
                throw(ArgumentError("multiple strands not supported in comparative mode (alifold)"))
            end
        end
        msa_strands = split.(msa, '&')
        if params == :RNA_Turner1999
            err = LibRNA.vrna_params_load_RNA_Turner1999()
        elseif params == :RNA_Turner2004
            err = LibRNA.vrna_params_load_RNA_Turner2004()
        elseif params == :RNA_Andronescu2007
            err = LibRNA.vrna_params_load_RNA_Andronescu2007()
        elseif params == :RNA_Langdon2018
            err = LibRNA.vrna_params_load_RNA_Langdon2018()
        else
            throw(ArgumentError("unknown energy parameters: $(params)"))
        end
        if err == 0
            error("Failed to load energy parameters $params_name")
        end
        temperature_nounit = ustrip(uconvert(u"°C", temperature))
        LibRNA.vrna_md_defaults_temperature(temperature_nounit)
        LibRNA.vrna_md_defaults_circ(Int(circular))
        LibRNA.vrna_md_defaults_dangles(dangles)
        LibRNA.vrna_md_defaults_min_loop_size(min_loop_size)
        LibRNA.vrna_md_defaults_uniq_ML(Int(uniq_ML))

        ptr = if length(msa) == 1
            LibRNA.vrna_fold_compound(first(msa), model_details, options)
        else
            LibRNA.vrna_fold_compound_comparative(msa, model_details, options)
        end
        ptr != C_NULL || error("vrna_fold_compound returned pointer == C_NULL")

        fc = new(ptr, UnsafePtr(ptr), msa, msa_strands)
        finalizer(fc) do x
            # TODO: do we have to call vrna_mx_mfe_free or
            #       vrna_mx_pf_free here ourselves?
            LibRNA.vrna_fold_compound_free(x.ptr)
        end
    end
end

Base.length(fc::FoldCompound) = sum(size(fc))

Base.size(fc::FoldCompound) = ntuple(i -> length(first(fc.msa_strands)[i]), fc.nstrands)

function Base.getproperty(fc::FoldCompound, sym::Symbol)
    if sym == :circular
        return Bool(fc.uptr.params[].model_details.circ[])
    elseif sym == :dangles
        return fc.uptr.params[].model_details.dangles[]
    elseif sym == :has_exp_matrices
        fc.uptr.exp_matrices[] != C_NULL
    elseif sym == :min_loop_size
        return fc.uptr.params[].model_details.min_loop_size[]
    elseif sym == :nstrands
        return length(first(fc.msa_strands))
    elseif sym == :params_name
        return unsafe_string(reinterpret(Ptr{UInt8}, pointer(fc.uptr.params[].param_file)))
    elseif sym == :temperature
        par_temperature = fc.uptr.params[].temperature[]
        md_temperature = fc.uptr.params[].model_details.temperature[]
        if par_temperature != md_temperature
            error("params temperature and model_details temperature don't agree")
        end
        return par_temperature * unit_temperature
    elseif sym == :uniq_ML
        return Bool(fc.uptr.params[].model_details.uniq_ML[])
    else
        # fallback
        getfield(fc, sym)
    end
end

Base.propertynames(fc::FoldCompound) =
    (fieldnames(typeof(fc))...,
     :circular, :dangles, :has_exp_matrices, :min_loop_size, :nstrands,
     :params_name, :temperature, :uniq_ML)

function Base.show(io::IO, mime::MIME"text/plain", fc::FoldCompound)
    strand = "$(fc.nstrands) strand" * (fc.nstrands > 1 ? "s" : "")
    nt = "$(length(fc)) nt$(fc.nstrands > 1 ? " total" : "")"
    circ = fc.circular ? " (circular)" : ""
    comparative = length(fc.msa) > 1 ? " [comparative]" : ""
    println(io, "FoldCompound, $strand, $nt$circ$comparative")
    println(io, "  params      : $(fc.params_name)")
    println(io, "  temperature : $(fc.temperature)")
    println(io, "  options     : circular=$(fc.circular), dangles=$(fc.dangles), min_loop_size=$(fc.min_loop_size), uniq_ML=$(fc.uniq_ML)")
    if length(fc.msa) == 1
        for (i,s) in enumerate(first(fc.msa_strands))
            println(io,   "  strand $i    : $(s)")
        end
    else
        println(io, "  MSA")
        for (i,strands) in enumerate(fc.msa_strands)
            println(io, "      ", join(strands, " & "))
        end
    end
end

"""
    Pairtable(structure)

A `Pairtable` represents a secondary structure given in dot-bracket
notation, e.g. `(((...)))`.  In dot-bracket notation, unpaired bases
are denoted by a dot `.` and base-pairs are denoted by matching
brackets `(` and `)`.
"""
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

basepairs(pt::Pairtable) = [(i, pt[i]) for i = 1:length(pt) if i < pt[i]]

# Notes
# - only works for unpseudoknotted structures
# - vrna_db_from_ptable does the same, but we would have to incur an
#   extra malloc/free to move the string to Julia, so we do it here
#   ourselves
Base.String(pt::Pairtable) =
    join(pt[i] == 0 ? '.' : pt[i] > i ? '(' : ')' for i = 1:length(pt))

function Base.show(io::IO, mime::MIME"text/plain", pt::Pairtable)
    print(io, String(pt))
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
        if i < 1 || val < 0 || i > len || val > len || i == val
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

"""
    bp_distance(structure1, structure2)

Base-pair distance between two secondary structures `structure1` and
`structure2`.
"""
function bp_distance(structure1::AbstractString, structure2::AbstractString)
    length(structure1) == length(structure2) ||
        throw(ArgumentError("input structures have different lengths"))
    return LibRNA.vrna_bp_distance(structure1, structure2)
end

"""
    tree_edit_dist(structure1, structure2; [hit_type])

Tree edit distance between two secondary structures `structure1` and
`structure2`.
"""
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

"""
    mean_bp_distance(fc)
    mean_bp_distance(sequence)

Base-pair distance of all possible pairs of secondary structures
weighted by their Boltzmann probabilities.
"""
function mean_bp_distance(fc::FoldCompound)
    fc.has_exp_matrices ||
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

"""
    ensemble_defect(fc, pt)
    ensemble_defect(sequence, pt)
    ensemble_defect(fc, structure)
    ensemble_defect(sequence, structure)

Ensemble defect between a sequence given as a `FoldCompound` or
`AbstractString` in dot-bracket notation and a secondary structure
given as a `Pairtable` or as an `AbstractString`.
"""
function ensemble_defect(fc::FoldCompound, pt::Pairtable)
    fc.has_exp_matrices ||
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

"""
    prob_of_structure(fc, structure)
    prob_of_structure(sequence, structure)

Boltzmann probability of a structure for a sequence.
"""
function prob_of_structure(fc::FoldCompound, structure::AbstractString)
    fc.has_exp_matrices ||
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

"""
    energy(fc, structure; [verbose])
    energy(sequence, structure; [verbose])

Calculate the free energy of folding for a given `sequence` and
`structure` in dot-bracket notation.
"""
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

"""
    mfe(fc)
    mfe(sequence)

Calculate the minimum free energy structure and energy.
"""
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

"""
    partfn(fc)
    partfn(sequence)

Calculate the partition function, returning a string representation of
secondary structure pairing probabilities and the ensemble free energy
`-RT log Z`.  The string consists of the letters `. , | { } ( )` which
denote bases that are essentially unpaired, weakly paired, strongly
paired without upstream/downstream preference, weakly
upstream/downstream paired, or strongly upstream/downstream paired.
"""
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

"""
    bpp(fc)
    bpp(sequence)

Calculate the matrix of base-pair probabilities.
"""
function bpp(fc::FoldCompound)
    # TODO: return upper triangular matrix type
    fc.has_exp_matrices ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    n = length(fc)
    p = zeros(n,n)
    f = unsafe_load(fc.ptr)
    index = f.iindx
    probs = unsafe_load(f.exp_matrices).probs
    for i = 1:n-1
        for j = i+1:n
            p[i,j] = unsafe_load(probs, unsafe_load(index, i + 1) - j + 1)
            p[j,i] = p[i,j]
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

"""
    pbacktrack(fc; [num_samples=1])
    pbacktrack(sequence; [num_samples=1])

Sample `num_samples` secondary structures according to their Boltzmann
probabilities.
"""
function pbacktrack(fc::FoldCompound;
                    num_samples::Integer=1,
                    options::Integer=LibRNA.VRNA_PBACKTRACK_DEFAULT)
    fc.has_exp_matrices ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    fc.uniq_ML ||
        throw(ArgumentError("must have fc.uniq_ML == true"))
    s = LibRNA.vrna_pbacktrack_num(fc.ptr, num_samples, options)
    samples = String[]
    s == C_NULL && return samples
    i = 1
    while unsafe_load(s, i) != C_NULL
        push!(samples, unsafe_string(unsafe_load(s, i)))
        i += 1
    end
    LibRNA_Helper.free_structure_list(s, num_samples)
    return samples
end

function pbacktrack(sequence::AbstractString; num_samples::Integer=1,
                    options::Integer=LibRNA.VRNA_PBACKTRACK_DEFAULT)
    fc = FoldCompound(sequence; uniq_ML=true)
    partfn(fc)
    return pbacktrack(fc; num_samples, options)
end

# maximum expected accuracy (MEA) structure

"""
    mea(fc); [gamma=1.0])
    mea(sequence); [gamma=1.0])

Maximum expected accuracy (MEA) structure.
"""
function mea(fc::FoldCompound; gamma=1.0)
    fc.has_exp_matrices ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    ptr_mea   = Ptr{Cfloat}(LibRNA.vrna_alloc(sizeof(Cfloat)))
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

"""
    centroid(fc)
    centroid(sequence)

Centroid structure of the ensemble, the secondary structure with the
smallest Boltzmann probability weighted base-pair distance to all
other structures.
"""
function centroid(fc::FoldCompound)
    fc.has_exp_matrices ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    ptr_dist = Ptr{Cdouble}(LibRNA.vrna_alloc(sizeof(Cdouble)))
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

"""
    subopt(fc; [delta], [sorted])
    subopt(sequence; [delta], [sorted])

All suboptimal structures with an energy not more than `delta` above
the minimum free energy structure.
"""
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
    while (si = unsafe_load(s, i)).structure != C_NULL
        push!(subs, (unsafe_string(si.structure), si.energy * en_unit))
        i += 1
    end
    LibRNA_Helper.free_subopt_solutions(s)
    return subs
end

function subopt(sequence::AbstractString; delta::Quantity, sorted::Bool=true)
    fc = FoldCompound(sequence; uniq_ML=true, options=LibRNA.VRNA_OPTION_MFE)
    return subopt(fc; delta, sorted)
end

"""
    subopt_zuker(fc)
    subopt_zuker(sequence)

All suboptimal structures according to the Zuker algorithm.
"""
function subopt_zuker(fc::FoldCompound)
    subs = Tuple{String,Quantity}[]
    s = LibRNA.vrna_subopt_zuker(fc.ptr)
    s == C_NULL && return subs
    i = 1
    while (si = unsafe_load(s, i)).structure != C_NULL
        push!(subs, (unsafe_string(si.structure), si.energy * en_unit))
        i += 1
    end
    LibRNA_Helper.free_subopt_solutions(s)
    return subs
end

function subopt_zuker(sequence::AbstractString)
    fc = FoldCompound(sequence; options=LibRNA.VRNA_OPTION_MFE)
    return subopt_zuker(fc)
end

# inverse folding / sequence design

"""
    inverse_fold(start, target)

Inverse folding (aka sequence design) for a given `target` structure,
starting from the sequence `start`.

The stopping criterion for the search is that the sequence has the
target structure as its minimum free energy structure.
"""
function inverse_fold(start::AbstractString, target::AbstractString)
    length(start) == length(target) ||
        throw(ArgumentError("start and target must have same size"))
    # TODO: there is no copy(::String), so this is a substitute
    #       needed because inverse_fold mutates its first argument
    design = join(collect(start))
    dist = LibRNA.inverse_fold(design, target)
    return design, dist
end

"""
    inverse_pf_fold(start, target)

Inverse folding (aka sequence design) for a given `target` structure,
starting from the sequence `start`.

The search in sequence space tries to maximise the probability of the
target structure.
"""
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

"""
    neighbors(fc, pt)

All neighbors of a secondary structure `pt` given as a `Pairtable`.
"""
function neighbors(fc::FoldCompound, pt::Pairtable;
                   options::Int=LibRNA.VRNA_MOVESET_DEFAULT)
    mptr = LibRNA.vrna_neighbors(fc.ptr, pt.ptr, options)
    moves = Vector{Tuple{Int,Int}}[]
    mptr == C_NULL && return moves
    i = 1
    while (mptr_i = unsafe_load(mptr, i)).pos_5 != 0 && mptr_i.pos_3 != 0
        mlist = [(mptr_i.pos_5, mptr_i.pos_3)]
        if mptr_i.next != C_NULL
            # this move is not an atomic move, so follow the next
            # pointer and collect all moves
            next = mptr_i.next
            j = 1
            while (next_j = unsafe_load(next, j)).pos_5 != 0 && next_j.pos_3 != 0
                if next_j.next != C_NULL
                    @warn "unexpected extra branching in neighbors()"
                end
                push!(mlist, (next_j.pos_5, next_j.pos_3))
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

"""
    plot_coords(structure; [plot_type])

Plots a secondary structure, returning the coordinates. The optional
`plot_type` parameter can be `:simple`, `:naview`, `:circular`,
`:turtle`, or `:puzzler`, with the default being `:naview`.
"""
function plot_coords(structure; plot_type::Symbol=:simple)
    _vrna_plot_coords(structure::AbstractString, cx, cy, type) =
        LibRNA.vrna_plot_coords(structure, cx, cy, type)
    _vrna_plot_coords(structure::Pairtable, cx, cy, type) =
        LibRNA.vrna_plot_coords_pt(structure.ptr, cx, cy, type)

    if length(structure) == 0
        return Float32[], Float32[]
    end
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
    ptr_cx = Ptr{Ptr{Cfloat}}(LibRNA.vrna_alloc(sizeof(Ptr{Cfloat})))
    ptr_cy = Ptr{Ptr{Cfloat}}(LibRNA.vrna_alloc(sizeof(Ptr{Cfloat})))
    n = _vrna_plot_coords(structure, ptr_cx, ptr_cy, type)
    if n == 0
        throw(ErrorException("error returned by vrna_plot_coords"))
    end
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

# heat capacity

"""
    heat_capacity(fc, Tmin, Tmax, [Tincrement=1u"°C"]; [mpoints=2])
    heat_capacity(sequence, Tmin, Tmax, [Tincrement=1u"°C"]; [mpoints=2])

Computes the specific heat of an RNA in a temperature range, given by
`Tmin`, `Tmax`, and `Tincrement`, from the partition function by
numeric differentiation.  The parameter `mpoints` determines how
smooth the curve should be: to calculate second derivatives, a
parabola is fit to `mpoints` + 1 data points, and increasing the
`mpoints` parameter produces a smoother curve.
"""
function heat_capacity(fc::FoldCompound, Tmin::Quantity, Tmax::Quantity,
                       Tincrement::Quantity=1u"°C"; mpoints::Integer=2)
    Tmin = ustrip(uconvert(u"°C", Tmin))
    Tmax = ustrip(uconvert(u"°C", Tmax))
    Tincrement = ustrip(uconvert(u"°C", Tincrement))
    type_T = typeof(1.0f0u"°C")
    type_c = typeof(1.0f0u"kcal/mol/K")
    ptr = LibRNA.vrna_heat_capacity(fc.ptr, Tmin, Tmax, Tincrement, mpoints)
    if ptr == C_NULL
        error("failure running vrna_heat_capacity")
    end
    hcs = Tuple{type_T, type_c}[]
    i = 1
    while (h = unsafe_load(ptr, i)).temperature >= Tmin
        T = h.temperature * u"°C"
        c = h.heat_capacity * u"kcal/mol/K"
        push!(hcs, (T, c))
        i += 1
    end
    Libc.free(ptr)
    return hcs
end

function heat_capacity(sequence::AbstractString, Tmin::Quantity, Tmax::Quantity,
                       Tincrement::Quantity=1.0u"°C"; mpoints::Integer=2)
    heat_capacity(FoldCompound(sequence), Tmin, Tmax, Tincrement; mpoints)
end

include("utils.jl")
include("plot_structure.jl")

end # module ViennaRNA
