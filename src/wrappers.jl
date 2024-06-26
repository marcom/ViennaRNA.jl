module Private
using Unitful: @u_str
import ViennaRNA.LibRNA

const unit_energy = 1.0u"kcal/mol"
const unit_energy_int = 0.01u"kcal/mol"
const unit_temperature = u"°C"

function _makedict(prefix::AbstractString)
    vs = filter(startswith(prefix), String.(names(LibRNA)))
    return Dict(Symbol(lowercase(replace(v, prefix => ""))) => getproperty(LibRNA, Symbol(v)) for v in vs)
end

const FOLDCOMPOUND_OPTIONS = _makedict("VRNA_OPTION_")
# const OPTIONS = Dict(
#     :default      => LibRNA.VRNA_OPTIONS_DEFAULT,
#     ...
# )

const FOLDCOMPOUND_TYPE = _makedict("VRNA_FC_TYPE_")
# const OPTIONS = Dict(
#     :single      => LibRNA.VRNA_FC_TYPE_SINGLE,
#     :comparative => LibRNA.VRNA_FC_TYPE_COMPARATIVE,
# )
const FOLDCOMPOUND_TYPE_C2JL = Dict(v => k for (k, v) in FOLDCOMPOUND_TYPE)

const PARAMS_LOADFNS =
    Dict(Symbol(replace(v, "vrna_params_load_" => "")) => getproperty(LibRNA, Symbol(v))
         for v in (String.(names(LibRNA)) |>
             vs -> filter!(startswith(r"vrna_params_load_(RNA|DNA)"), vs) |>
             vs -> filter!(s -> s ∉ ["vrna_params_load_RNA_misc_special_hairpins"], vs))
    )
# const PARAMS_LOADFNS = Dict(
#     :RNA_Turner1999     => LibRNA.vrna_params_load_RNA_Turner1999,
#     ...
# )

const PLOT_COORDS_PLOT_TYPE = _makedict("VRNA_PLOT_TYPE_")
# const PLOT_COORDS_PLOT_TYPE = Dict(
#     :default  => LibRNA.VRNA_PLOT_TYPE_DEFAULT,
#     ...
# )

const SAMPLE_STRUCTURES_OPTIONS = _makedict("VRNA_PBACKTRACK_")
# const SAMPLE_STRUCTURES_OPTIONS = Dict(
#     :default      => LibRNA.VRNA_PBACKTRACK_DEFAULT,
#     ...
# )

# data needed to auto-generate wrappers and tests for the
# vrna_sc_mod_* convenience functions for specific base modifications
const SC_MOD_PRESET_FUNCTIONS = [
    (short="m6A",            desc="N6-methyl-adenosine (m6A)",  unmod_base="A"),
    (short="pseudouridine",  desc="pseudouridine",              unmod_base="U"),
    (short="inosine",        desc="inosine",                    unmod_base="G"),
    (short="7DA",            desc="7-deaza-adenosine (7DA)",    unmod_base="A"),
    (short="purine",         desc="Purine (a.k.a. nebularine)", unmod_base="A"),
    (short="dihydrouridine", desc="dihydrouridine",             unmod_base="U"),
]

end # module Private
import .Private

# TODO: automate creation of params_load_(DNA|RNA)_* functions so it's
#       always in sync with ViennaRNA
"""
    params_load_defaults()

Load the default parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`. The default
energy parameter set is `:RNA_Turner2004`.
"""
function params_load_defaults()
    if LibRNA.vrna_params_load_defaults() == 0
        error("couldn't load params in params_load_defaults()")
    end
end

"""
    params_load_DNA_Mathews1999()

Load the `:DNA_Mathews1999` parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`.
"""
function params_load_DNA_Mathews1999()
    if LibRNA.vrna_params_load_DNA_Mathews1999() == 0
        error("couldn't load params in params_load_DNA_Mathews1999()")
    end
end

"""
    params_load_DNA_Mathews2004()

Load the `:DNA_Mathews2004` parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`.
"""
function params_load_DNA_Mathews2004()
    if LibRNA.vrna_params_load_DNA_Mathews2004() == 0
        error("couldn't load params in params_load_DNA_Mathews2004()")
    end
end

"""
    params_load_RNA_Andronescu2007()

Load the `:RNA_Andronescu2007` parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`.
"""
function params_load_RNA_Andronescu2007()
    if LibRNA.vrna_params_load_RNA_Andronescu2007() == 0
        error("couldn't load params in params_load_RNA_Andronescu2007()")
    end
end

"""
    params_load_RNA_Langdon2018()

Load the `:RNA_Langdon2018` parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`.
"""
function params_load_RNA_Langdon2018()
    if LibRNA.vrna_params_load_RNA_Langdon2018() == 0
        error("couldn't load params in params_load_RNA_Langdon2018()")
    end
end

"""
    params_load_RNA_Turner1999()

Load the `:RNA_Turner1999` parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`.
"""
function params_load_RNA_Turner1999()
    if LibRNA.vrna_params_load_RNA_Turner1999() == 0
        error("couldn't load params in params_load_RNA_Turner1999()")
    end
end

"""
    params_load_RNA_Turner2004()

Load the `:RNA_Turner2004` parameter set into the global variables holding the
energy parameters used when creating a `FoldCompound`.
"""
function params_load_RNA_Turner2004()
    if LibRNA.vrna_params_load_RNA_Turner2004() == 0
        error("couldn't load params in params_load_RNA_Turner2004()")
    end
end

function params_load_RNA_misc_special_hairpins()
    if LibRNA.vrna_params_load_RNA_misc_special_hairpins() == 0
        error("couldn't load params in params_load_RNA_misc_special_hairpins()")
    end
end

"""
    params_load_from_string(str, name)

Load the energy parameters from a string `str` and use `name` as the
name for the energy parmeter set.
"""
function params_load_from_string(str::AbstractString, name::AbstractString)
    if LibRNA.vrna_params_load_from_string(str, name, LibRNA.VRNA_PARAMETER_FORMAT_DEFAULT) == 0
        error("error loading energy parameters from string:\n$str")
    end
end

"""
    params_load(params::Symbol)

Load the energy parameter set specified by `params`. Options are
`$(keys(Private.PARAMS_LOADFNS))`.
"""
function params_load(params::Symbol)
    # TODO: support :defaults, :RNA_misc_special_hairpins
    if ! haskey(Private.PARAMS_LOADFNS, params)
        throw(ArgumentError("unknown energy parameter set $params, possible values are: " *
            "$(sort(collect(keys(Private.PARAMS_LOADFNS))))"))
    end
    loadparams = Private.PARAMS_LOADFNS[params]
    err = loadparams()
    if err == 0
        error("Failed to load energy parameters $params")
    end
    return nothing
end

"""
    params_load(filename)

Read energy parameters from a file `filename`.
"""
function params_load(filename::AbstractString)
    if LibRNA.vrna_params_load(filename, LibRNA.VRNA_PARAMETER_FORMAT_DEFAULT) == 0
        error("error loading energy parameters from file $filename")
    end
end

"""
    params_save(filename)

Write the current energy parameters to the file `filename`.
"""
function params_save(filename::AbstractString)
    if LibRNA.vrna_params_save(filename, LibRNA.VRNA_PARAMETER_FORMAT_DEFAULT) == 0
        error("error saving energy parameters to file $filename")
    end
end


"""
    FoldCompound(seq::AbstractString; [options, params, temperature], [model_details...])
    FoldCompound(msa::Vector{<:AbstractString}; [options, params, temperature], [model_details...])

A `FoldCompound` encapsulates nucleotide sequences, energy
parameters, and model details.

Input arguments:
- `seq`: nucleotide sequence, multiple strands are separated by an
  '&' character
- `msa`: multiple sequence alignment, for comparative folding
   (alifold). A vector of sequences which may contain multiple
   strands, denoted by '&', and gap '-' characters
- `options`: one or more of $(sort(collect(keys(Private.FOLDCOMPOUND_OPTIONS))))
- `temperature`: the temperature at which calculations are performed.
  Default is `37u"°C"`.

Model details (additional keyword arguments):
- `circular`: determines if the RNA strand is circular, i.e. the
  5'-end and 3'-end are covalently bonded. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_CIRC))`.
- `dangles`: how to treat dangling base pairs in multiloops and the
  exterior loop. Can be 0, 1, 2, or 3. See ViennaRNA docs for
  details. Default is `$(Int(LibRNA.VRNA_MODEL_DEFAULT_DANGLES))`.
- `gquadruplex`: allow G-quadruplexes in predictions. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_GQUAD))`.
- `log_ML`: use logarithmic energy model for multiloops. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_LOG_ML))`.
- `max_bp_span`: maximum number of bases over which a basepair can
  span. Default value is
  `$(Int(LibRNA.VRNA_MODEL_DEFAULT_WINDOW_SIZE))`.
- `min_loop_length`: the minimum size of a loop (without the closing
   base pair). Default is `$(Int(LibRNA.TURN))`.
- `no_GU_basepairs`: disallow G-U basepairs. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_NO_GU))`.
- `no_GU_closure`: disallow G-U basepairs as closing pairs for
  loops. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_NO_GU_CLOSURE))`.
- `no_lonely_pairs`: disallow isolated base pairs. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_NO_LP))`.
- `special_hairpins`: use special hairpin energies for certain tri-,
  tetra- and hexloops. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_SPECIAL_HP))`.
- `uniq_ML`: use unique decomposition for multiloops, needed for
  `sample_structures` and `subopt`. Default is
  `$(Bool(LibRNA.VRNA_MODEL_DEFAULT_UNIQ_ML))`.
- `window_size`: window size to be used for local calculations
  performed in a window moving over the sequence. This value is
  ignored unless the `:window` option is set in the `FoldCompound`
  `options`. The default value for `window_size` is
  `$(Int(LibRNA.VRNA_MODEL_DEFAULT_WINDOW_SIZE))`.

Global energy parameters used when constructing a `FoldCompound` can be changed with the following functions:
`ViennaRNA.params_load`, `ViennaRNA.params_load_defaults`,
`ViennaRNA.params_load_DNA_Mathews1999`, `ViennaRNA.params_load_DNA_Mathews2004`,
`ViennaRNA.params_load_RNA_Andronescu2007`, `ViennaRNA.params_load_RNA_Langdon2018`,
`ViennaRNA.params_load_RNA_Turner1999`, `ViennaRNA.params_load_RNA_Turner2004`.
"""
mutable struct FoldCompound
    ptr         :: Ptr{LibRNA.vrna_fc_s}
    uptr        :: UnsafePtr{LibRNA.vrna_fc_s}
    msa         :: Vector{String}
    msa_strands :: Vector{Vector{SubString{String}}}
    options     :: Vector{Symbol}

    FoldCompound(seq::AbstractString; kwargs...) = FoldCompound([seq]; kwargs...)
    function FoldCompound(msa::Vector{<:AbstractString};
                          model_details::Ptr{LibRNA.vrna_md_s}=Ptr{LibRNA.vrna_md_s}(C_NULL),
                          options::Vector{Symbol}=[:default],
                          temperature::Quantity=37.0u"°C",
                          # model_details options
                          circular::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_CIRC),
                          dangles::Int=Int(LibRNA.VRNA_MODEL_DEFAULT_DANGLES),
                          gquadruplex::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_GQUAD),
                          log_ML::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_LOG_ML),
                          max_bp_span::Int=Int(LibRNA.VRNA_MODEL_DEFAULT_MAX_BP_SPAN),
                          min_loop_size::Int=Int(LibRNA.TURN),
                          no_GU_basepairs::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_NO_GU),
                          no_GU_closure::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_NO_GU_CLOSURE),
                          no_lonely_pairs::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_NO_LP),
                          special_hairpins::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_SPECIAL_HP),
                          uniq_ML::Bool=Bool(LibRNA.VRNA_MODEL_DEFAULT_UNIQ_ML),
                          window_size::Int=Int(LibRNA.VRNA_MODEL_DEFAULT_WINDOW_SIZE))
        (dangles ∈ (0, 1, 2, 3)) || throw(ArgumentError("dangles must be 0, 1, 2, or 3"))
        min_loop_size ≥ 0 || throw(ArgumentError("min_loop_size must be ≥ 0"))
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
        if any(x -> !haskey(Private.FOLDCOMPOUND_OPTIONS, x), options)
            throw(ArgumentError("unknown FoldCompound options: $(setdiff(options, keys(Private.FOLDCOMPOUND_OPTIONS))),"
                                * " allowed values are: $(keys(Private.FOLDCOMPOUND_OPTIONS))"))
        end
        temperature_nounit = ustrip(uconvert(u"°C", temperature))
        LibRNA.vrna_md_defaults_temperature(temperature_nounit)
        LibRNA.vrna_md_defaults_circ(Int(circular))
        LibRNA.vrna_md_defaults_dangles(dangles)
        LibRNA.vrna_md_defaults_gquad(Int(gquadruplex))
        LibRNA.vrna_md_defaults_logML(Int(log_ML))
        LibRNA.vrna_md_defaults_max_bp_span(max_bp_span)
        LibRNA.vrna_md_defaults_min_loop_size(min_loop_size)
        LibRNA.vrna_md_defaults_noGU(Int(no_GU_basepairs))
        LibRNA.vrna_md_defaults_noGUclosure(Int(no_GU_closure))
        LibRNA.vrna_md_defaults_noLP(Int(no_lonely_pairs))
        LibRNA.vrna_md_defaults_special_hp(Int(special_hairpins))
        LibRNA.vrna_md_defaults_uniq_ML(Int(uniq_ML))
        LibRNA.vrna_md_defaults_window_size(window_size)

        options_int = UInt(0)
        for opt in options
            options_int |= Private.FOLDCOMPOUND_OPTIONS[opt]
        end
        ptr = if length(msa) == 1
            LibRNA.vrna_fold_compound(first(msa), model_details, options_int)
        else
            LibRNA.vrna_fold_compound_comparative(msa, model_details, options_int)
        end
        ptr != C_NULL || error("vrna_fold_compound returned pointer == C_NULL")

        fc = new(ptr, UnsafePtr(ptr), msa, msa_strands, options)
        if fc.min_loop_size != min_loop_size
            @warn "min_loop_size was specified as $min_loop_size, but got set to $(fc.min_loop_size) by ViennaRNA"
        end
        finalizer(fc) do x
            LibRNA.vrna_fold_compound_free(x.ptr)
        end
    end
end

# TODO: fc.uptr.length[] |> Int
Base.length(fc::FoldCompound) = sum(size(fc); init=0)

Base.size(fc::FoldCompound) = ntuple(i -> length(first(fc.msa_strands)[i]), fc.nstrands)

function Base.getproperty(fc::FoldCompound, sym::Symbol)
    # TODO: simplify this big if statement? but want to make sure the
    # Julia compiler can still constant-propagate through the desired
    # symbol (`sym` is a constant for most calls from the POV of the
    # caller of getproperty)

    # Notes
    # - must use getfield(fc, :uptr) here instead of fc.uptr,
    #   otherwise JET complains
    if sym === :circular
        return Bool(getfield(fc, :uptr).params[].model_details.circ[])
    elseif sym === :dangles
        return Int(getfield(fc, :uptr).params[].model_details.dangles[])
    elseif sym === :gquadruplex
        return Bool(getfield(fc, :uptr).params[].model_details.gquad[])
    elseif sym === :has_matrices
        return getfield(fc, :uptr).matrices[] != C_NULL
    elseif sym === :has_exp_matrices
        return getfield(fc, :uptr).exp_matrices[] != C_NULL
    elseif sym === :log_ML
        return Bool(getfield(fc, :uptr).params[].model_details.logML[])
    elseif sym === :max_bp_span
        return Int(getfield(fc, :uptr).params[].model_details.max_bp_span[])
    elseif sym === :min_loop_size
        return Int(getfield(fc, :uptr).params[].model_details.min_loop_size[])
    elseif sym === :no_GU_basepairs
        return Bool(getfield(fc, :uptr).params[].model_details.noGU[])
    elseif sym === :no_GU_closure
        return Bool(getfield(fc, :uptr).params[].model_details.noGUclosure[])
    elseif sym === :no_lonely_pairs
        return Bool(getfield(fc, :uptr).params[].model_details.noLP[])
    elseif sym === :nstrands
        return Int(getfield(fc, :uptr).strands[])
    elseif sym === :params_name
        return unsafe_string(reinterpret(Ptr{UInt8}, pointer(getfield(fc, :uptr).params[].param_file)))
    elseif sym === :sequence
        if getfield(fc, :uptr).sequence[] == C_NULL
            return nothing
        end
        return unsafe_string(pointer(getfield(fc, :uptr).sequence[]))
    elseif sym === :special_hairpins
        return Bool(getfield(fc, :uptr).params[].model_details.special_hp[])
    elseif sym === :temperature
        par_temperature = getfield(fc, :uptr).params[].temperature[]
        md_temperature = getfield(fc, :uptr).params[].model_details.temperature[]
        if par_temperature != md_temperature
            error("params temperature and model_details temperature don't agree: $par_temperature != $md_temperature")
        end
        return par_temperature * Private.unit_temperature
    elseif sym === :type
        return Private.FOLDCOMPOUND_TYPE_C2JL[getfield(fc, :uptr).type[]]
    elseif sym === :uniq_ML
        return Bool(getfield(fc, :uptr).params[].model_details.uniq_ML[])
    elseif sym === :window_size
        return Int(getfield(fc, :uptr).params[].model_details.window_size[])
    # matrices_{c,fML,fM1,f5,f3,fM2,Fc,FcH,FcI,FcM}
    elseif startswith(String(sym), "matrices_")
        fc.has_matrices || begin
            @warn "no information stored yet, run mfe() first (getfield(fc, :uptr).matrices[] == C_NULL)"
            return nothing
        end
        mat = getfield(fc, :uptr).matrices[]
        realsym = Symbol(last(split(string(sym), '_')))
        if realsym ∈ (:c, :fML, :fM1)
            jindx = getfield(fc, :uptr).jindx[]
            return unsafe_loadmat(fc, getproperty(mat, realsym)[];
                                  indexfn=(i,j)->jindx[j + 1] + i + 1)
        elseif realsym ∈ (:f5, :f3, :fM2)
            return unsafe_loadvec(fc, getproperty(mat, realsym)[])
        elseif realsym ∈ (:Fc, :FcH, :FcI, :FcM)
            return Int(getproperty(mat, realsym)[])
        else
            return getfield(fc, sym) # fallback
        end
    # exp_matrices_{q,qb,qm,qm1,probs,qm2,expMLbase,scale}
    elseif startswith(String(sym), "exp_matrices_")
        fc.has_exp_matrices || begin
            @warn "no information stored yet, run partfn() first (getfield(fc, :uptr).exp_matrices[] == C_NULL)"
            return nothing
        end
        expmat = getfield(fc, :uptr).exp_matrices[]
        realsym = Symbol(last(split(string(sym), '_')))
        if realsym ∈ (:q, :qb, :qm, :probs)
            iindx = getfield(fc, :uptr).iindx[]
            return unsafe_loadmat(fc, getproperty(expmat, realsym)[];
                                  indexfn=(i,j)->iindx[i + 1] - j + 1)
        elseif realsym ∈ (:qm1,)
            jindx = getfield(fc, :uptr).jindx[]
            return unsafe_loadmat(fc, getproperty(expmat, realsym)[];
                                  indexfn=(i,j)->jindx[j + 1] + i + 1)
        elseif realsym ∈ (:qm2, :expMLbase, :scale)
            return unsafe_loadvec(fc, getproperty(expmat, realsym)[])
        else
            return getfield(fc, sym) # fallback
        end
    else
        return getfield(fc, sym) # fallback
    end
end

Base.propertynames(fc::FoldCompound) =
    (fieldnames(typeof(fc))...,
     :circular, :dangles, :gquadruplex, :has_matrices,
     :has_exp_matrices, :log_ML, :max_bp_span, :min_loop_size,
     :no_GU_basepairs, :no_GU_closure, :no_lonely_pairs, :nstrands,
     :params_name, :sequence, :special_hairpins, :temperature, :type, :uniq_ML,
     :window_size,
     :matrices_c, :matrices_fML, :matrices_fM1,
     :matrices_f5, :matrices_f3, :matrices_fM2,
     :matrices_Fc, :matrices_FcH, :matrices_FcI, :matrices_FcM,
     :exp_matrices_q, :exp_matrices_qb, :exp_matrices_qm, :exp_matrices_qm1, :exp_matrices_probs,
     :exp_matrices_qm2, :exp_matrices_expMLbase, :exp_matrices_scale,
     )

function Base.show(io::IO, mime::MIME"text/plain", fc::FoldCompound)
    strand = "$(fc.nstrands) strand" * (fc.nstrands > 1 ? "s" : "")
    nt = "$(length(fc)) nt$(fc.nstrands > 1 ? " total" : "")"
    circ = fc.circular ? " (circular)" : ""
    println(io, "FoldCompound, $strand, $nt$circ")
    println(io, "  type          : $(fc.type)")
    println(io, "  options       : $(fc.options)")
    println(io, "  params        : $(fc.params_name)")
    println(io, "  temperature   : $(fc.temperature)")
    println(io, "  model_details : circular=$(fc.circular), dangles=$(fc.dangles), gquadruplex=$(fc.gquadruplex), log_ML=$(fc.log_ML),")
    println(io, "                  max_bp_span=$(fc.max_bp_span), min_loop_size=$(fc.min_loop_size), no_GU_basepairs=$(fc.no_GU_basepairs), no_GU_closure=$(fc.no_GU_closure),")
    println(io, "                  no_lonely_pairs=$(fc.no_lonely_pairs), special_hairpins=$(fc.special_hairpins), uniq_ML=$(fc.uniq_ML), window_size=$(fc.window_size)")
    if fc.type === :single
        for (i,s) in enumerate(first(fc.msa_strands))
            println(io,   "  strand $i      : $(s)")
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

# bpp: basepair probabilities

"""
    bpp(fc)
    bpp(sequence)

Calculate the matrix of base-pair probabilities.
"""
function bpp(fc::FoldCompound)
    # TODO: return symmetric (only store upper triangular part) matrix type
    fc.has_exp_matrices ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    p_ut = fc.exp_matrices_probs
    # p_ut is upper-triangular with empty diagonal, make it symmetric
    return p_ut + transpose(p_ut)
end

function bpp(sequence::AbstractString)
    fc = FoldCompound(sequence)
    partfn(fc)
    res = bpp(fc)
    finalize(fc)
    return res
end

# basepair distance

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
    res = centroid(fc)
    finalize(fc)
    return res
end

# energy: Gibbs free energy change of folding

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
    return en * Private.unit_energy
end

function energy(sequence::AbstractString,
                structure::AbstractString;
                verbose::Bool=false,
                verbosity_level::Integer=LibRNA.VRNA_VERBOSITY_DEFAULT)
    fc = FoldCompound(sequence; options=[:eval_only])
    res = energy(fc, structure; verbose, verbosity_level)
    finalize(fc)
    return res
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
    res = ensemble_defect(fc, Pairtable(structure))
    finalize(fc)
    return res
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
    fc = FoldCompound(sequence)
    res = heat_capacity(fc, Tmin, Tmax, Tincrement; mpoints)
    finalize(fc)
    return res
end

# initialise the random number generator seed value

"""
    init_rand_seed(seed::Integer)

Seed the random number generator used by ViennaRNA with the value
`seed`.
"""
init_rand_seed(seed::Integer) = LibRNA.vrna_init_rand_seed(seed)

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
    return design, RT_log_p * Private.unit_energy
end

# mea: maximum expected accuracy structure

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
    res = mea(fc; gamma)
    finalize(fc)
    return res
end

# mean all-pairs basepair distance, weighted by Boltzmann probabilities

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
    res = mean_bp_distance(fc)
    finalize(fc)
    return res
end

# mfe: minimum free energy structure

"""
    mfe(fc)
    mfe(sequence)

Calculate the minimum free energy structure and energy.
"""
function mfe(fc::FoldCompound)
    if :window in fc.options
        throw(ArgumentError("mfe() doesn't work with the `:window` option set for FoldCompound"))
    end
    seqlen = length(fc)
    cstr_structure = Ptr{Cchar}(LibRNA.vrna_alloc(seqlen + 1))
    en_mfe = LibRNA.vrna_mfe(fc.ptr, cstr_structure)
    structure = unsafe_string(cstr_structure)
    Libc.free(cstr_structure)
    return structure, en_mfe * Private.unit_energy
end

function mfe(sequence::AbstractString)
    fc = FoldCompound(sequence; options=[:mfe])
    res = mfe(fc)
    finalize(fc)
    return res
end

# mfe_window: sliding window over sequence

struct ResultWindowMFE
    seqidx :: UnitRange{Cint}
    dbn :: String
    en :: typeof(Private.unit_energy)
end

# callback function has to be defined before mfe_window
function _callback_mfe_window_save(startidx::Cint, endidx::Cint,
                                   dbn_ptr::Ptr{Cchar}, en::Cfloat,
                                   data::Ptr{Cvoid})
    # from ViennaRNA header: src/ViennaRNA/mfe_window.h
    #
    # @param start provides the first position of the hit (1-based, relative to entire sequence/alignment)
    # @param end provides the last position of the hit (1-based, relative to the entire sequence/alignment)
    # @param structure provides the (sub)structure in dot-bracket notation
    # @param en is the free energy of the structure hit in kcal/mol
    # @param data is some arbitrary data pointer passed through by the function executing the callback
    res = unsafe_pointer_to_objref(data)::Vector{ResultWindowMFE}
    dbn = unsafe_string(dbn_ptr)
    push!(res, ResultWindowMFE(startidx:endidx, dbn, en * Private.unit_energy))
    return nothing
end

"""
    mfe_window(fc::FoldCompound) -> Vector{ResultWindowMFE}
    mfe_window(seq::AbstractString; window_size) -> Vector{ResultWindowMFE}

Run a minimum free energy structure calculation in a sliding
window. The `FoldCompound` must have the `window_size` and correct
`options` set to be able to call `mfe_window` on it:
```
FoldCompound(seq; window_size, options=[:window])
```
"""
function mfe_window(fc::FoldCompound)
    :window in fc.options || throw(ArgumentError("fc must have :window set in options"))
    cb_fnptr = @cfunction(_callback_mfe_window_save, Cvoid,
                          (Cint, Cint, Ptr{Cchar}, Cfloat, Ptr{Cvoid}))
    res = ResultWindowMFE[]
    # TODO: is this GC.@preserve needed?
    GC.@preserve res ViennaRNA.LibRNA.vrna_mfe_window_cb(
        fc.ptr, cb_fnptr, pointer_from_objref(res)
    )
    return res
end

function mfe_window(seq::AbstractString; window_size::Int)
    fc = FoldCompound(seq; window_size, options=[:default, :window])
    res = mfe_window(fc)
    finalize(fc)
    return res
end

# mfe_window_channel: like mfe_window, but put results into a Channel

function _callback_mfe_window_channel(startidx::Cint, endidx::Cint,
                                      dbn_ptr::Ptr{Cchar}, en::Cfloat,
                                      data::Ptr{Cvoid})
    chan = unsafe_pointer_to_objref(data)::Channel
    dbn = unsafe_string(dbn_ptr)
    put!(chan, ResultWindowMFE(startidx:endidx, dbn, en * Private.unit_energy))
    return nothing
end

"""
    mfe_window_channel(fc::FoldCompound)
    mfe_window_channel(seq::AbstractString; window_size::Int)

Used to create a writer that puts mfe_window calculation results in a
`Channel` that can then be read from there and processed.

Usage:

```
chan = mfe_window_channel(seq; window_size)
take!(chan)
```
"""
function mfe_window_channel(fc::FoldCompound)
    :window in fc.options || throw(ArgumentError("fc must have :window set in options"))
    function mfe_window_channel_writer(chan::Channel{ResultWindowMFE})
        cb_fnptr = @cfunction(_callback_mfe_window_channel, Cvoid,
                              (Cint, Cint, Ptr{Cchar}, Cfloat, Ptr{Cvoid}))
        # TODO: is this GC.@preserve needed?
        GC.@preserve chan ViennaRNA.LibRNA.vrna_mfe_window_cb(
            fc.ptr, cb_fnptr, pointer_from_objref(chan)
        )
        return nothing
    end
    return Channel{ResultWindowMFE}(mfe_window_channel_writer, 1000; spawn=true)
end

function mfe_window_channel(seq::AbstractString; window_size::Int)
    fc = FoldCompound(seq; window_size, options=[:default, :window])
    return mfe_window_channel(fc)
end

# neighbors: move sets

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
    return structure, RTlogZ * Private.unit_energy
end

function partfn(sequence::AbstractString)
    fc = FoldCompound(sequence)
    res = partfn(fc)
    finalize(fc)
    return res
end

# plotting secondary structures

"""
    plot_coords(structure; [plot_type])

Plots a secondary structure, returning the coordinates. The optional
`plot_type` parameter can be one of
$(keys(Private.PLOT_COORDS_PLOT_TYPE)), with the default being
`:default`.
"""
function plot_coords(structure::Union{AbstractString,Pairtable};
                     plot_type::Symbol = :default)
    _vrna_plot_coords(structure::AbstractString, cx, cy, type) =
        LibRNA.vrna_plot_coords(structure, cx, cy, type)
    _vrna_plot_coords(structure::Pairtable, cx, cy, type) =
        LibRNA.vrna_plot_coords_pt(structure.ptr, cx, cy, type)
    has_zero_len_hairpin(structure::AbstractString) = contains(structure, "()")
    has_zero_len_hairpin(pt::Pairtable) = any(i -> pt[i] == i + 1, 1:length(pt)-1)

    if length(structure) == 0
        return Float32[], Float32[]
    end
    if !haskey(Private.PLOT_COORDS_PLOT_TYPE, plot_type)
        throw(ArgumentError("unknown plot_type $plot_type, options are: " *
            "$(keys(Private.PLOT_COORDS_PLOT_TYPE))"))
    end

    if plot_type ∈ (:default, :turtle, :puzzler) && has_zero_len_hairpin(structure)
        # LibRNA.vrna_plot_coords crashes/hangs in this case
        throw(ArgumentError("Cannot handle zero-length hairpins for plot_type $plot_type"))
    end

    type = Private.PLOT_COORDS_PLOT_TYPE[plot_type]
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
    res = prob_of_structure(fc, structure)
    finalize(fc)
    return res
end

# sample structures with stochastic backtracking

"""
    sample_structures(fc; [num_samples, options])
    sample_structures(sequence; [num_samples, options])

Sample `num_samples` secondary structures according to their Boltzmann
probabilities.  Possible values for `options` are
$(keys(Private.SAMPLE_STRUCTURES_OPTIONS)).
"""
function sample_structures(fc::FoldCompound;
                           num_samples::Integer=10,
                           options::Symbol=:default)
    num_samples ≥ 0 || throw(ArgumentError("num_samples must be ≥ 0"))
    fc.has_exp_matrices ||
        throw(ArgumentError("must call ViennaRNA.partfn(::FoldCompound) first"))
    fc.uniq_ML ||
        throw(ArgumentError("must have fc.uniq_ML == true"))
    if !haskey(Private.SAMPLE_STRUCTURES_OPTIONS, options)
        throw(ArgumentError("unknown option: $options, possible values are " *
            "$(keys(Private.SAMPLE_STRUCTURES_OPTIONS))"))
    end
    vrna_options = Private.SAMPLE_STRUCTURES_OPTIONS[options]
    s = LibRNA.vrna_pbacktrack_num(fc.ptr, num_samples, vrna_options)
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

function sample_structures(sequence::AbstractString;
                           num_samples::Integer=10,
                           options::Symbol=:default)
    fc = FoldCompound(sequence; uniq_ML=true)
    partfn(fc)
    res = try
        sample_structures(fc; num_samples, options)
    finally
        finalize(fc)
    end
    return res
end

# soft constraints
# modified bases via soft constraints functions (vrna_sc_mod_*)

# sc_mod_*! functions for preset (built-in) energy params of modified
# bases
for modfn_data in Private.SC_MOD_PRESET_FUNCTIONS
    short = modfn_data.short
    desc = modfn_data.desc
    func_name = Symbol("sc_mod_$(short)!")
    vrna_func_name = Symbol("vrna_sc_mod_$short")
    @eval begin
        @doc """
            $($func_name)(fc, modification_sites::Vector{<:Integer})

        Apply base modifications for $($desc) at the
        positions given by `modification_sites` via soft constraint callbacks.
        """
        function $func_name(fc::FoldCompound, modification_sites::Vector{T};
                            option_flags=LibRNA.VRNA_SC_MOD_DEFAULT) where {T<:Integer}
            # Note: ViennaRNA expects 1-based indexing for modification_sites
            nsites = length(modification_sites)
            sites = convert(Vector{Cuint}, modification_sites)
            push!(sites, 0)  # array needs to be terminated by 0
            nsites_mod = LibRNA.$vrna_func_name(fc.ptr, sites, option_flags)
            if nsites_mod != nsites
                throw(ArgumentError("Modified $nsites_mod positions but expected $nsites"))
            end
            return fc
        end
    end
end

# suboptimal structures

"""
    subopt(fc; [delta], [sorted])
    subopt(sequence; [delta], [sorted])

All suboptimal structures with an energy not more than `delta` above
the minimum free energy structure.
"""
function subopt(fc::FoldCompound; delta::Quantity, sorted::Bool=true)
    fc.uniq_ML == true || throw(ArgumentError("FoldCompound must use uniq_ML=true"))
    delta_nounit = uconvert(u"kcal/mol", delta) / Private.unit_energy_int
    delta_int = round(Int, delta_nounit)
    # TODO: warn if this float -> int conversion leads to loss of
    # precision, e.g. in the case that delta == 0.0001 kcal/mol where
    # we would have delta_int == 0
    subs = Tuple{String,Quantity}[]
    s = LibRNA.vrna_subopt(fc.ptr, delta_int, sorted, C_NULL)
    s == C_NULL && return subs
    i = 1
    while (si = unsafe_load(s, i)).structure != C_NULL
        push!(subs, (unsafe_string(si.structure), si.energy * Private.unit_energy))
        i += 1
    end
    LibRNA_Helper.free_subopt_solutions(s)
    return subs
end

function subopt(sequence::AbstractString; delta::Quantity, sorted::Bool=true)
    fc = FoldCompound(sequence; uniq_ML=true, options=[:mfe])
    res = subopt(fc; delta, sorted)
    finalize(fc)
    return res
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
        push!(subs, (unsafe_string(si.structure), si.energy * Private.unit_energy))
        i += 1
    end
    LibRNA_Helper.free_subopt_solutions(s)
    return subs
end

function subopt_zuker(sequence::AbstractString)
    fc = FoldCompound(sequence; options=[:mfe])
    res = subopt_zuker(fc)
    finalize(fc)
    return res
end

# tree edit distance

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
