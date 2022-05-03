# Julia interface to the ViennaRNA library

Unofficial Julia interface to the
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) library for RNA
secondary structure prediction and analysis.  Please cite the original
ViennaRNA publications if you use this library.

## Usage

```julia
using ViennaRNA, Unitful
# Unitful is needed to be able to specify units with @u_str,
# e.g. 4u"kcal/mol" or 37u"°C".  You can get the degree symbol ° by
# typing \degree and pressing the TAB key in the REPL or in an editor
# with Julia syntax support

# Notes
# The original C API functions can be found in the submodule `ViennaRNA.LibRNA`.
# Most functions can be called with a String containing the RNA
# sequence instead of a FoldCompound, e.g. `mfe("GGGAAACCC")`.

# sequence to be folded
# params, temperature, uniq_ML are optional keyword arguments
# - params sets the energy parameter set used, options are
#   :RNA_Turner1999, :RNA_Turner2004, :RNA_Andronescu2007,
#   :RNA_Langdon2018. The default is :RNA_Turner2004
# - temperature is used to rescale the free energies with the formula
#   ΔG = ΔH - TΔS (the energy parameter sets contain enthalpy and
#   entropy contributions). The default is 37u"°C"
# - uniq_ML=true (unique multiloop decomposition) is needed for some
#   functions, e.g. pbacktrack. The default is false
fc = FoldCompound("GGGGGAAAAACCCCCC";
		  params=:RNA_Turner2004,
		  temperature=37u"°C",
		  uniq_ML=true)

# minimum free energy structure (MFE) for a sequence
# please excuse the excess precision printed when displaying -9.4 kcal/mol
mfe(fc)                         # => ("(((((.....))))).", -9.399999618530273 kcal mol^-1)

# partition function
partfn(fc)                      # => ("(((((.....})))),", -9.81672180213034 kcal mol^-1)

# calculate energy for structure
energy(fc, "((((.......)))).")  # => -6.199999809265137 kcal mol^-1

# basepair probabilities
bpp(fc)                         # => 16×16 Matrix{Float64}

# Boltzmann probability of a structure
prob_of_structure(fc, "(((((.....))))).")  # => 0.5085737925408758

# ensemble defect
ensemble_defect(fc, "(((((.....))))).")  # => 0.33085374128228884

# probabilistic / stochastic backtrack, sample from Boltzmann ensemble of
# secondary structures
pbacktrack(fc)                  # => [ "((((......)))).." ]
pbacktrack(fc; num_samples=10)  # => 10-element Vector{String}

# all suboptimal structures with energies delta above the MFE
# structure
subopt(fc; delta=4u"kcal/mol")  # => Vector{Tuple{String, Quantity}}

# suboptimal structures with the method of Zuker
subopt_zuker(fc)                # => Vector{Tuple{String, Quantity}}

# move set to reach neighboring structures of a given structure
neighbors(fc, Pairtable("((.....))"))  # => Vector{Vector{Tuple{Int,Int}}}

# basepair distance between secondary structures
bp_distance("....", "(())")     # => 2

# tree edit distance between secondary structures
tree_edit_dist("(..)", "....")  # => 4.0f0

# mean basepair distance of all structures to each other,
# weighted by the structure's Boltzmann probabilities
mean_bp_distance(fc)            # => 5.266430215905888

# centroid structure of ensemble: structure with smallest sum of
# base-pair distances weighted by Boltzmann probabilities
centroid(fc)                    # => ("(((((.....))))).", 4.799131457924728)

# Maximum expected accuracy (MEA) structure. The gamma parameter
# trades off specificity (low gamma) and sensitivity (high gamma).
mea(fc; gamma=1.0)              # => ("(((((.....))))).", 10.706348f0)

# plot coordinates of a secondary structure, returns two arrays with
# x and y coordinates
plot_coords("(((...)))")        # => Tuple{Float32[], Float32[]}

# inverse folding / sequence design
inverse_fold("AAAAAAA", "((...))")    # => ("GCAAAGC", 2.0f0)
inverse_pf_fold("AAAAAAA", "((...))") # => ("GCCAAGC", 2.0244526863098145 kcal mol^-1)

# plot secondary structure
plot_structure("((...))"; sequence = "AAAAAAA")
plot_structure("((...))"; sequence = "AAAAAAA", targetdir = "2dstructure.pdf")
plot_structure("((...))"; sequence = "AAAAAAA", targetdir = "2dstructure.png")

```
