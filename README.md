# Julia interface to the ViennaRNA library

Unofficial Julia interface to the
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) library for RNA
secondary structure prediction and analysis.  Please cite the original
ViennaRNA publications if you use this library.

## Usage

```julia
using ViennaRNA

# Notes
# The original C API functions can be found in the submodule `ViennaRNA.LibRNA`.
# Most functions can be called with a String containing the RNA sequence,
# e.g. `ViennaRNA.mfe("GGGAAACCC")`.

# sequence to be folded
# uniq_ML=1 (unique multiloop decomposition) is needed for some functions,
# e.g. pbacktrack
fc = FoldCompound("GGGGGAAAAACCCCCC"; uniq_ML=1)

# minimum free energy structure (MFE) for a sequence
# please excuse the excess precision printed when displaying -9.4 kcal/mol
ViennaRNA.mfe(fc)                         # => ("(((((.....))))).", -9.399999618530273 kcal mol^-1)

# partition function
ViennaRNA.partfn(fc)                      # => ("(((((.....})))),", -9.81672191619873 kcal mol^-1)

# calculate energy for structure
ViennaRNA.energy(fc, "((((.......)))).")  # => -6.199999809265137 kcal mol^-1

# basepair probabilities
ViennaRNA.bpp(fc)                         # => 16×16 Matrix{Float64}

# Boltzmann probability of a structure
ViennaRNA.prob_of_structure(fc, "(((((.....))))).")  # => 0.5085737925408758

# ensemble defect
ViennaRNA.ensemble_defect(fc, "(((((.....))))).")  # => 0.3035942605397949

# probabilistic / stochastic backtrack, sample from Boltzmann ensemble of
# secondary structures
ViennaRNA.pbacktrack(fc)                  # => [ "((((......)))).." ]
ViennaRNA.pbacktrack(fc; num_samples=10)  # => 10-element Vector{String}

# all suboptimal structures with energies delta above the MFE
# structure
using Unitful    # to be able to write 4u"kcal/mol" with @u_str
ViennaRNA.subopt(fc; delta=4u"kcal/mol")  # => Vector{Tuple{String, Quantity}}

# basepair distance between secondary structures
ViennaRNA.bp_distance("....", "(())")     # => 2

# mean basepair distance of all structures to each other,
# weighted by the structure's Boltzmann probabilities
ViennaRNA.mean_bp_distance(fc)            # => 5.266430215905888

# centroid structure of ensemble: structure with smallest sum of
# base-pair distances weighted by Boltzmann probabilities
ViennaRNA.centroid(fc)                    # => ("(((((.....))))).", 4.799131457924728)

# Maximum expected accuracy (MEA) structure. The gamma parameter
# trades off specificity (low gamma) and sensitivity (high gamma).
ViennaRNA.mea(fc; gamma=1.0)              # => ("(((((.....))))).", 10.706348f0)

# plot coordinates of a secondary structure, returns two arrays with
# x and y coordinates
ViennaRNA.plot_coords("(((...)))")        # => Tuple{Float32[], Float32[]}

# inverse folding / sequence design
ViennaRNA.inverse_fold("AAAAAAA", "((...))")    # => ("GCAAAGC", 2.0f0)
ViennaRNA.inverse_pf_fold("AAAAAAA", "((...))") # => ("GCCAAGC", 2.0244526863098145 kcal mol^-1)
```
