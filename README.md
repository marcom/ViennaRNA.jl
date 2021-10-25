# Julia interface to the ViennaRNA library

Unofficial Julia interface to the
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) library for RNA
secondary structure prediction and analysis.

## Usage

```julia
using ViennaRNA

# sequence to be folded
# uniq_ML=1 (unique multiloop decomposition) is needed for some functions,
# e.g. pbacktrack
fc = FoldCompound("GGGGGAAAAACCCCCC"; uniq_ML=1)

# minimum free energy structure (MFE) for a sequence
ViennaRNA.mfe(fc)                         # => ("(((((.....))))).", -9.4f0)

# partition function
ViennaRNA.partfn(fc)                      # => ("(((((.....})))),", -9.816722f0)

# calculate energy for structure
ViennaRNA.energy(fc, "((((.......)))).")  # => -6.2f0

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

# suboptimal structures, with energies delta * 0.01 kcal/mol above the
# mfe structure
ViennaRNA.subopt(fc; delta=400, sorted=true)  # => Vector{Tuple{String, Float32}}

# basepair distance between secondary structures
ViennaRNA.bp_distance("....", "(())")     # => 2

# mean basepair distance of all structures to each other,
# weighted by the structure's Boltzmann probabilities
ViennaRNA.mean_bp_distance(fc)            # => 5.266430215905888

# plot coordinates of a secondary structure, returns two arrays with
# x and y coordinates
ViennaRNA.plot_coords("(((...)))")        # => Tuple{Float32[], Float32[]}

# inverse folding / sequence design
ViennaRNA.inverse_fold("AAAAAAA", "((...))")  # => ("GCAAAGC", 2.0f0)
```
