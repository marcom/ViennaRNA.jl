# Julia interface to the ViennaRNA library

Unofficial Julia interface to the
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) library for RNA
secondary structure prediction and analysis.  Please cite the original
ViennaRNA publications if you use this library.

## Usage

```julia
using ViennaRNA, Unitful
```

The `Unitful` library is needed to be able to specify units with
`@u_str`, e.g. `4u"kcal/mol"` or `37u"°C"`.  You can get the degree
symbol `°` by typing `\degree` and pressing the TAB key in the REPL or
in an editor with Julia syntax support

The original C API functions can be found in the submodule
`ViennaRNA.LibRNA`.  Most functions can be called with a String
containing the RNA sequence instead of a `FoldCompound`,
e.g. `mfe("GGGAAACCC")`.

### FoldCompound

A `FoldCompound` encapsulates nucleic acid strands and model details,
such as energy parameters, temperature, etc.

```julia
fc = FoldCompound("GGGGGAAAAACCCCCC";
                  params=:RNA_Turner2004,
                  temperature=37u"°C",
                  uniq_ML=true,
                  circular=false)
```

Optional keyword arguments

- `params` determines the energy parameter set used, options are
  `:RNA_Turner1999`, `:RNA_Turner2004`, `:RNA_Andronescu2007`,
  `:RNA_Langdon2018`. The default is `:RNA_Turner2004`

- `temperature` is used to rescale the free energies with the formula
  `ΔG = ΔH - TΔS` (the energy parameter sets contain enthalpy and
  entropy contributions). The default is `37u"°C"`

- `circular=false` determines if RNA strand is circular

- `uniq_ML=true` (unique multiloop decomposition) is needed for some
  functions, e.g. `pbacktrack`. The default is `false`


Mutiple strands can be given by separating them with an `&`, e.g.
`FoldCompound("GGGG&CCCC")`

### Minimum free energy structure (MFE)

```julia
# please excuse the excess precision printed when displaying -9.4 kcal/mol
mfe(fc)  # => ("(((((.....))))).", -9.399999618530273 kcal mol^-1)
```

### Partition function
```julia
partfn(fc)  # => ("(((((.....})))),", -9.81672180213034 kcal mol^-1)
```

### Free energy change of folding into a structure
```julia
energy(fc, "((((.......)))).")  # => -6.199999809265137 kcal mol^-1
```

### Basepair probabilities
```julia
bpp(fc)  # => 16×16 Matrix{Float64}
```

### Boltzmann probability of a structure
```julia
prob_of_structure(fc, "(((((.....))))).")  # => 0.5085737925408758
```

### Ensemble defect
```julia
ensemble_defect(fc, "(((((.....))))).")  # => 0.33085374128228884
```

### Probabilistic / stochastic backtrack

Sample from Boltzmann ensemble of secondary structures

```julia
pbacktrack(fc)                  # => [ "((((......)))).." ]
pbacktrack(fc; num_samples=10)  # => 10-element Vector{String}
```

### Suboptimal structures

All suboptimal structures with energies delta above the MFE structure

```julia
subopt(fc; delta=4u"kcal/mol")  # => Vector{Tuple{String, Quantity}}
```

Suboptimal structures with the method of Zuker

```julia
subopt_zuker(fc)  # => Vector{Tuple{String, Quantity}}
```

# Neighboring structures

Move set to reach neighboring structures of a given structure

```julia
neighbors(fc, Pairtable("((.....))"))  # => Vector{Vector{Tuple{Int,Int}}}
```

### Basepair distance between secondary structures

```julia
bp_distance("....", "(())")  # => 2
```

### Tree edit distance between secondary structures
```julia
tree_edit_dist("(..)", "....")  # => 4.0f0
```

### Mean basepair distance

Mean basepair distance of all structures to each other, weighted by
the structure's Boltzmann probabilities

```julia
mean_bp_distance(fc)  # => 5.266430215905888
```

### Centroid structure

Centroid structure of ensemble: structure with smallest sum of
base-pair distances weighted by Boltzmann probabilities:

```julia
centroid(fc)  # => ("(((((.....))))).", 4.799131457924728)
```

# Maximum expected accuracy (MEA) structure

The gamma parameter trades off specificity (low gamma) and sensitivity (high gamma).

```julia
mea(fc; gamma=1.0)  # => ("(((((.....))))).", 10.706348f0)
```

### Heat capacity calculation

```julia
# starting temperature, end temperature, temperature increment
heat_capacity(fc, 10u"°C", 60u"°C")  # => Vector{Tuple{Quantity,Quantity}}
heat_capacity(fc, 10u"°C", 60u"°C", 1u"°C"; mpoints=5)
```

### Plotting structures

```julia
# plot_structure (experimental), draws an image of a secondary structure
#
# There is also an experimental plot_structure_makie which looks a bit
# nicer but currently has a rather large time to first plot (half a
# minute). See issue #5. Subsequent plots are very fast.
dbn = "(((...)))"
seq = "GGGAAACCC"
plot_structure(dbn)
plot_structure(dbn; sequence=seq, base_colors=prob_of_basepairs(seq, dbn))
```

```julia
# plot coordinates of a secondary structure, returns two arrays with
# x and y coordinates
plot_coords("(((...)))")  # => Tuple{Float32[], Float32[]}
```

### Inverse folding / sequence design

```julia
inverse_fold("AAAAAAA", "((...))")     # => ("GCAAAGC", 2.0f0)
inverse_pf_fold("AAAAAAA", "((...))")  # => ("GCCAAGC", 2.0244526863098145 kcal mol^-1)
```
