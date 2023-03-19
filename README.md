# Julia interface to the ViennaRNA library

Unofficial Julia interface to the
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) library for RNA
secondary structure prediction and analysis.  Please cite the original
ViennaRNA publications if you use this library.

## Installation

Install ViennaRNA from the Julia package REPL, which can be accessed
by pressing `]` from the Julia REPL:

```
add ViennaRNA
```

## Usage

```julia
using ViennaRNA, Unitful
```

The `Unitful` library is needed to be able to specify units with
`@u_str`, e.g. `4u"kcal/mol"` or `37u"°C"`.  You can get the degree
symbol `°` by typing `\degree` and pressing the TAB key in the REPL or
in an editor with Julia syntax support.

The original C API functions can be found in the submodule
`ViennaRNA.LibRNA`.  Most functions can be called with a String
containing the RNA sequence instead of a `FoldCompound`,
e.g. `mfe("GGGAAACCC")`.

### FoldCompound

A `FoldCompound` encapsulates nucleic acid strands and model details,
such as energy parameters, temperature, etc.

```julia
fc = FoldCompound("GGGGGAAAAACCCCCC";
                  options=[:mfe, :pf],
                  temperature=37u"°C",
                  uniq_ML=true,
                  circular=false)
```

Important keyword arguments

- `options` is a subset of `[:default, :eval_only, :hybrid, :mfe,
  :pf, :window]`.

- `temperature` is used to rescale the free energies with the formula
  `ΔG = ΔH - TΔS` (the energy parameter sets contain enthalpy and
  entropy contributions). The default is `37u"°C"`

Model details (additional keyword arguments):
- `circular`: determines if the RNA strand is circular, i.e. the
  5'-end and 3'-end are covalently bonded. Default is `false`.
- `dangles`: how to treat dangling base pairs in multiloops and the
  exterior loop. Can be 0, 1, 2, or 3. See ViennaRNA docs for
  details. Default is `2`.
- `gquadruplex`: allow G-quadruplexes in predictions. Default is
  `false`.
- `log_ML`: use logarithmic energy model for multiloops. Default is
  `false`.
- `max_bp_span`: maximum number of bases over which a basepair can
  span. Default value is `-1` (which means unlimited).
- `min_loop_length`: the minimum size of a loop (without the closing
   base pair). Default is `3`.
- `no_GU_basepairs`: disallow G-U basepairs. Default is `false`.
- `no_GU_closure`: disallow G-U basepairs as closing pairs for
  loops. Default is `false`.
- `no_lonely_pairs`: disallow isolated base pairs. Default is `false`.
- `special_hairpins`: use special hairpin energies for certain tri-,
  tetra- and hexloops. Default is `true`.
- `uniq_ML`: use unique decomposition for multiloops, needed for
  `sample_structures` and `subopt`. Default is `false`.
- `window_size`: window size to be used for local calculations
  performed in a window moving over the sequence. This value is
  ignored unless the `:window` option is set in the `FoldCompound`
  `options`. The default value for `window_size` is `-1`.

#### Changing the energy parameter set

ViennaRNA stores energy parameters in global variables after loading
them from a file.  Each time a new `FoldCompound` is created, the
parameters are copied from the global variables and saved inside the
`FoldCompound`.

The global variables storing energy parameters can be changed by
calling a function specific to each parameter set, or via a Symbol
with `ViennaRNA.params_load(:RNA_Turner1999)`. Subsequent calls to
`FoldCompound` will use the new parameters and store a copy of the
parameters in the newly created `FoldCompound`.

The default energy set loaded on startup is `:RNA_Turner2004`.

```julia
ViennaRNA.params_load_defaults()  # default is :RNA_Turner2004
ViennaRNA.params_load_DNA_Mathews1999()
ViennaRNA.params_load_DNA_Mathews2004()
ViennaRNA.params_load_RNA_Andronescu2007()
ViennaRNA.params_load_RNA_Langdon2018()
ViennaRNA.params_load_RNA_Turner1999()
ViennaRNA.params_load_RNA_Turner2004()
ViennaRNA.params_load_RNA_misc_special_hairpins()

ViennaRNA.params_load(:DNA_Mathews2004)
# Options are
#   :DNA_Mathews1999, :DNA_Mathews2004,
#   :RNA_Andronescu2007, :RNA_Langdon2018,
#   :RNA_Turner1999, :RNA_Turner2004
```

#### Multiple strands

Mutiple strands can be given by separating them with an `&`, e.g.
`FoldCompound("GGGG&CCCC")`.

#### Comparative folding with an MSA (alifold)

Pass multiple sequences to the FoldCompound constructor for comparative mode (alifold):
`FoldCompound(["GG-GAAAACCCC", "GCCGAAA-CGGC"])`.

It is currently not possible to have multiple strands in alifold mode.

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

### Sample structures (probabilistic / stochastic backtrack)

Sample from Boltzmann ensemble of secondary structures

```julia
sample_structures(fc)                         # => [ "((((......)))).." ]
sample_structures(fc; options=:nonredundant,
                      num_samples=20)         # => 20-element Vector{String}
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

### Sliding window prediction of MFE

`mfe_window` saves all the results in a `Vector`.

```julia
seq = "G"^50 * "A"^4 * "C"^50
mfe_window(seq; window_size=30)
fc = FoldCompound(seq; options=[:default, :window], window_size=30)
mfe_window(fc)  # => Vector{ResultWindowMFE}
```

`mfe_window_channel` returns a `Channel` that can be used to
iteratively process the results.

```julia
seq = "G"^50 * "A"^4 * "C"^50
chan = mfe_window_channel(seq; window_size=30)
take!(chan)
fc = FoldCompound(seq; options=[:default, :window], window_size=30)
chan = mfe_window_channel(fc)
take!(chan)
```

### Neighboring structures

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

### Maximum expected accuracy (MEA) structure

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
# plot coordinates of a secondary structure, returns two arrays with
# x and y coordinates
plot_coords("(((...)))")  # => Tuple{Float32[], Float32[]}
```

See [PlotRNA.jl](https://github.com/marcom/PlotRNA.jl) for more
secondary structure plotting functionality.

### Inverse folding / sequence design

```julia
inverse_fold("AAAAAAA", "((...))")     a# => ("GCAAAGC", 2.0f0)
inverse_pf_fold("AAAAAAA", "((...))")  # => ("GCCAAGC", 2.0244526863098145 kcal mol^-1)
```

### Seeding the random number generator

```julia
ViennaRNA.init_rand_seed(42)
```

## Reducing memory usage

When creating many `FoldCompound`s, running `finalize` manually will
avoid excessive memory buildup.

```julia
for i = 1:100_000
    fc = FoldCompound("ACGU")
    # do something with fc
    finalize(fc)
end
```


## Related Julia packages

- [RNAstructure.jl](https://github.com/marcom/RNAstructure.jl)
- [LinearFold.jl](https://github.com/marcom/LinearFold.jl)
- [PlotRNA.jl](https://github.com/marcom/PlotRNA.jl)
