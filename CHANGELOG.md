# CHANGELOG

## 0.10.0

- Energy parameter sets are now changed by calling functions like
  `ViennaRNA.params_load_RNA_Turner1999()` or
  `ViennaRNA.params_load(:RNA_Turner1999)` that set global
  variables. This more closely mirrors the ViennaRNA API and avoids
  the costly reloading of energy parameters inside the
  `FoldCompound()` constructor.  The global variables are copied by
  ViennaRNA to a newly created `FoldCompound` when the
  `FoldCompound()` constructor is called.  Previously, the kwarg
  `param` when creating a `FoldCompound` was used to set the parameter
  set.


## 0.9.0

- `FoldCompound`: new properties `max_bp_span` and `window_size` that
  can be set in constructor

- added `mfe_window()` to calculate mfe substructures for a
  sliding window

- added `mfe_window_channel()` to process results from `mfe_window`
  through a `Channel`

- added `ViennaRNA.init_rand_seed(seedval)` to seed the random
  number generator used by ViennaRNA

- rename `pbacktrack` to `sample_structures`, the same name used in
  RNAstructure.jl and LinearFold.jl

- `sample_structures` (formerly known as `pbacktrack`)
  - kwarg `num_samples` now defaults to 10 (previously: 1)
  - kwarg `options` must now be a `Symbol`, currently either `:default`
    or `:nonredundant`

- `plot_coords`
  - restrict input types to `Union{AbstractString,Pairtable}`
  - new `:default` plot type that uses ViennaRNA's default plot type,
    which currently is 'puzzler'. This means the default plot type of
    plot_coords has changed from `:simple` to `:puzzler`

- regenerated lib/LibRNA.jl with gen/generator.jl
  - now with ViennaRNA_jll 2.5.1 and not 2.5.0 (this was a mistake in
    0.8.x, but it still worked as the differences in the headers were
    small)
  - always keep the version of ViennaRNA_jll used in gen/ in sync with
    the parent dir by using stacked pkg environments
  - include ViennaRNA `mfe_window.h` and `part_func_window.h` headers
