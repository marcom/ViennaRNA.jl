# CHANGELOG

## 0.9.0

- rename `pbacktrack` to `sample_structures`, the same name used
  in RNAstructure.jl and LinearFold.jl

- plot_coords: restrict input types, new `:default` plot type that
  uses ViennaRNA's default plot type

- regenerated lib/LibRNA.jl with gen/generator.jl
    - now with ViennaRNA_jll 2.5.1 and not 2.5.0 (this was a mistake in 0.8.x,
      but it still worked as the differences in the headers were
      small)
    - always keep the version of ViennaRNA_jll used in gen/ in sync
      with the parent dir by using stacked pkg environments

