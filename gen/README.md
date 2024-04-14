# ViennaRNA libRNA wrapper generator

## How to generate the wrappers

1. `cd` to this directory
2. run `julia --project generator.jl`

## How to upgrade Clang.jl

1. `cd` to this directory
2. run `julia --project` and then in the Julia REPL, run `pkg> up`

## Notes on patches to headers

Some headers of ViennaRNA have to be patched to be able to be parsed
by Clang.jl.  These patches can be found are under the directory
`patches-for-headers/`.  The patches are applied to a temporary copy
of the header files that is deleted after parsing by the
`generator.jl` script that generates the `../lib/LibRNA.jl` file.
This only needs to be done once on a developer's machine when the
header files change.
