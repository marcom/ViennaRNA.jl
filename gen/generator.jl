using Clang.Generators
using ViennaRNA_jll

cd(@__DIR__)

include_dir = normpath(ViennaRNA_jll.artifact_dir, "include")
viennarna_dir = joinpath(include_dir, "ViennaRNA")

options = load_options(joinpath(@__DIR__, "generator.toml"))

args = get_default_args()
append!(args, [
    "-I$include_dir",
    "-DVRNA_DISABLE_C11_FEATURES",
])

accept_headers = [
    "equilibrium_probs.h",
    "eval.h",
    "fold_compound.h",
    "inverse.h",
    "landscape/neighbor.h",
    "MEA.h",
    "mfe.h",
    "params/io.h",
    "part_func.h",
    "plotting/layouts.h",
    "subopt.h",
    "subopt_zuker.h",
    "treedist.h",
    "utils/basic.h",
]
headers = [joinpath(viennarna_dir, header) for header in accept_headers]
# there is also an experimental `detect_headers` function for auto-detecting top-level headers in the directory
# headers = detect_headers(clang_dir, args)

# create context
ctx = create_context(headers, args, options)

# run generator
build!(ctx)
