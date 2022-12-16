module ViennaRNA

import Base
using Unitful: @u_str, Quantity, uconvert, ustrip
using UnsafePointers: UnsafePtr

# types
export FoldCompound, Pairtable

# wrapper functions
export bp_distance, bpp, centroid, energy, ensemble_defect,
    heat_capacity, inverse_fold, inverse_pf_fold, mea,
    mean_bp_distance, mfe, mfe_window, neighbors, partfn, plot_coords,
    prob_of_structure, sample_structures, subopt, subopt_zuker,
    tree_edit_dist

# additional utility functions
export basepairs, prob_of_basepairs

# generated C bidings for libRNA from ViennaRNA
# see src/gen/generator.jl
include("../lib/LibRNA.jl")
import .LibRNA

# helper functions for LibRNA
include("librna-helpers.jl")
import .LibRNA_Helper

include("wrappers.jl")
include("extras.jl")
include("utils.jl")

end # module ViennaRNA
