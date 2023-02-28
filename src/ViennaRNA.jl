module ViennaRNA

import Base
using Unitful: @u_str, Quantity, uconvert, ustrip
using UnsafePointers: UnsafePtr

# types
export FoldCompound, Pairtable

# wrapper functions
export bp_distance, bpp, centroid, energy, ensemble_defect,
    heat_capacity, inverse_fold, inverse_pf_fold, mea,
    mean_bp_distance, mfe, mfe_window, mfe_window_channel, neighbors,
    partfn, plot_coords, prob_of_structure, sample_structures, subopt,
    subopt_zuker, tree_edit_dist

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

function __init__()
    # otherwise the name of the default params is not set
    # i.e., we would have FoldCompound("GAC").params_name == ""
    # time cost is approx. 5ms
    params_load_defaults()
end

end # module ViennaRNA
