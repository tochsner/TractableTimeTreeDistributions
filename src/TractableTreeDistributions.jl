module TractableTreeDistributions

using StatsBase
using DataStructures
using Memoization
using Distributions

# utils

include("utils/trees.jl")

export load_trees, load_trees_from_newick, write_tree

# Clades

include("clades/clade.jl")
include("clades/clade_split.jl")
include("clades/cladify_tree.jl")

export AbstractClade, Leaf, Clade, Split, TreeWithClades, cladify_tree
export CCD1, fit, log_density, most_likely_tree

export construct_tree

include("utils/construct_tree.jl")

# Distributions

include("distributions/distribution.jl")
include("distributions/short_branch_dist.jl")
include("distributions/height_ratio_dist.jl")
include("distributions/independent_dist.jl")
include("distributions/ccd1.jl")
include("distributions/tractable_time_tree_dist.jl")

export CladifiedTree, TractableTimeTreeDist
export LastDivergenceBranchDist, transform_last_div, transform_branches
export HeightRatioDist, transform_height, transform_ratios, invert_height, invert_ratios, transform_short_branches, invert_short_branches, ShorterBranchDist
export IndependentDist, sample_tree

end # module TractableTreeDistributions
