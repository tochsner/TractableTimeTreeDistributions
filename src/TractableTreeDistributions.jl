module TractableTimeTreeDistributions

using StatsBase
using DataStructures
using Memoization
using Distributions
using SpecialFunctions

# utils

include("utils/trees.jl")

export load_trees, load_trees_from_newick, write_tree

# Clades

include("clades/clade.jl")
include("clades/clade_split.jl")
include("clades/cladify_tree.jl")

export AbstractClade, Leaf, Clade, Split, TreeWithClades, cladify_tree
export CCD1, CCD0, fit, log_density, point_estimate

export construct_tree

include("utils/construct_tree.jl")

include("utils/rf_distance.jl")
include("utils/tree_ess.jl")

# Distributions

include("distributions/distribution.jl")
include("distributions/short_branch_dist.jl")
include("distributions/height_ratio_dist.jl")
include("distributions/independent_dist.jl")
include("distributions/dirichlet_dist.jl")
include("distributions/ccd.jl")
include("distributions/ccd0.jl")
include("distributions/ccd1.jl")
include("distributions/tractable_time_tree_dist.jl")
include("distributions/last_divergence_branch_dist.jl")
include("distributions/mrca.jl")
include("distributions/numerical_distributions_utils.jl")

export CladifiedTree, TractableTimeTreeDist, TreeDirichletDist
export LastDivergenceBranchDist, transform_last_divergence, transform_branches, invert_last_divergence_branches
export HeightRatioDist, transform_height, transform_ratios, invert_height, invert_ratios, transform_short_branches, invert_short_branches, ShorterBranchDist
export IndependentDist, sample_tree, MultivariateDist, get_entropy

export mrca_tree

export readable_name
export robinson_foulds_distance, get_ess

end # module TractableTimeTreeDistributions
