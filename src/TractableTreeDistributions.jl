module TractableTreeDistributions

using StatsBase
using DataStructures
using Memoization
using Distributions

# utils

include("utils/trees.jl")

export load_trees, load_trees_from_newick

# CCD

include("ccd/clade.jl")
include("ccd/clade_split.jl")
include("ccd/cladify_tree.jl")
include("ccd/ccd1.jl")

export AbstractClade, Leaf, Clade, Split, TreeWithClades, cladify_tree
export CCD1, fit, sample, log_density, get_most_likely_tree

export construct_tree

include("utils/construct_tree.jl")

# Transformations

include("transformations/parameterized_tree.jl")
include("transformations/transform.jl")
include("transformations/height_ratio_transform.jl")
include("transformations/independent_dist.jl")

export ParameterizedTree
export HeightRatioTransform, transform_height, transform_ratios, invert_height, invert_ratios
export IndependentDist
export transform, invert, fit!, sample, log_density

end # module TractableTreeDistributions
