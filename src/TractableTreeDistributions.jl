module TractableTreeDistributions

# utils

include("utils/trees.jl")

export load_trees, load_trees_from_newick

# models

include("models/clade.jl")
include("models/clade_split.jl")
include("models/cladify_tree.jl")
include("models/ccd1.jl")

export Clade, CladeSplit, TreeWithClades, cladify_tree
export CCD1, fit, sample, get_log_probability, get_most_likely_tree

export construct_tree

include("utils/construct_tree.jl")

end # module TractableTreeDistributions
