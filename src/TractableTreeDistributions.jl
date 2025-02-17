module TractableTreeDistributions

# utils

include("utils/trees.jl")

export load_trees

# models

include("models/clade.jl")
include("models/ccd1.jl")

export Clade
export CCD1, fit, sample, get_probability

end # module TractableTreeDistributions
