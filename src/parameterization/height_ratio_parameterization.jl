struct HeightRatioParameterization <: Parameterization end

function get_embeddings(::HeightRatioParameterization, trees::Vector{Vector{AbstractClade}})::Vector{Dict{AbstractClade,Float64}}

end

function get_trees(::HeightRatioParameterization, embeddings::Vector{Dict{AbstractClade,Float64}})::Vector{Vector{AbstractClade}}

end