struct IndependentGamma <: Distribution end

struct EstimatedIndependentGamma <: EstimatedDistribution
    shapes::Dict{AbstractClade,Float64}
    scales::Dict{AbstractClade,Float64}
end

function estimate(distribution::IndependentGamma, embeddings::Vector{Dict{AbstractClade,Float64}})::EstimatedIndependentGamma
end

function log_density(estimated_distribution::EstimatedIndependentGamma, embeddings::Vector{Dict{AbstractClade,Float64}})::Float64
end

function sample(estimated_distribution::EstimatedIndependentGamma)::Vector{AbstractClade}

end
