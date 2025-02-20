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


height_ratio_transform(
    independent_normal |> logtransform, 
    independent_beta
)

height_ratio_transform(
    independent_gamma,
    independent_beta
)

height_ratio_transform(
    independent_normal |> logtransform,
    multivariate_normal |> logit_transform
)

height_ratio_transform(
    independent_gamma,
    multivariate_normal |> logit_transform
)

independent_evolution_transform(
    independent_normal |> logtransform, 
    multivariate_normal |> logtransform
)

independent_evolution_transform(
    independent_gamma,
    multivariate_normal |> logtransform
)

independent_evolution_transform(
    independent_normal |> logtransform, 
    independent_gamma
)

independent_evolution_transform(
    independent_gamma,
    independent_gamma
)

independent_normal() |> logtransform |> short_branch_transform
multivariate_normal() |> logtransform |> short_branch_transform
independent_gamma() |> short_branch_transform

- fit
- sample
- density

fit(transform, trees::Vector{ParameterizedTree})
sample(transform)::ParameterizedTree
density(transform, tree::ParameterizedTree)

transform(transform, tree::ParameterizedTree)::ParameterizedTree
invert(transform, tree::ParameterizedTree)::ParameterizedTree

create()