struct IndependentDist{D} <: Transform
    distributions::Dict{Clade,D}
    IndependentDist{D}() where D = new(Dict())
end

function fit!(transformation::IndependentDist{D}, trees::Vector{ParameterizedTree}) where D
    empty!(transformation.distributions)

    observations = DefaultDict{Clade,Vector{Float64}}([])

    for tree in trees
        for (clade, parameter) in tree.parameters
            push!(observations[clade], parameter)
        end
    end

    global_distribution = fit_mle(D, reduce(vcat, values(observations)))

    for (clade, clade_observations) in observations
        if length(clade_observations) < 3
            transformation.distributions[clade] = global_distribution
        else
            transformation.distributions[clade] = fit_mle(D, clade_observations)
        end
    end
end

function sample(transformation::IndependentDist{D}, tree::CladifiedTree)::ParameterizedTree where D
    parameters = Dict(
        clade => rand(dist) for (clade, dist) in transformation.distributions
    )
    return ParameterizedTree(parameters, tree)
end

function log_density(transformation::IndependentDist{D}, tree::ParameterizedTree) where D
    sum(logpdf(transformation.distributions[clade], param) for (clade, param) in tree.parameters)
end