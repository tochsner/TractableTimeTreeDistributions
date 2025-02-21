struct IndependentDist{D} <: AbstractDistribution
    distributions::Dict{Clade,D}
end

function IndependentDist{D}(trees::Vector{ParameterizedTree}) where D
    observations = DefaultDict{Clade,Vector{Float64}}([])

    for tree in trees
        for (clade, parameter) in tree.parameters
            push!(observations[clade], parameter)
        end
    end

    global_distribution = fit_mle(D, reduce(vcat, values(observations)))

    distributions::Dict{Clade,D} = Dict()
    for (clade, clade_observations) in observations
        if length(clade_observations) < 3
            distributions[clade] = global_distribution
        else
            distributions[clade] = fit_mle(D, clade_observations)
        end
    end

    return IndependentDist{D}(distributions)
end

function sample_tree(distribution::IndependentDist{D}, tree::Union{CladifiedTree,ParameterizedTree})::Union{CladifiedTree,ParameterizedTree} where D
    parameters = Dict(
        clade => rand(dist) for (clade, dist) in distribution.distributions
    )
    return ParameterizedTree(parameters, tree)
end

function log_density(distribution::IndependentDist{D}, tree::ParameterizedTree) where D
    sum(logpdf(distribution.distributions[clade], param) for (clade, param) in tree.parameters)
end