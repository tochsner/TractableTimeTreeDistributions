struct IndependentDist{D} <: AbstractDistribution
    distributions::Dict{Clade,D}
end

function IndependentDist{D}(trees::Vector{CladifiedTree}) where D
    observations = DefaultDict{Clade,Vector{Float64}}([])

    for tree in trees
        for (clade, parameter) in tree.parameters
            push!(observations[clade], parameter)
        end
    end

    global_distribution = fit_mle(D, reduce(vcat, values(observations)))

    distributions::Dict{Clade,D} = Dict()
    for (clade, clade_observations) in observations
        if length(clade_observations) < 5
            distributions[clade] = global_distribution
        else
            try
                distributions[clade] = fit_mle(D, clade_observations)
            catch
                @info "MLE not working"
                distributions[clade] = fit(D, clade_observations)
            end
        end
    end

    return IndependentDist{D}(distributions)
end

function readable_name(distribution::Type{IndependentDist{D}}) where D
    readable_name(D)
end

function sample_tree(distribution::IndependentDist{D}, tree::Union{CladifiedTree,CladifiedTree})::Union{CladifiedTree,CladifiedTree} where D
    parameters = Dict(
        clade => rand(dist) for (clade, dist) in distribution.distributions
    )
    return CladifiedTree(parameters, tree)
end

function log_density(distribution::IndependentDist{D}, tree::CladifiedTree) where D
    sum(logpdf(distribution.distributions[clade], param) for (clade, param) in tree.parameters)
end