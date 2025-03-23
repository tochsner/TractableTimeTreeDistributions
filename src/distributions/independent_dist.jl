struct IndependentDist{D} <: AbstractDistribution
    distributions::Dict{Clade,D}
end

function IndependentDist{D}(trees::Vector{CladifiedTree}) where {D}
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
            continue
        end

        try
            distributions[clade] = fit_mle(D, clade_observations)
            continue
        catch
        end

        @warn "Numerical MLE threw an error, we fall back to a heruristic fit"

        try
            distributions[clade] = fit(D, clade_observations)
            continue
        catch
        end

        @warn "Numerical heuristic threw an error as well, we fall back to the global distribution"
        distributions[clade] = global_distribution
    end

    return IndependentDist{D}(distributions)
end

readable_name(::Type{IndependentDist{D}}) where {D} = readable_name(D)

function sample_tree(distribution::IndependentDist{D}, tree::CladifiedTree)::CladifiedTree where {D}
    parameters = Dict(clade => rand(dist) for (clade, dist) in distribution.distributions)
    return CladifiedTree(parameters, tree)
end

function point_estimate(distribution::IndependentDist{D}, tree::CladifiedTree)::CladifiedTree where {D}
    parameters = Dict(
        clade => mean(distribution.distributions[clade]) for
        clade in keys(tree.splits) if haskey(distribution.distributions, clade)
    )
    return CladifiedTree(parameters, tree)
end

function log_density(distribution::IndependentDist{D}, tree::CladifiedTree) where {D}
    sum(
        haskey(distribution.distributions, clade) ? logpdf(distribution.distributions[clade], param) : 0.0 for
        (clade, param) in tree.parameters
    )
end
