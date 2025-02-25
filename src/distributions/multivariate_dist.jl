struct MultivariateDist{D} <: AbstractDistribution
    clade_indices::Dict{Clade,Int}
    μ::Vector{Float64}
    Σ::Matrix{Float64}
end

function MultivariateDist{D}(trees::Vector{CladifiedTree}) where {D}
    all_clades = Set()
    for tree in trees
        union!(all_clades, keys(tree.parameters))
    end

    clade_indices = Dict(clade => i for (i, clade) in enumerate(all_clades))

    observations = Matrix{Float64}(undef, length(trees), length(all_clades))
    for (i, tree) in enumerate(trees)
        for (clade, parameter) in tree.parameters
            observations[i, clade_indices[clade]] = parameter
        end
    end

    @info "Fitting distribution"

    complete_distribution = fit_mle(D, observations)

    @info "Fitted distribution"

    return MultivariateDist(clade_indices, complete_distribution.μ, complete_distribution.Σ)
end

function sample_tree(distribution::MultivariateDist{D}, tree::Union{CladifiedTree,CladifiedTree})::Union{CladifiedTree,CladifiedTree} where {D}
    clade_indices = [distribution.clade_indices[clade] for clade in keys(tree.parameters)]
    parameters = rand(D(distribution.μ[clade_indices], distribution.Σ[clade_indices, clade_indices]))
    return CladifiedTree(Dict(clade => param for (clade, param) in zip(keys(tree.parameters), parameters)), tree)
end

function log_density(distribution::MultivariateDist{D}, tree::CladifiedTree) where {D}
    clade_indices = [distribution.clade_indices[clade] for clade in keys(tree.parameters)]
    log_density = logpdf(D(distribution.μ[clade_indices], distribution.Σ[clade_indices, clade_indices]), [param for (clade, param) in tree.parameters])
    return log_density
end