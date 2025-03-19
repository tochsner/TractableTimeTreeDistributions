struct TreeDirichletDist <: AbstractDistribution
    K::Float64
    alphas::Dict{Clade,Float64}
end

function TreeDirichletDist(trees::Vector{CladifiedTree})
    K = estimate_k(trees)
    
    observations = DefaultDict{Clade,Vector{Float64}}([])
    observed_parents = DefaultDict{Clade,Vector{Clade}}([])

    for tree in trees
        for (clade, parameter) in tree.parameters
            push!(observations[clade], parameter)
        end

        for split in values(tree.splits)
            if !is_leaf(split.clade1)
                push!(observed_parents[split.clade1], split.parent)
            end
            if !is_leaf(split.clade2)
                push!(observed_parents[split.clade2], split.parent)
            end
        end
    end

    alphas = estimate_clade_alphas(observations, observed_parents, K, trees[1].root)

    return TreeDirichletDist(K, alphas)
end

function estimate_k(trees::Vector{CladifiedTree})
    total_fractions = DefaultDict{Clade, Vector{Float64}}([])

    for tree in trees
        tree_height = tree.heights[tree.root]

        for split in values(tree.splits)
            if !is_leaf(split.clade1)
                push!(total_fractions[split.parent], (tree.heights[split.parent] - tree.heights[split.clade1]) / tree_height)
            end
            if !is_leaf(split.clade2)
                push!(total_fractions[split.parent], (tree.heights[split.parent] - tree.heights[split.clade2]) / tree_height)
            end
        end
    end

    sum_K = 0.0
    normalization = 0.0

    for fractions in values(total_fractions)
        if length(fractions) < 5
            continue
        end
        
        beta = try
            fit_mle(Beta, fractions)
        catch
            fit(Beta, fractions)
        end

        K = beta.α + beta.β
        n = length(fractions)

        sum_K += n*K
        normalization += n
    end

    return sum_K / normalization
end

function estimate_clade_alphas(observations, observed_parents, K, root)
    alphas = Dict{Clade,Float64}(root => K)

    for _ in observations
        for (clade, observations) in observations
            if haskey(alphas, clade)
                continue
            end

            parent_alphas = [
                alphas[parent] for parent in observed_parents[clade] if haskey(alphas, parent)
            ]
            if length(parent_alphas) != length(observed_parents[clade])
                continue
            end

            alpha = estimate_clade_alpha(observations, parent_alphas)
            alphas[clade] = alpha
        end
    end

    return alphas
end

function estimate_clade_alpha(observations, parent_alphas)
    n = length(observations)
    logF = sum(log.(observations))
    logOneMinusF = sum(log.(1 .- observations))

    max_alpha = minimum(parent_alphas)

    if length(observations) < 5
        return 0.5 * max_alpha
    end

    try
        return bisection(1e-5, max_alpha - 1e-5; tol=1e-5, maxiter=500) do alpha
            n * digamma(alpha) - sum(digamma.(parent_alphas .- alpha)) + logF - logOneMinusF
        end
    catch
        @warn "Error estimating dirichlet parameters"
        return 0.5 * max_alpha
    end
end

function readable_name(distribution::Type{TreeDirichletDist})
    "TreeDirichletDist"
end

function collect_beta_distributions(distribution::TreeDirichletDist, tree::CladifiedTree)
    distributions = Dict{Clade,Beta}()
    collect_beta_distributions(distribution, tree, tree.root, distributions)
    return distributions
end
function collect_beta_distributions(distribution::TreeDirichletDist, tree::CladifiedTree, parent::Clade, distributions::Dict{Clade,Beta})
    child1 = tree.splits[parent].clade1
    if !is_leaf(child1) && haskey(distribution.alphas, child1)
        distributions[child1] = Beta(
            max(0.05 * distribution.alphas[parent], distribution.alphas[parent] - distribution.alphas[child1]),
            distribution.alphas[child1]
        )
        collect_beta_distributions(distribution, tree, child1, distributions)
    end

    child2 = tree.splits[parent].clade2
    if !is_leaf(child2) && haskey(distribution.alphas, child2)
        distributions[child2] = Beta(
            max(0.05 * distribution.alphas[parent], distribution.alphas[parent] - distribution.alphas[child2]),
            distribution.alphas[child2]
        )
        collect_beta_distributions(distribution, tree, child2, distributions)
    end
end
function collect_beta_distributions(distribution::TreeDirichletDist, tree::CladifiedTree, clade::Leaf, distributions::Dict{Clade,Beta})
end

function sample_tree(distribution::TreeDirichletDist, tree::CladifiedTree)::CladifiedTree
    distributions = collect_beta_distributions(distribution, tree)
    parameters = Dict(
        clade => rand(dist) for (clade, dist) in distributions
    )
    return CladifiedTree(parameters, tree)
end

function point_estimate(distribution::TreeDirichletDist, tree::CladifiedTree)::CladifiedTree
    distributions = collect_beta_distributions(distribution, tree)
    parameters = Dict(
        clade => mean(distributions[clade])
        for clade in keys(tree.splits)
        if haskey(distributions, clade)
    )
    return CladifiedTree(parameters, tree)
end

function log_density(distribution::TreeDirichletDist, tree::CladifiedTree)
    distributions = collect_beta_distributions(distribution, tree)
    sum(
        haskey(distributions, clade) ?
        logpdf(distributions[clade], param)
        : 0.0 for (clade, param) in tree.parameters
    )
end