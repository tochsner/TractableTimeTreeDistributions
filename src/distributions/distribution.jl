abstract type AbstractDistribution end

log_density(distribution::AbstractDistribution, tree::Tree) = log_density(distribution, cladify_tree(tree))
readable_name(distribution) = nameof(distribution)
