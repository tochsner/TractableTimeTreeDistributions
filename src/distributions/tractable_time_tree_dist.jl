struct TractableTimeTreeDist{Topopolgy,Times} <: AbstractDistribution
    topology::Topopolgy
    times::Times
end

function TractableTimeTreeDist{Topopolgy,Times}(trees::Vector{Tree}) where {Topopolgy,Times}
    TractableTimeTreeDist{Topopolgy,Times}(
        Topopolgy(trees),
        Times(cladify_tree.(trees))
    )
end

function log_density(distribution::TractableTimeTreeDist{Topopolgy,Times}, tree::Tree) where {Topopolgy,Times}
    cladified_tree = cladify_tree(tree)
    return log_density(distribution.topology, cladified_tree) + log_density(distribution.times, cladified_tree)
end

function sample_tree(distribution::TractableTimeTreeDist{Topopolgy,Times})::CladifiedTree where {Topopolgy,Times}
    sample_tree(
        distribution.times,
        sample_tree(distribution.topology)
    )
end