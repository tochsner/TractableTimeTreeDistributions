struct TractableTimeTreeDist{Topology,Times} <: AbstractDistribution
    topology::Topology
    times::Times
end

TractableTimeTreeDist{Topology,Times}(trees::Vector{Tree}) where {Topology,Times} =
    TractableTimeTreeDist{Topology,Times}(Topology(trees), Times(cladify_tree.(trees)))

readable_name(::Type{TractableTimeTreeDist{Topology,Times}}) where {Topology,Times} = "$(readable_name(Times))"

log_density(distribution::TractableTimeTreeDist{Topology,Times}, tree::Tree) where {Topology,Times} =
    log_density(distribution, cladify_tree(tree))
function log_density(
    distribution::TractableTimeTreeDist{Topology,Times},
    cladified_tree::CladifiedTree,
) where {Topology,Times}
    log_density(distribution.topology, cladified_tree) + log_density(distribution.times, cladified_tree)
end

function point_estimate(distribution::TractableTimeTreeDist{Topology,Times}) where {Topology,Times}
    point_estimate(distribution.times, point_estimate(distribution.topology))
end

function sample_tree(distribution::TractableTimeTreeDist{Topology,Times})::CladifiedTree where {Topology,Times}
    sample_tree(distribution.times, sample_tree(distribution.topology))
end
