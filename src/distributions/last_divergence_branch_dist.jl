struct LastDivergenceBranchDist{D,B} <: AbstractDistribution
    last_div::AbstractDistribution
    branches::AbstractDistribution
end

LastDivergenceBranchDist{D,B}(trees::Vector{CladifiedTree}) where {D,B} = LastDivergenceBranchDist{D,B}(
    D(transform_last_divergence.(trees)),
    B(transform_branches.(trees))
)

readable_name(::Type{LastDivergenceBranchDist{D,B}}) where {D,B} = "Last Divergence ($(readable_name(D))), Branches ($(readable_name(B)))"

# high-level functions

function sample_tree(distribution::LastDivergenceBranchDist, tree::CladifiedTree)::CladifiedTree
    sampled_tree_with_branches = sample_tree(distribution.branches, tree)
    sampled_tree_with_last_div = sample_tree(distribution.last_div, tree)
    sampled_tree = CladifiedTree(
        merge(sampled_tree_with_last_div.parameters, sampled_tree_with_branches.parameters),
        tree
    )
    return invert_last_divergence_branches(sampled_tree)
end

function point_estimate(distribution::LastDivergenceBranchDist, tree::CladifiedTree)
    map_tree_with_branches = point_estimate(distribution.branches, tree)
    map_tree_with_last_div = point_estimate(distribution.last_div, tree)
    map_tree = CladifiedTree(
        merge(map_tree_with_last_div.parameters, map_tree_with_branches.parameters),
        tree
    )
    return invert_last_divergence_branches(map_tree)
end

log_density(distribution::LastDivergenceBranchDist, tree::CladifiedTree) = (
    log_density(distribution.last_div, transform_last_divergence(tree)) +
    log_density(distribution.branches, transform_branches(tree))
)

# transformations

function transform_last_divergence(tree::CladifiedTree)::CladifiedTree
    last_divergence = minimum(param for (clade, param) in tree.parameters if size(clade) == 2)
    return CladifiedTree(Dict(tree.root => last_divergence), tree)
end
function transform_branches(tree::CladifiedTree)::CladifiedTree
    parameters = Dict{Clade,Float64}()

    for split in values(tree.splits)
        if !is_leaf(split.clade1)
            parameters[split.clade1] = tree.parameters[split.parent] - tree.parameters[split.clade1]
        end

        if !is_leaf(split.clade2)
            parameters[split.clade2] = tree.parameters[split.parent] - tree.parameters[split.clade2]
        end
    end

    return CladifiedTree(parameters, tree)
end

function invert_last_divergence_branches(tree::CladifiedTree)::CladifiedTree
    root_split = tree.splits[tree.root]
    set_height_from_branches!(tree, root_split.clade1, tree.parameters[tree.root])
    set_height_from_branches!(tree, root_split.clade2, tree.parameters[tree.root])

    current_last_divergence = minimum(param for (clade, param) in tree.parameters if size(clade) == 2)
    actual_last_divergence = tree.parameters[tree.root]
    new_parameters = Dict(
        clade => old_param - current_last_divergence + actual_last_divergence
        for (clade, old_param) in tree.parameters if !is_leaf(clade)
    )

    return CladifiedTree(new_parameters, tree)
end

function set_height_from_branches!(tree::CladifiedTree, clade::Clade, parent_height::Float64,)
    tree.parameters[clade] = parent_height - tree.parameters[clade]

    split = tree.splits[clade]
    set_height_from_branches!(tree, split.clade1, tree.parameters[clade])
    set_height_from_branches!(tree, split.clade2, tree.parameters[clade])
end
function set_height_from_branches!(parameterized_tree::CladifiedTree, clade::Leaf, parent_height::Float64,) end