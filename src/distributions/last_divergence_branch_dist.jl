struct LastDivergenceBranchDist{D,B} <: AbstractDistribution
    last_div::AbstractDistribution
    branches::AbstractDistribution
end

function LastDivergenceBranchDist{D,B}(trees::Vector{CladifiedTree}) where {D,B}
    LastDivergenceBranchDist{D,B}(
        H(transform_last_div.(trees)),
        R(transform_branches.(trees))
    )
end

function sample_tree(distribution::LastDivergenceBranchDist, tree::CladifiedTree)::CladifiedTree
    sampled_tree_with_branches = sample_tree(distribution.branches, tree)
    sampled_tree_with_last_div = sample_tree(distribution.last_div, tree)
    sampled_tree = CladifiedTree(
        merge(sampled_tree_with_last_div.parameters, sampled_tree_with_branches.parameters),
        tree
    )
    return (invert_branches ∘ invert_last_div)(sampled_tree)
end

function log_density(distribution::LastDivergenceBranchDist, tree::CladifiedTree)
    (
        log_density(distribution.last_div, transform_last_div(tree)) +
        log_density(distribution.branches, transform_branches(tree))
    )
end

function transform_last_div(tree::CladifiedTree)::CladifiedTree
    CladifiedTree(Dict(tree.root => tree.parameters[tree.root]), tree)
end
function transform_branches(tree::CladifiedTree)::CladifiedTree
    parameters = Dict{Clade,Float64}()

    for split in values(tree.splits)
        if !is_leaf(split.clade1)
            parameters[split.clade1] = 1.0 - tree.parameters[split.clade1] / tree.parameters[split.parent]
        end

        if !is_leaf(split.clade2)
            parameters[split.clade2] = 1.0 - tree.parameters[split.clade2] / tree.parameters[split.parent]
        end
    end

    return CladifiedTree(parameters, tree)
end

function invert_last_div(tree::CladifiedTree)::CladifiedTree
    tree
end

function invert_branches(tree::CladifiedTree)::CladifiedTree
    root_split = tree.splits[tree.root]
    branches!(tree, root_split.clade1, tree.parameters[tree.root])
    branches!(tree, root_split.clade2, tree.parameters[tree.root])
    return tree
end
