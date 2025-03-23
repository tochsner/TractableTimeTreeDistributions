struct ShorterBranchDist{B} <: AbstractDistribution
    branches::AbstractDistribution
end

ShorterBranchDist{B}(trees::Vector{CladifiedTree}) where {B} = ShorterBranchDist{B}(B(transform_short_branches.(trees)))

readable_name(::Type{ShorterBranchDist{B}}) where {B} = "Shortest Branch ($(readable_name(B)))"

# high-level functions

sample_tree(distribution::ShorterBranchDist, tree::CladifiedTree)::CladifiedTree =
    sample_tree(distribution.branches, tree) |> invert_short_branches
point_estimate(distribution::ShorterBranchDist, tree::CladifiedTree) =
    point_estimate(distribution.branches, tree) |> invert_short_branches
log_density(distribution::ShorterBranchDist, tree::CladifiedTree) =
    log_density(distribution.branches, transform_short_branches(tree))

# transformations

transform_short_branches(tree::CladifiedTree)::CladifiedTree = CladifiedTree(
    Dict(
        s.parent =>
            tree.parameters[s.parent] - max(get(tree.parameters, s.clade1, 0.0), get(tree.parameters, s.clade2, 0.0)) for
        s in values(tree.splits)
    ),
    tree,
)


function invert_short_branches(tree::CladifiedTree)::CladifiedTree
    root_split = tree.splits[tree.root]
    set_height_from_short_branches!(tree, root_split)
    return tree
end

function set_height_from_short_branches!(tree::CladifiedTree, split::Split)
    if !is_leaf(split.clade1)
        set_height_from_short_branches!(tree, tree.splits[split.clade1])
    end

    if !is_leaf(split.clade2)
        set_height_from_short_branches!(tree, tree.splits[split.clade2])
    end

    short_branch_length = tree.parameters[split.parent]
    clade1_height = get(tree.parameters, split.clade1, 0.0)
    clade2_height = get(tree.parameters, split.clade2, 0.0)

    tree.parameters[split.parent] = short_branch_length + max(clade1_height, clade2_height)
end
