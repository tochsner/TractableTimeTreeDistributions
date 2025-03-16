struct HeightRatioDist{H,R} <: AbstractDistribution
    height::AbstractDistribution
    ratios::AbstractDistribution
end

function HeightRatioDist{H,R}(trees::Vector{CladifiedTree}) where {H,R}
    HeightRatioDist{H,R}(
        H(transform_height.(trees)),
        R(transform_ratios.(trees))
    )
end

function readable_name(distribution::Type{HeightRatioDist{H, R}}) where {H, R}
    "Height ($(readable_name(H))), Ratios ($(readable_name(R)))"
end

function most_likely_tree(distribution::HeightRatioDist, tree::CladifiedTree)::CladifiedTree
    map_tree_with_height = most_likely_tree(distribution.height, tree)
    map_tree_with_ratios = most_likely_tree(distribution.ratios, tree)
    map_tree_with_height_and_ratios = CladifiedTree(
        merge(map_tree_with_height.parameters, map_tree_with_ratios.parameters),
        map_tree_with_height
    )
    return (invert_ratios ∘ invert_height)(map_tree_with_height_and_ratios)
end

function sample_tree(distribution::HeightRatioDist, tree::CladifiedTree)::CladifiedTree
    sampled_tree_with_height = sample_tree(distribution.height, tree)
    sampled_tree_with_ratios = sample_tree(distribution.ratios, tree)
    sampled_tree_with_height_and_ratios = CladifiedTree(
        merge(sampled_tree_with_height.parameters, sampled_tree_with_ratios.parameters),
        sampled_tree_with_height
    )
    return (invert_ratios ∘ invert_height)(sampled_tree_with_height_and_ratios)
end

function log_density(distribution::HeightRatioDist, tree::CladifiedTree)
    (
        log_density(distribution.height, transform_height(tree)) +
        log_density(distribution.ratios, transform_ratios(tree)) +
        log_abs_det_jacobian(distribution, tree)
    )
end

function transform_height(tree::CladifiedTree)::CladifiedTree
    CladifiedTree(Dict(tree.root => tree.parameters[tree.root]), tree)
end
function transform_ratios(tree::CladifiedTree)::CladifiedTree
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
function log_abs_det_jacobian(distribution::HeightRatioDist, tree::CladifiedTree)
    log_abs_det_jacobian = 0.0

    for split in values(tree.splits)
        if !is_leaf(split.clade1)
            log_abs_det_jacobian -= log(tree.parameters[split.parent])
        end

        if !is_leaf(split.clade2)
            log_abs_det_jacobian -= log(tree.parameters[split.parent])
        end
    end

    return log_abs_det_jacobian
end

function invert_height(tree::CladifiedTree)::CladifiedTree
    tree
end

function invert_ratios(tree::CladifiedTree)::CladifiedTree
    root_split = tree.splits[tree.root]
    set_height_from_ratios!(tree, root_split.clade1, tree.parameters[tree.root])
    set_height_from_ratios!(tree, root_split.clade2, tree.parameters[tree.root])
    return tree
end

function set_height_from_ratios!(tree::CladifiedTree, clade::Clade, parent_height::Float64)
    tree.parameters[clade] = parent_height * (1.0 - tree.parameters[clade])

    split = tree.splits[clade]
    set_height_from_ratios!(tree, split.clade1, tree.parameters[clade])
    set_height_from_ratios!(tree, split.clade2, tree.parameters[clade])
end
function set_height_from_ratios!(parameterized_tree::CladifiedTree, clade::Leaf, parent_height::Float64) end