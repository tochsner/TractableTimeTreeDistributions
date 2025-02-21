struct HeightRatioDist{H, R} <: AbstractDistribution
    height::AbstractDistribution
    ratios::AbstractDistribution
end

function HeightRatioDist{H, R}(trees::Vector{CladifiedTree}) where {H, R}
    HeightRatioDist{H, R}(
        H(transform_height.(trees)),
        R(transform_ratios.(trees))
    )
end

function sample_tree(distribution::HeightRatioDist, tree::CladifiedTree)::CladifiedTree
    sampled_tree = sample_tree(distribution.height, tree)
    sampled_tree = sample_tree(distribution.ratios, sampled_tree)
    return (invert_ratios ∘ invert_height)(sampled_tree) |> to_cladified_tree
end

function log_density(distribution::HeightRatioDist, tree::CladifiedTree)
    (
        log_density(distribution.height, transform_height(tree)) + 
        log_density(distribution.ratios, transform_ratios(tree)) +
        log_abs_det_jacobian(distribution, tree)
    )
end

function transform_height(tree::CladifiedTree)::ParameterizedTree
    ParameterizedTree(Dict(tree.root => tree.root.height), tree)
end
function transform_ratios(tree::CladifiedTree)::ParameterizedTree
    parameters = Dict{Clade, Float64}()

    for split in values(tree.splits)
        if !is_leaf(split.clade1)
            parameters[split.clade1] = 1.0 - split.clade1.height / split.parent.height
        end

        if !is_leaf(split.clade2)
            parameters[split.clade2] = 1.0 - split.clade2.height / split.parent.height
        end
    end

    return ParameterizedTree(parameters, tree)
end
function log_abs_det_jacobian(distribution::HeightRatioDist, tree::CladifiedTree)
    log_abs_det_jacobian = 0.0

    for split in values(tree.splits)
        if !is_leaf(split.clade1)
            log_abs_det_jacobian -= log(split.parent.height)
        end
        
        if !is_leaf(split.clade2)
            log_abs_det_jacobian -= log(split.parent.height)
        end
    end

    return log_abs_det_jacobian
end

function invert_height(tree::ParameterizedTree)::ParameterizedTree
    tree.root.height = tree.parameters[tree.root]
    return tree
end

function invert_ratios(tree::ParameterizedTree)::ParameterizedTree
    root_split = tree.splits[tree.root]
    set_height_from_ratios!(tree , root_split.clade1, tree.root.height)
    set_height_from_ratios!(tree , root_split.clade2, tree.root.height)
    return tree
end

function set_height_from_ratios!(parameterized_tree::ParameterizedTree, clade::Clade, parent_height::Float64)
    clade.height = parent_height * (1.0 - parameterized_tree.parameters[clade])

    split = parameterized_tree.splits[clade]
    set_height_from_ratios!(parameterized_tree, split.clade1, clade.height)
    set_height_from_ratios!(parameterized_tree, split.clade2, clade.height)
end
function set_height_from_ratios!(parameterized_tree::ParameterizedTree, clade::Leaf, parent_height::Float64) end