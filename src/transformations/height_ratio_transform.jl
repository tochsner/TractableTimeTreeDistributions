struct HeightRatioTransform <: Transform
    height::Transform
    ratios::Transform
end

function fit(transformation::HeightRatioTransform, trees::Vector{CladifiedTree})
    fit(transformation.height, transform_height.(trees))
    fit(transformation.ratios, transform_ratios.(trees))
end

function sample(transformation::HeightRatioTransform, tree::CladifiedTree)::ParameterizedTree
    sampled_tree = sample(transformation.height, tree)
    sampled_tree = sample(transformation.ratios, sampled_tree)
    return (invert_height ∘ invert_ratios)(sampled_tree)
end

function log_density(transformation::HeightRatioTransform, tree::CladifiedTree)
    (
        log_density(transformation.height, transform_height(tree)) + 
        log_density(transformation.ratios, transform_ratios(tree))
    )
end

function transform(::HeightRatioTransform, tree::CladifiedTree)::ParameterizedTree
    tree |> transform_height |> transform_ratios
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

function log_density(transform::HeightRatioTransform, tree::ParameterizedTree)
    log_density(transform.wrapped, transform(transform, tree)) # TODO multiply with jacobian
end
