struct HeightRatioParameterization <: Parameterization end

function parameterize(::HeightRatioParameterization, tree::CladifiedTree)::ParameterizedTree
    parameters = Dict{Clade, Float64}()

    for split in values(tree.splits)
        if !is_leaf(split.clade1)
            parameters[split.clade1] = 1.0 - split.clade1.height / split.parent.height
        end

        if !is_leaf(split.clade2)
            parameters[split.clade2] = 1.0 - split.clade2.height / split.parent.height
        end
    end

    parameters[tree.root] = tree.root.height

    return ParameterizedTree(parameters, tree)
end

function set_heights!(::HeightRatioParameterization, parameterized_tree::ParameterizedTree)::CladifiedTree
    root = parameterized_tree.root

    root_split = parameterized_tree.splits[root]
    root.height = parameterized_tree.parameters[root]
    set_height!(HeightRatioParameterization(), parameterized_tree , root_split.clade1, root.height)
    set_height!(HeightRatioParameterization(), parameterized_tree , root_split.clade2, root.height)

    return CladifiedTree(parameterized_tree.tip_names, root, parameterized_tree.splits)
end

function set_height!(::HeightRatioParameterization, parameterized_tree::ParameterizedTree, clade::Clade, parent_height::Float64)
    clade.height = parent_height * (1.0 - parameterized_tree.parameters[clade])

    split = parameterized_tree.splits[clade]
    set_height!(HeightRatioParameterization(), parameterized_tree, split.clade1, clade.height)
    set_height!(HeightRatioParameterization(), parameterized_tree, split.clade2, clade.height)
end

function set_height!(::HeightRatioParameterization, parameterized_tree::ParameterizedTree, clade::Leaf, parent_height::Float64) end