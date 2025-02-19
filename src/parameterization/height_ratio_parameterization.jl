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

function set_heights(::HeightRatioParameterization, parameterized_trees::Vector{ParameterizedTree})::Vector{CladifiedTree}

end