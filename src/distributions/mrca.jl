function mrca_tree(tree::CladifiedTree, observed_trees::Vector{CladifiedTree})
    parameters = Dict{Clade,Float64}()

    for clade in keys(tree.splits)
        total_heights = 0.0

        for tree in observed_trees
            if haskey(tree.parameters, clade)
                total_heights += tree.parameters[clade]
            else
                mrca_height = Inf
                for (potential_mrca, height) in tree.parameters
                    if clade in potential_mrca
                        mrca_height = min(mrca_height, height)
                    end
                end
                total_heights += mrca_height
            end
        end
        parameters[clade] = total_heights / length(observed_trees)
    end

    return CladifiedTree(parameters, tree)
end