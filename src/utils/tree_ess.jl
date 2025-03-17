using MCMCDiagnosticTools

function get_ess(distances_to_map_tree::Vector)
    ess(distances_to_map_tree, split_chains=1)
end

function get_ess(cladified_trees::Vector{CladifiedTree})
    
    ccd = CCD1(cladified_trees)
    map_tree = point_estimate(ccd)
    distances_to_map_tree = robinson_foulds_distance.(cladified_trees, Ref(map_tree))
    distances_to_map_tree_ess = get_ess(distances_to_map_tree)

    total_branch_lengths = [sum(values(tree.parameters)) for tree in cladified_trees]
    total_branch_lengths_ess = get_ess(total_branch_lengths)

    if distances_to_map_tree_ess == 0
        # all topologies are identical
        return total_branch_lengths_ess
    end
    if total_branch_lengths_ess == 0
        # all topologies have the same total branch lengths
        return distances_to_map_tree_ess
    end

    return min(distances_to_map_tree_ess, total_branch_lengths_ess)
end