using MCMCDiagnosticTools

function get_ess(distances_to_map_tree::Vector)
    ess(distances_to_map_tree, split_chains=1)
end

function get_ess(cladified_trees::Vector{CladifiedTree})
    total_branch_lengths = [sum(values(tree.parameters)) for tree in cladified_trees]
    
    try
        total_branch_lengths_ess = get_ess(total_branch_lengths)
    catch
        total_branch_lengths_ess = Inf
    end
    
    ccd = CCD1(cladified_trees)
    map_tree = point_estimate(ccd)
    distances_to_map_tree = robinson_foulds_distance.(cladified_trees, Ref(map_tree))
    
    try
        distances_to_map_tree_ess = get_ess(distances_to_map_tree)
    catch
        distances_to_map_tree_ess = Inf
    end
    
    return min(distances_to_map_tree_ess, total_branch_lengths_ess)
end