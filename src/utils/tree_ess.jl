using MCMCDiagnosticTools

function get_ess(distances_to_map_tree::Vector{Int64})
    ess(distances_to_map_tree, split_chains=1)
end

function get_ess(cladified_trees::Vector{CladifiedTree})
    ccd = CCD1(cladified_trees)
    map_tree = point_estimate(ccd)
    distances_to_map_tree = robinson_foulds_distance.(cladified_trees, Ref(map_tree))
    return min(
        get_ess(distances_to_map_tree), 
        length(cladified_trees)
    )
end