using PhyloNetworks

const Tree = HybridNetwork

function load_trees(mcmc_path::String)::Vector{Tree}
    readnexus_treeblock(mcmc_path)
end

function load_trees_from_newick(newick::String)::Vector{Tree}
    [readnewick(newick)]
end

function load_trees_from_newick(newicks::Vector{String})::Vector{Tree}
    [readnewick(newick) for newick in newicks]
end

function get_leaf_index(tree::Tree, leaf::String)
    return findfirst(x -> x == leaf, sort(tiplabels(tree)))
end