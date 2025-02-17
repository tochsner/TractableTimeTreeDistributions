using PhyloNetworks

const Tree = HybridNetwork

function load_trees(mcmc_path::String)::Vector{Tree}
    readnexus_treeblock(mcmc_path)
end

function load_trees_from_newick(mcmc_path::String)::Vector{Tree}
    [readnewick(mcmc_path)]
end
