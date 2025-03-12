using PhyloNetworks

const Tree = HybridNetwork

function load_trees(mcmc_path::String)::Vector{Tree}
    trees = readnexus_treeblock(mcmc_path)

    if length(trees) == 0
        open(mcmc_path) do io
            trees = map(readnewick, io |> eachline)
        end
    end

    return trees
end

function load_trees_from_newick(newick::String)::Vector{Tree}
    [readnewick(newick)]
end

function load_trees_from_newick(newicks::Vector{String})::Vector{Tree}
    [readnewick(newick) for newick in newicks]
end

function get_tip_names(tree::Tree)::Vector{String}
    tiplabels(tree) |> sort
end

function get_leaf_index_mapping(tree::Tree)
    Dict(leaf => i for (i, leaf) in get_tip_names(tree) |> enumerate)
end

function get_leaf_index(tree::Tree, leaf::String)
    findfirst(x -> x == leaf, sort(tiplabels(tree)))
end