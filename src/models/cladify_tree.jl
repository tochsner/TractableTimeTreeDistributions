struct TreeWithClades
    tree::Tree
    clades::Vector{Clade}
    splits::Vector{CladeSplit}
end

function cladify_tree(tree::Tree)
    tree_with_clades = TreeWithClades(tree, [], [])
    cladify_node(getroot(tree), tree_with_clades)
    return tree_with_clades
end

function cladify_node(node, tree_with_clades::TreeWithClades)
    if isleaf(node)
        taxa_index = findfirst(x -> x == node.name, tiplabels(tree_with_clades.tree))
        clade = Clade(taxa_index, tree_with_clades.tree.numtaxa)
        push!(tree_with_clades.clades, clade)
        return clade
    end

    children = getchildren(node)
    clade1 = cladify_node(children[1], tree_with_clades)
    clade2 = cladify_node(children[2], tree_with_clades) 

    split = CladeSplit(clade1, clade2)
    push!(tree_with_clades.splits, split)

    combined_clade = union(clade1, clade2)
    push!(tree_with_clades.clades, combined_clade)

    return combined_clade
end