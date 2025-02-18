struct CladifiedTree
    tree::Tree
    tip_indices::Dict{String,Int}
    clades::Set{AbstractClade}
    splits::Set{CladeSplit}
end

function cladify_tree(tree::Tree)::CladifiedTree
    tip_indices = get_leaf_index_mapping(tree)
    tree_with_clades = CladifiedTree(tree, tip_indices, Set([]), Set([]))
    cladify_node!(getroot(tree), tree_with_clades)
    return tree_with_clades
end

function cladify_node!(node, tree_with_clades::CladifiedTree)
    if isleaf(node)
        taxa_index = tree_with_clades.tip_indices[node.name]
        clade = Leaf(taxa_index, tree_with_clades.tree.numtaxa)
        return clade
    end

    children = getchildren(node)
    clade1 = cladify_node!(children[1], tree_with_clades)
    clade2 = cladify_node!(children[2], tree_with_clades)

    combined_clade = union(clade1, clade2)
    push!(tree_with_clades.clades, combined_clade)

    split = CladeSplit(clade1, clade2, combined_clade)
    push!(tree_with_clades.splits, split)

    return combined_clade
end