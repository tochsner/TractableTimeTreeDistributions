struct CladifiedTree
    tree::Tree
    tip_indices::Dict{String,Int}
    clades::Set{AbstractClade}
    splits::Set{CladeSplit}
end

function cladify_tree(tree::Tree)::CladifiedTree
    tip_indices = get_leaf_index_mapping(tree)

    tree_with_clades = CladifiedTree(tree, tip_indices, Set([]), Set([]))
    
    clade_visitor = x -> push!(tree_with_clades.clades, x)
    split_visitor = x -> push!(tree_with_clades.splits, x)
    
    cladify_tree(tree, clade_visitor, split_visitor)
    
    return tree_with_clades
end

function cladify_tree(tree::Tree, clade_visitor, split_visitor)
    tip_indices = get_leaf_index_mapping(tree)
    num_taxa = tree.numtaxa
    cladify_node(getroot(tree), clade_visitor, split_visitor, tip_indices, num_taxa)
end

function cladify_node(
    node, 
    clade_visitor, 
    split_visitor, 
    tip_indices::Dict{String,Int}, 
    num_taxa::Int
)
    if isleaf(node)
        taxa_index = tip_indices[node.name]
        leaf = Leaf(taxa_index, num_taxa, node.name)
        clade_visitor(leaf)
        return leaf
    end

    children = getchildren(node)
    clade1 = cladify_node(children[1], clade_visitor, split_visitor, tip_indices, num_taxa)
    clade2 = cladify_node(children[2], clade_visitor, split_visitor, tip_indices, num_taxa)

    combined_clade = union(clade1, clade2)
    clade_visitor(combined_clade)

    split = CladeSplit(clade1, clade2, combined_clade)
    split_visitor(split)

    return combined_clade
end