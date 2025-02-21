struct CladifiedTree
    tip_names::Vector{String}
    root::AbstractClade
    splits::Dict{Clade,Split}
end

function cladify_tree(tree::Tree)::CladifiedTree
    tip_names = get_tip_names(tree)
    num_tips = length(tip_names)

    root_clade = Clade(1:num_tips, num_tips)
    splits = Dict()

    function clade_visitor(clade::AbstractClade)
        if is_root(clade)
            root_clade = clade
        end
    end

    function split_visitor(split::Split)
        splits[split.parent] = split
    end

    cladify_tree(tree, clade_visitor, split_visitor)

    return CladifiedTree(tip_names, root_clade, splits)
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

    combined_clade_height = clade1.height + getparentedge(children[1]).length
    combined_clade = union(clade1, clade2, combined_clade_height)
    clade_visitor(combined_clade)

    split = Split(clade1, clade2, combined_clade)
    split_visitor(split)

    return combined_clade
end