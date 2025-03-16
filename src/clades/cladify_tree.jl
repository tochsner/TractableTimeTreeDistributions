struct CladifiedTree
    tip_names::Vector{String}
    root::AbstractClade
    splits::Dict{Clade,Split}
    parameters::Dict{Clade,Float64}
end

function CladifiedTree(new_parameters::Dict{Clade,Float64}, cladified_tree::CladifiedTree)
    CladifiedTree(cladified_tree.tip_names, cladified_tree.root, cladified_tree.splits, new_parameters)
end

function cladify_tree(tree::Tree)::CladifiedTree
    tip_names = get_tip_names(tree)
    num_tips = length(tip_names)

    root_clade = Clade(1:num_tips, num_tips)
    splits = Dict{Clade,Split}()
    heights = Dict{Clade,Float64}()

    function clade_visitor(clade::AbstractClade, height::Float64)
        if is_root(clade)
            root_clade = clade
        end
        if !is_leaf(clade)
            heights[clade] = height
        end
    end

    function split_visitor(split::Split)
        splits[split.parent] = split
    end

    cladify_tree(tree, clade_visitor, split_visitor)

    return CladifiedTree(tip_names, root_clade, splits, heights)
end

function cladify_tree(tree::Tree, clade_visitor, split_visitor)
    tip_indices = get_leaf_index_mapping(tree)
    num_taxa = tree.numtaxa
    cladify_node(getroot(tree), clade_visitor, split_visitor, tip_indices, num_taxa, Dict{AbstractClade,Float64}())
end

function cladify_node(
    node,
    clade_visitor,
    split_visitor,
    tip_indices::Dict{String,Int},
    num_taxa::Int,
    heights::Dict{AbstractClade,Float64}
)
    if isleaf(node)
        taxa_index = tip_indices[node.name]
        leaf = Leaf(taxa_index, num_taxa, node.name)
        heights[leaf] = 0.0
        clade_visitor(leaf, 0.0)
        return leaf
    end

    children = getchildren(node)
    clade1 = cladify_node(children[1], clade_visitor, split_visitor, tip_indices, num_taxa, heights)
    clade2 = cladify_node(children[2], clade_visitor, split_visitor, tip_indices, num_taxa, heights)

    combined_clade = union(clade1, clade2)
    
    combined_clade_height = heights[clade1] + getparentedge(children[1]).length
    heights[combined_clade] = combined_clade_height
    
    clade_visitor(combined_clade, combined_clade_height)

    split = Split(clade1, clade2, combined_clade)
    split_visitor(split)

    return combined_clade
end