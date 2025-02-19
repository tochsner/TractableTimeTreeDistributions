import Phylo

function construct_tree(cladified_tree::CladifiedTree)::Phylo.RootedTree 
    tree = Phylo.RootedTree()
    root_clade = cladified_tree.root
    construct_tree!(root_clade, cladified_tree, tree)
    return tree
end

function construct_tree!(clade::Clade, cladified_tree, tree)
    split = cladified_tree.splits[clade]

    child_node1 = construct_tree!(split.clade1, cladified_tree, tree)
    child_node2 = construct_tree!(split.clade2, cladified_tree, tree)
    
    clade_node = Phylo.createnode!(tree, missing)

    Phylo.createbranch!(tree, clade_node, child_node1, clade.height - split.clade1.height)
    Phylo.createbranch!(tree, clade_node, child_node2, clade.height - split.clade2.height)

    return clade_node
end

function construct_tree!(clade::Leaf, cladified_tree, tree)
    Phylo.createnode!(tree, clade.name)
end