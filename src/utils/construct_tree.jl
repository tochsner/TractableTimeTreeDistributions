import Phylo

function construct_tree(cladified_tree::CladifiedTree)::Phylo.RootedTree
    tree = Phylo.RootedTree()
    root_clade = cladified_tree.root
    construct_tree!(root_clade, cladified_tree, tree)
    return tree
end
function construct_tree!(clade::Clade, cladified_tree::CladifiedTree, tree::Phylo.RootedTree)
    split = cladified_tree.splits[clade]

    child_node1 = construct_tree!(split.clade1, cladified_tree, tree)
    child_node2 = construct_tree!(split.clade2, cladified_tree, tree)

    clade_node = Phylo.createnode!(tree, missing)

    branch1 = cladified_tree.parameters[clade] - (is_leaf(split.clade1) ? 0.0 : cladified_tree.parameters[split.clade1])
    branch2 = cladified_tree.parameters[clade] - (is_leaf(split.clade2) ? 0.0 : cladified_tree.parameters[split.clade2])

    Phylo.createbranch!(tree, clade_node, child_node1, branch1)
    Phylo.createbranch!(tree, clade_node, child_node2, branch2)

    return clade_node
end
construct_tree!(clade::Leaf, cladified_tree, tree) = Phylo.createnode!(tree, clade.name)

write_tree(path::String, cladified_tree::CladifiedTree) =
    Phylo.write(path, construct_tree(cladified_tree), format = Phylo.Nexus())
write_tree(io::IOStream, cladified_tree::CladifiedTree) =
    Phylo.write(io, construct_tree(cladified_tree), format = Phylo.Nexus())
write_tree(io::IOStream, labels::Vector{String}, cladified_trees::Vector{CladifiedTree}) = Phylo.write(
    io,
    Phylo.TreeSet(Dict(l => construct_tree(cladified_tree) for (l, cladified_tree) in zip(labels, cladified_trees))),
)
