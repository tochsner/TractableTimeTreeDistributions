import Phylo

function construct_tree(clades::Set{Clade}, leaves::Set{Leaf})::Phylo.RootedTree 
    tree = Phylo.RootedTree()

    clade_to_node = Dict()

    for leaf in leaves
        clade_to_node[leaf] = Phylo.createnode!(tree, leaf.name)
    end

    for _ in clades
        for clade in clades
            if is_leaf(clade)
                clade_to_node[clade] = Phylo.createnode!(tree, clade.name)
            end
        end
    end

    return tree
end