import Phylo

function construct_tree(clades::Set{AbstractClade})::Phylo.RootedTree 
    tree = Phylo.RootedTree()

    clade_to_node = Dict()

    sorted_clades = sort(collect(clades), by = size, rev = true)

    for clade in sorted_clades
        if is_leaf(clade)
            clade_to_node[clade] = Phylo.createnode!(tree, clade.name)
        else
            clade_to_node[clade] = Phylo.createnode!(tree, missing)
        end

        if !is_root(clade)
            potential_parents = clade_to_node |> keys |> filter(x -> clade in x && x != clade)
            parent = argmin(size, potential_parents)

            Phylo.createbranch!(tree, clade_to_node[parent], clade_to_node[clade], parent.height - clade.height)
        end
    end

    return tree
end