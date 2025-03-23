function robinson_foulds_distance(tree1::CladifiedTree, tree2::CladifiedTree)
    clades1 = values(tree1.splits)
    clades2 = values(tree2.splits)
    common_clades = intersect(clades1, clades2)
    return length(clades1) - length(common_clades)
end
