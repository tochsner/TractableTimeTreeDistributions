abstract type CCD end

struct CCD1 <: CCD
    cladifies_trees::Vector{CladifiedTree}
    num_taxa::Int
    num_trees::Int

    root_clade::Clade
    clades::Set{Clade}
    splits::Set{CladeSplit}

    num_occurrences::Dict{Union{Clade,CladeSplit},Int}
end

function CCD1(trees::Vector{Tree})
    num_taxa = length(tiplabels(trees[1]))
    num_trees = length(trees)
    cladified_trees = map(cladify_tree, trees)

    root_clade = Clade(1:num_taxa, num_taxa)
    clades = Set()
    splits = Set()
    num_occurrences = Dict()

    for tree in cladified_trees
        for clade in tree.clades
            push!(clades, clade)
            num_occurrences[clade] = get(num_occurrences, clade, 0) + 1
        end

        for split in tree.splits
            push!(splits, split)
            num_occurrences[split] = get(num_occurrences, split, 0) + 1
        end
    end

    return CCD1(cladified_trees, num_taxa, num_trees, root_clade, clades, splits, num_occurrences)
end

function get_log_probability(ccd::CCD1, tree::Tree)
    cladified_tree = cladify_tree(tree)
    return sum(get_log_probability(ccd, split) for split in cladified_tree.splits)
end

function get_log_probability(ccd::CCD1, split::CladeSplit)
    return log(ccd.num_occurrences[split] / ccd.num_occurrences[split.parent])
end

function get_most_likely_tree(ccd::CCD1)
end

function sample(ccd::CCD1)
end
