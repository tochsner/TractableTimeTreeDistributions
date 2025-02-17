abstract type CCD end

struct CCD1 <: CCD
    cladifies_trees::Vector{CladifiedTree}
    num_taxa::Int
    num_trees::Int

    root_clade::Clade
    clades::Set{Clade}
    splits::Set{CladeSplit}
    splits_per_clade::Dict{Clade,Set{CladeSplit}}

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
    splits_per_clade = Dict()

    for tree in cladified_trees
        for clade in tree.clades
            push!(clades, clade)
            num_occurrences[clade] = get(num_occurrences, clade, 0) + 1
        end

        for split in tree.splits
            push!(splits, split)
            splits_per_clade[split.parent] = push!(get(splits_per_clade, split.parent, Set()), split)
            num_occurrences[split] = get(num_occurrences, split, 0) + 1
        end
    end

    return CCD1(cladified_trees, num_taxa, num_trees, root_clade, clades, splits, splits_per_clade, num_occurrences)
end

# get log probability

function get_log_probability(ccd::CCD1, tree::Tree)
    cladified_tree = cladify_tree(tree)
    return sum(get_log_probability(ccd, split) for split in cladified_tree.splits)
end

function get_log_probability(ccd::CCD1, split::CladeSplit)
    return log(ccd.num_occurrences[split] / ccd.num_occurrences[split.parent])
end

# get most likely tree

function get_most_likely_tree(ccd::CCD1)
    most_likely_splits::Vector{CladeSplit} = []
    collect_most_likely_splits!(ccd, ccd.root_clade, most_likely_splits)
    return construct_tree_from_splits(most_likely_splits)
end

function collect_most_likely_splits!(ccd::CCD1, current_clade::Clade, most_likely_splits::Vector{CladeSplit})
    if isLeaf(current_clade)
        return
    end
    
    current_splits = ccd.splits_per_clade[current_clade]
    most_likely_split = argmax(split -> get_max_log_ccp(ccd, split), current_splits)
    push!(most_likely_splits, most_likely_split)

    collect_most_likely_splits!(ccd, most_likely_split.clade1, most_likely_splits)
    collect_most_likely_splits!(ccd, most_likely_split.clade2, most_likely_splits)
end

function get_max_log_ccp(ccd::CCD1, split::CladeSplit)
    return get_log_probability(ccd, split) + get_max_log_ccp(ccd, split.clade1) + get_max_log_ccp(ccd, split.clade2)
end

function get_max_log_ccp(ccd::CCD1, clade::Clade)
    if isLeaf(clade)
        return 0.0
    else
        return maximum(get_max_log_ccp(ccd, split) for split in ccd.splits_per_clade[clade])
    end
end

# get most likely tree

function sample(ccd::CCD1)
end
