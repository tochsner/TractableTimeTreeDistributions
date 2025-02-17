abstract type CCD end

struct CCD1 <: CCD
    cladifies_trees::Vector{CladifiedTree}
    num_taxa::Int
    num_trees::Int

    root_clade::AbstractClade
    clades::Set{AbstractClade}
    splits::Set{CladeSplit}
    splits_per_clade::Dict{AbstractClade,Set{CladeSplit}}

    num_occurrences::Dict{Union{AbstractClade,CladeSplit},Int}
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
    log(ccd.num_occurrences[split] / ccd.num_occurrences[split.parent])
end

# get most likely tree

function get_most_likely_tree(ccd::CCD1)
    most_likely_clades::Vector{Clade} = []
    collect_most_likely_clades!(ccd, ccd.root_clade, most_likely_clades)
    return most_likely_clades
end

function collect_most_likely_clades!(ccd::CCD1, current_clade::Clade, most_likely_clades::Vector{Clade})
    current_splits = ccd.splits_per_clade[current_clade]

    most_likely_split = argmax(split -> get_max_log_ccp(ccd, split), current_splits)
    push!(most_likely_clades, most_likely_split.parent)

    collect_most_likely_clades!(ccd, most_likely_split.clade1, most_likely_clades)
    collect_most_likely_clades!(ccd, most_likely_split.clade2, most_likely_clades)
end

function collect_most_likely_clades!(ccd::CCD1, current_clade::Leaf, most_likely_clades::Vector{Clade}) end

function get_max_log_ccp(ccd::CCD1, split::CladeSplit)
    get_log_probability(ccd, split) + get_max_log_ccp(ccd, split.clade1) + get_max_log_ccp(ccd, split.clade2)
end

function get_max_log_ccp(ccd::CCD1, clade::Leaf)
    0.0
end

function get_max_log_ccp(ccd::CCD1, clade::Clade)
    maximum(get_max_log_ccp(ccd, split) for split in ccd.splits_per_clade[clade])
end

# sample tree

function sample(ccd::CCD1)
    sampled_clades::Vector{Clade} = []
    collect_sampled_clades!(ccd, ccd.root_clade, sampled_clades)
    return sampled_clades
end

function collect_sampled_clades!(ccd::CCD1, current_clade::Clade, sampled_clades::Vector{Clade})
    current_splits = collect(ccd.splits_per_clade[current_clade])
    weights = AnalyticWeights([get_max_log_ccp(ccd, split) for split in current_splits])

    sampled_split = StatsBase.sample(current_splits, weights)

    push!(sampled_clades, sampled_split.parent)
    collect_sampled_clades!(ccd, sampled_split.clade1, sampled_clades)
    collect_sampled_clades!(ccd, sampled_split.clade2, sampled_clades)
end

function collect_sampled_clades!(ccd::CCD1, current_clade::Leaf, sampled_clades::Vector{Clade}) end
