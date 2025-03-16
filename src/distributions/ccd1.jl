struct CCD1 <: AbstractDistribution
    num_taxa::Int
    num_trees::Int

    tip_names::Vector{String}
    root_clade::AbstractClade
    clades::Set{AbstractClade}
    splits::Set{Split}
    splits_per_clade::Dict{AbstractClade,Set{Split}}

    num_clade_occurrences::Dict{AbstractClade,Int}
    num_split_occurrences::Dict{Split,Int}
end

function CCD1(trees::Vector{Tree})
    num_taxa = length(tiplabels(trees[1]))
    tip_names = get_tip_names(trees[1])
    num_trees = length(trees)
    root_clade = Clade(1:num_taxa, num_taxa)

    clades = Set()
    num_clade_occurrences = DefaultDict{AbstractClade,Int64}(0)

    function clade_visitor(clade::AbstractClade, height::Float64)
        push!(clades, clade)
        num_clade_occurrences[clade] += 1
    end

    splits = Set()
    num_split_occurrences = DefaultDict{Split,Int64}(0)
    splits_per_clade = DefaultDict{Clade,Set{Split}}(Set)

    function split_visitor(split::Split)
        push!(splits, split)
        push!(splits_per_clade[split.parent], split)

        if !(split.clade1 in split.parent) || !(split.clade2 in split.parent)
            throw("Not a child")
        end

        num_split_occurrences[split] += 1
    end

    for tree in trees
        cladify_tree(tree, clade_visitor, split_visitor)
    end

    return CCD1(num_taxa, num_trees, tip_names, root_clade, clades, splits, splits_per_clade, num_clade_occurrences, num_split_occurrences)
end

function CCD1(cladified_trees::Vector{CladifiedTree})
    tip_names = cladified_trees[1].tip_names
    num_taxa = length(tip_names)
    num_trees = length(cladified_trees)
    root_clade = Clade(1:num_taxa, num_taxa)

    clades = Set()
    num_clade_occurrences = DefaultDict{AbstractClade,Int64}(0)

    splits = Set()
    num_split_occurrences = DefaultDict{Split,Int64}(0)
    splits_per_clade = DefaultDict{Clade,Set{Split}}(Set)

    for tree in cladified_trees
        push!(clades, tree.root)
        num_clade_occurrences[tree.root] += 1

        for split in values(tree.splits)
            push!(splits, split)
            push!(splits_per_clade[split.parent], split)
            num_split_occurrences[split] += 1

            push!(clades, split.clade1)
            num_clade_occurrences[split.clade1] += 1

            push!(clades, split.clade2)
            num_clade_occurrences[split.clade2] += 1
        end
    end

    return CCD1(num_taxa, num_trees, tip_names, root_clade, clades, splits, splits_per_clade, num_clade_occurrences, num_split_occurrences)
end

# get log probability

function log_density(ccd::CCD1, cladified_tree::CladifiedTree)
    sum(log_density(ccd, split) for split in values(cladified_tree.splits))
end

function log_density(ccd::CCD1, split::Split)
    split_occurrences = get(ccd.num_split_occurrences, split, 0.0)

    if split_occurrences == 0.0
        -Inf
    else
        log(split_occurrences / ccd.num_clade_occurrences[split.parent])
    end
end

# get most likely tree

function point_estimate(ccd::CCD1)
    most_likely_splits::Dict{Clade,Split} = Dict()
    collect_most_likely_splits!(ccd, ccd.root_clade, most_likely_splits)
    return CladifiedTree(ccd.tip_names, ccd.root_clade, most_likely_splits, Dict())
end

function collect_most_likely_splits!(ccd::CCD1, current_clade::Clade, most_likely_splits::Dict{Clade,Split})
    current_splits = ccd.splits_per_clade[current_clade]

    most_likely_split = argmax(split -> get_max_log_ccp(ccd, split), current_splits)
    most_likely_splits[most_likely_split.parent] = most_likely_split

    collect_most_likely_splits!(ccd, most_likely_split.clade1, most_likely_splits)
    collect_most_likely_splits!(ccd, most_likely_split.clade2, most_likely_splits)
end

function collect_most_likely_splits!(ccd::CCD1, current_clade::Leaf, most_likely_splits::Dict{Clade,Split})
end

function get_max_log_ccp(ccd::CCD1, split::Split)
    log_density(ccd, split) + get_max_log_ccp(ccd, split.clade1) + get_max_log_ccp(ccd, split.clade2)
end

@memoize function get_max_log_ccp(ccd::CCD1, clade::Clade)
    maximum(get_max_log_ccp(ccd, split) for split in ccd.splits_per_clade[clade])
end

function get_max_log_ccp(ccd::CCD1, clade::Leaf)
    0.0
end

# sample tree

function sample_tree(ccd::CCD1)
    sampled_splits::Dict{Clade,Split} = Dict()
    collect_sampled_splits!(ccd, ccd.root_clade, sampled_splits)
    return CladifiedTree(ccd.tip_names, ccd.root_clade, sampled_splits, Dict())
end

function collect_sampled_splits!(ccd::CCD1, current_clade::Clade, sampled_splits::Dict{Clade,Split})
    current_splits = collect(ccd.splits_per_clade[current_clade])
    weights = AnalyticWeights([exp(log_density(ccd, split)) for split in current_splits])

    sampled_split = StatsBase.sample(current_splits, weights)
    sampled_splits[sampled_split.parent] = sampled_split

    collect_sampled_splits!(ccd, sampled_split.clade1, sampled_splits)
    collect_sampled_splits!(ccd, sampled_split.clade2, sampled_splits)
end

function collect_sampled_splits!(ccd::CCD1, current_clade::Leaf, sampled_splits::Dict{Clade,Split}) end

function readable_name(ccd::Type{CCD1})
    "CCD1"
end