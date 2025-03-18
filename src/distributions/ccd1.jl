struct CCD1 <: CCD
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

function log_density(ccd::CCD1, split::Split)
    split_occurrences = get(ccd.num_split_occurrences, split, 0.0)

    if split_occurrences == 0.0
        -Inf
    else
        log(split_occurrences / ccd.num_clade_occurrences[split.parent])
    end
end

# utils

function readable_name(ccd::Type{CCD1})
    "CCD1"
end