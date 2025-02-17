abstract type CCD end

struct CCD1 <: CCD
    cladifies_trees::Vector{CladifiedTree}
    num_taxa::Int

    root_clade::Clade
    clades::Set{Clade}
    splits::Set{CladeSplit}

    num_occurrences::Dict{Union{Clade, CladeSplit}, Int}
end

function CCD1(trees::Vector{Tree})
    num_taxa = length(tiplabels(trees[1]))
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

    return CCD1(cladified_trees, num_taxa, root_clade, clades, splits, num_occurrences)
end

function get_probability(ccd::CCD1, tree::Tree)
end

function sample(ccd::CCD1)
end
