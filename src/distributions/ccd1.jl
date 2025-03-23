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

CCD1(trees::Vector{Tree}) = cladify_tree.(trees) |> CCD1

function CCD1(cladified_trees::Vector{CladifiedTree})
    tip_names = cladified_trees[begin].tip_names
    num_taxa = length(tip_names)
    num_trees = length(cladified_trees)
    root_clade = Clade(1:num_taxa, num_taxa)

    clades = Set{AbstractClade}()
    num_clade_occurrences = DefaultDict{AbstractClade,Int64}(0)

    splits = Set{Split}()
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

    return CCD1(
        num_taxa,
        num_trees,
        tip_names,
        root_clade,
        clades,
        splits,
        splits_per_clade,
        num_clade_occurrences,
        num_split_occurrences,
    )
end

readable_name(ccd::Type{CCD1}) = "CCD1"

function log_density(ccd::CCD1, split::Split)
    split_occurrences = get(ccd.num_split_occurrences, split, 0.0)

    if split_occurrences == 0.0
        -Inf
    else
        log(split_occurrences / ccd.num_clade_occurrences[split.parent])
    end
end
