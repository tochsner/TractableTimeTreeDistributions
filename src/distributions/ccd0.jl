struct CCD0 <: CCD
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

function CCD0(cladified_trees::Vector{CladifiedTree})
    tip_names = cladified_trees[begin].tip_names
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

    ccd = CCD0(num_taxa, num_trees, tip_names, root_clade, clades, splits, splits_per_clade, num_clade_occurrences, num_split_occurrences)
    expand_splits!(ccd)

    return ccd
end

function expand_splits!(ccd::CCD0)
    clade_per_size = DefaultDict{Int,Set{AbstractClade}}(Set)
    for clade in ccd.clades
        push!(clade_per_size[size(clade)], clade)
    end

    for parent in ccd.clades
        if size(parent) <= 2
            # we ignore leaves and cherries
            continue
        end

        for size in 1:(size(parent)÷2)
            for potential_child in clade_per_size[size]
                if !(potential_child in parent)
                    continue
                end

                potential_sibling_bits = parent.bits .⊻ potential_child.bits
                potential_sibling = Clade(potential_sibling_bits)

                if potential_sibling in ccd.clades
                    split = Split(potential_child, potential_sibling, parent)
                    push!(ccd.splits, split)
                    push!(ccd.splits_per_clade[split.parent], split)
                end
            end
        end
    end
end

readable_name(ccd::Type{CCD0})  ="CCD0"

function log_density(ccd::CCD0, split::Split)
    if min(
        get(ccd.num_clade_occurrences, split.clade1, 0),
        get(ccd.num_clade_occurrences, split.clade2, 0),
        get(ccd.num_clade_occurrences, split.parent, 0),
    ) == 0
        return -Inf
    end

    return (
        log_clade_normalization(ccd, split.clade1) +
        log_clade_normalization(ccd, split.clade2) -
        log_clade_normalization(ccd, split.parent) +
        log_clade_credibility(ccd, split.parent)
    )
end

@memoize function log_clade_normalization(ccd::CCD0, clade::Clade)
    if is_cherry(clade)
        return log_clade_credibility(ccd, clade)
    end

    split_log_probabilities = [
        log_clade_normalization(ccd, split.clade1) + log_clade_normalization(ccd, split.clade2)
        for split in ccd.splits_per_clade[clade]
    ]

    # we use the log-sum-exp trick to compute the normalization constant in a numerically
    # more stable
    # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
    max_log_clade_credibility = maximum(split_log_probabilities)
    log_normalization = max_log_clade_credibility + log(sum(exp.(split_log_probabilities .- max_log_clade_credibility)))

    return log_clade_credibility(ccd, clade) + log_normalization
end
log_clade_normalization(ccd::CCD0, clade::Leaf) = 0.0

log_clade_credibility(ccd::CCD0, clade::Clade) = log(get(ccd.num_clade_occurrences, clade, 0) / ccd.num_trees)
