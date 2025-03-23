abstract type CCD <: AbstractDistribution end

log_density(ccd::CCD, cladified_tree::CladifiedTree) =
    sum(log_density(ccd, split) for split in values(cladified_tree.splits))

function point_estimate(ccd::CCD)
    most_likely_splits = Dict{Clade,Split}()
    collect_most_likely_splits!(ccd, ccd.root_clade, most_likely_splits)
    return CladifiedTree(
        ccd.tip_names,
        ccd.root_clade,
        most_likely_splits,
        Dict{Clade,Float64}(),
        Dict{Clade,Float64}(),
    )
end

function collect_most_likely_splits!(ccd::CCD, current_clade::Clade, most_likely_splits::Dict{Clade,Split})
    current_splits = ccd.splits_per_clade[current_clade]

    most_likely_split = argmax(split -> get_max_log_ccp(ccd, split), current_splits)
    most_likely_splits[most_likely_split.parent] = most_likely_split

    collect_most_likely_splits!(ccd, most_likely_split.clade1, most_likely_splits)
    collect_most_likely_splits!(ccd, most_likely_split.clade2, most_likely_splits)
end
function collect_most_likely_splits!(ccd::CCD, current_clade::Leaf, most_likely_splits::Dict{Clade,Split}) end

get_max_log_ccp(ccd::CCD, split::Split) =
    log_density(ccd, split) + get_max_log_ccp(ccd, split.clade1) + get_max_log_ccp(ccd, split.clade2)
@memoize get_max_log_ccp(ccd::CCD, clade::Clade) =
    maximum(get_max_log_ccp(ccd, split) for split in ccd.splits_per_clade[clade])
get_max_log_ccp(ccd::CCD, clade::Leaf) = 0.0

function sample_tree(ccd::CCD)
    sampled_splits::Dict{Clade,Split} = Dict()
    collect_sampled_splits!(ccd, ccd.root_clade, sampled_splits)
    return CladifiedTree(ccd.tip_names, ccd.root_clade, sampled_splits, Dict(), Dict())
end
function collect_sampled_splits!(ccd::CCD, current_clade::Clade, sampled_splits::Dict{Clade,Split})
    current_splits = collect(ccd.splits_per_clade[current_clade])
    weights = AnalyticWeights([exp(log_density(ccd, split)) for split in current_splits])

    sampled_split = StatsBase.sample(current_splits, weights)
    sampled_splits[sampled_split.parent] = sampled_split

    collect_sampled_splits!(ccd, sampled_split.clade1, sampled_splits)
    collect_sampled_splits!(ccd, sampled_split.clade2, sampled_splits)
end
function collect_sampled_splits!(ccd::CCD, current_clade::Leaf, sampled_splits::Dict{Clade,Split}) end


TractableTimeTreeDistributions.entropy(ccd::CCD) = entropy(ccd, ccd.root_clade)
@memoize TractableTimeTreeDistributions.entropy(ccd::CCD, clade::Clade) = sum(
    -exp(log_density(ccd, split)) * (log_density(ccd, split) - entropy(ccd, split.clade1) - entropy(ccd, split.clade2)) for
    split in ccd.splits_per_clade[clade]
)
TractableTimeTreeDistributions.entropy(::CCD, ::Leaf) = 0.0
