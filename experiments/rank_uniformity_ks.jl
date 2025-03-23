using TractableTimeTreeDistributions
using Distributions
using Logging
using Plots
using StatsBase
using HypothesisTests

theme(:wong)

function get_ref_probabilities(distribution, ref_trees)
    exp.(log_density.(Ref(distribution), ref_trees))
end

function get_sample_probabilities(distribution)
    sampled_trees = [sample_tree(distribution) for _ = 1:10000]
    sample_probabilities = log_density.(Ref(distribution), sampled_trees)
    return exp.(sample_probabilities)
end

function get_ks_statistic(ref_probabilities, sample_probabilities)
    test = ApproximateTwoSampleKSTest(ref_probabilities, sample_probabilities)
    n = test.n_x * test.n_y / (test.n_x + test.n_y)
    ks_statistic = sqrt(n) * test.δ
    return ks_statistic
end

function plot_rank_uniformity(distribution_constructors, tree_file)
    @info "Load reference trees"
    ref_trees = load_trees(tree_file)

    @info "Get distributions"
    distributions = [distribution_constructor(ref_trees) for distribution_constructor in distribution_constructors]

    @info "Get probabilities"
    ref_probabilities = get_ref_probabilities.(distributions, Ref(ref_trees))
    sample_probabilities = get_sample_probabilities.(distributions)

    @info "Get KS statistics"
    ks_statistics = get_ks_statistic.(ref_probabilities, sample_probabilities)
    n = length(ks_statistics)
    palette = cgrad(:Set2_8)
    plot(size = (750, 500), xticks = false, xticklabels = false, legend = :outerbottom)
    bar!((1:n)', ks_statistics', label = permutedims(readable_name.(distribution_constructors)), color = palette[1:n]')

    ylabel!("KS Statistic")
end

distributions = [
    TractableTimeTreeDist{CCD1,HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}}},
    TractableTimeTreeDist{CCD1,HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}},
    TractableTimeTreeDist{CCD1,ShorterBranchDist{IndependentDist{LogNormal}}},
    TractableTimeTreeDist{CCD1,ShorterBranchDist{IndependentDist{Gamma}}},
    TractableTimeTreeDist{CCD1,LastDivergenceBranchDist{IndependentDist{LogNormal},IndependentDist{LogNormal}}},
    TractableTimeTreeDist{CCD1,LastDivergenceBranchDist{IndependentDist{Gamma},IndependentDist{Gamma}}},
]

plot_rank_uniformity(distributions, "/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-10_1.trees")
