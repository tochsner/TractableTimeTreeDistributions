using TractableTreeDistributions
using Distributions
using Logging
# using Plots
using StatsBase
using StatsPlots

theme(:wong)

function get_ref_probabilities(distribution, ref_trees)
    log_density.(Ref(distribution), ref_trees)
end

function get_sample_probabilities(distribution)
    sampled_trees = [sample_tree(distribution) for _ in 1:5000]
    sample_probabilities = log_density.(Ref(distribution), sampled_trees)
    return sample_probabilities
end

function plot_rank_uniformity(distribution_constructors, tree_file)
    @info "Load reference trees"
    ref_trees = load_trees(tree_file)

    @info "Get distributions"
    distributions = [distribution_constructor(ref_trees) for distribution_constructor in distribution_constructors]

    @info "Get probabilities"
    ref_probabilities = get_ref_probabilities.(distributions, Ref(ref_trees))
    sample_probabilities = get_sample_probabilities.(distributions)

    @info "Plot ECDFs"
    plot(
        [qqplot(exp.(ref), exp.(sample)) for (ref, sample, distribution_constructor) in zip(ref_probabilities, sample_probabilities, distribution_constructors)]...
    )
end

distribution_constructors = [
    TractableTimeTreeDist{
        CCD1,
        HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}}
    },
    TractableTimeTreeDist{
        CCD1,
        HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}
    },
    TractableTimeTreeDist{
        CCD1,
        ShorterBranchDist{IndependentDist{LogNormal}}
    },
    TractableTimeTreeDist{
        CCD1,   
        ShorterBranchDist{IndependentDist{Gamma}}
    },
    TractableTimeTreeDist{
        CCD1,
        LastDivergenceBranchDist{IndependentDist{LogNormal}, IndependentDist{LogNormal}}
    },
    TractableTimeTreeDist{
        CCD1,
        LastDivergenceBranchDist{IndependentDist{Gamma}, IndependentDist{Gamma}}
    }
]

plot_rank_uniformity(
    distribution_constructors,
    "/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-10_1.trees"
)