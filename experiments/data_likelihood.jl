using TractableTreeDistributions
using Distributions
using Logging
using Plots
using StatsBase

theme(:wong)

function get_data_likelihood(distribution_constructor, ref_trees)
    distributions = distribution_constructor(ref_trees)
    ref_probabilities = log_density.(Ref(distributions), ref_trees)
    return sum(ref_probabilities)
end

function plot_get_data_likelihood(distributions, tree_file)
    @info "Load reference trees"
    ref_trees = load_trees(tree_file)

    @info "Get data likelihoods"
    data_likelihoods = get_data_likelihood.(distributions, Ref(ref_trees))

    @info "Plot data likelihoods"
    n = length(data_likelihoods)
    palette = cgrad(:Set2_8);
    plot(
        size=(500, 400),
        xticks=false,
        xticklabels=false,
        legend=:outerbottom,
    )
    bar!(
        (1:n)', 
        data_likelihoods', 
        label=permutedims(readable_name.(distributions)),
        color=palette[1:n]',
    )
    ylabel!("Log Data Likelihood")
end

distributions = [
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
        ShorterBranchDist{IndependentDist{Weibull}}
    },
    TractableTimeTreeDist{
        CCD1,
        LastDivergenceBranchDist{IndependentDist{LogNormal},IndependentDist{LogNormal}}
    },
    TractableTimeTreeDist{
        CCD1,
        LastDivergenceBranchDist{IndependentDist{LogNormal},IndependentDist{Gamma}}
    },
    TractableTimeTreeDist{
        CCD1,
        LastDivergenceBranchDist{IndependentDist{LogNormal},IndependentDist{Weibull}}
    }
]

plot_get_data_likelihood(
    distributions,
    "/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-50_1.trees"
)