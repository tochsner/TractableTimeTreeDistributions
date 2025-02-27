using TractableTreeDistributions
using Distributions
using Logging
using Plots
using StatsBase

theme(:wong)

function get_ecdf(distribution_constructor, ref_trees)
    @info "Fit distributions"
    distributions = distribution_constructor(ref_trees)

    @info "Sample trees"
    sampled_trees = [sample_tree(distributions) for _ in 1:5000]

    @info "Get probabilities of reference trees"
    ref_probabilities = log_density.(Ref(distributions), ref_trees)

    @info "Get probabilities of sampled trees"
    sample_probabilities = log_density.(Ref(distributions), sampled_trees)

    @info "Get ECDF"
    ecdf_func = ecdf(sample_probabilities)
    sorted_ecdf = ecdf_func(ref_probabilities) |> sort

    return sorted_ecdf
end

function plot_rank_uniformity(distribution_constructors, tree_file)
    @info "Load reference trees"
    ref_trees = load_trees(tree_file)

    @info "Get ECDFs"
    ecdfs = get_ecdf.(distribution_constructors, Ref(ref_trees))

    @info "Plot ECDFs"
    n = length(ecdfs[1])
    x_ticks = (0:0.25 * n:n, 0:0.25:1)
    y_ticks = 0:0.25:1
    plot(
        xticks=x_ticks, 
        yticks=y_ticks,
        size=(500, 500),
        legend=:bottomright
    )
    plot!([0, length(ecdfs[1])], [0, 1.0]; c = :black, lw = 1, label=nothing)
    plot!(ecdfs, label=permutedims(readable_name.(distribution_constructors)))
    
    xlabel!("Model CDF")
    ylabel!("MCMC CDF")
end

distributions = [
    CCD1
    # TractableTimeTreeDist{
    #     CCD1,
    #     HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}}
    # },
    # TractableTimeTreeDist{
    #     CCD1,
    #     HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}
    # },
    # TractableTimeTreeDist{
    #     CCD1,
    #     ShorterBranchDist{IndependentDist{LogNormal}}
    # },
    # TractableTimeTreeDist{
    #     CCD1,   
    #     ShorterBranchDist{IndependentDist{Gamma}}
    # },
    # TractableTimeTreeDist{
    #     CCD1,
    #     LastDivergenceBranchDist{IndependentDist{LogNormal}, IndependentDist{LogNormal}}
    # },
    # TractableTimeTreeDist{
    #     CCD1,
    #     LastDivergenceBranchDist{IndependentDist{Gamma}, IndependentDist{Gamma}}
    # },
]

plot_rank_uniformity(
    distributions,
    "/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-20_1.trees"
)