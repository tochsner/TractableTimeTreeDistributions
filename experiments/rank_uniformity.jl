using TractableTreeDistributions
using Distributions
using Logging
using Plots
using StatsBase

theme(:wong)

function get_ecdf(distribution_constructor, ccd, ref_trees)
    distribution = TractableTimeTreeDist(
        ccd, distribution_constructor(cladify_tree.(ref_trees))
    )

    ref_probabilities = log_density.(Ref(distribution), ref_trees)
    replace!(ref_probabilities, NaN => -Inf)
    
    sampled_trees = [sample_tree(distribution) for _ in 1:5000]
    sample_probabilities = log_density.(Ref(distribution), sampled_trees)

    ecdf_func = ecdf(sample_probabilities)
    sorted_ecdf = (1.0 .- ecdf_func(ref_probabilities)) |> sort

    return sorted_ecdf
end

function plot_rank_uniformity(distribution_constructors, tree_file)
    @info "Load reference trees"
    ref_trees = load_trees(tree_file)

    @info "Fit CCD"
    ccd = CCD1(ref_trees)

    @info "Get ECDFs"
    ecdfs = get_ecdf.(distribution_constructors, Ref(ccd), Ref(ref_trees))

    @info "Plot ECDFs"
    n = length(ecdfs[1])
    plot([0, length(ecdfs[1])], [0, 1.0]; c = :black, lw = 1, label=nothing)
    plot!(ecdfs, label=permutedims(readable_name.(distribution_constructors)))
    plot!(
        xticks=0:0.25 * n:n, 
        yticks=0:0.25:1,
        size=(750, 500),
        legend=:topright,
        yflip = true,
        yformatter = y -> "$(100 * y) %",
        xformatter = x -> "$(100 * x / n) %"
    )
    xlabel!("Fraction of Reference Trees in Credible Set")
    ylabel!("Credible Set")
end

distribution_constructors = [
    # HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}}
    # HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}
    # ShorterBranchDist{IndependentDist{LogNormal}}
    # ShorterBranchDist{IndependentDist{Gamma}}
    # LastDivergenceBranchDist{IndependentDist{LogNormal}, IndependentDist{LogNormal}}
    # LastDivergenceBranchDist{IndependentDist{Gamma}, IndependentDist{Gamma}}
    HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}},
    LastDivergenceBranchDist{IndependentDist{Gamma},IndependentDist{Gamma}}
]

plot_rank_uniformity(
    distribution_constructors,
    "/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-50_40.trees"
)