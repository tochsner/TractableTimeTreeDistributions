using TractableTreeDistributions
using Distributions
using Logging
using StatsBase
using HypothesisTests

trees_file = "/Users/tobiaochsner/Documents/Thesis/TractableTreeDistributions/test/ref_trees.trees"
output_dir = "/Users/tobiaochsner/Documents/Thesis/TractableTreeDistributions"

if length(ARGS) == 2
    trees_file = ARGS[1]
    output_dir = ARGS[2]
end

distributions = [
    HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}}
    HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}
    ShorterBranchDist{IndependentDist{LogNormal}}
    ShorterBranchDist{IndependentDist{Gamma}}
    ShorterBranchDist{IndependentDist{Weibull}}
    LastDivergenceBranchDist{IndependentDist{LogNormal}, IndependentDist{LogNormal}}
    LastDivergenceBranchDist{IndependentDist{Gamma}, IndependentDist{Gamma}}
    LastDivergenceBranchDist{IndependentDist{Weibull}, IndependentDist{Weibull}}
]
num_samples = 10_000

@info "Load and prepare reference trees"

trees = load_trees(trees_file)
cladified_trees = cladify_tree.(trees)

@info "Calculate tree ESS and subsample trees down to ESS"

tree_ess = get_ess(cladified_trees)
trees_subsampled = sample(cladified_trees, floor(Int, tree_ess), replace=false)

@info "Split into train and val set for validation"

half = length(trees_subsampled) ÷ 2
trees_train = trees_subsampled[1:half]
trees_val = trees_subsampled[half+1:end]

@info "Fit distributions on train set"

ccd_train = CCD1(trees_train)
distributions_train = [
    distribution(trees_train)
    for distribution in distributions
]

@info "Validate distributions"

log_data_likelihoods_val = []

ad_test_statistics_val = []
ad_test_p_values_val = []

credible_sets_val = []

for distribution in distributions_train
    log_densities_val = log_density.(Ref(distribution), trees_val) .+ log_density.(Ref(ccd_train), trees_val)
    log_densities_val = filter(x -> -Inf < x, log_densities_val)
    
    log_data_likelihood_val = sum(log_densities_val)
    push!(log_data_likelihoods_val, log_data_likelihood_val)

    samples = [sample_tree(distribution, sample_tree(ccd_train)) for _ in 1:num_samples]
    log_densities_samples = log_density.(Ref(distribution), samples) .+ log_density.(Ref(ccd_train), samples)
    log_densities_samples = filter(x -> -Inf < x, log_densities_samples)

    ad_test = KSampleADTest(log_densities_val, log_densities_samples)
    ad_test_statistic = ad_test.A²k
    ad_test_p_value = pvalue(ad_test)
    push!(ad_test_statistics_val, ad_test_statistic)
    push!(ad_test_p_values_val, ad_test_p_value)

    ecdf_func = ecdf(log_densities_samples)
    credible_sets = (1.0 .- ecdf_func(log_densities_val))
    push!(credible_sets_val, credible_sets)
end

@info "Get point estimates based on all ESS trees"

ccd_subsampled = CCD1(trees_subsampled)
distributions_subsampled = [
    distribution(trees_subsampled)
    for distribution in distributions
]

ccd_map_tree = point_estimate(ccd_subsampled)
point_estimates = point_estimate.(distributions_subsampled, Ref(ccd_map_tree))
mrca_point_estimate = mrca_tree(ccd_map_tree, trees_subsampled)

@info "Store results on disk"

base_file_name = basename(trees_file) |> splitext |> first

statistics_file = joinpath(output_dir, "$(base_file_name)_stats.log")
open(statistics_file, "w") do io
    println(io, "distribution;metric;value")
    println(io, "-;tree_ess;$(tree_ess)")

    for i in eachindex(distributions)
        println(io, "$(readable_name(distributions[i]));log_data_likelihood;$(log_data_likelihoods_val[i])")
        println(io, "$(readable_name(distributions[i]));anderson_darling_test_statistics;$(ad_test_statistics_val[i])")
        println(io, "$(readable_name(distributions[i]));anderson_darling_p_value;$(ad_test_p_values_val[i])")
    end
end

credible_sets_file = joinpath(output_dir, "$(base_file_name)_credible_sets.log")
open(credible_sets_file, "w") do io
    println(io, "distribution;sample_credible_set")

    for i in eachindex(distributions)
        for credible_set in credible_sets_val[i]
            println(io, "$(readable_name(distributions[i]));$(credible_set)")
        end
    end
end

point_estimate_file = joinpath(output_dir, "$(base_file_name)_point_estimate.trees")
open(point_estimate_file, "w") do io
    write_tree(
        io, 
        push!(map(x -> "'$(x)'", readable_name.(distributions)), "'MRCA'"), 
        push!(point_estimates, mrca_point_estimate)
    )
end