using TractableTreeDistributions
using Distributions
using Logging

function plot_rank_uniformity(distribution_constructor, tree_file)
    @info "Load trees"
    trees = load_trees(tree_file)

    @info "Fit distribution"
    distribution = distribution_constructor(trees)

    @info "Sample trees"
    sampled_trees = [sample_tree(distribution) for _ in 1:1]

    write_tree("sampled_trees.trees", sampled_trees[1])

    @info "Get probabilities of sampled trees"
    sample_probabilities = log_density.(Ref(distribution), sampled_trees)
end

plot_rank_uniformity(
    TractableTimeTreeDist{
        CCD1,
        HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}
    },
    "/Users/tobiaochsner/Documents/Thesis/TractableTreeDistributions/test/ref_trees.trees"
)