import REPL
using REPL.TerminalMenus

function run_interactive_workflow(trees_file::String)
    distributions = [
        HeightRatioDist{IndependentDist{LogNormal},IndependentDist{LogitNormal}}
        HeightRatioDist{IndependentDist{LogNormal},IndependentDist{Beta}}
        ShorterBranchDist{IndependentDist{LogNormal}}
        ShorterBranchDist{IndependentDist{Gamma}}
        ShorterBranchDist{IndependentDist{Weibull}}
        LastDivergenceBranchDist{IndependentDist{LogNormal},IndependentDist{LogNormal}}
        LastDivergenceBranchDist{IndependentDist{Gamma},IndependentDist{Gamma}}
        LastDivergenceBranchDist{IndependentDist{Weibull},IndependentDist{Weibull}}
    ]
    num_samples = 10_000
    train_fraction = 0.75
    burn_in_fraction = 0.1

    @info "Load and prepare MCMC trees"

    trees = load_trees(trees_file)
    trees = trees[ceil(Int, length(trees) * burn_in_fraction):end]
    cladified_trees = cladify_tree.(trees)

    @info "Calculate tree ESS and subsample trees down to ESS"

    tree_ess = get_ess(cladified_trees)
    trees_subsampled = sample(cladified_trees, floor(Int, tree_ess); replace = false, ordered = true)

    @info "ESS is $(tree_ess)"

    num_train_trees = floor(Int, length(trees_subsampled) * train_fraction)
    trees_train = trees_subsampled[1:num_train_trees]
    trees_val = trees_subsampled[num_train_trees+1:end]

    @info "Fit distributions on train set"

    ccd_train = CCD1(trees_train)
    distributions_train = [distribution(trees_train) for distribution in distributions]

    @info "Validate distributions on val set"

    log_ccd_densities_val = log_density.(Ref(ccd_train), trees_val)

    log_data_likelihoods_val = []
    for (i, distribution) in enumerate(distributions_train)
        log_densities_val = log_density.(Ref(distribution), trees_val) .+ log_ccd_densities_val
        finite_log_densities_val = filter(x -> -Inf < x, log_densities_val)

        log_data_likelihood_val = 0 < length(finite_log_densities_val) ? finite_log_densities_val |> sum : -Inf
        push!(log_data_likelihoods_val, log_data_likelihood_val)
    end

    distribution_names = readable_name.(distributions)
    max_distribution_name_length = maximum(length.(distribution_names))

    println("-"^length("Log Data Likelihoods"))
    println("Log Data Likelihoods")
    println("-"^length("Log Data Likelihoods"))

    for (i, log_data_likelihood) in sort(enumerate(log_data_likelihoods_val) |> collect, by = x -> x[2], rev = true)
        @info "$(rpad(distribution_names[i] * ":", max_distribution_name_length))\t$(log_data_likelihood)"
    end

    println("-"^length("Model Selection"))
    println("Model Selection")
    println("-"^length("Model Selection"))

    best_distribution = distributions[argmax(log_data_likelihoods_val)]
    println("We recommend to use the $(readable_name(best_distribution)) distribution.")
    println("Do you want to use it? [y/n]")

    answer = readline()
    if answer == "y"
        selected_distribution = best_distribution
    else
        menu = RadioMenu(distribution_names)
        choice = request("Select one of the following distributions:", menu)
        selected_distribution = distributions[choice]
    end

    @info "You selected the $(readable_name(selected_distribution)) distribution."

    println("-"^length("Tractable Operations"))
    println("Tractable Operations")
    println("-"^length("Tractable Operations"))

    operations_menu = RadioMenu(
        [
            "Create a credible-region plot"
            "Create a log-likelihood plot"
            "Create a point estimate"
            "Get the likelihood of a given tree"
            "Test if a given tree is in the smallest 95%-credible region"
        ],
    )
    operations_choice = request("Select one or more tasks to perform using the selected distribution:", operations_menu)
end
