abstract type AbstractDistribution end

log_density(distribution::AbstractDistribution, tree::Tree) = log_density(distribution, cladify_tree(tree))
readable_name(distribution) = nameof(distribution)

function average_parameter_correlation(cladified_trees::Vector{CladifiedTree}; n = 100)
    subsample = sample(cladified_trees, n)

    all_clades = Set{Clade}()
    for tree in subsample
        union!(all_clades, keys(tree.heights))
    end

    clade_pairs = Dict{Clade,Dict{Clade,Vector{Tuple{Float64,Float64}}}}()
    for tree in subsample
        for clade1 in keys(tree.heights)
            if !haskey(clade_pairs, clade1)
                clade_pairs[clade1] = Dict{Clade,Vector{Tuple{Float64,Float64}}}()
            end

            for clade2 in keys(tree.heights)
                if !haskey(clade_pairs[clade1], clade2)
                    clade_pairs[clade1][clade2] = Vector{Tuple{Float64,Float64}}()
                end

                push!(clade_pairs[clade1][clade2], (tree.heights[clade1], tree.heights[clade2]))
            end
        end
    end
end
