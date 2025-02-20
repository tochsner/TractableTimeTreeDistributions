struct IndependentDist{D} <: Transform
    distributions::Dict{Clade,D}
    IndependentDist() = new(Dict())
end

function fit(transformation::IndependentDist{D}, trees::Vector{ParameterizedTree})
    empty!(transformation.distributions)

    observations = DefaultDict{Clade,Vector[Float64]}([])

    for tree in trees
        for (clade, parameter) in keys(tree.parameters)
            push!(observations[clade], parameter)
        end
    end

    for (clade, observations) in observations
        transformation.distributions[clade] = fit_mle(D, observations)
    end
end

function sample(transformation::IndependentDist{D}, tree::CladifiedTree)::ParameterizedTree
    parameters = {
        clade => rand(dist) for (clade, dist) in transformation.distributions
    }
    return ParameterizedTree(parameters, tree)
end

function log_density(transformation::IndependentDist{D}, tree::ParameterizedTree)
    sum(logpdf(get(transformation.distributions, clade, 0), param) for (clade, param) in tree.parameters)
end