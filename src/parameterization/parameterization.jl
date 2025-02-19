abstract type Parameterization end

function parameterize(parameterization::Parameterization, trees::Vector{CladifiedTree})::Vector{ParameterizedTree}
    [parameterize(parameterization, tree) for tree in trees]
end