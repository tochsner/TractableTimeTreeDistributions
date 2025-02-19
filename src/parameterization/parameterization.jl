abstract type Parameterization end

function parameterize(parameterization::Parameterization, trees::Vector{CladifiedTree})::Vector{ParameterizedTree}
    [parameterize(parameterization, tree) for tree in trees]
end

function set_heights!(parameterization::Parameterization, parameterized_trees::Vector{ParameterizedTree})::Vector{CladifiedTree}
    [set_heights!(parameterization, tree) for tree in parameterized_trees]
end