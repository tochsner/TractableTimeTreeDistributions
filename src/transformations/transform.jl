abstract type Transform end

function transform(transformation::Transform, trees)
    [transform(transformation, tree) for tree in trees]
end

function transform(transformation::Transform, tree::ParameterizedTree) end

function invert(transformation::Transform, trees)
    [invert(transformation, tree) for tree in trees]
end

function invert(transformation::Transform, tree::ParameterizedTree) end

function fit(transformation::Transform, trees::Vector{ParameterizedTree})
    fit(transformation.wrapped, transform(transformation, trees))
end

function sample(transformation::Transform, tree::CladifiedTree)::ParameterizedTree
    invert(transformation, sample(transformation.wrapped, tree))
end

function log_density(transformation::Transform, tree::ParameterizedTree)
    log_density(transformation.wrapped, transform(transformation, tree))
end