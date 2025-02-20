struct LogTransform <: Transform
    wrapped::Transform
end

function transform(transformation::Transform, tree::ParameterizedTree)
    transform(transformation.wrapped, ParameterizedTree(log.(tree.parameters), tree))
end

function invert(transformation::Transform, tree::ParameterizedTree)
    invert(transformation.wrapped, ParameterizedTree(exp.(tree.parameters), tree))
end

function log_density(transformation::Transform, tree::ParameterizedTree)
    log_density(transformation.wrapped, transform(transformation, tree)) - sum(tree.parameters)
end