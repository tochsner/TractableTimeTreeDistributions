struct LogTransform <: Transform
    wrapped::Transform
end

function transform(transformation::LogTransform, tree::ParameterizedTree)
    transform(transformation.wrapped, ParameterizedTree(log.(tree.parameters), tree))
end

function invert(transformation::LogTransform, tree::ParameterizedTree)
    invert(transformation.wrapped, ParameterizedTree(exp.(tree.parameters), tree))
end

function log_abs_det_jacobian(transformation::LogTransform, tree::ParameterizedTree)
    sum(tree.parameters)
end