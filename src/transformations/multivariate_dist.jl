struct MultivariateDist{D} <: Transform
    distributions::Dict{Clade,D}
    MultivariateDist() = new(Dict())
end

function fit(transformation::MultivariateDist{D}, trees::Vector{ParameterizedTree})
    # todo
end

function sample(transformation::MultivariateDist{D}, tree::CladifiedTree)::ParameterizedTree
    # todo
end

function log_density(transformation::MultivariateDist{D}, tree::ParameterizedTree)
    # todo
end