struct ParameterizedTree
    tip_names::Vector{String}
    root::AbstractClade
    splits::Dict{Clade,Split}
    parameters::Dict{Clade,Float64}
end

function ParameterizedTree(parameters::Dict{Clade,Float64}, cladified_tree::CladifiedTree)
    ParameterizedTree(cladified_tree.tip_names, cladified_tree.root, cladified_tree.splits, parameters)
end

function to_cladified_tree(parameterized_tree::ParameterizedTree)::CladifiedTree
    CladifiedTree(parameterized_tree.tip_names, parameterized_tree.root, parameterized_tree.splits)
end