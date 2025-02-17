abstract type CCD end

struct CCD1 <: CCD
    trees::Vector
end

function CCD1(trees::Vector{Tree})
    println("Fitting...")
end

function get_probability(ccd::CCD1, tree::Tree)
end

function sample(ccd::CCD1)
end
