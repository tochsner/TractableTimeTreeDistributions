abstract type AbstractClade end

mutable struct Clade <: AbstractClade
    bits::BitVector
    hash::UInt
    height::Float64

    function Clade(bits::BitVector, height::Float64=0.0)
        new(bits, hash(bits), height)
    end
end

struct Leaf <: AbstractClade
    bits::BitVector
    name::String
    hash::UInt
    height::Float64

    function Leaf(bits::BitVector, name::String, height::Float64=0.0)
        new(bits, name, hash(bits), height)
    end
end

function Leaf(leaf::Int, num_taxa::Int, name::String)
    bits = BitVector(undef, num_taxa)
    bits[leaf] = true
    return Leaf(bits, name)
end


function Clade(leaves, num_taxa::Int)
    bits = BitVector(undef, num_taxa)

    for leaf in leaves
        bits[leaf] = true
    end

    return Clade(bits)
end

function Base.:(==)(clade1::AbstractClade, clade2::AbstractClade)
    return clade1.bits == clade2.bits
end

function Base.hash(clade::AbstractClade, h::UInt)
    clade.hash
end

function Base.in(leaf::Int, clade::AbstractClade)
    clade.bits[leaf]
end

function Base.in(child_clade::AbstractClade, clade::AbstractClade)
    clade.bits .& child_clade.bits == child_clade.bits
end

function Base.union(clade1::AbstractClade, clade2::AbstractClade)
    Clade(clade1.bits .| clade2.bits)
end

function Base.union(clade1::AbstractClade, clade2::AbstractClade, height::Float64)
    Clade(clade1.bits .| clade2.bits, height)
end

function is_leaf(clade::AbstractClade)
    sum(clade.bits) == 1
end

function is_root(clade::AbstractClade)
    sum(clade.bits) == clade.bits.len
end

function size(clade::AbstractClade)
    sum(clade.bits)
end

function Base.show(io::IO, clade::AbstractClade)
    print(io, "(")
    for (i, bit) in enumerate(clade.bits)
        if bit
            print(io, i)
            print(io, ",")
        end
    end
    print(io, ")")
    print(io, ":")
    print(io, clade.height)
end
