abstract type AbstractClade end

struct Clade <: AbstractClade
    bits::BitVector
end

struct Leaf <: AbstractClade
    bits::BitVector
end

function Leaf(leaf::Int, num_taxa::Int)
    bits = BitVector(undef, num_taxa)
    bits[leaf] = true
    return Leaf(bits)
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
    hash(clade.bits, h)
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

function isLeaf(clade::AbstractClade)
    sum(clade.bits) == 1
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
end
