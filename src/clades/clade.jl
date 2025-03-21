abstract type AbstractClade end

struct Clade <: AbstractClade
    bits::BitVector
    hash::UInt

    function Clade(bits::BitVector)
        new(bits, hash(bits))
    end
end
function Clade(leaves, num_taxa::Int)
    bits = BitVector(undef, num_taxa)

    for leaf in leaves
        bits[leaf] = true
    end

    return Clade(bits)
end

struct Leaf <: AbstractClade
    bits::BitVector
    name::String
    hash::UInt

    function Leaf(bits::BitVector, name::String)
        new(bits, name, hash(bits))
    end
end
function Leaf(leaf::Int, num_taxa::Int, name::String)
    bits = falses(num_taxa)
    bits[leaf] = true
    Leaf(bits, name)
end

Base.:(==)(clade1::AbstractClade, clade2::AbstractClade) = clade1.hash == clade2.hash && clade1.bits == clade2.bits
Base.hash(clade::AbstractClade, h::UInt) = clade.hash

Base.in(leaf::Int, clade::AbstractClade) = clade.bits[leaf]
Base.in(child_clade::AbstractClade, clade::AbstractClade) = (clade.bits .& child_clade.bits) == child_clade.bits

Base.union(clade1::AbstractClade, clade2::AbstractClade) = Clade(clade1.bits .| clade2.bits)

size(clade::AbstractClade) = sum(clade.bits)
is_leaf(clade::AbstractClade) = size(clade) == 1
is_cherry(clade::AbstractClade) = size(clade) == 2
is_root(clade::AbstractClade) = size(clade) == clade.bits.len

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
