struct Clade
    bits::BitVector
end

function Clade(leaf::Int, num_taxa::Int)
    bits = BitVector(undef, num_taxa)
    bits[leaf] = true
    return Clade(bits)
end


function Clade(leaves, num_taxa::Int)
    bits = BitVector(undef, num_taxa)
    
    for leaf in leaves
        bits[leaf] = true
    end

    return Clade(bits)
end

function Base.:(==)(clade1::Clade, clade2::Clade)
    return clade1.bits == clade2.bits
end

function Base.hash(clade::Clade, h::UInt)
    return hash(clade.bits, h)
end

function Base.in(leaf::Int, clade::Clade)
    clade.bits[leaf]
end

function Base.in(child_clade::Clade, clade::Clade)
    clade.bits .& child_clade.bits == child_clade.bits
end

function Base.union(clade1::Clade, clade2::Clade)
    return Clade(clade1.bits .| clade2.bits)
end
