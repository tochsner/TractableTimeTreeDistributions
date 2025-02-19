struct CladeSplit
    clade1::AbstractClade
    clade2::AbstractClade
    parent::AbstractClade
    
    hash::UInt
    
    function CladeSplit(clade1::AbstractClade, clade2::AbstractClade, parent::AbstractClade)
        if findfirst(x -> x, clade1.bits) < findfirst(x -> x, clade2.bits)
            return new(clade1, clade2, parent, hash(clade1, hash(clade2)))
        else
            return new(clade2, clade1, parent, hash(clade2, hash(clade1)))
        end
    end
end

function CladeSplit(clade1::AbstractClade, clade2::AbstractClade)
    CladeSplit(clade1, clade2, union(clade1, clade2))
end

function Base.:(==)(clade_split1::CladeSplit, clade_split2::CladeSplit)
    clade_split1.clade1 == clade_split2.clade1 && clade_split1.clade2 == clade_split2.clade2
end

function Base.hash(clade_split::CladeSplit, h::UInt)
    clade_split.hash
end

function Base.in(clade::AbstractClade, clade_split::CladeSplit)
    clade_split.clade1 == clade || clade_split.clade2 == clade
end

function construct_tree_from_splits(splits::Vector{CladeSplit})
    splits
end