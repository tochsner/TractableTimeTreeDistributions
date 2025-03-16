struct Split
    clade1::AbstractClade
    clade2::AbstractClade
    parent::Clade
    
    hash::UInt
    
    function Split(clade1::AbstractClade, clade2::AbstractClade, parent::Clade)
        if findfirst(x -> x, clade1.bits) < findfirst(x -> x, clade2.bits)
            return new(clade1, clade2, parent, hash(clade1) + hash(clade2))
        else
            return new(clade2, clade1, parent, hash(clade1) + hash(clade2))
        end
    end
end

function Split(clade1::AbstractClade, clade2::AbstractClade)
    Split(clade1, clade2, union(clade1, clade2))
end

function Base.:(==)(split1::Split, split2::Split)
    split1.hash == split2.hash && split1.clade1.hash == split2.clade1.hash && split1.clade2.hash == split2.clade2.hash && split1.clade1 == split2.clade1 && split1.clade2 == split2.clade2
end

function Base.hash(split::Split, h::UInt)
    split.hash
end

function Base.in(clade::AbstractClade, split::Split)
    split.clade1 == clade || split.clade2 == clade
end

function construct_tree_from_splits(splits::Vector{Split})
    splits
end