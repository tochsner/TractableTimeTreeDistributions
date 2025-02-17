struct CladeSplit
    clade1::Clade
    clade2::Clade

    function CladeSplit(clade1::Clade, clade2::Clade)
        if findfirst(x -> x, clade1.bits) < findfirst(x -> x, clade2.bits)
            return new(clade1, clade2)
        else
            return new(clade2, clade1)
        end
    end
end

function Base.:(==)(clade_split1::CladeSplit, clade_split2::CladeSplit)
    return clade_split1.clade1 == clade_split2.clade1 && clade_split1.clade2 == clade_split2.clade2
end

function Base.hash(clade_split::CladeSplit, h::UInt)
    return hash(clade_split.clade1, hash(clade_split.clade2, h))
end

function Base.in(clade::Clade, clade_split::CladeSplit)
    return clade_split.clade1 == clade || clade_split.clade2 == clade
end