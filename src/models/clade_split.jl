struct CladeSplit
    clade1::Clade
    clade2::Clade
end

function Base.:(==)(clade_split1::CladeSplit, clade_split2::CladeSplit)
    return clade_split1.clade1 == clade_split2.clade1 && clade_split1.clade2 == clade_split2.clade2 ||
           clade_split1.clade1 == clade_split2.clade2 && clade_split1.clade2 == clade_split2.clade1
end

function Base.in(clade::Clade, clade_split::CladeSplit)
    return clade_split.clade1 == clade || clade_split.clade2 == clade
end