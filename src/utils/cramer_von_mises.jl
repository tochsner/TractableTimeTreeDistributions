function cramer_von_mises(x::AbstractVector{<:Real})
    n = length(x)
    x = sort(x)
    (1 / 12n + sum((x .- (2(1:n) .- 1) ./ 2n) .^ 2)) / n
end
