using Distributions

"""
A wrapper for the fit_mle function from Distributions.jl with
    a better initial guess for the Weibull distribution.

    This is based on the paper "Approximate Formula of Coefficient
    of Variation for Weibull Distribution" by S. Tanaka & M. Ichikawa.
"""
function Distributions.fit_mle(::Type{<:Weibull}, x::AbstractArray{<:Real})
    coefficient_of_variation = std(x) / mean(x)
    alpha0 = coefficient_of_variation^(-1.075)
    fit_mle(Weibull, x, alpha0=alpha0)
end

"""
A numerical bisection method to find the mode of a univariate
LogitNormal distribution.
"""
function Distributions.mode(d::LogitNormal; tol::Real=1e-16, maxiter::Int=1000)
    logit(x) = log(x) - log1p(-x)
    f(x) = d.σ^2 * (2x - 1) + d.μ - logit(x)

    extrema = 1e-10

    left, right = extrema, 1.0 - extrema

    Δx = right - left
    ϵ = Δx / left

    N = 0

    while ϵ > tol && N < maxiter
        middle = (left + right) / 2

        fleft = f(left)
        fmiddle = f(middle)

        if sign(fleft) == sign(fmiddle)
            left = middle
        else
            right = middle
        end

        Δx = abs(right - left)
        ϵ = Δx / left

        N += 1
    end

    if pdf(d, left) < pdf(d, extrema) || pdf(d, left) < pdf(d, 1.0 - extrema)
        throw(ArgumentError("LogitNormal distribution is bivariate"))
    else
        # distribution is univariate
        return left
    end
end

"""
A simple Monte Carlo estimator for the mean of a LogitNormal distribution.
"""
function Distributions.mean(d::LogitNormal; num_samples::Int=10_000)
    return mean(rand(d, num_samples))
end