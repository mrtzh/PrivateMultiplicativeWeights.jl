"""
    maximum_error(mw::MWState)

Compute maximum error of synthetic data on the query set. Result is not
differentially private.
"""
function maximum_error(mw::MWState)
    norm(evaluate(mw.queries, mw.synthetic) - mw.real_answers, Inf)
end

"""
    mean_squared_error(mw::MWState)

Compute mean_squared_error of synthetic data on the query set. Result is not
differentially private.
"""
function mean_squared_error(mw::MWState)

    errors = evaluate(mw.queries, mw.synthetic) - mw.real_answers
    (norm(errors)^2)/length(errors)
end

function total_squared_error(mw::MWState)

    errors = evaluate(mw.queries, mw.synthetic) - mw.real_answers
    norm(errors)^2
end

"""
    kl_divergence(a, b)

Compute the base-2 KL divergence `D_KL(a || b)` between two nonnegative vectors.
The vectors are expected to be normalized distributions. Returns `Inf` when
`a` assigns positive mass where `b` assigns zero mass.
"""
function kl_divergence(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    length(a) == length(b) || throw(DimensionMismatch("distributions must have the same length"))
    all(isfinite, a) && all(isfinite, b) || throw(ArgumentError("distributions must be finite"))
    all(>=(0), a) && all(>=(0), b) || throw(ArgumentError("distributions must be nonnegative"))

    r = 0.0
    for i in eachindex(a, b)
        @inbounds ai = a[i]
        @inbounds bi = b[i]
        if ai > 0
            bi > 0 || return Inf
            r += ai * log2(ai / bi)
        end
    end
    r
end

"""
    kl_divergence_error(mw::MWState)

Compute `D_KL(mw.real.weights || mw.synthetic.weights)` for an explicit
histogram MWEM result. This diagnostic uses exact source information and is not
differentially private.
"""
function kl_divergence_error(mw::MWState)

    #error = kl_divergence(mw.synthetic.weights, mw.real.weights)
    error = kl_divergence(mw.real.weights, mw.synthetic.weights)
    error
end
