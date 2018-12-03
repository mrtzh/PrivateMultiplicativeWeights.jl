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
    (norm(errors)^2)/length(errors)
end

"""
   kl_divergence(a, b)

Compute KL-divergence between two normalized histograms a and b.
"""
function kl_divergence(a::Array{Float64, 1}, b::Array{Float64, 1})
    r = 0.0
    for i = 1 : length(a)
        @inbounds ai = a[i]
        @inbounds bi = b[i]
        if ai > 0 && bi > 0
            r += ai * log(2, ai / bi)
        end
    end
    r
end

function kl_divergence_error(mw::MWState)

    #error = kl_divergence(mw.synthetic.weights, mw.real.weights)
    error = kl_divergence(mw.real.weights, mw.synthetic.weights)
    error
end