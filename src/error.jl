# measure error on query set achieved by MWEM instance

function maximum_error(mw::MWState)
    maximum(abs(evaluate(mw.queries,mw.synthetic) - mw.real_answers))
end

function mean_squared_error(mw::MWState)
    errors = evaluate(mw.queries,mw.synthetic) - mw.real_answers
    norm(errors)^2/length(errors)
end

function kl_divergence(a::Array{Float64,1}, b::Array{Float64,1})
    r = 0.0
    for i = 1 : length(a)
        @inbounds ai = a[i]
        @inbounds bi = b[i]
        if ai > 0
            r += ai * log(2, ai / bi)
        end
    end
    r
end
