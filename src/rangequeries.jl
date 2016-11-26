"""
    RangeQueries

Range queries are type of histogram query corresponding to all possible
intervals over the domain that start with the first element.
"""
type RangeQueries <: Queries
    domain::Int
end

"""
    get(queries::RangeQueries, i::Int)

Return -1/+1 indicator vector of interval [1:i].
"""
function get(queries::RangeQueries, i::Int)
    HistogramQuery(vcat(ones(i), -1.0*ones(queries.domain-i)))
end

"""
   evaluate(queries::RangeQueries, h::Histogram)

Evaluate all range queries on the given histogram.
"""
function evaluate(queries::RangeQueries, h::Histogram)
    @assert queries.domain == length(h.weights)
    answers = Array(Float64, queries.domain)
    running_sum = -1.0
    @simd for j=1:queries.domain
        @inbounds running_sum += 2.0 * h.weights[j]
        @inbounds answers[j] = running_sum
    end
    answers
end
