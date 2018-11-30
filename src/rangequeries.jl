"""
    RangeQuery

Range query is a type of histogram query corresponding to an
intervals over the domain. i -> j (1 <= i <= j <= domain).
"""
struct RangeQuery <: Query
    interval::Interval
    domain::Int # the domain size
end

"""
    RangeQueriesRangeQueries

Range queries are type of histogram query corresponding to
intervals over the domain.
"""
struct RangeQueries <: Queries
    intervals::Array{Interval, 1}
    domain::Int # the domain size
end

"""
    SeriesRangeQueriesRangeQueries

Range queries are type of histogram query corresponding to all possible
intervals over the domain that start with the first element.
"""
struct SeriesRangeQueries <: Queries
    domain::Int # the domain size
end

"""
    get(interval::Interval, domain::Int64)

Return 0/1 indicator vector for the interval in the size of the domain.
"""
function get_query_vector_from_interval(interval::Interval, domain::Int64)
    query_vector = zeros(domain)
    start_index, end_index = interval
    query_vector[start_index : end_index] .= 1
    query_vector
end

"""
    get(query::RangeQuery)

Return 0/1 indicator vector of interval [j:k] of the query. (query.interval = (j,k))
"""
function get(query::RangeQuery)
    HistogramQuery(get_query_vector_from_interval(query.interval, query.domain))
end

"""
    get(queries::RangeQueries, i::Int)

Return 0/1 indicator vector of interval [j:k] of the i'th range query (queries[i] = (j,k)).
"""
function get(queries::RangeQueries, i::Int)
    HistogramQuery(get_query_vector_from_interval(queries.intervals[i], queries.domain))
end

"""
    get(queries::RangeQueries, i::Int)

Return 0/1 indicator vector of interval [1:i].
"""
function get(queries::SeriesRangeQueries, i::Int)
    HistogramQuery(vcat(ones(i), zeros(queries.domain-i)))#-1.0*ones(queries.domain-i)))
end

"""
   evaluate(query::RangeQuery, h::Histogram)

Evaluate the query on the given histogram.
"""
function evaluate(query::RangeQuery, h::Histogram)
    @assert query.domain == length(h.weights)
    answer = evaluate(get(query), h)
    answer
end

"""
   evaluate(queries::RangeQueries, h::Histogram)

Evaluate all the given range queries on the given histogram.
"""
function evaluate(queries::RangeQueries, h::Histogram)
    @assert queries.domain == length(h.weights)
    num_queries = length(queries.intervals)
    answers = Array{Float64}(undef, num_queries)
    @simd for j=1:num_queries
        @inbounds answers[j] = evaluate(get(queries, j), h)
    end
    answers
end

"""
   evaluate(queries::RangeQueries, h::Histogram)

Evaluate all range queries on the given histogram.
"""
function evaluate(queries::SeriesRangeQueries, h::Histogram)
    @assert queries.domain == length(h.weights)
    answers = zeros(Float64, queries.domain)
    running_sum = 0
    @simd for j=1:queries.domain
        @inbounds running_sum += h.weights[j]
        @inbounds answers[j] = running_sum
    end
    answers
end

function queriesMatrix(queries::SeriesRangeQueries)

    query_matrix = zeros(queries.domain, queries.domain)
    @simd for i=1:queries.domain
        @simd for j=1:i
            @inbounds query_matrix[i,j] = 1
        end
    end
    query_matrix
end
