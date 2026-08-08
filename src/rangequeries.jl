"""
    RangeQuery

Range query is a histogram query corresponding to an interval over the domain,
`i:j`, where `1 <= i <= j <= domain`.
"""
struct RangeQuery <: Query
    interval::Interval
    domain::Int # the domain size

    function RangeQuery(interval::Tuple{<:Integer, <:Integer}, domain::Integer)
        domain > 0 || throw(ArgumentError("domain must be positive"))
        normalized = (Int(interval[1]), Int(interval[2]))
        1 <= normalized[1] <= normalized[2] <= domain ||
            throw(ArgumentError("interval must satisfy 1 <= start <= stop <= domain"))
        new(normalized, Int(domain))
    end
end

"""
    RangeQueries

Range queries are type of histogram query corresponding to
intervals over the domain.
"""
struct RangeQueries <: Queries
    intervals::Vector{Interval}
    domain::Int # the domain size

    function RangeQueries(intervals::AbstractVector{<:Tuple{<:Integer, <:Integer}}, domain::Integer)
        domain > 0 || throw(ArgumentError("domain must be positive"))
        isempty(intervals) && throw(ArgumentError("at least one interval is required"))
        normalized = Interval[(Int(interval[1]), Int(interval[2])) for interval in intervals]
        all(interval -> 1 <= interval[1] <= interval[2] <= domain, normalized) ||
            throw(ArgumentError("every interval must satisfy 1 <= start <= stop <= domain"))
        new(normalized, Int(domain))
    end
end

"""
    SeriesRangeQueries

Range queries are type of histogram query corresponding to all possible
intervals over the domain that start with the first element.
"""
struct SeriesRangeQueries <: Queries
    domain::Int # the domain size

    function SeriesRangeQueries(domain::Integer)
        domain > 0 || throw(ArgumentError("domain must be positive"))
        new(Int(domain))
    end
end

"""
    get_query_vector_from_interval(interval, domain)

Return 0/1 indicator vector for the interval in the size of the domain.
"""
function get_query_vector_from_interval(interval::Tuple{<:Integer, <:Integer}, domain::Integer)
    domain > 0 || throw(ArgumentError("domain must be positive"))
    1 <= interval[1] <= interval[2] <= domain ||
        throw(ArgumentError("interval must satisfy 1 <= start <= stop <= domain"))
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
    query.domain == length(h.weights) || throw(DimensionMismatch("query and histogram dimensions differ"))
    answer = evaluate(get(query), h)
    answer
end

"""
   evaluate(queries::RangeQueries, h::Histogram)

Evaluate all the given range queries on the given histogram.
"""
function evaluate(queries::RangeQueries, h::Histogram)
    queries.domain == length(h.weights) || throw(DimensionMismatch("query and histogram dimensions differ"))
    num_queries = length(queries.intervals)
    answers = Vector{Float64}(undef, num_queries)
    @simd for j=1:num_queries
        @inbounds answers[j] = evaluate(get(queries, j), h)
    end
    answers
end

"""
   evaluate(queries::SeriesRangeQueries, h::Histogram)

Evaluate all range queries on the given histogram.
"""
function evaluate(queries::SeriesRangeQueries, h::Histogram)
    queries.domain == length(h.weights) || throw(DimensionMismatch("query and histogram dimensions differ"))
    answers = zeros(Float64, queries.domain)
    running_sum = 0
    @simd for j=1:queries.domain
        @inbounds running_sum += h.weights[j]
        @inbounds answers[j] = running_sum
    end
    answers
end

"""
    queriesMatrix(queries::SeriesRangeQueries)

Return the `N × N` linear-query matrix equivalent to `queries`, where column
`i` is the indicator vector for prefix `1:i`.
"""
function queriesMatrix(queries::SeriesRangeQueries)

    query_matrix = zeros(queries.domain, queries.domain)
    @simd for i=1:queries.domain
        @simd for j=1:i
            @inbounds query_matrix[j,i] = 1
        end
    end
    query_matrix
end
