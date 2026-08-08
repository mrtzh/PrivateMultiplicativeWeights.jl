"""
    Histogram(weights, num_samples=0)

Histogram represents the data as a vector where each coordinate corresponds to
one element of the data universe. Weights must be finite and nonnegative, but
need not be normalized. Pass the positive source record count as `num_samples`
when using a directly constructed histogram as private input to [`mwem`](@ref).
The constructor copies and converts `weights` to `Vector{Float64}`.
"""
mutable struct Histogram <: Data
    weights::Vector{Float64}
    num_samples::Int

    function Histogram(weights::AbstractVector{<:Real}, num_samples::Integer=0)
        isempty(weights) && throw(ArgumentError("a histogram must contain at least one weight"))
        all(isfinite, weights) || throw(ArgumentError("histogram weights must be finite"))
        all(>=(0), weights) || throw(ArgumentError("histogram weights must be nonnegative"))
        num_samples >= 0 || throw(ArgumentError("num_samples must be nonnegative"))
        new(Vector{Float64}(weights), Int(num_samples))
    end
end

struct HistogramQuery <: Query
    weights::Vector{Float64}

    function HistogramQuery(weights::AbstractVector{<:Real})
        isempty(weights) && throw(ArgumentError("a histogram query must contain at least one weight"))
        all(isfinite, weights) || throw(ArgumentError("query weights must be finite"))
        new(Vector{Float64}(weights))
    end
end

"""
    HistogramQueries(query_matrix)

Construct a linear workload over histograms. For a histogram with `N` bins and
`k` queries, `query_matrix` must have shape `N × k`; each column contains one
query's weights.
"""
struct HistogramQueries <: Queries
    queries::Matrix{Float64}

    function HistogramQueries(queries::AbstractMatrix{<:Real})
        size(queries, 1) > 0 || throw(ArgumentError("a query matrix must have at least one row"))
        size(queries, 2) > 0 || throw(ArgumentError("a query matrix must have at least one column"))
        all(isfinite, queries) || throw(ArgumentError("query weights must be finite"))
        new(Matrix{Float64}(queries))
    end
end

function get(queries::HistogramQueries, i::QueryIndex)
    HistogramQuery(queries.queries[:, i])
end

function evaluate(query::HistogramQuery, h::Histogram)
    length(query.weights) == length(h.weights) || throw(DimensionMismatch("query and histogram dimensions differ"))
    dot(query.weights, h.weights)
end

function evaluate(queries::HistogramQueries, h::Histogram)
    size(queries.queries, 1) == length(h.weights) || throw(DimensionMismatch("query and histogram dimensions differ"))
    queries.queries' * h.weights
end

function normalize!(h::Histogram)
    total = sum(h.weights)
    isfinite(total) && total > 0 || throw(ArgumentError("histogram weights must have a finite, positive sum"))
    h.weights ./= total
    h
end

# if error is > 0, increase weight on positive elements of queries[query] and
# decrease weight on negative elements. Magnitude of error determines step size.
function update!(q::HistogramQuery, h::Histogram, error::Float64)
    length(q.weights) == length(h.weights) || throw(DimensionMismatch("query and histogram dimensions differ"))
    isfinite(error) || throw(ArgumentError("update error must be finite"))
    logweights = log.(h.weights) .+ error .* q.weights ./ 2.0
    logweights .-= maximum(logweights)
    h.weights .= exp.(logweights)
    normalize!(h)
end

function initialize(queries::Queries, data::Histogram, ps::MWParameters)
    data.num_samples > 0 || throw(ArgumentError("num_samples must be positive when running MWEM on a histogram"))
    data = normalize!(Histogram(data.weights, data.num_samples))
    histogram_length = length(data.weights)
    num_samples = data.num_samples
    effective_epsilon = ps.epsilon
    if ps.noisy_init
        # Noisy init incurs an additional `epsilon` privacy cost
        weights = zeros(Float64, histogram_length)
        noise = rand(Laplace(0.0, 1/(num_samples*ps.epsilon*ps.init_budget)), histogram_length)
        @simd for i = 1:histogram_length
             @inbounds weights[i] =
                 max(data.weights[i] + noise[i] - 1.0/(exp(1)*num_samples*ps.epsilon), 0.0)
        end
        total = sum(weights)
        if total > 0
            weights ./= total
        else
            fill!(weights, 1/histogram_length)
        end
        synthetic = Histogram(0.5 .* weights .+ 0.5/histogram_length)
        effective_epsilon *= 1-ps.init_budget
    else
        weights = ones(Float64, histogram_length)
        synthetic = Histogram(weights/histogram_length)
    end
    real_answers = evaluate(queries, data)
    scale = 2.0/(effective_epsilon*num_samples)
    MWState(data, synthetic, queries, real_answers, Dict{Int, Float64}(), scale)
end

function initialize(queries::HistogramQueries, data::Tabular, ps::MWParameters)
    initialize(queries, Histogram(data), ps)
end


"""
    Histogram(table)

Create histogram representation from tabular data.
"""
function Histogram(table::Tabular)
    d, n = size(table.data)
    n > 0 || throw(ArgumentError("cannot build a histogram from an empty table"))
    d < 8*sizeof(Int) || throw(ArgumentError("the table has too many binary attributes for this platform"))
    all(x -> x == 0 || x == 1, table.data) || throw(ArgumentError("histogram conversion requires binary tabular data"))
    histogram = zeros(Float64,2^d)
    for i = 1:n
        num = 0
        x = vec(table.data[:, i])
        for j = 1:d
            num += Int(x[d-j+1]) * 2^(j-1)
        end
        histogram[num+1] += 1.0
    end
    normalize!(Histogram(histogram, n))
end

"""
    Tabular(histogram, n)

Create tabular data by sampling n times weighted by histogram.
"""
function Tabular(histogram::Histogram, n::Integer)
    n >= 0 || throw(ArgumentError("the number of samples must be nonnegative"))
    n = Int(n)
    N = length(histogram.weights)
    ispow2(N) || throw(ArgumentError("histogram length must be a power of two"))
    total = sum(histogram.weights)
    isfinite(total) && total > 0 || throw(ArgumentError("histogram weights must have a finite, positive sum"))
    d = trailing_zeros(N)
    idx = wsample(0:N-1, histogram.weights ./ total, n)
    data_matrix = zeros(d, n)
    for i = 1:n
        data_matrix[:, i] = reverse(digits(idx[i], base = 2, pad = d))
    end
    Tabular(data_matrix)
end
