# data and queries are represented as vectors in the histogram space

type Histogram <: Data
    weights::Vector
end

type HistogramQuery <: Query
    weights::Vector
end

type HistogramQueries <: Queries
    queries::Matrix
end

function get(queries::HistogramQueries, i::QueryIndex)
    HistogramQuery(queries.queries[:, i])
end

function evaluate(query::HistogramQuery, h::Histogram)
    dot(query.weights, h.weights)
end

function evaluate(queries::HistogramQueries, h::Histogram)
    queries.queries' * vec
end

function normalize!(h::Histogram)
    h.weights /= sum(h.weights)
    h
end

# if error is > 0, we want to increase weight on positive elements of queries[query]
# and decrease weight on negative elements. Magnitude of error determines step size.
function update!(q::HistogramQuery, h::Histogram, error::Float64)
    @simd for j = 1:length(h.weights)
        @inbounds h.weights[j] *= exp(error * q.weights[j] / 2.0)
    end
end

function histogram_initialize(queries::Queries, table::Tabular, parameters)
    d, n  = size(table.data)
    epsilon, iterations, repetitions, smart = parameters
    N = 2^d
    real = Histogram(table)
    if smart
        # spend half of epsilon on histogram initialization
        weights = Array(Float64, N)
        noise = rand(Laplace(0.0, 1.0/(n*epsilon)), N)
        @simd for i = 1:N
             @inbounds weights[i] = max(real.weights[i]+noise[i]-1.0/(e*n*epsilon), 0.0)
        end
        weights /= sum(weights)
        synthetic = Histogram(0.5 * weights + 0.5/N)
        epsilon = 0.5*epsilon
    else
        synthetic = Histogram(ones(N)/N)
    end
    real_answers = evaluate(queries, real)
    scale = 2*iterations/(epsilon*n)
    mwstate = MWState(real, synthetic, queries, real_answers, Dict{Int, Float}(), scale, repetitions)
    mwstate
end

function initialize(queries::HistogramQueries, table::Tabular, parameters)
    histogram_initialize(queries, table, parameters)
end

# convert 0/1 data matrix to its histogram representation
function Histogram(table::Tabular)
    d, n = size(table.data)
    histogram = zeros(2^d)
    for i = 1:n
        num = 0
        x = vec(table.data[:, i])
        for i = 1:d
            num += x[d-i+1] * 2^(i-1)
        end
        histogram[num+1] += 1.0
    end
    normalize!(Histogram(histogram))
end

# convert histogram to 0/1 data matrix
function Tabular(histogram::Histogram, n::Int)
    N = length(histogram.weights)
    d = int(log(2, length(histogram.weights)))
    idx = wsample([0:N-1], histogram.weights, n)
    data_matrix = zeros(d, n)
    for i = 1:n
        data_matrix[:, i] = reverse(digits(idx[i], 2, d))
    end
    Tabular(data_matrix)
end
