"""
    Histogram

Histogram represents the data as a vector where each coordinate corresponds to
one element of the data universe. This is the default representation for MWEM.
"""
type Histogram <: Data
    weights::Array{Float64, 1}
    num_samples::Int64
end

function Histogram(weights::Array{Float64, 1})
    Histogram(weights, 0)
end

type HistogramQuery <: Query
    weights::Array{Float64, 1}
end

type HistogramQueries <: Queries
    queries::Array{Float64, 2}
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

# if error is > 0, increase weight on positive elements of queries[query] and 
# decrease weight on negative elements. Magnitude of error determines step size.
function update!(q::HistogramQuery, h::Histogram, error::Float64)
    @simd for j = 1:length(h.weights)
        @inbounds h.weights[j] *= exp(error * q.weights[j] / 2.0)
    end
    normalize!(h)
end

function initialize(queries::Queries, data::Histogram, parameters)
    data = normalize!(data)
    histogram_length = length(data.weights)
    num_samples = data.num_samples
    epsilon, iterations, repetitions, noisy_init = parameters
    if noisy_init
        # Noisy init incurs an additional `epsilon` privacy cost
        weights = Array(Float64, histogram_length)
        noise = rand(Laplace(0.0, 1.0/(epsilon*num_samples)), histogram_length)
        @simd for i = 1:histogram_length
             @inbounds weights[i] = 
                 max(data.weights[i] + noise[i] - 1.0/(e*n*epsilon), 0.0)
        end
        weights /= sum(weights)
        synthetic = Histogram(0.5 * weights + 0.5/histogram_length)
    else
        synthetic = Histogram(ones(histogram_length)/histogram_length)
    end
    real_answers = evaluate(queries, data)
    scale = 2*iterations/(epsilon*num_samples)
    MWState(data, synthetic, queries, real_answers, Dict{Int, Float64}(),
                                                    scale, repetitions)
end

function initialize(queries::Queries, data::Tabular, parameters)
    initialize(queries, Histogram(data), parameters)
end


"""
    Histogram(table)

Create histogram representation from tabular data.
"""
function Histogram(table::Tabular)
    d, n = size(table.data)
    histogram = zeros(2^d)
    for i = 1:n
        num = 0
        x = vec(table.data[:, i])
        for i = 1:d
            num += convert(Int64, x[d-i+1]) * 2^(i-1)
        end
        histogram[num+1] += 1.0
    end
    normalize!(Histogram(histogram, n))
end

"""
    Tabular(histogram, n)

Create tabular data by sampling n times weighted by histogram.
"""
function Tabular(histogram::Histogram, n::Int)
    N = length(histogram.weights)
    d = convert(Int64, log(2, N))
    idx = wsample(collect(0:N-1), histogram.weights, n)
    data_matrix = zeros(d, n)
    for i = 1:n
        data_matrix[:, i] = reverse(digits(idx[i], 2, d))
    end
    Tabular(data_matrix)
end
