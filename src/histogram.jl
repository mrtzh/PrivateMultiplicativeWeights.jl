"""
    Histogram

Histogram represents the data as a vector where each coordinate corresponds to
one element of the data universe. This is the default representation for MWEM.
"""
mutable struct Histogram <: Data
    weights::Array{Float64, 1}
    num_samples::Int64
end

function Histogram(weights::Array{Float64, 1})
    Histogram(weights, 0)
end

struct HistogramQuery <: Query
    weights::Array{Float64, 1}
end

struct HistogramQueries <: Queries
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

    Threads.@threads for j = 1:length(h.weights)# @threads can help to speed up
        @inbounds h.weights[j] *= exp(error * q.weights[j] / 2.0)
    end
    normalize!(h)
end

function initialize(queries::Queries, data::Histogram, ps::MWParameters)
    data = normalize!(data)
    histogram_length = length(data.weights)
    num_samples = data.num_samples
    if ps.noisy_init
        # Noisy init incurs an additional `epsilon` privacy cost
        weights = zeros(Float64, histogram_length)
        noise = rand(Laplace(0.0, 1/(num_samples*ps.epsilon*ps.init_budget)), histogram_length)
        @simd for i = 1:histogram_length
             @inbounds weights[i] =
                 max(data.weights[i] + noise[i] - 1.0/(e*num_samples*ps.epsilon), 0.0)
        end
        weights /= sum(weights)
        synthetic = Histogram(0.5 * weights + 0.5/histogram_length)
        ps.epsilon = (1-ps.init_budget)*ps.epsilon
    else
        weights = ones(Float64, histogram_length)
        synthetic = Histogram(weights/histogram_length)
    end
    real_answers = evaluate(queries, data)
    scale = 2.0/(ps.epsilon*num_samples)
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
    histogram = zeros(Float64,2^d)
    for i = 1:n
        num = 0
        x = vec(table.data[:, i])
        for j = 1:d
            num += convert(Int64, x[d-j+1]) * 2^(j-1)
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
        #println(idx[i], 2, d)
        data_matrix[:, i] = reverse(digits(idx[i], base = 2, pad = d))
    end
    Tabular(data_matrix)
end
