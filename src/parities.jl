"""
    Parities

Implementation of parity queries using Fast Hadamard Walsh Transform and Gosper
iteration.

Releasing privacy prserving marginals for contingency table.
The measurments in this case are not the marginals, but the fourier coefficients of the contingency table,
that are sufficient for recovering the marginals.

Details can be found in:
"B. Barak, K. Chaudhuri, C. Dwork, S. Kale,F. McSherry, and K. Talwar. 
Privacy, accuracy, and consistency too: a holistic solution to contingency table release. In PODS, 2007."
"""
struct Parities <: Queries
    dimension::Int
    order::Int
    idx::Array{Int, 1}
end

"""
    hadamard_basis_vector(index, dimension)

Returns hadamard basis vector given by `index`.
"""
function hadamard_basis_vector(index::Int, dimension::Int)
    hadamard = zeros(Float64, 1 << dimension)
    hadamard[1] = 1.0
    for i = 0:dimension-1
        sign = (index & (1 << i)) > 0 ? -1.0 : 1.0
        @simd for j = 1:(1 << i)
            @inbounds hadamard[j + (1 << i)] = sign * hadamard[j]
        end
    end
    hadamard
end

"""
return the downward closure of all the incides of the k-way marginals, of dimension d.

# Example:
    k = 2, d = 4
    k-way indices= [0011, 0101, 0110, 1001, 1010, 1100]
    downward_closure(0011) = [0000, 0001, 0010, 0011]
"""
function complete_way_marginals_indices(k::Int64, d::Int64)
    
    idx = Set(zeros(Int64, 0))
    
    for s = gosper(d, k)
        for alpha in BinaryItr(d, s)
            temp = alpha + 1
            push!(idx, temp)
        end
    end

    sort(collect(idx))
end


function Parities(dimension, order)
    Parities(dimension, order, complete_way_marginals_indices(order, dimension))
end

function get(queries::Parities, i::Int)
    HistogramQuery(hadamard_basis_vector(queries.idx[i]-1, queries.dimension))
end

"""
evaluate the measurements of the linear queries defined by 
queries.idx (the indices needed for calculating the marginals of the contingency table)
"""
function evaluate(queries::Parities, h::Histogram)
    2^queries.dimension * fwht_natural(h.weights)[queries.idx]
end

function fourierCoefficients(queries::Parities, h::Histogram)
    2^queries.dimension * fwht_natural(h.weights)
end

"""
Calculate the marginal C_beta : R^(2^d) -> R^(2^(norm1(beta))), for beta in {0,1}^d, 
of the fourier basis vector in the location alpha
"""
function fourierMarginal(dimension::Int64, beta::Int64, alpha::Int64)

    @assert beta & alpha == alpha

    beta_norm1, beta_indices = norm_1_with_indices(beta)
    alpha_norm1, alpha_indices = norm_1_with_indices(alpha)

    marginal = zeros(Int64, 2^beta_norm1)

    inverse_beta = setdiff(collect(1:dimension), beta_indices)

    for (gamma_index, gamma) in enumerate(BinaryItr(dimension, beta))
        marginal[gamma_index] = 2^length(inverse_beta)*((-1)^(norm_1(alpha&gamma)))
    end

    marginal
end

"""
Calculate the marginal C_beta : R^(2^d) -> R^(2^(norm1(beta))), for beta in {0,1}^d, 
of the contingencyTable (represent as histogram).
"""
function calculateMarginal(queries::Parities, h::Histogram, beta::Int64)

    beta_norm1, beta_indices = norm_1_with_indices(beta)

    marginal = zeros(Int64, 2^beta_norm1)

    coefficients = fourierCoefficients(queries, h)

    for alpha in BinaryItr(queries.dimension, beta)
        marginal += fourierMarginal(queries.dimension, beta, alpha)*coefficients[alpha+1]
    end

    marginal/(2^queries.dimension)
end

function Marginals(queries::Parities, h::Histogram)

    marginals = []

    for beta in queries.idx
        push!(marginals, calculateMarginal(queries, h, beta))
    end

    marginals
end

function initialize(queries::Parities, data::Tabular, ps::MWParameters)
    initialize(queries, Histogram(data), ps)
end

struct FactorParity <: FactorHistogramQuery
    attributes::Array{Int, 1}
end

"""
    FactorParities

Implementation of parity queries over factored histograms.
"""
struct FactorParities <: FactorHistogramQueries
    queries::Array{FactorParity, 1}
end

attributes(q::FactorParity) = q.attributes
default_value(q::FactorParity) = 0.0
get(qs::FactorParities, i::QueryIndex) = qs.queries[i]

function restrict(q::FactorParity, attributes::Array{Int, 1})
    idx = 0
    d = length(attributes)
    for a in q.attributes
        i = findfirst(x -> x==a, attributes)
        idx += 2^(d-i)
    end
    HistogramQuery(hadamard_basis_vector(idx, d))
end

function FactorParities(dimension, order)
    qs = FactorParity[]
    for k = 1:order
        append!(qs, [FactorParity(s) for s = subsets(collect(1:dimension), k)])
    end
    FactorParities(qs)
end

function parity(x::Array{Float64, 1}, attributes::Array{Int64, 1})
    (-1.0)^sum(x[attributes])
end

function evaluate(q::FactorParity, table::Tabular)
    s = 0.0
    n = size(table.data)[2]
    a = attributes(q)
    @simd for i = 1:n
        @inbounds s += parity(table.data[:, i], a)
    end
    s/n
end

function evaluate(qs::FactorParities, table::Tabular)
    evals = zeros(Float64, length(qs.queries))
    @simd for i in 1:length(qs.queries)
        @inbounds evals[i] = evaluate(qs.queries[i], table)
    end
    evals
end

function evaluate(qs::FactorParities, d::FactorHistogram)
    [evaluate(q, d) for q in qs.queries]
end
