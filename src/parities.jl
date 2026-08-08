"""
    Parities

Implementation of parity queries using Fast Hadamard Walsh Transform and Gosper
iteration.

Releasing privacy-preserving marginals for a contingency table.
The measurements in this case are not the marginals, but the Fourier coefficients of the contingency table,
that are sufficient for recovering the marginals.

Details can be found in:
"B. Barak, K. Chaudhuri, C. Dwork, S. Kale,F. McSherry, and K. Talwar. 
Privacy, accuracy, and consistency too: a holistic solution to contingency table release. In PODS, 2007."
"""
struct Parities <: Queries
    dimension::Int
    order::Int
    idx::Vector{Int}
end

"""
    hadamard_basis_vector(index, dimension)

Returns hadamard basis vector given by `index`.
"""
function hadamard_basis_vector(index::Integer, dimension::Integer)
    0 <= dimension <= 62 || throw(ArgumentError("dimension must be between 0 and 62"))
    0 <= index < (1 << dimension) || throw(ArgumentError("index is outside the Hadamard basis"))
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
Return the downward closure of all indices of the `k`-way marginals in dimension `d`.

# Example:
    k = 2, d = 4
    k-way indices= [0011, 0101, 0110, 1001, 1010, 1100]
    downward_closure(0011) = [0000, 0001, 0010, 0011]
"""
function complete_way_marginals_indices(k::Integer, d::Integer)
    1 <= k <= d <= 62 || throw(ArgumentError("order and dimension must satisfy 1 <= order <= dimension <= 62"))

    idx = Set{Int}()
    
    for s = gosper(d, k)
        for alpha in BinaryItr(d, s)
            temp = alpha + 1
            push!(idx, temp)
        end
    end

    sort(collect(idx))
end


function Parities(dimension::Integer, order::Integer)
    1 <= order <= dimension <= 62 || throw(ArgumentError("order and dimension must satisfy 1 <= order <= dimension <= 62"))
    Parities(Int(dimension), Int(order), complete_way_marginals_indices(order, dimension))
end

function get(queries::Parities, i::Int)
    HistogramQuery(hadamard_basis_vector(queries.idx[i]-1, queries.dimension))
end

"""
evaluate the measurements of the linear queries defined by 
queries.idx (the indices needed for calculating the marginals of the contingency table)
"""
function evaluate(queries::Parities, h::Histogram)
    length(h.weights) == 1 << queries.dimension || throw(DimensionMismatch("histogram length must equal 2^dimension"))
    2^queries.dimension * fwht_natural(h.weights)[queries.idx]
end

"""
    fourierCoefficients(queries::Parities, histogram)

Return all Walsh–Hadamard/Fourier coefficients of `histogram` using the
dimension specified by `queries`.
"""
function fourierCoefficients(queries::Parities, h::Histogram)
    length(h.weights) == 1 << queries.dimension || throw(DimensionMismatch("histogram length must equal 2^dimension"))
    2^queries.dimension * fwht_natural(h.weights)
end

"""
Calculate the marginal C_beta : R^(2^d) -> R^(2^(norm1(beta))), for beta in {0,1}^d, 
of the fourier basis vector in the location alpha
"""
function fourierMarginal(dimension::Integer, beta::Integer, alpha::Integer)
    0 <= dimension <= 62 || throw(ArgumentError("dimension must be between 0 and 62"))
    0 <= beta < (1 << dimension) || throw(ArgumentError("beta is outside the data domain"))
    0 <= alpha < (1 << dimension) || throw(ArgumentError("alpha is outside the data domain"))
    beta & alpha == alpha || throw(ArgumentError("alpha must be a subset of beta"))

    beta_norm1, beta_indices = norm_1_with_indices(beta)

    marginal = zeros(Float64, 2^beta_norm1)

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
function calculateMarginal(queries::Parities, h::Histogram, beta::Integer)
    0 <= beta < (1 << queries.dimension) || throw(ArgumentError("beta is outside the data domain"))

    beta_norm1, _ = norm_1_with_indices(beta)

    marginal = zeros(Float64, 2^beta_norm1)

    coefficients = fourierCoefficients(queries, h)

    for alpha in BinaryItr(queries.dimension, beta)
        marginal += fourierMarginal(queries.dimension, beta, alpha)*coefficients[alpha+1]
    end

    marginal/(2^queries.dimension)
end

"""
    Marginals(queries::Parities, histogram)

Reconstruct every exact-order marginal requested by `queries`. The returned
vector contains one probability vector for each size-`queries.order` subset of
attributes, in Gosper-iterator order.
"""
function Marginals(queries::Parities, h::Histogram)
    [calculateMarginal(queries, h, beta) for beta in gosper(queries.dimension, queries.order)]
end

function initialize(queries::Parities, data::Tabular, ps::MWParameters)
    initialize(queries, Histogram(data), ps)
end

struct FactorParity <: FactorHistogramQuery
    attributes::Vector{Int}

    function FactorParity(attributes::AbstractVector{<:Integer})
        isempty(attributes) && throw(ArgumentError("a factor parity must contain at least one attribute"))
        all(>(0), attributes) || throw(ArgumentError("attribute indices must be positive"))
        allunique(attributes) || throw(ArgumentError("attribute indices must be unique"))
        new(Int.(attributes))
    end
end

"""
    FactorParities

Implementation of parity queries over factored histograms.
"""
struct FactorParities <: FactorHistogramQueries
    queries::Vector{FactorParity}
end

attributes(q::FactorParity) = q.attributes
default_value(q::FactorParity) = 0.0
get(qs::FactorParities, i::QueryIndex) = qs.queries[i]

function restrict(q::FactorParity, attributes::Vector{Int})
    idx = 0
    d = length(attributes)
    for a in q.attributes
        i = findfirst(x -> x==a, attributes)
        idx += 2^(d-i)
    end
    HistogramQuery(hadamard_basis_vector(idx, d))
end

function FactorParities(dimension::Integer, order::Integer)
    1 <= order <= dimension || throw(ArgumentError("order and dimension must satisfy 1 <= order <= dimension"))
    qs = FactorParity[]
    for k = 1:order
        append!(qs, [FactorParity(s) for s = subsets(collect(1:dimension), k)])
    end
    FactorParities(qs)
end

function parity(x::AbstractVector{<:Real}, attributes::AbstractVector{<:Integer})
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

function evaluate(q::FactorParity, d::FactorHistogram)
    parent_indices = unique(d.lookup[a] for a in q.attributes)
    mapreduce(*, parent_indices; init=1.0) do parent_index
        factor = d.factors[parent_index]
        local_attributes = [a for a in q.attributes if d.lookup[a] == parent_index]
        evaluate(restrict(FactorParity(local_attributes), factor.attributes), factor.histogram)
    end
end
