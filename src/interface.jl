abstract type Data end
abstract type Query end
abstract type Queries end


const QueryIndex = Int
const Interval = Tuple{Int, Int}


struct MWState
    real::Data
    synthetic::Data
    queries::Queries
    real_answers::Vector{Float64}
    measurements::Dict{Int, Float64}
    scale::Float64
end


struct MWParameters
    epsilon::Float64
    iterations::Int
    repetitions::Int
    noisy_init::Bool
    verbose::Bool
    init_budget::Float64
    noisy_max_budget::Float64
end


"""
    Tabular(data)

Store a finite numeric matrix with attributes in rows and records in columns.
Built-in histogram conversion and parity workloads require every value to be
binary (`0` or `1`). The input is copied and converted to `Matrix{Float64}`.
"""
struct Tabular <: Data
    data::Matrix{Float64}

    function Tabular(data::AbstractMatrix{<:Real})
        all(isfinite, data) || throw(ArgumentError("tabular data must contain only finite values"))
        new(Matrix{Float64}(data))
    end
end


"""
    get(qs::Queries, qindex::QueryIndex)

Returns Query corresponding to given QueryIndex.
"""
function get(qs::Queries, qindex::QueryIndex)
    throw("`get` not implemented for `$(typeof(qs))`.")
end


"""
    evaluate(q::Query, d::Data)

Returns float value of Query evaluated on Data.
"""
function evaluate(q::Query, d::Data)
    throw("`evaluate` not implemented for `$(typeof(q))`, `$(typeof(d))`.")
end


"""
    evaluate(qs::Queries, d::Data)

Returns vector of floats corresponding to Queries evaluated on Data.
"""
function evaluate(qs::Queries, d::Data)
    throw("`evaluate` not implemented for `$(typeof(qs))`, `$(typeof(d))`.")
end


"""
    initialize(qs::Queries, d::Data, ps::MWParameters)

Returns MWState initialization for given queries, data, and parameters.
"""
function initialize(qs::Queries, d::Data, ps::MWParameters)
    throw("`initialize` not implemented for `$(typeof(qs))`, `$(typeof(d))`.")
end


"""
    update!(q::Query, d::Data, error::Float64)

Performs multiplicative weights update of synthetic data `d` for given query `q`
and error value `error`.
"""
function update!(q::Query, d::Data, error::Float64)
    throw("`update!` not implemented for `$(typeof(q))`, `$(typeof(d))`.")
end


"""
   normalize!(d::Data)

Normalizes data.
"""
function normalize!(d::Data)
    throw("`normalize!` not implemented for `$(typeof(d))`.")
end
