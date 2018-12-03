abstract type Data end
abstract type Query end
abstract type Queries end


const QueryIndex = Int
const Interval = Tuple{Int64, Int64}


struct MWState
    real::Data
    synthetic::Data
    queries::Queries
    real_answers::Array{Float64, 1}
    measurements::Dict{Int, Float64}
    scale::Float64
end


struct MWParameters
    epsilon::Float64
    iterations::Int64
    repetitions::Int64
    noisy_init::Bool
    verbose::Bool
    init_budget::Float64
    noisy_max_budget::Float64
end


struct Tabular <: Data
    data::Array{Float64, 2}
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
