abstract Data
abstract Query
abstract Queries

typealias QueryIndex Int

type MWState
    real::Data
    synthetic::Data
    queries::Queries
    real_answers::Array{Float64, 1}
    measurements::Dict{Int, Float64}
    scale::Float64
    repetitions::Int
end

type Tabular <: Data
    data::Array{Float64, 2}
end

function get(qs::Queries, ::QueryIndex)
    throw("`get` not implemented for `$(typeof(qs))`.")
end

function evaluate(q::Query, d::Data)
    throw("`evaluate` not implemented for `$(typeof(q))`, `$(typeof(d))`.")
end

function evaluate(qs::Queries, d::Data)
    throw("`evaluate` not implemented for `$(typeof(qs))`, `$(typeof(d))`.")
end

function initialize(qs::Queries, d::Data, parameters)
    throw("`initialize` not implemented for `$(typeof(qs))`, `$(typeof(d))`.")
end

function update!(q::Query, d::Data, error::Float64)
    throw("`update!` not implemented for `$(typeof(q))`, `$(typeof(d))`.")
end

function normalize!(d::Data)
    throw("`normalize!` not implemented for `$(typeof(d))`.")
end
