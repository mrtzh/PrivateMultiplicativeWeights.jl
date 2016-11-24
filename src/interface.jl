typealias Float Float64
typealias Vector Array{Float, 1}
typealias Matrix Array{Float, 2}

abstract Data
abstract Query
abstract Queries

typealias QueryIndex Int

type MWState
    real::Data
    synthetic::Data
    queries::Queries
    real_answers::Vector
    measurements::Dict{Int, Float}
    scale::Float
    repetitions::Int
end

type Tabular <: Data
    data::Matrix
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
