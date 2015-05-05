module PrivateMultiplicativeWeights

using
    Hadamard,
    Distributions.Laplace,
    Distributions.wsample

export
    mwem,
    Tabular,
    Histogram,
    HistogramQueries,
    Parities,
    FactorParities,
    maximum_error,
    mean_squared_error

import
    Base: start, next, done, eltype, length

include("interface.jl")
include("histogram.jl")
include("factors.jl")
include("iterators.jl")
include("parities.jl")
include("error.jl")
include("algorithm.jl")


end # module
