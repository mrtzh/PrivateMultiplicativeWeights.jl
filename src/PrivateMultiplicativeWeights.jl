module PrivateMultiplicativeWeights

using
    Distributions: Laplace, wsample
using
    Printf,
    Hadamard,
    LinearAlgebra,
    Random,
    IterTools,
    Statistics,
    Distributed

export
    mwem,
    MWParameters,
    Tabular,
    Histogram,
    HistogramQueries,
    SeriesRangeQueries,
    RangeQueries,
    RangeQuery,
    Parities,
    FactorParities,
    Interval,
    gosper,
    BinaryItr,
    
    maximum_error,
    kl_divergence_error,
    mean_squared_error,
    queriesMatrix,

    evaluate,
    get,

    kl_divergence,
    complete_way_marginals_indices,
    fourierCoefficients,
    calculateMarginal,
    Marginals,

    get_query_vector_from_interval,
    hadamard_basis_vector

import
    Base: eltype, length, iterate

include("interface.jl")
include("histogram.jl")
include("rangequeries.jl")
include("factors.jl")
include("gosper.jl")
include("parities.jl")
include("error.jl")
include("mwem.jl")

"""
    PrivateMultiplicativeWeights

A simple and practical algorithm for differentially private data analysis.
"""
PrivateMultiplicativeWeights

end # module
