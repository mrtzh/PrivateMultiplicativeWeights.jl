abstract type FactorHistogramQuery <: Query end

abstract type FactorHistogramQueries <: Queries end

mutable struct Factor
    attributes::Array{Int, 1}
    histogram::Histogram
end

"""
    FactorHistogram

FactorHistogram represents the data as a product of histograms on disjoint
attributes. This corresponds to a distribution where variables in different
factors are independent.
"""
struct FactorHistogram <: Data
    factors::Array{Factor, 1}
    lookup::Dict{Int, Int}
end

function attributes(q::FactorHistogramQuery)
    throw("`attributes` not implemented for `$(typeof(q))`.")
end

function default_value(q::FactorHistogramQuery)
    throw("`default_value` not implemented for `$(typeof(q))`.")
end

function restrict(q::FactorHistogramQuery, attributes::Array{Int, 1})
    throw("`restrict` not implemented for `$(typeof(q))`.")
end

function get_parent(d::FactorHistogram, attributes::Array{Int, 1})
    factor_idxs = [d.lookup[i] for i in attributes]
    all(i -> i == factor_idxs[1], factor_idxs) ? factor_idxs[1] : 0
end

function merge_factors(f::Factor, g::Factor)
    Factor(append!(f.attributes, g.attributes),
           Histogram(kron(f.histogram.weights, g.histogram.weights)))
end

function merge_components!(q::FactorHistogramQuery, d::FactorHistogram)
    parent_index = d.lookup[attributes(q)[1]]
    for i in 2:length(attributes(q))
        current_index = d.lookup[attributes(q)[i]]
        if current_index != parent_index
            d.factors[parent_index] = merge_factors(d.factors[parent_index],
                                                    d.factors[current_index])
            for a in d.factors[current_index].attributes
                d.lookup[a] = parent_index
            end
            d.factors[current_index].attributes = Int[]
            d.factors[current_index].histogram = Histogram([1.0])
        end
    end
end

function evaluate(q::FactorHistogramQuery, d::FactorHistogram)
    parent_index = get_parent(d, attributes(q))
    if parent_index > 0
        c = d.factors[parent_index]
        qc = restrict(q, c.attributes)
        evaluate(qc, c.histogram)
    else
        default_value(q)
    end
end

function update!(q::FactorHistogramQuery, d::FactorHistogram, error::Float64)
    parent_index = get_parent(d, attributes(q))
    if parent_index > 0
        c = d.factors[parent_index]
        qc = restrict(q, c.attributes)
        update!(qc, c.histogram, error)
        normalize!(c.histogram)
    else
        merge_components!(q, d)
        update!(q, d, error)
    end
end

# normalize happens inside update!
normalize!(h::FactorHistogram) = h

function initialize(queries::FactorHistogramQueries, table::Tabular, ps::MWParameters)
    d, n = size(table.data)
    components = [Factor([i], Histogram([0.5, 0.5])) for i = 1:d]
    lookup = Dict{Int, Int}()
    for i = 1:d
        lookup[i] = i
    end
    synthetic = FactorHistogram(components, lookup)
    real_answers = evaluate(queries, table)
    scale = 2.0/(ps.epsilon*n)
    MWState(table, synthetic, queries, real_answers, Dict{Int, Float64}(), scale)
end
