"""
    noisy_max(mwstate::MWState)

Select index of query with largest error after noise addition.
"""
function noisy_max(mwstate::MWState)
    diffs = mwstate.real_answers - evaluate(mwstate.queries, mwstate.synthetic)
    # do not select previously measured queries
    diffs[collect(keys(mwstate.measurements))] = 0.0
    indmax(abs(diffs) + rand(Laplace(0.0, mwstate.scale), length(diffs)))
end


"""
    update!(mwstate::MWState)

Perform multiplicative weights update.
"""
function update!(mwstate::MWState)
    # repeatedly update queries that were already measured before
    for i = 1:mwstate.repetitions
        for qindex in keys(mwstate.measurements)
            query = get(mwstate.queries, qindex)
            error = mwstate.measurements[qindex] 
                        - evaluate(query, mwstate.synthetic)
            update!(query, mwstate.synthetic, error)
        end
    end
    normalize!(mwstate.synthetic)
    mwstate
end

"""
    mwem

Private Multiplicative Weights (MWEM) repeatedly selects largest error query
and performs multiplicative weights update.

# Arguments
* `queries::Queries`: Set of queries
* `data::Data`: Input data
* `epsilon::Float=1.0`: Each iterations is `epsilon`-differentially private
* `repetitions::Integer=10`: Repeatedly update with previously selected queries.
* `verbose::Boolean=false`: Print output
* `noisy_init::Boolean=false`: Initialize with noisy histogram if true.
"""
function mwem(queries::Queries, 
              data::Data;
              epsilon=1.0, 
              iterations=10, 
              repetitions=10, 
              verbose=false, 
              noisy_init=false)

    # Initialization
    parameters = (epsilon, iterations, repetitions, noisy_init)
    time = @elapsed mwstate = initialize(queries, data, parameters)
    verbose ? print("Iter.\t Max error\t\t time (sec)\n") : nothing
    verbose ? print("0\t", maximum_error(mwstate), "\t", time, "\n") : nothing

    # Iterations
    for t = 1:iterations
        time = @elapsed begin
            i = noisy_max(mwstate)
            mwstate.measurements[i] = mwstate.real_answers[i]
                                          + rand(Laplace(0.0, mwstate.scale))
            update!(mwstate)
        end
        verbose ? print(t, "\t", maximum_error(mwstate), "\t", time, "\n") : nothing
    end

    mwstate
end
