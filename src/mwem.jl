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
    update!(mwstate::MWState, qindex::QueryIndex)

Perform multiplicative weights update with query given by `qindex`.
"""
function update!(mwstate::MWState, qindex::QueryIndex)
    query = get(mwstate.queries, qindex)
    error = (mwstate.measurements[qindex] - evaluate(query, mwstate.synthetic))
    update!(query, mwstate.synthetic, error)
end

"""
    mwem(queries, data[, epsilon, repetitions, verbose, noisy_init])

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

    if verbose
        error = maximum_error(mwstate)
        print("Iter.\t Max error\t time (sec)\n")
        @printf("0\t %.3f\t\t %.3f\n", error, time)
    end

    # Iterations
    for t = 1:iterations
        time = @elapsed begin
            qindex = noisy_max(mwstate)
            mwstate.measurements[qindex] = 
              mwstate.real_answers[qindex] + rand(Laplace(0.0, mwstate.scale))
            update!(mwstate, qindex)

            # repeatedly update on previously measured queries in random order
            for i = 1:mwstate.repetitions
                for qindex in shuffle(collect(keys(mwstate.measurements)))
                    update!(mwstate, qindex)
                end
            end
           #update!(mwstate)
        end

        if verbose
            error = maximum_error(mwstate)
            @printf("%d\t %.3f\t\t %.3f\n", t, error, time)
        end
    end

    mwstate
end
