"""
    MWParameters

Returns MWParameters with default settings where not specified.

# Arguments
* `epsilon::Float=1.0`: Each iterations is `epsilon`-differentially private
* `iterations::Integer=10`: Number of iterations
* `repetitions::Integer=10`: Repeatedly update with previously selected queries.
* `noisy_init::Boolean=false`: Initialize with noisy histogram if true.
* `verbose::Boolean=false`: Print timing and error information (not private)
"""
function MWParameters(; epsilon=1.0,
                        iterations=10,
                        repetitions=10,
                        noisy_init=false,
                        verbose=false)
    MWParameters(epsilon, iterations, repetitions, noisy_init, verbose)
end


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
    mwem(queries::Queries, data::Data, ps::MWParameters)

Private Multiplicative Weights (MWEM) repeatedly selects largest error query
and performs multiplicative weights update.

"""
function mwem(queries::Queries, data::Data, ps=MWParameters())

    # Initialization
    time = @elapsed mwstate = initialize(queries, data, ps)

    if ps.verbose
        error = mean_squared_error(mwstate)
        print("Iter.\t Mean sq err\t time (sec)\n")
        @printf("0\t %.3f\t\t %.3f\n", error, time)
    end

    # Iterations
    for t = 1:ps.iterations
        time = @elapsed begin
            # select query via noisy max
            qindex = noisy_max(mwstate)
            mwstate.measurements[qindex] = 
              mwstate.real_answers[qindex] + rand(Laplace(0.0, mwstate.scale))

            # update synthetic data approximation
            update!(mwstate, qindex)

            # repeatedly update on previously measured queries in random order
            for i = 1:ps.repetitions
                for qindex in shuffle(collect(keys(mwstate.measurements)))
                    update!(mwstate, qindex)
                end
            end
        end

        if ps.verbose
            error = mean_squared_error(mwstate)
            @printf("%d\t %.3f\t\t %.3f\n", t, error, time)
        end
    end

    mwstate
end
