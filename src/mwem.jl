"""
    MWParameters

Returns MWParameters with default settings where not specified.

# Arguments
* `epsilon::Float=1.0`: Each iterations is `epsilon`-differentially private
* `iterations::Integer=10`: Number of iterations
* `repetitions::Integer=10`: Repeatedly update with previously selected queries.
* `noisy_init::Boolean=false`: Initialize with noisy histogram if true.
* `verbose::Boolean=false`: Print timing and error information (not private)
* `init_budget=0.05: fraction of the epsilon for the noisy initialization
* `noisy_max_budget: fraction of the badget from every step for for the noisy max
"""
function MWParameters(; epsilon=1.0,
                        iterations=10,
                        repetitions=10,
                        noisy_init=false,
                        verbose=false,
                        init_budget=0.05,
                        noisy_max_budget=0.5
                        )
    MWParameters(epsilon, iterations, repetitions, noisy_init, verbose, init_budget, noisy_max_budget)
end


"""
    noisy_max(mwstate::MWState)

Select index of query with largest error after noise addition.
"""
function noisy_max(mwstate::MWState, noisy_max_budget::Float64)
    diffs = mwstate.real_answers - evaluate(mwstate.queries, mwstate.synthetic)
    diffs[collect(keys(mwstate.measurements))] .= 0.0
    real_index = argmax(abs.(diffs) + rand(Laplace(0.0, mwstate.scale/noisy_max_budget), length(diffs)))
    real_index
end


"""
    update!(mwstate::MWState, qindex::QueryIndex)

Perform multiplicative weights update with query given by `qindex`.
"""
function update!(mwstate::MWState, qindex::QueryIndex)
    query = get(mwstate.queries, qindex)
    error = mwstate.measurements[qindex] - evaluate(query, mwstate.synthetic)
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
            qindex = noisy_max(mwstate, ps.noisy_max_budget)
            mwstate.measurements[qindex] =
              mwstate.real_answers[qindex] + rand(Laplace(0.0, mwstate.scale/(1-ps.noisy_max_budget)))

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
            error_mean = mean_squared_error(mwstate)
            error_max = maximum_error(mwstate)
            @printf("itr: %d\t mean_error: %.3f\t\t %.3f\n", t, error_mean, time)
            @printf("itr: %d\t max_error: %.3f\t\t %.3f\n", t, error_max, time)
        end
    end

    mwstate
end
