"""
    MWParameters(; epsilon=1.0, iterations=10, repetitions=10,
                   noisy_init=false, verbose=false, init_budget=0.05,
                   noisy_max_budget=0.5)

Construct validated settings for [`mwem`](@ref).

# Keywords
- `epsilon`: positive privacy parameter allocated to initialization and each
  private iteration.
- `iterations`: number of private query-selection and measurement rounds.
- `repetitions`: additional post-processing passes over released measurements.
- `noisy_init`: use private data-dependent initialization for explicit
  histograms. Factored histograms do not support this option.
- `verbose`: print exact error diagnostics; this output is not private.
- `init_budget`: fraction of `epsilon` assigned to noisy initialization.
- `noisy_max_budget`: fraction of each iteration assigned to query selection.
"""
function MWParameters(; epsilon::Real=1.0,
                        iterations::Integer=10,
                        repetitions::Integer=10,
                        noisy_init::Bool=false,
                        verbose::Bool=false,
                        init_budget::Real=0.05,
                        noisy_max_budget::Real=0.5)
    isfinite(epsilon) && epsilon > 0 || throw(ArgumentError("epsilon must be finite and positive"))
    iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
    repetitions >= 0 || throw(ArgumentError("repetitions must be nonnegative"))
    isfinite(init_budget) && 0 <= init_budget < 1 || throw(ArgumentError("init_budget must be in [0, 1)"))
    noisy_init && init_budget == 0 && throw(ArgumentError("init_budget must be positive when noisy_init is enabled"))
    isfinite(noisy_max_budget) && 0 < noisy_max_budget < 1 || throw(ArgumentError("noisy_max_budget must be in (0, 1)"))

    MWParameters(Float64(epsilon), Int(iterations), Int(repetitions), noisy_init,
                 verbose, Float64(init_budget), Float64(noisy_max_budget))
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
    mwem(queries::Queries, data::Data, parameters::MWParameters=MWParameters())

Private Multiplicative Weights (MWEM) repeatedly selects largest error query
and performs multiplicative-weights updates. Returns an `MWState`; its
`synthetic` field contains the released approximation and `measurements`
contains the selected noisy answers.

`parameters.epsilon` applies to each iteration. This function does not account
for privacy composition across iterations or invocations. When `verbose=true`,
the printed error values use exact source answers and are not differentially
private.
"""
function mwem(queries::Queries, data::Data, ps::MWParameters=MWParameters())

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
