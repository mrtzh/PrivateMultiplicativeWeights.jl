using PrivateMultiplicativeWeights
using Hadamard
using LinearAlgebra
using Random
using Test

const PMW = PrivateMultiplicativeWeights

@testset "Histogram and tabular representations" begin
    @test Histogram(Tabular(zeros(Int, 1, 10))).weights == [1.0, 0.0]
    @test Histogram(Tabular(ones(Int, 1, 10))).weights == [0.0, 1.0]
    @test Histogram(Tabular([1 0; 0 1])).weights == [0.0, 0.5, 0.5, 0.0]
    @test Histogram(Tabular([1 0 1; 0 1 0])).weights == [0.0, 1.0, 2.0, 0.0] ./ 3
    @test Histogram(Tabular([1 0 1 1; 0 1 0 1])).weights == [0.0, 1.0, 2.0, 1.0] ./ 4
    @test Histogram(Tabular([1 0 1 0 1; 0 1 0 0 1])).weights == [1.0, 1.0, 2.0, 1.0] ./ 5

    Random.seed!(1234)
    @test Tabular(Histogram([1.0, 0.0, 0.0, 0.0]), 10).data == zeros(2, 10)
    @test Tabular(Histogram([0.0, 0.0, 0.0, 1.0]), 10).data == ones(2, 10)
    @test size(Tabular(Histogram([0.5, 0.5]), 0).data) == (1, 0)

    @test_throws ArgumentError Histogram(Float64[])
    @test_throws ArgumentError Histogram([-1.0, 2.0])
    @test_throws ArgumentError PMW.normalize!(Histogram([0.0, 0.0]))
    @test_throws ArgumentError Histogram(Tabular(zeros(2, 0)))
    @test_throws ArgumentError Histogram(Tabular([0 2]))
    @test_throws ArgumentError Tabular(Histogram([0.2, 0.3, 0.5]), 2)
    @test_throws ArgumentError Tabular(Histogram([0.5, 0.5]), -1)
end

@testset "Histogram and range queries" begin
    h = Histogram([0.1, 0.2, 0.3, 0.4])
    matrix_queries = HistogramQueries([
        1.0 0.0 1.0
        0.0 1.0 1.0
        0.0 0.0 0.0
        0.0 0.0 0.0
    ])
    @test evaluate(matrix_queries, h) ≈ [0.1, 0.2, 0.3]

    prefixes = SeriesRangeQueries(4)
    @test evaluate(prefixes, h) ≈ [0.1, 0.3, 0.6, 1.0]
    @test evaluate(HistogramQueries(queriesMatrix(prefixes)), h) ≈ evaluate(prefixes, h)
    @test evaluate(RangeQuery((2, 3), 4), h) ≈ 0.5
    @test evaluate(RangeQueries([(1, 1), (2, 4)], 4), h) ≈ [0.1, 0.9]
    @test get_query_vector_from_interval((2, 3), 4) == [0.0, 1.0, 1.0, 0.0]

    @test_throws DimensionMismatch evaluate(PMW.HistogramQuery([1.0]), h)
    @test_throws DimensionMismatch evaluate(HistogramQueries(ones(3, 2)), h)
    @test_throws DimensionMismatch evaluate(SeriesRangeQueries(3), h)
    @test_throws ArgumentError SeriesRangeQueries(0)
    @test_throws ArgumentError RangeQuery((0, 2), 4)
    @test_throws ArgumentError RangeQueries([(3, 2)], 4)
end

@testset "Combinatorial iterators" begin
    @test collect(gosper(3, 2)) == [3, 5, 6]
    @test length(gosper(5, 2)) == 10
    @test collect(BinaryItr(4, [1, 3])) == [0, 1, 4, 5]
    @test collect(BinaryItr(4, 5, 2)) == [2, 3, 6, 7]
    @test_throws ArgumentError gosper(3, 0)
    @test_throws ArgumentError BinaryItr(3, [1, 1])
    @test_throws ArgumentError BinaryItr(3, 8)
    @test_throws ArgumentError BinaryItr(3, [1], 1)
end

@testset "Hadamard parities and marginals" begin
    for dimension in 0:8, index in 0:(1 << dimension)-1
        @test PMW.hadamard_basis_vector(index, dimension) == Hadamard.hadamard(1 << dimension)[:, index + 1]
    end

    @test complete_way_marginals_indices(2, 3) == collect(1:7)
    @test_throws ArgumentError Parities(2, 3)

    correlated = Histogram([0.5, 0.0, 0.0, 0.5])
    one_way = Marginals(Parities(2, 1), correlated)
    @test length(one_way) == 2
    @test all(marginal -> marginal ≈ [0.5, 0.5], one_way)
    @test Marginals(Parities(2, 2), correlated) == [[0.5, 0.0, 0.0, 0.5]]
    @test_throws DimensionMismatch evaluate(Parities(3, 2), correlated)
end

@testset "Factored histograms" begin
    @test_throws ArgumentError PMW.FactorParity(Int[])
    @test_throws ArgumentError PMW.FactorParity([1, 1])
    factors = [
        PMW.Factor([1], Histogram([0.8, 0.2])),
        PMW.Factor([2], Histogram([0.25, 0.75])),
    ]
    factored = PMW.FactorHistogram(factors, Dict(1 => 1, 2 => 2))
    query = PMW.FactorParity([1, 2])
    @test evaluate(query, factored) ≈ (0.8 - 0.2) * (0.25 - 0.75)

    PMW.update!(query, factored, 0.5)
    @test PMW.get_parent(factored, [1, 2]) > 0
    @test sum(factored.factors[1].histogram.weights) ≈ 1.0
    @test_throws ArgumentError mwem(
        FactorParities(2, 1), Tabular([0 1; 1 0]),
        MWParameters(iterations=0, noisy_init=true),
    )
end

@testset "MWEM parameters and updates" begin
    @test MWParameters().epsilon == 1.0
    @test MWParameters(epsilon=2, iterations=0).iterations == 0
    @test_throws ArgumentError MWParameters(epsilon=0)
    @test_throws ArgumentError MWParameters(iterations=-1)
    @test_throws ArgumentError MWParameters(repetitions=-1)
    @test_throws ArgumentError MWParameters(noisy_max_budget=0)
    @test_throws ArgumentError MWParameters(noisy_max_budget=1)
    @test_throws ArgumentError MWParameters(noisy_init=true, init_budget=0)

    histogram = Histogram([70.0, 30.0], 100)
    original_weights = copy(histogram.weights)
    state = mwem(HistogramQueries(Matrix{Float64}(I, 2, 2)), histogram,
                 MWParameters(epsilon=10, iterations=3, repetitions=2))
    @test histogram.weights == original_weights
    @test sum(state.synthetic.weights) ≈ 1.0
    @test all(isfinite, state.synthetic.weights)

    stable = Histogram([0.5, 0.5])
    PMW.update!(PMW.HistogramQuery([0.0, 1.0]), stable, 10_000.0)
    @test stable.weights ≈ [0.0, 1.0]

    Random.seed!(1234)
    settings = MWParameters(epsilon=2, iterations=0, noisy_init=true, init_budget=0.25)
    initialized = mwem(SeriesRangeQueries(2), Histogram([50.0, 50.0], 100), settings)
    @test initialized.scale ≈ 2 / (100 * 2 * (1 - 0.25))
    @test sum(initialized.synthetic.weights) ≈ 1.0
    @test all(isfinite, initialized.synthetic.weights)
    @test settings.epsilon == 2

    @test_throws ArgumentError mwem(SeriesRangeQueries(2), Histogram([0.5, 0.5]))
end

@testset "End-to-end workloads" begin
    Random.seed!(1234)
    range_data = Histogram([0.0, 1.0, 0.0], 1000)
    @test maximum_error(mwem(SeriesRangeQueries(3), range_data)) < 0.1

    Random.seed!(1234)
    table = Tabular([0 0 1 1; 0 0 1 1; 0 1 0 1])
    settings = MWParameters(epsilon=100.0)
    @test maximum_error(mwem(Parities(3, 2), table, settings)) < 0.1
    @test maximum_error(mwem(FactorParities(3, 2), table, settings)) < 0.1
end

@testset "Error metrics" begin
    state = mwem(SeriesRangeQueries(2), Histogram([3.0, 1.0], 4),
                 MWParameters(iterations=0))
    @test maximum_error(state) == 0.25
    @test mean_squared_error(state) == 0.03125
    @test kl_divergence_error(state) > 0
    @test kl_divergence([0.5, 0.5], [0.5, 0.5]) ≈ 0.0
    @test kl_divergence([1.0, 0.0], [0.5, 0.5]) ≈ 1.0
    @test isinf(kl_divergence([1.0, 0.0], [0.0, 1.0]))
    @test_throws DimensionMismatch kl_divergence([1.0], [0.5, 0.5])
    @test_throws ArgumentError kl_divergence([-1.0, 2.0], [0.5, 0.5])
end

@testset "Method hygiene" begin
    @test isempty(Test.detect_ambiguities(PMW; recursive=true))
    if isdefined(Docs, :undocumented_names)
        @test isempty(Docs.undocumented_names(PMW; private=false))
    else
        @test_skip "Docs.undocumented_names requires a newer Julia version"
    end
end
