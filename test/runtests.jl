using PrivateMultiplicativeWeights
using Hadamard
using Test
using Random

# test histogram transform
@test Histogram(Tabular(convert(Array{Int64, 2}, zeros(1,10)))).weights == [1.0, 0.0]
@test Histogram(Tabular(convert(Array{Int64, 2}, ones(1,10)))).weights == [0.0, 1.0]
@test Histogram(Tabular([1 0; 0 1])).weights == [0.0, 0.5, 0.5, 0]
@test Histogram(Tabular([1 0 1; 0 1 0])).weights == [0.0, 1.0, 2.0, 0]/3.0
@test Histogram(Tabular([1 0 1 1 ; 0 1 0 1])).weights == [0.0, 1.0, 2.0, 1.0]/4.0
@test Histogram(Tabular([1 0 1 0 1; 0 1 0 0 1])).weights == [1.0, 1.0, 2.0, 1.0]/5.0
@test Tabular(Histogram([1.0,0.0,0.0,0.0]), 10).data == zeros(2, 10)
@test Tabular(Histogram([0.0,0.0,0.0,1.0]), 10).data == ones(2, 10)

# test our hadamard basis vectors agree with Hadamard module
for j = 0:10
    for i = j:10
        @test PrivateMultiplicativeWeights.hadamard_basis_vector(j,i) == Hadamard.hadamard(2^i)[:,j+1]
    end
end

# test range queries
Random.seed!(1234)
data = Histogram([0.0, 1.0, 0.0], 1000)
@test maximum_error(mwem(SeriesRangeQueries(3), data)) < 0.1

# test parities
Random.seed!(1234)
data = Tabular([0 0 1 1; 0 0 1 1; 0 1 0 1])
@test maximum_error(mwem(Parities(3, 2), data, MWParameters(epsilon=100.0))) < 0.1
@test maximum_error(mwem(FactorParities(3, 2), data, MWParameters(epsilon=100.0))) < 0.1
