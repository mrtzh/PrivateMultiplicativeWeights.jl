using PrivateMultiplicativeWeights

function marginals(d=20, order=3, n=1000)
    data_matrix = rand(0:1, d, n)
    data_matrix[3, :] = data_matrix[1, :] .* data_matrix[2, :]
    mwem(Parities(d, order), Tabular(data_matrix), MWParameters(verbose=true))
end

function factored_marginals(d=20, order=3, n=1000)
    data_matrix = rand(0:1, d, n)
    data_matrix[3, :] = data_matrix[1, :] .* data_matrix[2, :]
    mwem(FactorParities(d, order), Tabular(data_matrix), MWParameters(verbose=true))
end

function range_queries()
    histogram = vcat(zeros(100), ones(100), zeros(100))
    data = Histogram(histogram, 100)
    mwem(SeriesRangeQueries(300), data, MWParameters(repetitions=100, verbose=true))
end

print("Elapsed time: ", @elapsed marginals())
print("\n\n")

print("Elapsed time: ", @elapsed factored_marginals())
print("\n\n")

print("Elapsed time: ", @elapsed range_queries())
print("\n\n")
