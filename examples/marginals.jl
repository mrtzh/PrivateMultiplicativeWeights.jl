using PrivateMultiplicativeWeights

function marginals(d=20, order=3, n=1000)
    data_matrix = rand(0:1, d, n)
    data_matrix[3, :] = data_matrix[1, :] .* data_matrix[2, :]
    mwem(Parities(d, order), Tabular(data_matrix))
end

function factored_marginals(d=20, order=3, n=1000)
    data_matrix = rand(0:1, d, n)
    data_matrix[3, :] = data_matrix[1, :] .* data_matrix[2, :]
    mwem(FactorParities(d, order), Tabular(data_matrix))
end
