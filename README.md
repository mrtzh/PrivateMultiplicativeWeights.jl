# PrivateMultiplicativeWeights.jl

A simple and practical algorithm for differentially private data release.

MIT Licensed. See `LICENSE.md`.

## Installation

Open a Julia prompt and call: `Pkg.clone("https://github.com/mrtzh/PrivateMultiplicativeWeights.jl.git")`

## Main Features

* Differentially private synthetic data preserving lower order marginals of an input data set
* Optimized in-memory implementation for small number of data attributes
* Scalable heuristic for large number of data attributes
* Easy-to-use interfaces for custom query sets and data representations

## Usage
For illustration, we create a random data set with hidden correlations. Columns correspond to data points.
```
d, n = 20, 1000
data_matrix = rand(0:1,d,n)
data_matrix[3,:] = data_matrix[1,:] .* data_matrix[2,:]
```

### Histograms

We can run MWEM to produce synthetic data accurate for 1st, 2nd, 3rd order marginals of the source data.
```
mw = mwem(Parities(d,3),Tabular(data_matrix))
```
This will convert the data to its explicit histogram representation of size 2^d and may not be useful when d is large. See section below on factored histograms for a scalable alternative.

We can tweak various parameters:
```
mw = mwem(Parities(d,3),Tabular(data_matrix),epsilon=0.5,iterations=10,repetitions=10,smart=false)
```
Parameters:

| Name | Default | Description |
| ---- | ------- | ----------- |
| `epsilon` | `1.0` |  privacy parameter |
| `iterations` | `10` | number of iterations |
| `repetitions`| `5` | controls how often MWEM cycles through previously measured queries per iteration |
| `smart` | `false` | smart histogram initalization through noise addition. When `smart` is set to false, the initialization is uniform. |
| `verbose` | `false` | print timing and error statistics per iteration (information is not differentially private)

We can convert synthetic data in histogram representation to matrix tabular representation.
```
table = Tabular(mw.synthetic)
```
Compute error achieved by MWEM:
```
maximum_error(mw), mean_squared_error(mw)
```
Note that these statistics are *not* differentially private.

#### Custom query sets

You can define custom query workloads by using `HistogramQueries(query_matrix)` instead of `Parities(d,3)`. Here `query matrix` is an `N x k` matrix specifying the query set in its Histogram representation, `N` is the histogram length and `k` is the `k` is the number of queries.

To build query sets with your own implicit representations, create a sub-type of `Query` and `Queries`, respectively. See `interface.jl`.

### Factored histograms
When the histogram representation is too large, try using factored histograms. Factored histograms maintain a product distribution over clusters of attributes of the data. Each component is represented using a single histogram. Components are merged as it becomes necessary. This often allows to scale up MWEM by orders of magnitude.
```
d, n = 100, 1000
data_matrix = rand(0:1,d,n)
data_matrix[3,:] = data_matrix[1,:] .* data_matrix[2,:]
mw = mwem(FactorParities(d,3),Tabular(data_matrix))
```

Also see `examples.jl`.

## Citing this package

The MWEM algorithm was presented in the following paper:
```
@inproceedings{HLM12,
	author = "Moritz Hardt and Katrina Ligett and Frank McSherry",
	title = "A simple and practical algorithm for differentially-private data release",
    title = {Proc.\ 26th Neural Information Processing Systems (NIPS)},
    booktitle = {Proc.\ 26th Neural Information Processing Systems (NIPS)},
    year = {2012},
}
```
