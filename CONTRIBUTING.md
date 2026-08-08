# Contributing

Contributions are welcome. Bug fixes should include a regression test, and
changes to public behavior should update both the README and changelog.

## Development setup

Install Julia 1.10 or a newer stable Julia release, clone the repository, and
instantiate the project:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

`Manifest.toml` is intentionally ignored because this is a library package. Do
not include a locally generated manifest in a commit.

## Running tests

Run the package suite from the repository root:

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Before submitting a change, also run the stricter configuration used during
maintenance:

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test(julia_args=["--startup-file=no", "--depwarn=error", "--check-bounds=yes"])'
```

The tests cover:

- histogram and binary-table conversion;
- histogram, range, parity, and factored workloads;
- noisy initialization and multiplicative-weights updates;
- marginal reconstruction and error metrics;
- invalid-input behavior and method ambiguity checks;
- the DataFrames accessibility helper.

CI runs the suite on the oldest compatible Julia minor and the latest stable
Julia 1.x release across Linux, macOS, and Windows.

## Code guidelines

- Keep public constructors permissive about concrete numeric container types,
  but validate invariants at API boundaries.
- Preserve the `d × n` tabular convention: attributes are rows and records are
  columns.
- Do not expose exact answers or source-dependent diagnostics as private output.
- Avoid mutating caller-owned inputs unless the function name ends in `!` and
  the behavior is documented.
- Prefer deterministic unit tests. Seed `Random` for stochastic MWEM tests and
  assert meaningful tolerances rather than exact random draws.
- Keep explicit-histogram tests small; their memory use grows as `2^d`.

## Adding a workload

Custom workloads subtype `PrivateMultiplicativeWeights.Query` and `Queries`.
Depending on the data representation, implement:

- `get(queries, index)`;
- `evaluate(query, data)`;
- `evaluate(queries, data)`;
- `initialize(queries, data, parameters)`;
- `update!(query, synthetic, error)`;
- `normalize!(synthetic)`.

Use `src/parities.jl` for an implicit workload and `src/histogram.jl` for a
linear histogram workload as references.

## Pull-request checklist

- Tests pass on a supported Julia release.
- New behavior has focused tests.
- Public API or privacy changes are documented.
- `git diff --check` reports no whitespace errors.
- No generated manifest, coverage files, notebooks checkpoints, or local data
  are included.
