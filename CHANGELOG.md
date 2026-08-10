# Changelog

All notable changes to PrivateMultiplicativeWeights.jl are documented here.
The project follows semantic versioning.

## Unreleased

### Fixed

- Modernized the histogram notebook for Julia 1.10 and newer, replaced its
  obsolete Julia 0.5 and PyPlot setup, corrected endpoint binning, and updated
  it to the current MWEM API.

### Removed

- Removed the undocumented legacy DataFrames accessibility helper and its
  test-only DataFrames dependency.

## [0.2.0] - 2026-08-08

This release revives the package for current Julia versions and repairs several
previously untested execution paths.

### Added

- Support for Julia 1.10 LTS through the latest stable Julia 1.x release.
- Input validation for privacy budgets, query dimensions, intervals,
  histograms, binary tables, sample counts, and combinatorial iterators.
- Regression coverage for histogram conversion, matrix and range workloads,
  Fourier marginals, factored histograms, noisy initialization, numerical
  stability, and error metrics.
- Cross-platform GitHub Actions CI for the minimum and latest supported Julia
  versions.

### Fixed

- Declared `Random` as a runtime dependency so the package precompiles in a
  clean environment.
- Repaired `HistogramQueries` evaluation, which previously multiplied by the
  `vec` function rather than the histogram weights.
- Corrected the orientation of `queriesMatrix(SeriesRangeQueries(...))`.
- Repaired noisy histogram initialization, including the undefined Euler
  constant, scalar broadcasting, immutable parameter mutation, zero-mass
  fallback, and remaining-budget calculation.
- Made multiplicative-weights updates stable for very large errors by updating
  in log space.
- Corrected `Marginals` to return the requested exact-order marginals.
- Corrected parity evaluation across independent factors.
- Corrected KL divergence when the reference distribution has zero support and
  fixed total squared error calculation.

### Changed

- Julia 1.10 is now the minimum supported version.
- Direct histogram input to `mwem` requires a positive `num_samples` value so
  privacy noise can be calibrated correctly.
- `mwem` normalizes a copy of histogram input instead of mutating caller-owned
  weights.
- Invalid inputs now raise descriptive `ArgumentError` or `DimensionMismatch`
  exceptions instead of assertions, `NaN` results, or downstream bounds errors.
- Factored workloads now reject `noisy_init=true`; older versions silently
  ignored the option.
- Obsolete Travis CI configuration was replaced by GitHub Actions.

## [0.1.2] - 2020-12-17

- Last release before the Julia 1.10+ modernization.
