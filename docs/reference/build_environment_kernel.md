# Build an Environmental Similarity Kernel

Build an Environmental Similarity Kernel

## Usage

``` r
build_environment_kernel(
  environmental_covariates,
  environment_col = "environment",
  covariate_cols = NULL,
  standardise = TRUE
)
```

## Arguments

- environmental_covariates:

  Data frame containing one row per environment.

- environment_col:

  Column containing unique environment IDs.

- covariate_cols:

  Numeric environmental covariates. Defaults to every numeric column
  other than `environment_col`.

- standardise:

  Logical. Centre and scale covariates before calculating the linear
  kernel.

## Value

A positive-semidefinite environment-by-environment kernel.
