# Construct an Economic or Desired-Gain Selection Index

Constructs a Smith-Hazel economic index or a Pesek-Baker desired-gain
index from named genetic and phenotypic covariance matrices. Trait units
and directions are checked explicitly before values are combined.

## Usage

``` r
build_selection_index(
  trait_values,
  genetic_cov,
  phenotypic_cov,
  economic_weights = NULL,
  desired_gains = NULL,
  directions,
  units = NULL,
  method = NULL,
  n_select = max(1L, floor(0.1 * nrow(trait_values))),
  quadratic_weights = NULL,
  dgsi_control = list(),
  selection_intensity = 1
)
```

## Arguments

- trait_values:

  Numeric matrix or data frame with candidates in rows and traits in
  columns. Row names are required.

- genetic_cov:

  Named genetic covariance matrix.

- phenotypic_cov:

  Named positive-definite phenotypic covariance matrix.

- economic_weights:

  Named numeric vector of economic weights. Supply exactly one of
  `economic_weights` and `desired_gains`.

- desired_gains:

  Named numeric vector of desired gains.

- directions:

  Named character vector containing `"increase"` or `"decrease"` for
  every trait.

- units:

  Optional named character vector of trait units.

- method:

  Index method. The classical methods are \`"smith_hazel"\` and
  \`"pesek_baker"\`. \`"dgsi"\` and \`"qgsi"\` are delegated to
  DesiredGainR. If \`NULL\`, the method is inferred from the supplied
  classical objective for backward compatibility.

- n_select:

  Number of candidates used to calibrate a DGSI response.

- quadratic_weights:

  Symmetric quadratic and cross-product weight matrix required by QGSI.

- dgsi_control:

  Named list of additional arguments passed to
  [`DesiredGainR::run_dgsi()`](https://rdrr.io/pkg/DesiredGainR/man/run_dgsi.html).

- selection_intensity:

  Positive selection intensity used to report the expected response
  vector. Default `1`.

## Value

A `hapblockr_result` containing index coefficients, candidate scores,
covariance diagnostics, expected responses, and index accuracy.
