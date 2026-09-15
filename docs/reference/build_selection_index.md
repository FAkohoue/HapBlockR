# Construct an Economic or Desired-Gain Selection Index

Constructs a Smith-Hazel economic index or a Pesek-Baker desired-gain
index, or delegates an optimised Desired-Gain Selection Index (DGSI) or
Quadratic Genomic Selection Index (QGSI) to DesiredGainR. Trait
directions are required; declared units are optional metadata.

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
  n_select = NULL,
  quadratic_weights = NULL,
  dgsi_control = list(),
  qgsi_control = list(),
  selection_intensity = NULL
)
```

## Arguments

- trait_values:

  Numeric matrix or data frame with candidates in rows and traits in
  columns. Row names are required. DGSI uses internally modelled trait
  values; QGSI uses genomic estimated breeding values. The trait name
  `"id"` is reserved by the delegated interface.

- genetic_cov:

  Named genetic covariance matrix in original trait units. Used in the
  classical and DGSI coefficient calculations. For QGSI, this matrix is
  retained as upstream context; it is not substituted for Gamma.

- phenotypic_cov:

  Named positive-definite phenotypic or index-variable covariance matrix
  in original trait units. Used by classical methods and DGSI; retained
  as upstream context for QGSI.

- economic_weights:

  Named numeric vector of economic weights. Supply exactly one of
  `economic_weights` and `desired_gains`. Values must be non-negative
  importance magnitudes after trait orientation; `directions`, rather
  than a negative weight, declares decrease.

- desired_gains:

  Named numeric vector of non-negative desired-gain magnitudes in the
  favourable-direction trait space. For Pesek-Baker these are in
  original trait units. For DGSI these are in candidate standard
  deviations, regardless of `dgsi_control$scale_traits`; divide an
  original-unit target by its candidate standard deviation before
  passing it.

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

  Number of candidates selected by DGSI or QGSI. DGSI may retain fewer
  when fewer candidates satisfy its eligibility thresholds. When omitted
  for either delegated method, the default is 10 percent of candidates,
  with a minimum of one. It does not apply to classical indices.

- quadratic_weights:

  Symmetric quadratic and cross-product weight matrix required by QGSI.
  It must describe the favourable-direction trait space obtained after
  applying `directions` and any requested scaling. QGSI linear economic
  weights must also refer to that analysis scale.

- dgsi_control:

  Named list of additional arguments passed to
  [`DesiredGainR::run_dgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_dgsi.html).
  Names must match exactly; partial names and structural overrides are
  rejected. Reference and validation data contain the same original-unit
  trait columns as `trait_values`.

- qgsi_control:

  Named list of additional arguments passed to
  [`DesiredGainR::run_qgsi()`](https://FAkohoue.github.io/DesiredGainR/reference/run_qgsi.html).
  Structural arguments supplied by HapBlockR, including the candidate
  data, weights, directions and selection count, cannot be overridden.
  Names must match exactly. Supply `Gamma` in original trait units, or
  `relationship_matrix` with row and column names matching the reference
  IDs. Reference data use an `id` column and the original-unit trait
  columns. Without an explicit reference, candidate values form the
  reference population.

- selection_intensity:

  Optional positive selection intensity used to report the expected
  response vector for the classical Smith-Hazel and Pesek-Baker methods;
  the default for those methods is one. DGSI and QGSI derive normal
  selection intensity from the proportion actually selected under their
  requested rule, so this argument does not apply to delegated methods.

## Value

A `hapblockr_result` containing candidate scores, objective information,
covariance diagnostics, and method-appropriate response summaries.
Classical and DGSI methods return linear coefficients. QGSI instead
returns linear weights, the quadratic-weight matrix, and
candidate-specific contributions because it has no single global
coefficient vector. The unmodified delegated fit is in `engine_result`
and can be passed to DesiredGainR comparison tools. DGSI's
`coefficients$coefficient` remains on the engine's favourable-direction
analysis scale. Use the named `coefficients_original_units` vector to
combine original-unit marker, haplotype or block effects. Candidate
scores equal
`trait_values %*% coefficients_original_units + score_intercept`; the
intercept accounts for reference centring.
