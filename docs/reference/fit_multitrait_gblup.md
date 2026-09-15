# Fit a Multivariate Genomic BLUP

Fits a multivariate GBLUP by restricted maximum likelihood. Without a
full sampling covariance, the model uses \\V = Sigma_g \otimes K +
Sigma_e \otimes R\\, where \\R\\ carries relative precision. With a full
sampling covariance \\S\\, the default is \\V = Sigma_g \otimes K + S\\;
an explicit option adds \\Sigma_e \otimes I\\. Missing phenotype cells
are allowed; all candidates must be represented in the genomic
relationship matrix.

## Usage

``` r
fit_multitrait_gblup(
  phenotypes,
  K,
  genetic_cov = NULL,
  residual_cov = NULL,
  sampling_covariance_mode = c("sampling_only", "sampling_plus_residual"),
  estimate_covariances = is.null(genetic_cov) && is.null(residual_cov),
  maxit = 200L,
  reltol = 1e-08,
  min_reliability = 0.3
)
```

## Arguments

- phenotypes:

  Numeric candidate-by-trait matrix with dimnames, or a multi-trait
  result from
  [`prepare_breeding_targets`](https://FAkohoue.github.io/HapBlockR/reference/prepare_breeding_targets.md).

- K:

  Named genomic relationship matrix.

- genetic_cov:

  Optional named genetic covariance matrix. Supply with `residual_cov`
  to fit with fixed covariance components.

- residual_cov:

  Optional named residual covariance matrix.

- sampling_covariance_mode:

  Treatment of a full sampling covariance supplied through
  [`prepare_breeding_targets`](https://FAkohoue.github.io/HapBlockR/reference/prepare_breeding_targets.md).
  With `"sampling_only"` (default), the known matrix is the complete
  record-error covariance and no additional residual covariance is
  fitted. With `"sampling_plus_residual"`, an additional \\Sigma_e
  \otimes I\\ nugget is fitted. The option has no effect when no full
  sampling covariance is present.

- estimate_covariances:

  Logical. Estimate covariance matrices by REML.

- maxit:

  Maximum optimiser iterations.

- reltol:

  Relative optimiser tolerance.

- min_reliability:

  Minimum reliability for a recommendation.

## Value

A `hapblockr_result` with trait covariance diagnostics, multivariate
predictions, fixed-effect-adjusted PEV, reliability, fitted means, the
REML log-likelihood, and sampling-covariance provenance. Compare REML
likelihoods only between models with the same fixed-effect design.
