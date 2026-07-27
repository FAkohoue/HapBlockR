# Fit a Reaction-Norm Genomic GxE Model

Fits a genomic main-effect plus reaction-norm model with record
covariance \\sigma_g^2 K_g + sigma\_{ge}^2(K_g \circ K_e) + sigma_e^2
I\\. Environment means are fitted as fixed effects by default.

## Usage

``` r
fit_gxe_gblup(
  data,
  K,
  K_environment = NULL,
  id_col = "id",
  environment_col = "environment",
  phenotype_col = "phenotype",
  precision_col = NULL,
  environment_fixed = TRUE,
  min_reliability = 0.3,
  maxit = 200L
)
```

## Arguments

- data:

  Data frame with genotype, environment, and phenotype columns.

- K:

  Named genomic relationship matrix.

- K_environment:

  Optional named environment kernel. The identity matrix gives
  independent environment-specific deviations.

- id_col, environment_col, phenotype_col:

  Column names.

- precision_col:

  Optional positive precision-weight column. Prepared breeding targets
  supply their normalised precision automatically.

- environment_fixed:

  Fit an independent fixed mean per environment.

- min_reliability:

  Minimum reliability for recommendation.

- maxit:

  Maximum REML optimiser iterations.

## Value

A `hapblockr_result` with variance components and environment-specific
genomic predictions.
