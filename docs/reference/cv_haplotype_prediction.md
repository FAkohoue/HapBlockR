# K-Fold Cross-Validation for Haplotype-Based Genomic Prediction

Estimates the predictive ability of the haplotype GBLUP model via k-fold
cross-validation. In each fold, a subset of individuals is masked from
the phenotype and predicted from the haplotype GRM; Pearson correlation
between predicted and observed BLUEs is returned as the predictive
ability (PA). Runs per trait when multiple traits are supplied.

## Usage

``` r
cv_haplotype_prediction(
  geno_matrix,
  snp_info,
  blocks,
  blues,
  k = 5L,
  n_rep = 1L,
  top_n = NULL,
  min_freq = 0.05,
  min_snps = 3L,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  seed = 42L,
  validation = c("random", "grouped", "forward"),
  groups = NULL,
  time = NULL,
  on_fit_error = c("error", "record"),
  verbose = TRUE
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), MAF-filtered dosage.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md).

- blues:

  Pre-adjusted phenotype means. Accepts the same four formats as
  [`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md):
  named numeric vector, single-trait data frame, multi-trait data frame,
  or named list.

- k:

  Integer. Number of folds. Default `5L`.

- n_rep:

  Integer. Number of CV replications (each with a different random fold
  assignment). Default `1L`.

- top_n:

  Integer or `NULL`. Maximum haplotype alleles per block passed to
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md).
  Default `NULL` (all alleles above `min_freq`).

- min_freq:

  Numeric. Minimum haplotype allele frequency. Default `0.05`.

- min_snps:

  Integer. Minimum SNPs per block for haplotype extraction. Default
  `3L`.

- id_col:

  Character. Name of the individual ID column when `blues` is a data
  frame. Default `"id"`.

- blue_col:

  Character. Name of the BLUE column for single-trait data frames.
  Default `"blue"`.

- blue_cols:

  Character vector. Trait column names for multi-trait data frames.
  Default `NULL` (auto-detect all numeric non-ID columns).

- seed:

  Integer. RNG seed for reproducible fold assignment. Default `42L`.

- validation:

  Validation design. `"random"` assigns individuals to balanced folds;
  `"grouped"` keeps all individuals in the same group in one test fold;
  `"forward"` trains only on earlier time points. Default `"random"`.

- groups:

  Named vector mapping individual identifiers to groups. Required when
  `validation = "grouped"`.

- time:

  Named numeric or ordered vector mapping individual identifiers to
  breeding cycles, years, or other ordered time points. Required when
  `validation = "forward"`.

- on_fit_error:

  Behaviour when
  [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
  fails. `"error"` stops with trait, repetition, and fold context;
  `"record"` retains the failed fold with missing predictions and an
  explicit error message. Default `"error"`.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

A named list of class `HapBlockR_cv`:

- `pa_summary`:

  Data frame: `trait`, `rep`, `fold`, `n_train`, `n_test`, `PA` (Pearson
  r), `RMSE`, fit status, and any fit error.

- `pa_pooled`:

  One row per trait and repetition, calculated from the complete pooled
  out-of-fold predictions: predictive ability, root mean squared error,
  mean absolute error, bias, and calibration slope.

- `pa_mean`:

  Mean pooled predictive ability and root mean squared error per trait
  across replications.

- `gebv_all`:

  Out-of-fold predictions for every tested individual, trait, and
  repetition. `gebv` is the predicted phenotype on the BLUE scale;
  `breeding_value` is the centred random genetic deviation returned by
  `kin.blup()`.

- `k`:

  Number of folds used.

- `n_rep`:

  Number of replications.

- `validation`:

  Validation design used.

## See also

[`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "HapBlockR")
cv <- cv_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info    = ldx_snp_info,
  blocks      = ldx_blocks,
  blues       = ldx_blues,
  k           = 5L,
  id_col      = "id",
  verbose     = FALSE
)
cv$pa_mean
#>   trait          PA      RMSE PA_sd RMSE_sd
#> 1   RES  0.32163303 0.9464001    NA      NA
#> 2   YLD -0.06058831 1.0222290    NA      NA
# }
```
