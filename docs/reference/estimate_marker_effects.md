# Estimate Per-SNP Marker Effects via GBLUP, rrBLUP (SNP-BLUP), or Bayesian Regression

Unified "Step B" marker-effect estimator for the haplotype-stacking
pipeline (see
[`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)),
mirroring the flexibility HapSelect's own `create_marker_effects_file()`
offers (rrBLUP / BGLR / Sommer / ASReml-R – "any model that returns a
vector of marker effects aligned to the same map"). Returns a named
vector of per-SNP additive effects on the centred-dosage (\\2p\\) scale
used throughout HapBlockR, ready for
[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md).

- `"gblup"`:

  (default) Fits
  [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
  against a supplied genomic/haplotype relationship matrix `G`, then
  backsolves per-SNP effects from the resulting GEBV via
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  (VanRaden 2008; Tong et al. 2025), using every individual present in
  `geno_matrix`/`G`. This is the same computation
  [`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
  has always used internally, except that function's own bespoke "gblup"
  code path restricts backsolving to the phenotyped subset for backward
  compatibility; calling `estimate_marker_effects(method = "gblup")`
  directly backsolves over the full population instead. Requires `G`.

- `"rrblup"`:

  Fits
  [`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)
  directly as a marker-effect (SNP-BLUP / ridge-regression BLUP) model,
  \\y = \mu + Zu + e\\ with \\u \sim N(0, \sigma_u^2 I)\\ (Meuwissen et
  al. 2001; Endelman 2011), returning `u` directly – no relationship
  matrix or backsolving step needed. Under the standard
  equal-per-SNP-variance assumption this is mathematically equivalent to
  VanRaden GBLUP-then-backsolve, computed the other way around.

- `"bayesb"`, `"bayesc"`, `"bayesa"`:

  Fit [`BGLR::BGLR()`](https://rdrr.io/pkg/BGLR/man/BGLR.html) with the
  corresponding variable-selection prior (Perez & de los Campos 2014)
  directly on the genotype matrix, taking the posterior mean marker
  effects as `alpha`. `"bayesa"` is BayesA (marker-specific variances
  from a scaled-inverse-chi-squared prior, no point mass at zero);
  `"bayesb"`/`"bayesc"` add a point mass at zero (variable selection).
  Note there is no `"bayesr"` option: BGLR does not implement the
  four-component mixture "BayesR" method (Erbe et al. 2012) – an earlier
  version of this function incorrectly assumed it did, which fails at
  runtime with
  [`BGLR::BGLR()`](https://rdrr.io/pkg/BGLR/man/BGLR.html)'s own "model
  BayesR not implemented" error; this was caught by running the test
  suite and `"bayesr"` was replaced with the real BGLR model `"bayesa"`.
  Requires the BGLR package (`Suggests`). BGLR runs an MCMC sampler and
  writes temporary diagnostic files (removed automatically unless
  `bglr_dir` is supplied); this is markedly slower than
  `"gblup"`/`"rrblup"` and, being stochastic, is only exactly
  reproducible when `seed` is set.

All four methods return `gebv` on the same scale (\\\hat g = M
\hat\alpha\\, the centred genotype matrix times the estimated marker
effects), so results are directly comparable across `method` choices –
e.g. to sanity-check one method's block ranking against another's.

## Usage

``` r
estimate_marker_effects(
  geno_matrix,
  y,
  method = c("gblup", "rrblup", "bayesb", "bayesc", "bayesa"),
  G = NULL,
  n_iter = 6000L,
  burn_in = 1000L,
  bglr_dir = NULL,
  seed = NULL,
  verbose = FALSE,
  ploidy = 2L
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), values 0/1/2/NA. Row names are
  genotype IDs; column names are SNP IDs.

- y:

  Named numeric vector of phenotypes (BLUEs/BLUPs/adjusted means). Names
  must match `rownames(geno_matrix)` for at least some individuals. `NA`
  values and individuals absent from `geno_matrix` are dropped before
  fitting (allowed for all methods; `"gblup"` additionally supports true
  missing-phenotype prediction candidates via `kin.blup()`'s own
  handling, but the estimation step here always fits on the non-missing
  subset).

- method:

  One of `"gblup"` (default), `"rrblup"`, `"bayesb"`, `"bayesc"`,
  `"bayesa"`.

- G:

  Genomic/haplotype relationship matrix (n x n, dimnames = genotype
  IDs). Required for `method = "gblup"`; ignored otherwise.

- n_iter, burn_in:

  Integer. MCMC iterations / burn-in for the Bayesian methods. Defaults
  `6000`/`1000`. Ignored for `"gblup"`/`"rrblup"`.

- bglr_dir:

  Directory for BGLR's temporary diagnostic files. Default `NULL`: a
  per-call [`tempfile()`](https://rdrr.io/r/base/tempfile.html)
  directory is created and removed automatically when the fit completes.
  Ignored for `"gblup"`/`"rrblup"`.

- seed:

  Integer or `NULL`. Random seed set before the Bayesian sampler runs
  (via [`set.seed()`](https://rdrr.io/r/base/Random.html)), for
  reproducibility. Ignored for `"gblup"`/`"rrblup"`. Default `NULL`.

- verbose:

  Logical. For Bayesian methods, passed through to `BGLR()`'s own
  (otherwise very chatty) console output. Default `FALSE`.

- ploidy:

  Integer \>= 2. Ploidy level of `geno_matrix`'s dosage encoding.
  Default `2L` (diploid). See
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  for the generalised centering used. Passed through to
  [`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  for `method = "gblup"`; for `"rrblup"`/Bayesian methods it only
  affects the `gebv` centering (\\M\hat\alpha\\), since those solvers
  fit directly on the raw dosage matrix and are otherwise
  ploidy-agnostic. **In practice this centering step is
  ploidy-invariant**: the centring term is \\\text{ploidy} \cdot \hat
  p\\, and \\\hat p = \text{colMeans(geno\\matrix)} / \text{ploidy}\\,
  so \\\text{ploidy} \cdot \hat p\\ always simplifies to
  `colMeans(geno_matrix)` regardless of `ploidy` (as long as the
  internal safety clamp on \\\hat p\\ does not engage, which it will not
  for realistic dosage ranges). So changing `ploidy` will *not* change
  `"rrblup"`/Bayesian `gebv` in practice – this is a mathematical
  property of mean-centering, not a bug. `ploidy`'s real, discriminating
  effect elsewhere in the package is on *scaling* terms (e.g.
  [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)'s
  \\\text{ploidy} \cdot \sum \hat p (1-\hat p)\\ denominator).

## Value

A named list:

- `snp_effects`:

  Named numeric vector, length `ncol(geno_matrix)`, one effect per SNP
  (0 for any SNP a variable-selection method shrank exactly to zero).

- `gebv`:

  Named numeric vector (one per individual in `geno_matrix`), \\M
  \hat\alpha\\ – the additive genomic merit implied by `snp_effects`, on
  a consistent scale across methods.

- `method`:

  Character, echoes the `method` argument.

- `fit`:

  The raw fitted model object (`kin.blup()`, `mixed.solve()`, or
  `BGLR()` output) for diagnostics.

## References

Meuwissen THE, Hayes BJ, Goddard ME (2001). Prediction of total genetic
value using genome-wide dense marker maps. *Genetics*
**157**(4):1819-1829.

Endelman JB (2011). Ridge regression and other kernels for genomic
selection with R package rrBLUP. *Plant Genome* **4**:250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

Perez P, de los Campos G (2014). Genome-wide regression and prediction
with the BGLR statistical package. *Genetics* **198**(2):483-495.
[doi:10.1534/genetics.114.164442](https://doi.org/10.1534/genetics.114.164442)

VanRaden PM (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science* **91**(11):4414-4423.

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

## See also

[`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md),
[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
me_gblup  <- estimate_marker_effects(geno, y, method = "gblup",  G = G)
me_rr     <- estimate_marker_effects(geno, y, method = "rrblup")
me_bayesb <- estimate_marker_effects(geno, y, method = "bayesb", seed = 1)
loc <- compute_local_gebv(geno, snp_info, blocks, me_rr$snp_effects)
} # }
```
