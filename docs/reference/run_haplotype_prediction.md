# Haplotype Prediction and Block Importance from Pre-Adjusted Phenotypes

Runs the complete Tong et al. (2025) haplotype stacking pipeline using
pre-adjusted phenotype values (BLUEs, BLUPs, or adjusted entry means).
Accepts either a single trait or multiple traits simultaneously.

When a single trait is supplied, GBLUP is fitted via
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) and
block importance is ranked by `Var(local GEBV)` for that trait.

When multiple traits are supplied, a single trait-agnostic GRM is
computed once and shared across all traits.
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) is
fitted per trait using this shared GRM. Block importance is summarised
across traits with per-trait columns plus cross-trait aggregates.

## Usage

``` r
run_haplotype_prediction(
  geno_matrix,
  snp_info,
  blocks,
  blues,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  importance_rule = c("any", "all", "mean"),
  top_n = NULL,
  min_freq = 0.01,
  min_snps = 3L,
  bend = TRUE,
  marker_effect_method = c("gblup", "rrblup", "bayesb", "bayesc", "bayesa"),
  include_dominance = FALSE,
  complete_decomposition = TRUE,
  importance_threshold = 0.9,
  n_iter = 6000L,
  burn_in = 1000L,
  seed = NULL,
  min_reliability = 0.3,
  verbose = TRUE,
  ploidy = 2L
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), values 0/1/2/NA. Row names must
  be genotype IDs.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  Block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md).

- blues:

  Pre-adjusted phenotype values. See section above for accepted formats.

- id_col:

  Name of the ID column when `blues` is a data frame. Default `"id"`.

- blue_col:

  Name of the BLUE column for single-trait data frames. Default
  `"blue"`.

- blue_cols:

  Character vector of trait column names for multi-trait data frames.
  Default `NULL` – all numeric non-ID columns are used.

- importance_rule:

  How to set the combined `important` flag for multi-trait results:

  - `"any"` (default): TRUE if important for \>= 1 trait.

  - `"all"`: TRUE only if important for all traits.

  - `"mean"`: TRUE if `var_scaled_mean` \>= 0.9.

  Ignored for single-trait runs (single-trait `important` is always
  scaled variance \>= 0.9).

- top_n:

  Integer or `NULL`. Max haplotype alleles per block. Default `NULL`.

- min_freq:

  Minimum haplotype allele frequency. Default `0.01`.

- min_snps:

  Minimum SNPs per block. Default `3L`. Only affects the haplotype
  feature matrix used for the GRM/diversity (`hap_matrix`,
  `haplotypes`); the local GEBV decomposition is governed separately by
  `complete_decomposition`.

- bend:

  Logical. Add ridge to GRM diagonal. Default `TRUE`.

- marker_effect_method:

  One of `"gblup"` (default), `"rrblup"`, `"bayesb"`, `"bayesc"`,
  `"bayesa"`. See the *Marker-effect estimation methods* section.

- include_dominance:

  Logical, default `FALSE`. Fit a joint additive + dominance model
  instead of additive-only. Only valid with
  `marker_effect_method = "gblup"`; requires BGLR. See the *Additive +
  dominance GBLUP* section.

- complete_decomposition:

  Logical. Default `TRUE`. Passed to
  [`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md):
  guarantees every SNP with a known effect is attributed to some block
  (its own singleton pseudo-block if not covered by `blocks`), so
  `rowSums(local_gebv)` reconstructs each individual's full backsolved
  GEBV – matching HapSelect's completeness guarantee for `localGEBV`.
  Set `FALSE` to restore the previous behaviour (SNPs outside every
  block window are dropped).

- importance_threshold:

  Numeric in (0,1\]. Scaled-variance cutoff for the `important` flag
  (single-trait, and multi-trait `importance_rule = "mean"`). Default
  `0.9` (unchanged from previous releases). See also
  [`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  for count/percentage/cumulative-variance selection instead of a fixed
  cutoff.

- n_iter, burn_in:

  Integer. MCMC iterations / burn-in, only used when
  `marker_effect_method` is `"bayesb"`/`"bayesc"`/`"bayesa"`. Defaults
  `6000`/`1000`.

- seed:

  Integer or `NULL`. Random seed for the Bayesian methods (ignored
  otherwise). Default `NULL`.

- min_reliability:

  Numeric in \[0, 1\]. Minimum GEBV reliability for an individual to be
  marked `recommendable` in `gebv_uncertainty`. Default `0.30`. Methods
  that do not expose prediction error variance return missing
  reliability and are not promoted as recommendations.

- verbose:

  Logical. Print progress. Default `TRUE`.

- ploidy:

  Integer \>= 2. Ploidy level of `geno_matrix`'s dosage encoding.
  Default `2L` (diploid, unchanged from previous releases). Passed to
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)/
  [`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)/[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
  so marker effects and local GEBVs use ploidy-generalised VanRaden
  centring (\\\text{ploidy} \cdot p\\ instead of \\2p\\). Does **not**
  affect the shared GRM (`G`), which is always built from the
  haplotype-allele feature matrix – itself derived from HapBlockR's
  diploid-only `hap1`/`hap2` phased representation, so its own dosage
  columns are always 0/1/2 regardless of this argument. Set this only
  when `geno_matrix` itself holds raw 0..ploidy dosage calls for a
  higher-ploidy species.

## Value

Named list. For single-trait runs, contains:

- `blocks`, `diversity`, `hap_matrix`, `haplotypes`, `G`:

  Core pipeline outputs.

- `n_blocks`, `n_hap_columns`:

  Summary counts.

- `n_traits`:

  Integer: 1 for single-trait.

- `traits`:

  Character vector of trait names.

- `solver_used`:

  Character: always `"rrBLUP"` for the additive-only default, or
  `"BGLR (dual-kernel RKHS: additive + dominance)"` when
  `include_dominance = TRUE`.

- `include_dominance`:

  Logical, echoes the argument.

- `gebv`:

  Named numeric vector of (additive) GEBV.

- `gebv_uncertainty`:

  Data frame containing individual GEBV, prediction error variance,
  reliability, and the minimum-reliability recommendation gate.

- `dominance_deviation`:

  Named numeric vector of dominance deviations, or `NULL` when
  `include_dominance = FALSE`.

- `total_genetic_value`:

  `gebv + dominance_deviation`, or `NULL` when
  `include_dominance = FALSE`.

- `snp_effects`:

  Per-SNP additive effects.

- `local_gebv`:

  Matrix (individuals x blocks).

- `block_importance`:

  Data frame with `block_id`, coordinates, `var_local_gebv`,
  `var_scaled`, `important`, and `singleton` (TRUE for SNPs not covered
  by any row of `blocks`, added as their own pseudo-block under
  `complete_decomposition = TRUE`; see
  [`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)).

- `G_dominance`:

  Dominance relationship matrix, or `NULL` when
  `include_dominance = FALSE`.

- `n_train`, `n_predict`:

  Training/prediction counts.

For multi-trait runs, the same structure but `gebv`, `snp_effects`,
`local_gebv` are named lists (one per trait), and `block_importance`
contains additional per-trait and cross-trait columns as described in
the *Block importance* section. `block_importance_list` is added as a
named list of single-trait block importance data frames.

## Input format for `blues`

Accepts any of three formats:

- A **named numeric vector** – single trait, names are genotype IDs:
  `c(G001 = 4.2, G002 = 3.8)`. `id_col`/`blue_col` are ignored.

- A **data frame with one trait column** – single trait, `id_col` and
  `blue_col` specify the columns.

- A **data frame with multiple trait columns** – multi-trait, `id_col`
  names the ID column, `blue_cols` names the trait columns (or `NULL` to
  use all numeric non-ID columns).

- A **named list of named numeric vectors** – multi-trait with
  potentially different individuals per trait:
  `list(YLD = c(G001=4.2, ...), DIS = c(G001=0.3, ...))`.

## Multi-trait GBLUP solver strategy

All traits are fitted using
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) in
a per-trait loop sharing the same GRM. This produces numerically
identical results to
[`sommer::mmes()`](https://rdrr.io/pkg/sommer/man/mmes.html) for
single-trait models (verified empirically to 4 decimal places), while
avoiding the multi-trait [`cbind()`](https://rdrr.io/r/base/cbind.html)
formula which fails in sommer \<= 4.4.5 (Armadillo fixed-size matrix
error in the C++ coefficient matrix construction). Because all traits
use the same GRM, cross-trait block importance values are directly
comparable. `solver_used` is `"rrBLUP"` for the default
`marker_effect_method = "gblup"`; see the next section for the other
methods.

## Marker-effect estimation methods

`marker_effect_method` selects "Step B" of the pipeline – how per-SNP
effects are derived before summing them within each block (Step C). All
four options ultimately call
[`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
and produce a `snp_effects` vector on the same centred-dosage scale, so
`local_gebv`/`block_importance` are comparable across methods:

- `"gblup"`:

  (default)
  [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
  against the haplotype-block GRM, then backsolved to per-SNP effects
  ([`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md);
  VanRaden 2008; Tong et al. 2025). This is the method the pipeline has
  always used.

- `"rrblup"`:

  [`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)
  fit directly as a marker-effect (SNP-BLUP / ridge-regression) model –
  no relationship matrix or backsolving step (Meuwissen et al. 2001;
  Endelman 2011).

- `"bayesb"`, `"bayesc"`, `"bayesa"`:

  Bayesian variable-selection regression via
  [`BGLR::BGLR()`](https://rdrr.io/pkg/BGLR/man/BGLR.html) (Perez & de
  los Campos 2014). Requires the BGLR package. Slower (MCMC) and
  stochastic unless `seed` is set. (No `"bayesr"` option – BGLR does not
  implement the four-component-mixture "BayesR" method; see
  [`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)'s
  documentation.)

Mirrors the flexibility of HapSelect's own
`create_marker_effects_file()` (rrBLUP / BGLR / Sommer / ASReml-R).

## Additive + dominance GBLUP (`include_dominance`)

By default (`include_dominance = FALSE`), this pipeline – like
[`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
– is additive-only. Set `include_dominance = TRUE` (only with
`marker_effect_method = "gblup"`) to fit a joint additive + dominance
model instead:
[`compute_dominance_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md)
builds a dominance relationship matrix \\G_D\\ (Vitezica, Varona &
Legarra 2013) from the raw SNP dosage matrix, paired with the existing
haplotype-block additive GRM \\G_A\\ as two random-effect kernels in a
single
`BGLR::BGLR(ETA = list(list(K = G_A, model = "RKHS"), list(K = G_D, model = "RKHS")))`
call –
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
(the default GBLUP solver) only accepts a single relationship matrix and
cannot partition additive/dominance variance. Relevant for any trait
suspected of heterosis/overdominance, which
[`estimate_diplotype_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
can already flag per block via its d/a ratio. Requires the BGLR package.
Adds `dominance_deviation`, `total_genetic_value`
(`gebv + dominance_deviation`), and `G_dominance` to the return value
(all `NULL` when `include_dominance = FALSE`, so the return list's field
names are stable regardless of this argument). **Scope limit:**
per-SNP/per-block decomposition (`snp_effects`, `local_gebv`,
`block_importance`) remains additive-only even under
`include_dominance = TRUE` –
[`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
is an additive-only formula, and a block-level dominance decomposition
is not a settled methodology (see
[`compute_dominance_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md)'s
own documentation).

## Block importance – single trait

`Var(local GEBV)` is computed per block after backsolving per-SNP
effects from GEBV. Blocks are scaled to \[0,1\]; those with scaled
variance \>= 0.9 are flagged `important`. This follows Tong et al.
(2024).

## Block importance – multi-trait

Per-trait columns `var_scaled_<trait>` and `important_<trait>` are added
for every trait. Cross-trait aggregates:

- `var_scaled_mean`:

  Mean scaled variance across all traits. Primary ranking criterion –
  blocks consistently important across traits are more robust stacking
  candidates.

- `n_traits_important`:

  Count of traits for which this block is flagged important.

- `important_any`:

  TRUE if important for at least one trait.

- `important_all`:

  TRUE if important for all traits.

- `important`:

  Controlled by `importance_rule`.

## References

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection. *Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)

Endelman JB (2011). Ridge regression and other kernels for genomic
selection with R package rrBLUP. *Plant Genome* **4**:250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

Covarrubias-Pazaran G (2016). Genome-assisted prediction of quantitative
traits using the R package sommer. *PLOS ONE* **11**:e0156744.
[doi:10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744)

Vitezica ZG, Varona L & Legarra A (2013). On the additive and dominant
variance and covariance of individuals within the genomic selection
scope. *Genetics* **195**:1223-1230.
[doi:10.1534/genetics.113.155176](https://doi.org/10.1534/genetics.113.155176)

Su G, Christensen OF, Ostersen T, Henryon M & Lund MS (2012). Estimating
additive and non-additive genetic variances and predicting genetic
merits using genome-wide dense single nucleotide polymorphism markers.
*PLoS ONE* **7**:e45293.
[doi:10.1371/journal.pone.0045293](https://doi.org/10.1371/journal.pone.0045293)

## Examples

``` r
if (FALSE) { # \dontrun{
library(HapBlockR)
be     <- read_geno("mydata.vcf.gz")
blocks <- run_Big_LD_all_chr(be, CLQcut = 0.70)
geno   <- read_chunk(be, seq_len(be$n_snps))
rownames(geno) <- be$sample_ids
colnames(geno) <- be$snp_info$SNP

# -- Single trait: named vector ---------------------------------------------
blues_vec <- c(G001 = 4.21, G002 = 3.87, G003 = 5.14)
res <- run_haplotype_prediction(geno, be$snp_info, blocks, blues = blues_vec)
res$block_importance[res$block_importance$important, ]
sort(res$gebv, decreasing = TRUE)

# -- Single trait: data frame -----------------------------------------------
blues_df <- read.csv("blues.csv")   # columns: Genotype, YLD_BLUE
res <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                 blues    = blues_df,
                                 id_col   = "Genotype",
                                 blue_col = "YLD_BLUE")

# -- Single trait: additive + dominance (heterosis/overdominance) -----------
res_ad <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                    blues = blues_vec,
                                    include_dominance = TRUE)
res_ad$dominance_deviation
res_ad$total_genetic_value          # gebv + dominance_deviation
sort(res_ad$total_genetic_value, decreasing = TRUE)

# -- Multi-trait: data frame ------------------------------------------------
blues_mt <- read.csv("blues_mt.csv")  # columns: id, YLD, DIS, PHT
res_mt <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                    blues           = blues_mt,
                                    id_col          = "id",
                                    blue_cols       = c("YLD","DIS","PHT"),
                                    importance_rule = "any")
res_mt$n_traits        # 3
res_mt$solver_used     # "rrBLUP"
res_mt$block_importance[
  res_mt$block_importance$important_any,
  c("block_id","var_scaled_YLD","var_scaled_DIS","var_scaled_mean",
    "n_traits_important")]

# -- Multi-trait: named list (different individuals per trait) --------------
res_mt2 <- run_haplotype_prediction(geno, be$snp_info, blocks,
  blues = list(
    YLD = c(G001 = 4.2, G002 = 3.8),
    DIS = c(G001 = 0.3, G003 = 0.7)
  ))
} # }
```
