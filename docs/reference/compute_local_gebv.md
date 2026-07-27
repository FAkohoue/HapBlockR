# Compute Local Haplotype GEBV per Block (Tong et al. 2025)

For each LD block, computes the local GEBV (haplotype effect) for every
individual by summing the per-SNP additive effects of the alleles they
carry within the block:

\$\$\text{local GEBV}\_f = \sum\_{t \in f} (x_t - 2 p_t)\\\alpha_t\$\$

where \\x_t\\ is the allele dosage (0/1/2) at SNP \\t\\, \\p_t\\ is the
SNP's own population allele frequency, and \\\alpha_t\\ is its additive
effect. Centring at \\2p_t\\ (rather than at the heterozygote midpoint,
dosage = 1) matches the convention used by
[`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
to derive `alpha` in the first place (\\\text{GEBV} = M\alpha\\ with
\\M\\ centred at \\2p\\), so that summing `local_gebv` across every
block for one individual reconstructs that individual's overall
backsolved GEBV.

Blocks are ranked by `Var(local GEBV)` – blocks with high variance
contribute strongly to trait differences among individuals and likely
harbour causal loci (Tong et al. 2025). This ranking is unaffected by
the choice of centring point (variance is shift-invariant), but the
absolute `local_gebv` values are only meaningful – i.e. only summable
back to the genome-wide GEBV – under the \\2p_t\\ convention used here.

## Usage

``` r
compute_local_gebv(
  geno_matrix,
  snp_info,
  blocks,
  snp_effects,
  snp_effect_se = NULL,
  scale = TRUE,
  complete_decomposition = TRUE,
  importance_threshold = 0.9,
  ploidy = 2L
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), values 0/1/2/NA.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  Block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md).

- snp_effects:

  Named numeric vector of per-SNP additive effects from
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md),
  [`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md),
  or from a marker model directly.

- snp_effect_se:

  Optional named non-negative numeric vector of per-SNP effect standard
  errors. When supplied, local GEBV standard errors are propagated under
  an independent-marker-error approximation.

- scale:

  Logical. If `TRUE` (default), scale `Var(local GEBV)` to \[0,1\] so
  blocks are comparable across traits and datasets.

- complete_decomposition:

  Logical. If `TRUE` (default), append a singleton pseudo-block for
  every SNP with a known effect that isn't covered by any row of
  `blocks`, so the decomposition is complete (see Complete decomposition
  section). If `FALSE`, restores the previous behaviour where such SNPs
  are silently excluded.

- importance_threshold:

  Numeric in (0,1\]. Scaled-variance cutoff for the `important` flag.
  Default `0.9` (unchanged from previous releases, where this was
  hardcoded). Raise it (e.g. `0.95`) for a more conservative set of
  "important" blocks, or lower it to flag more. See also
  [`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  for count/percentage/ cumulative-variance selection instead of a fixed
  cutoff.

- ploidy:

  Integer \>= 2. Ploidy level of `geno_matrix`'s dosage encoding.
  Default `2L` (diploid, unchanged from previous releases). Generalises
  the \\2p_t\\ centring above to \\\text{ploidy} \cdot p_t\\; must match
  the `ploidy` used to derive `snp_effects` (via
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  or
  [`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md))
  for local GEBVs to sum correctly to the genome-wide GEBV.

## Value

A list with two elements:

- `local_gebv`:

  Numeric matrix (individuals x blocks) of per-block local GEBV values.
  Includes singleton pseudo-block columns when
  `complete_decomposition = TRUE` and any exist.

- `local_gebv_se`:

  Numeric matrix of propagated local GEBV standard errors, or `NULL`
  when `snp_effect_se` is not supplied.

- `block_importance`:

  Data frame with one row per block: block_id, CHR, start_bp, end_bp,
  n_snps, var_local_gebv, var_scaled, important (logical: scaled
  variance \>= `importance_threshold`), and `singleton` (logical: TRUE
  for SNPs never covered by an LD block, added under
  `complete_decomposition = TRUE`).

## Complete decomposition (`complete_decomposition = TRUE`)

`blocks` (typically from
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md))
may not cover every SNP: isolated markers that never cliqued with a
neighbour are dropped from the block table entirely unless
`singleton_as_block = TRUE` was used at detection time. Any SNP absent
from every block window would otherwise silently vanish from the
decomposition – its share of the genome-wide GEBV is real (it has a
non-zero entry in `snp_effects`) but never attributed to any block, so
`sum(local_gebv[i, ])` would fall short of individual `i`'s true
backsolved GEBV.

When `complete_decomposition = TRUE` (the default), every SNP present in
both `colnames(geno_matrix)` and `names(snp_effects)` that isn't already
covered by a row of `blocks` is appended as its own single-SNP
"singleton block" (`block_id` prefixed `"singleton_"`, `n_snps = 1`,
`singleton = TRUE` in `block_importance`). This mirrors HapSelect's
`def_blocks()` behaviour, where an isolated marker is always retained as
its own one-marker haploblock rather than dropped, and guarantees
`rowSums(local_gebv)` reconstructs the individual's full backsolved GEBV
exactly (up to any SNPs missing from `snp_info` entirely, which cannot
be spatially placed and are reported via a warning). Set to `FALSE` to
restore the previous behaviour (only SNPs inside an existing block row
contribute).

**This guarantee depends on `geno_matrix` being the full, unfiltered,
genome-wide genotype matrix** – the same one used to derive
`snp_effects` – not a matrix already subset to block-member SNPs (e.g. a
haplotype feature matrix). If it has been pre-filtered, any SNP dropped
upstream was never a column here and cannot be recovered; the output
will look complete without necessarily being so. A warning fires when
this can be detected: if `snp_info`'s own coordinates place SNPs outside
every block window, but none of those SNPs are present in
`colnames(geno_matrix)` at all, `geno_matrix` is almost certainly
pre-filtered.

## References

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

## Examples

``` r
if (FALSE) { # \dontrun{
snp_fx  <- backsolve_snp_effects(my_geno, gebv)
loc     <- compute_local_gebv(my_geno, snp_info, blocks, snp_fx)
# Blocks with scaled variance >= 0.9 are most important
head(loc$block_importance[loc$block_importance$important, ])
# Every individual's local GEBVs now sum to their backsolved GEBV:
recon_gebv <- rowSums(loc$local_gebv)
all.equal(unname(recon_gebv[names(gebv)]), unname(gebv))
} # }
```
