# Finlay-Wilkinson Stability Analysis of Haplotype Effects Across Environments

Estimates the stability of each haplotype block's effect across multiple
environments using Finlay-Wilkinson (1963) regression. For each block,
the per-environment GEBV contribution is regressed on the environment
mean (the environmental index). A regression slope b_i = 1 indicates
average stability; b_i \> 1 = above-average response (exploits good
environments); b_i \< 1 = below-average response (robust across
environments).

## Usage

``` r
run_haplotype_stability(
  geno_matrix,
  snp_info,
  blocks,
  blues_list,
  top_n = NULL,
  min_freq = 0.05,
  min_snps = 3L,
  verbose = TRUE
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs).

- snp_info:

  Data frame with `SNP`, `CHR`, `POS`.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md).

- blues_list:

  Named list of named numeric vectors, one per environment:
  `list(env1 = c(id1=val,...), env2 = ...)`.

- top_n:

  Integer or `NULL`. Max haplotype alleles per block. Default `NULL`.

- min_freq:

  Numeric. Minimum haplotype allele frequency. Default `0.05`.

- min_snps:

  Integer. Minimum SNPs per block. Default `3L`.

- verbose:

  Logical. Default `TRUE`.

## Value

Data frame with one row per block, sorted by `CHR` and `start_bp`.

- `block_id`, `CHR`, `start_bp`, `end_bp`:

  Block coordinates.

- `b`:

  Finlay-Wilkinson slope. b=1: average stability; b\>1: exploits good
  environments; b\<1: robust across environments.

- `b_se`:

  Standard error of b.

- `r2_fw`:

  R-squared of the FW regression.

- `s2d`:

  Deviation mean square (non-linear instability).

- `stable`:

  Logical; TRUE when b is not significantly different from 1 at
  alpha=0.05.

## References

Finlay KW, Wilkinson GN (1963). The analysis of adaptation in a
plant-breeding programme. *Australian Journal of Agricultural Research*
**14**(6):742-754.

## See also

[`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md),
[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "HapBlockR")
# Simulate two environments by splitting the BLUEs
b1 <- setNames(ldx_blues$YLD + rnorm(nrow(ldx_blues), 0, 0.1),
               ldx_blues$id)
b2 <- setNames(ldx_blues$YLD + rnorm(nrow(ldx_blues), 0.5, 0.15),
               ldx_blues$id)
stab <- run_haplotype_stability(
  geno_matrix = ldx_geno,
  snp_info    = ldx_snp_info,
  blocks      = ldx_blocks,
  blues_list  = list(env1 = b1, env2 = b2),
  verbose     = FALSE
)
head(stab[order(stab$b), ])
#>                block_id CHR start.bp end.bp b b_se r2_fw s2d p_b1 stable
#> 1    block_1_1000_25027   1     1000  25027 0  NaN     1  NA  NaN  FALSE
#> 2 block_1_155368_179371   1   155368 179371 0  NaN     1  NA  NaN  FALSE
#> 3   block_1_81064_99022   1    81064  99022 0  NaN     1  NA  NaN  FALSE
#> 4    block_2_1000_30023   2     1000  30023 0  NaN     1  NA  NaN  FALSE
#> 5 block_2_161515_180473   2   161515 180473 0  NaN     1  NA  NaN  FALSE
#> 6  block_2_86236_105290   2    86236 105290 0  NaN     1  NA  NaN  FALSE
# }
```
