# Compute Haplotype-Based Genomic Relationship Matrix

Computes the additive genomic relationship matrix (GRM) from a haplotype
feature matrix using the VanRaden (2008) method extended to
multi-allelic haplotype blocks. The resulting G matrix can be used
directly in GBLUP, rrBLUP, ASReml-R, or any software that accepts a
realized relationship matrix.

## Usage

``` r
compute_haplotype_grm(hap_matrix, bend = FALSE, phased = NULL, ploidy = 2L)
```

## Arguments

- hap_matrix:

  Numeric matrix (individuals x haplotype alleles) from
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md).

- bend:

  Logical. If `TRUE`, add a small constant to the diagonal to ensure
  positive-definiteness: `diag(G) + 0.001`. Useful when passing G to
  mixed model solvers. Default `FALSE`.

- phased:

  Logical or `NULL` (default). Whether the haplotype feature matrix uses
  `additive_012` encoding for phased data (dosage values 0/1/2, dose
  scale = 2) or unphased data (presence/absence 0/1, dose scale = 1).
  `NULL` auto-detects from column maxima: any column with a value \> 1
  implies phased encoding. For guaranteed correctness when all haplotype
  alleles lack homozygous carriers (max=1 even in phased data), supply
  `phased = TRUE` or `FALSE` explicitly.

- ploidy:

  Integer \>= 2. Dose scale to use for `phased = TRUE` (or auto-detected
  phased) columns. Default `2L` (diploid, i.e. `additive_012` 0/1/2
  dosage – unchanged behaviour from previous releases). Set to e.g. `4L`
  if `hap_matrix` instead holds 0..4 autotetraploid dosage counts.
  Ignored for unphased (presence/absence 0/1) columns, which always use
  dose scale 1 regardless of `ploidy`. Note: HapBlockR's own
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md)
  always produces 0/1/2 dosage for phased data (its `hap1`/`hap2`
  representation is diploid-only), so this parameter is only meaningful
  when `hap_matrix` was constructed externally with higher-ploidy
  dosage.

## Value

Symmetric n x n numeric matrix (individuals x individuals). Row and
column names match `rownames(hap_matrix)`.

## Details

The GRM is computed as: \$\$G = \frac{ZZ^\top}{2\sum_j p_j(1-p_j)}\$\$
where \\Z\\ is the centred haplotype dosage matrix (columns centred by
\\2p_j\\) and \\p_j\\ is the frequency of haplotype allele \\j\\. This
matches the standard G matrix of Weber et al. (2023) and Difabachew et
al. (2023), ensuring the relationship scale is compatible with
conventional SNP-based GRMs.

Missing dosage values (`NA`) are mean-imputed per column before
centering.

## References

VanRaden PM (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science* **91**(11):4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

Weber SE et al. (2023). Haplotype blocks for genomic prediction: a
comparative evaluation in multiple crop datasets. *Frontiers in Plant
Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)

## Examples

``` r
data(ldx_geno, ldx_snp_info, ldx_blocks, package = "HapBlockR")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
feat <- build_haplotype_feature_matrix(haps, top_n = 5)$matrix
G    <- compute_haplotype_grm(feat)
dim(G)
#> [1] 120 120
round(range(diag(G)), 3)  # diagonal ~= 1 for typical populations
#> [1] 0.162 0.797
```
