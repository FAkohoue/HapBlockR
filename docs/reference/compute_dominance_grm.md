# Genomic Dominance Relationship Matrix (Vitezica et al. 2013)

Computes a genomic dominance relationship matrix `D` from a diploid,
biallelic SNP genotype matrix (dosage-coded 0/1/2), using the natural
and orthogonal parametrization of Vitezica, Varona & Legarra (2013) –
the same parametrization implemented by
`AGHmatrix::Gmatrix(method = "Dominance")` and
[`sommer::D.mat()`](https://rdrr.io/pkg/sommer/man/D.mat.html). This is
a separate, complementary relationship matrix to the additive GRM from
[`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md):
fit together (\\G_A + G_D\\ as two random effects, e.g. via
[`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html) or ASReml)
to partition additive and dominance genetic variance – relevant for any
trait suspected of heterosis/overdominance, which
[`estimate_diplotype_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
can already flag per block via its d/a ratio.

## Usage

``` r
compute_dominance_grm(geno_matrix, ploidy = 2L, bend = TRUE)
```

## Arguments

- geno_matrix:

  Numeric matrix, individuals x SNPs, dosage-coded 0/1/2 (rows must have
  names = individual IDs).

- ploidy:

  Integer, must be `2` (diploid). Included for API consistency with
  other functions in this package; dominance coding for higher ploidy is
  not a settled area and is not implemented here.

- bend:

  Logical, default `TRUE`. Add a small ridge (`0.001`) to the diagonal
  for numerical stability in downstream mixed-model fitting, matching
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)'s
  `bend` convention.

## Value

A numeric matrix (dominance relationship matrix `D`), dimnamed by
individual ID.

## Details

For SNP \\j\\ with reference-allele frequency \\p_j\\ (from dosage mean
/ 2) and \\q_j = 1-p_j\\, the dominance incidence coding is \$\$S\_{ij}
= \begin{cases} -2q_j^2 & \text{dose} = 0 \\ 2p_jq_j & \text{dose} = 1
\\ -2p_j^2 & \text{dose} = 2 \end{cases}\$\$ and \\D = SS' / \sum_j
(2p_jq_j)^2\\. Dosage values are rounded to the nearest genotype class
(0/1/2) before coding, with a warning if any value is far from an
integer call (e.g. imputed dosage probabilities rather than hard
genotype calls) – allele frequencies themselves use the original,
unrounded dosage for a less biased estimate.

## What this does and does not do

This function computes `D` only – it does not itself fit a dual-kernel
(\\G_A + G_D\\) mixed model.
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html),
the solver
[`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
uses for its default GBLUP path, only accepts a single relationship
matrix and cannot fit \\G_A\\ and \\G_D\\ simultaneously. Pair the
output of this function with the output of
[`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
(or
[`prepare_gblup_inputs`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md)`()$G`)
and fit both as random effects in
[`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html), ASReml, or
another dual-kernel-capable solver – mirroring how
[`prepare_gblup_inputs`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md)
already hands off a single `G` matrix for external GBLUP fitting rather
than fitting the model itself.

Diploid only (`ploidy = 2`): the dominance coding below is specific to
biallelic diploid loci. Computed from the raw per-SNP genotype matrix
(not the haplotype-block feature matrix from
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md)),
since the dominance formula requires biallelic loci and haplotype-block
"alleles" can be multi-allelic (dominance coding for
multi-allelic/haplotype-block loci is not a settled area and is
deliberately not attempted here).

## References

Vitezica, Z.G., Varona, L. & Legarra, A. (2013). On the additive and
dominant variance and covariance of individuals within the genomic
selection scope. *Genetics*, 195, 1223-1230.

Su, G., Christensen, O.F., Ostersen, T., Henryon, M. & Lund, M.S.
(2012). Estimating additive and non-additive genetic variances and
predicting genetic merits using genome-wide dense single nucleotide
polymorphism markers. *PLoS ONE*, 7, e45293.

## See also

[`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md),
[`prepare_gblup_inputs`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md),
[`estimate_diplotype_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
