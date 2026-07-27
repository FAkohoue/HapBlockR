# Assess Phasing Accuracy against a Truth Set

Compares estimated phased, diploid, biallelic genotypes with a phased
truth set. Variant, sample, chromosome, position, and available allele
identities are checked before calculation. Haplotype labels may be
globally swapped within each sample and chromosome. A switch error is
counted when the estimated orientation relative to truth changes between
consecutive informative heterozygous variants on the same chromosome.

## Usage

``` r
assess_phasing_accuracy(
  truth,
  estimate,
  max_switch_error_rate = 1,
  min_dosage_accuracy = 0,
  min_allele_concordance = 0,
  min_call_rate = 0,
  require_switch_information = TRUE,
  strict = TRUE
)
```

## Arguments

- truth:

  A phased VCF path or a list containing `hap1`, `hap2`, and preferably
  `snp_info`, as returned by
  [`read_phased_vcf`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md).

- estimate:

  Estimated phased data in the same form as `truth`.

- max_switch_error_rate:

  Maximum acceptable switch-error rate. Default `1`.

- min_dosage_accuracy:

  Minimum acceptable exact dosage accuracy among truth genotypes for
  which an estimated call is available. Default `0`.

- min_allele_concordance:

  Minimum acceptable phased-allele concordance after the best global
  haplotype-label orientation is selected within each sample and
  chromosome. Default `0`.

- min_call_rate:

  Minimum acceptable estimated call rate among observed truth genotypes.
  Default `0`.

- require_switch_information:

  Logical. Require at least one transition between informative
  heterozygous variants. Default `TRUE`.

- strict:

  Logical. Stop when a quality gate fails. Default `TRUE`.

## Value

A `"HapBlockR_phasing_accuracy"` and `"hapblockr_result"` object
containing overall metrics, per-sample metrics, per-chromosome metrics,
quality gates, input hashes, and immutable identifiers.

## Examples

``` r
truth <- list(
  hap1 = matrix(c(0, 0, 1), ncol = 1,
                dimnames = list(paste0("s", 1:3), "P1")),
  hap2 = matrix(c(1, 1, 0), ncol = 1,
                dimnames = list(paste0("s", 1:3), "P1")),
  snp_info = data.frame(
    SNP = paste0("s", 1:3), CHR = "1", POS = 1:3,
    REF = "A", ALT = "G"
  )
)
accuracy <- assess_phasing_accuracy(truth, truth)
accuracy$metrics$switch_error_rate
#> [1] 0
validate(accuracy)
#>                 check passed                  detail
#> 1     required_fields   TRUE                        
#> 2      schema_version   TRUE                   1.0.0
#> 3              method   TRUE assess_phasing_accuracy
#> 4   sample_ids_unique   TRUE  1 sample identifier(s)
#> 5  variant_ids_unique   TRUE 3 variant identifier(s)
#> 6 input_hashes_sha256   TRUE        2 input hash(es)
#> 7       quality_gates   TRUE       5 quality gate(s)
#> 8   validation_status   TRUE                  passed
```
