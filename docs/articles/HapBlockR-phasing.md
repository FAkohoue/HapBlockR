# Statistical Phasing with External Beagle 5.x

## Scientific distinction

An unphased dosage states how many alternate alleles an individual
carries. Phasing adds the chromosome assignment needed to identify
alleles transmitted together. Use validated phase for analyses of
gametic haplotypes, cis relationships or recombination; use dosage when
allele count is the relevant representation.

HapBlockR integrates Beagle 5.x as an external scientific dependency.
The user supplies the official JAR explicitly, preserving control of its
version, licence and execution.

## Configure Beagle

Download Beagle and the corresponding source from the official site,
review the GNU General Public License terms, and verify the downloaded
file:

<https://faculty.washington.edu/browning/beagle/beagle.html>

Configure one explicit path:

``` r
options(HapBlockR.beagle_jar = "/absolute/path/to/beagle.jar")

# Alternatively:
Sys.setenv(HAPBLOCKR_BEAGLE_JAR = "/absolute/path/to/beagle.jar")
```

The `beagle_jar` argument has highest priority, followed by the package
option, environment variable, and `beagle.jar` beside `out_prefix`.

Java 8 or later is required. Record `java -version` and the Beagle JAR
SHA-256 checksum with the analysis.

## Run and validate phasing

``` r
phased <- phase_with_beagle(
  input_vcf = "genotypes.vcf.gz",
  out_prefix = "results/genotypes_phased",
  nthreads = 2,
  seed = 42,
  ref_panel = "reference.vcf.gz",
  map_file = "genetic.map",
  min_genotype_concordance = 0.99,
  min_imputation_rate = 0.95,
  return_details = TRUE
)

phased$output_vcf
phased$quality_control
phased$provenance
```

Only one or two Beagle threads are accepted. HapBlockR verifies:

- the reported Beagle major version;
- input and output sample identity and order;
- chromosome, position, variant ID, REF, ALT, and variant order;
- dosage concordance at initially observed genotypes; and
- the imputation rate at initially missing genotypes.

Execution stops below user-defined thresholds. A structured provenance
file is written beside the phased VCF.

## Executable phased-VCF specification

This small deterministic example verifies the parser without requiring
Beagle during vignette construction.

``` r
vcf <- tempfile(fileext = ".vcf")
writeLines(
  c(
    "##fileformat=VCFv4.2",
    paste(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
      "FILTER", "INFO", "FORMAT", "P01", "P02",
      sep = "\t"
    ),
    paste("1", "100", "s1", "A", "G", ".", "PASS", ".",
          "GT", "0|1", "1|1", sep = "\t"),
    paste("1", "200", "s2", "C", "T", ".", "PASS", ".",
          "GT", "0|0", "1|0", sep = "\t")
  ),
  vcf
)

parsed <- read_phased_vcf(vcf)
parsed$hap1
#>    P01 P02
#> s1   0   1
#> s2   0   1
parsed$hap2
#>    P01 P02
#> s1   1   1
#> s2   0   0
parsed$dosage
#>    P01 P02
#> s1   1   2
#> s2   0   1
stopifnot(
  identical(parsed$dosage, parsed$hap1 + parsed$hap2),
  identical(parsed$sample_ids, c("P01", "P02"))
)
```

## End-to-end pipeline

`run_ldx_pipeline(phase = TRUE)` accepts only VCF or VCF.gz input. It
passes the explicitly configured Beagle path, aligns phased output by
immutable variant and allele identity, and can cache phased matrices in
a file-backed `bigmemory` backend.

``` r
result <- run_ldx_pipeline(
  geno_source = "genotypes.vcf.gz",
  out_dir = "results",
  out_blocks = "results/blocks.csv",
  out_diversity = "results/diversity.csv",
  out_hap_matrix = "results/haplotypes.csv",
  phase = TRUE,
  beagle_jar = "/absolute/path/to/beagle.jar",
  beagle_threads = 2,
  beagle_seed = 42,
  use_bigmemory = TRUE
)
```

## Interpretation and integration testing

Phased output remains model-dependent. Compare against a truth set when
possible and report switch error, genotype concordance, dosage accuracy,
imputation accuracy, and stratified performance by allele frequency.

The repository’s required Beagle integration workflow downloads a pinned
Beagle 5.x release from the official source, records its checksum and
Java version, and prevents the phasing tests from silently skipping when
that lane is configured.

``` r
packageVersion("HapBlockR")
#> [1] '0.3.12.9000'
sessionInfo()
#> R version 4.5.0 (2025-04-11 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/Bogota
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] HapBlockR_0.3.12.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.51        rlang_1.3.0       xfun_0.57        
#>  [5] otel_0.2.0        rrBLUP_4.6.3      textshaping_1.0.5 jsonlite_2.0.0   
#>  [9] data.table_1.18.4 htmltools_0.5.9   ragg_1.5.2        sass_0.4.10      
#> [13] rmarkdown_2.32    evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
#> [17] yaml_2.3.12       lifecycle_1.0.5   compiler_4.5.0    igraph_2.3.1     
#> [21] fs_2.1.0          pkgconfig_2.0.3   htmlwidgets_1.6.4 Rcpp_1.1.1-1.1   
#> [25] rstudioapi_0.18.0 systemfonts_1.3.2 digest_0.6.39     R6_2.6.1         
#> [29] parallel_4.5.0    magrittr_2.0.5    bslib_0.12.0      tools_4.5.0      
#> [33] pkgdown_2.2.0     cachem_1.1.0      desc_1.4.3
```
