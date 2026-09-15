# Large-Data Backends and Reproducible Performance Evidence

## Memory model

HapBlockR provides dense, subset-oriented and file-backed access
pathways. The appropriate memory model depends on the input format and
downstream analysis.

| Pathway | Access model | Best use and memory implication |
|----|----|----|
| GDS | Subset-oriented through SNPRelate | Downstream methods may still allocate dense matrices |
| PLINK BED | Operating-system memory mapping through BEDMatrix | Access pattern and downstream copies determine memory |
| `bigmemory` | File-backed matrix | Some algorithms materialise requested blocks |
| Text dosage or HapMap | Parsed text objects | Can require substantial in-memory representation |
| Ordinary VCF | Parser and optional GDS cache | Conversion, caching, and downstream objects have separate costs |

Peak memory depends on sample count, variant count, storage type, window
size, algorithm, missingness, threads, and copying by dependencies.
Measure the actual workflow on representative data.

## Executable correctness fixture

The small included data verify equality between one and two compiled
threads. This is a correctness test, not a performance claim.

``` r
data(ldx_geno)
x <- ldx_geno[seq_len(min(40L, nrow(ldx_geno))),
              seq_len(min(30L, ncol(ldx_geno))), drop = FALSE]

r_one <- compute_r2(x, n_threads = 1)
r_two <- compute_r2(x, n_threads = 2)

stopifnot(isTRUE(all.equal(r_one, r_two, tolerance = 1e-12)))
dim(r_one)
#> [1] 30 30
```

HapBlockR permits at most two compiled threads. Avoid nested parallel
plans that multiply process-level and compiled threads.

## File-backed workflows

Convert and validate data once, retain source and cache hashes, and
preserve sample and variant order.

``` r
library(SNPRelate)

SNPRelate::snpgdsVCF2GDS(
  vcf.fn = "genotypes.vcf.gz",
  out.fn = "genotypes.gds",
  method = "biallelic.only"
)

gds <- read_geno("genotypes.gds")
```

``` r
bed <- read_geno("genotypes.bed")
```

``` r
backed <- read_geno_bigmemory(
  "genotypes.vcf.gz",
  backingpath = "cache/genotypes"
)
```

HapBlockR cache manifests use SHA-256 source identity and parsing
parameters. A cache is reused only when that identity matches; otherwise
it is refused or explicitly rebuilt.

## Sparse and scalable LD options

Physical-distance restriction and community detection can reduce work:

``` r
blocks <- run_Big_LD_all_chr(
  geno_source = "genotypes.gds",
  CLQmode = "Leiden",
  max_bp_distance = 500000,
  subSegmSize = 1500,
  n_threads = 2
)
```

These options define the search space or segmentation method. Record
their biological rationale and compare them with a trusted scalar or
smaller exact analysis on a representative subset.

## Benchmark record

A publishable benchmark should retain:

``` r
benchmark_schema <- data.frame(
  field = c(
    "dataset_id", "input_sha256", "n_samples", "n_variants",
    "missing_rate", "method", "parameters", "threads", "R_version",
    "HapBlockR_version", "compiler", "hardware", "wall_seconds",
    "cpu_seconds", "peak_resident_mb", "replicate", "correctness_passed"
  ),
  required = TRUE
)
benchmark_schema
#>                 field required
#> 1          dataset_id     TRUE
#> 2        input_sha256     TRUE
#> 3           n_samples     TRUE
#> 4          n_variants     TRUE
#> 5        missing_rate     TRUE
#> 6              method     TRUE
#> 7          parameters     TRUE
#> 8             threads     TRUE
#> 9           R_version     TRUE
#> 10  HapBlockR_version     TRUE
#> 11           compiler     TRUE
#> 12           hardware     TRUE
#> 13       wall_seconds     TRUE
#> 14        cpu_seconds     TRUE
#> 15   peak_resident_mb     TRUE
#> 16          replicate     TRUE
#> 17 correctness_passed     TRUE
```

Report medians and variability from repeated runs, and distinguish
measured results from extrapolations. Promote a performance result after
numerical equality, sample identity and allele identity have passed.

## Scheduled integration

Large public or simulated fixtures belong in scheduled CI, with:

- fixed data checksums and package versions;
- hardware and compiler records;
- one- versus two-thread equality;
- a trusted scalar comparison;
- peak resident memory;
- performance non-regression thresholds; and
- retained benchmark artefacts.

The examples above remain non-executing only where they require external
files or optional backends; the vignette’s correctness and schema
examples are executed on every build.

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
