# From Genotypes to a Validated Breeding Decision

## Purpose

This executable vignette is a small specification of the primary
workflow. It uses included deterministic data, prints current result
schemas, and validates the final decision object. File-backed,
external-tool, and large-data routes are covered in their own
integration workflows.

## Load and inspect the data

``` r
data(ldx_geno)
data(ldx_snp_info)
data(ldx_blocks)
data(ldx_blues)

dim(ldx_geno)
#> [1] 120 230
head(ldx_snp_info)
#>      SNP CHR  POS REF ALT
#> 1 rs1001   1 1000   G   C
#> 2 rs1002   1 2192   C   G
#> 3 rs1003   1 3253   T   A
#> 4 rs1004   1 4132   A   G
#> 5 rs1005   1 5188   G   A
#> 6 rs1006   1 6314   C   T
head(ldx_blocks)
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp n_snps
#> 1     1  25     rs1001   rs1025     1000  25027   1     24028     25
#> 2    31  50     rs1031   rs1050    81064  99022   1     17959     20
#> 3    56  80     rs1056   rs1080   155368 179371   1     24004     25
#> 4    81 110     rs2001   rs2030     1000  30023   2     29024     30
#> 5   116 135     rs2036   rs2055    86236 105290   2     19055     20
#> 6   141 160     rs2061   rs2080   161515 180473   2     18959     20
```

The genotype matrix contains individuals in rows and variants in
columns. Row and column identifiers are part of the scientific identity
and must be unique and stable.

``` r
stopifnot(
  !is.null(rownames(ldx_geno)),
  !is.null(colnames(ldx_geno)),
  !anyDuplicated(rownames(ldx_geno)),
  !anyDuplicated(colnames(ldx_geno)),
  all(ldx_snp_info$SNP %in% colnames(ldx_geno))
)
```

For an external file, use
[`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
and retain its import and cache manifest. Memory behaviour depends on
the format and selected backend.

``` r
imported <- read_geno(
  "genotypes.vcf.gz",
  multiallelic = "split",
  verbose = TRUE
)
```

## Haplotype prediction

``` r
prediction <- run_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info = ldx_snp_info,
  blocks = ldx_blocks,
  blues = setNames(ldx_blues$YLD, ldx_blues$id),
  seed = 42,
  min_reliability = 0.30,
  verbose = FALSE
)

names(prediction)
#>  [1] "blocks"              "diversity"           "hap_matrix"         
#>  [4] "haplotypes"          "snp_info_filtered"   "n_blocks"           
#>  [7] "n_hap_columns"       "n_traits"            "traits"             
#> [10] "solver_used"         "include_dominance"   "gebv"               
#> [13] "gebv_uncertainty"    "dominance_deviation" "total_genetic_value"
#> [16] "snp_effects"         "local_gebv"          "block_importance"   
#> [19] "G"                   "G_dominance"         "n_train"            
#> [22] "n_predict"           "target_provenance"
head(prediction$block_importance)
#>                             block_id CHR start_bp end_bp n_snps var_local_gebv
#> 1                 block_1_1000_25027   1     1000  25027     25   2.434197e-02
#> singleton_1_75860  singleton_1_75860   1    75860  75860      1   1.614996e-04
#> singleton_1_77008  singleton_1_77008   1    77008  77008      1   8.681093e-05
#> singleton_1_77873  singleton_1_77873   1    77873  77873      1   4.061862e-03
#> singleton_1_79041  singleton_1_79041   1    79041  79041      1   3.675717e-04
#> singleton_1_80183  singleton_1_80183   1    80183  80183      1   2.677999e-05
#>                   mean_local_gebv_se singleton   var_scaled important
#> 1                                 NA     FALSE 0.1059798127     FALSE
#> singleton_1_75860                 NA      TRUE 0.0007031353     FALSE
#> singleton_1_77008                 NA      TRUE 0.0003779565     FALSE
#> singleton_1_77873                 NA      TRUE 0.0176844912     FALSE
#> singleton_1_79041                 NA      TRUE 0.0016003297     FALSE
#> singleton_1_80183                 NA      TRUE 0.0001165944     FALSE
head(prediction$gebv_uncertainty)
#>       id         gebv PEV reliability min_reliability recommendable
#> 1 ind001 -0.021344120  NA          NA             0.3         FALSE
#> 2 ind002  0.075066051  NA          NA             0.3         FALSE
#> 3 ind003 -0.052453052  NA          NA             0.3         FALSE
#> 4 ind004  0.005175617  NA          NA             0.3         FALSE
#> 5 ind005  0.230417687  NA          NA             0.3         FALSE
#> 6 ind006  0.127726805  NA          NA             0.3         FALSE
```

Reliability is `1 - PEV / Vg` where the fitted method provides
prediction error variance. When a model does not estimate reliability,
the result records it as unavailable so the breeder can apply the
programme’s evidence policy.

## Leakage-aware validation

``` r
cv <- cv_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info = ldx_snp_info,
  blocks = ldx_blocks,
  blues = setNames(ldx_blues$YLD, ldx_blues$id),
  k = 3,
  seed = 42,
  validation = "random",
  verbose = FALSE
)

cv$pa_pooled
#>   trait rep   n         PA     RMSE       MAE        bias calibration_slope
#> 1 trait   1 120 -0.0842114 1.100226 0.8970519 -0.02631208         -0.214659
head(cv$gebv_all)
#>       id trait rep fold observed       gebv breeding_value fitted_mean status
#> 1 ind001 trait   1    1  -0.5175 0.01652875  -1.998683e-10  0.01652875     ok
#> 2 ind004 trait   1    1  -1.1162 0.01652875   8.750502e-10  0.01652875     ok
#> 3 ind006 trait   1    1   0.9307 0.01652875   2.190268e-10  0.01652875     ok
#> 4 ind012 trait   1    1  -0.8182 0.01652875  -1.079175e-09  0.01652875     ok
#> 5 ind014 trait   1    1  -0.8170 0.01652875   1.054357e-10  0.01652875     ok
#> 6 ind025 trait   1    1   1.0970 0.01652875   3.842155e-10  0.01652875     ok
#>   error
#> 1  <NA>
#> 2  <NA>
#> 3  <NA>
#> 4  <NA>
#> 5  <NA>
#> 6  <NA>
validate(cv)
#>                 check passed                    detail
#> 1     required_fields   TRUE                          
#> 2      schema_version   TRUE                     1.0.0
#> 3              method   TRUE   cv_haplotype_prediction
#> 4   sample_ids_unique   TRUE  120 sample identifier(s)
#> 5  variant_ids_unique   TRUE 230 variant identifier(s)
#> 6 input_hashes_sha256   TRUE          3 input hash(es)
#> 7       quality_gates   TRUE         3 quality gate(s)
#> 8   validation_status   TRUE                    passed
```

Random folds estimate interpolation within a reasonably exchangeable
population and are used for this demonstration. Grouped validation
estimates transfer across families, populations or sites, while forward
validation estimates prediction into later cycles, years or
environments.

``` r
cv_grouped <- cv_haplotype_prediction(
  geno_matrix, snp_info, blocks, blues,
  k = 5,
  validation = "grouped",
  groups = family_id,
  seed = 42
)
```

## Parent selection

``` r
top_blocks <- select_top_blocks(
  prediction$block_importance,
  n = min(10L, nrow(prediction$block_importance))
)

value_matrix <- prediction$local_gebv[
  , top_blocks$block_id, drop = FALSE
]

selection <- select_parents_ga(
  value_matrix = value_matrix,
  n_founders = 8,
  seed = 42,
  n_reps = 3,
  verbose = FALSE
)

selection$selected
#> [1] "ind029" "ind032" "ind035" "ind044" "ind049" "ind083" "ind095" "ind104"
selection$stability
#> $n_reps
#> [1] 3
#> 
#> $fitness_values
#> [1] 2.877236 2.884246 2.845721
#> 
#> $fitness_range
#> [1] 2.845721 2.884246
#> 
#> $best_rep
#> [1] 2
#> 
#> $converged
#> [1] TRUE TRUE TRUE
#> 
#> $feasible
#> [1] TRUE TRUE TRUE
#> 
#> $selection_freq
#>    ind032    ind035    ind049    ind095    ind104    ind029    ind044    ind083 
#> 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.6666667 0.6666667 0.6666667 
#>    ind063    ind115    ind111 
#> 0.3333333 0.3333333 0.3333333 
#> 
#> $mean_relationship
#> [1] NA NA NA
#> 
#> $mean_merit
#> [1] NA NA NA
#> 
#> $seeds
#> [1] 42 43 44
#> 
#> $termination_reason
#> [1] "fitness_plateau" "fitness_plateau" "fitness_plateau"
validate(selection)
#>                 check passed                   detail
#> 1     required_fields   TRUE                         
#> 2      schema_version   TRUE                    1.0.0
#> 3              method   TRUE        select_parents_ga
#> 4   sample_ids_unique   TRUE 120 sample identifier(s)
#> 5  variant_ids_unique   TRUE 10 variant identifier(s)
#> 6 input_hashes_sha256   TRUE         4 input hash(es)
#> 7       quality_gates   TRUE        3 quality gate(s)
#> 8   validation_status   TRUE                   passed
```

The validation status, exact founder count, warnings, exclusions, seed,
input hashes, and stability evidence should be retained with the
decision.

## Result contract

``` r
names(selection$result_contract)
#>  [1] "schema_version"    "method"            "created_at_utc"   
#>  [4] "call"              "parameters"        "random_seed"      
#>  [7] "identifiers"       "input_hashes"      "transformations"  
#> [10] "software"          "quality_gates"     "warnings"         
#> [13] "fallbacks"         "excluded_records"  "decision_table"   
#> [16] "uncertainty"       "validation_status"
summary(selection)
#> $method
#> [1] "select_parents_ga"
#> 
#> $schema_version
#> [1] "1.0.0"
#> 
#> $validation_status
#> [1] "passed"
#> 
#> $quality_gates
#>                         gate passed detail
#> 1        exact_founder_count   TRUE   <NA>
#> 2 hard_constraints_satisfied   TRUE   <NA>
#> 3   all_repetitions_feasible   TRUE   <NA>
#> 
#> $warnings
#> character(0)
#> 
#> $fallbacks
#> character(0)
#> 
#> $n_decisions
#> [1] 8
#> 
#> $n_uncertainty_records
#> [1] 3
as.data.frame(selection)
#>       id selection_frequency
#> 1 ind029           0.6666667
#> 2 ind032           1.0000000
#> 3 ind035           1.0000000
#> 4 ind044           0.6666667
#> 5 ind049           1.0000000
#> 6 ind083           0.6666667
#> 7 ind095           1.0000000
#> 8 ind104           1.0000000
```

Promote a result after
[`validate()`](https://FAkohoue.github.io/HapBlockR/reference/validate.md)
reports passing gates. Retain fallbacks, excluded candidates and
uncertainty with the breeding decision.

## Session information

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
#>  [1] crayon_1.5.3      cli_3.6.6         knitr_1.51        rlang_1.3.0      
#>  [5] xfun_0.57         otel_0.2.0        rrBLUP_4.6.3      textshaping_1.0.5
#>  [9] jsonlite_2.0.0    data.table_1.18.4 htmltools_0.5.9   ragg_1.5.2       
#> [13] sass_0.4.10       rmarkdown_2.32    evaluate_1.0.5    jquerylib_0.1.4  
#> [17] fastmap_1.2.0     foreach_1.5.2     yaml_2.3.12       lifecycle_1.0.5  
#> [21] compiler_4.5.0    codetools_0.2-20  igraph_2.3.1      fs_2.1.0         
#> [25] pkgconfig_2.0.3   htmlwidgets_1.6.4 Rcpp_1.1.1-1.1    rstudioapi_0.18.0
#> [29] systemfonts_1.3.2 digest_0.6.39     R6_2.6.1          parallel_4.5.0   
#> [33] magrittr_2.0.5    bslib_0.12.0      tools_4.5.0       iterators_1.0.14 
#> [37] GA_3.2.5          pkgdown_2.2.0     cachem_1.1.0      desc_1.4.3
```
