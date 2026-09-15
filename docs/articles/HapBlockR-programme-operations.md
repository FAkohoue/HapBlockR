# Programme-Level Models, Feasibility, and Data Exchange

## Externally analysed targets

HapBlockR does not analyse raw plot records or fit replicate, block, row
and column effects. Supply genotype-level estimates from the approved
external analysis, with uncertainty. Units may be retained but are
optional.

``` r
external_values <- data.frame(
  id = paste0("G", 1:4),
  trait = "yield",
  value = c(5.2, 4.8, 5.0, 4.6),
  SE = c(0.20, 0.30, 0.25, 0.35)
)

targets <- prepare_breeding_targets(
  external_values,
  input_type = "BLUE",
  se_col = "SE"
)
targets$targets[
  , c("id", "value", "model_value", "precision_weight")
]
#>   id value model_value precision_weight
#> 1 G1   5.2         5.2        1.6590798
#> 2 G2   4.8         4.8        0.7373688
#> 3 G3   5.0         5.0        1.0618111
#> 4 G4   4.6         4.6        0.5417403
```

Accepted types are `adjusted_mean`, `BLUE`, `BLUP_identity`, `PBLUP`,
`BV`, `GCA` and `TGV`;
[`breeding_target_types()`](https://FAkohoue.github.io/HapBlockR/reference/breeding_target_types.md)
gives their complete definitions. Random predictions are deregressed
from reliability or prediction error variance. Pedigree BLUPs require
the named numerator relationship matrix `A`. External genomic
predictions and external selection-index values are rejected because
HapBlockR estimates genomic effects and the final index internally.

## Multi-trait index

Improvement directions are model inputs. Trait units are optional
metadata.

``` r
values <- cbind(
  yield = c(5.2, 4.8, 5.0, 4.6),
  disease = c(2.0, 1.0, 2.5, 0.8)
)
rownames(values) <- paste0("G", 1:4)
G_trait <- matrix(
  c(1.0, -0.2, -0.2, 0.7), 2, 2,
  dimnames = list(colnames(values), colnames(values))
)
P_trait <- matrix(
  c(1.5, -0.1, -0.1, 1.1), 2, 2,
  dimnames = list(colnames(values), colnames(values))
)

index <- build_selection_index(
  values,
  genetic_cov = G_trait,
  phenotypic_cov = P_trait,
  economic_weights = c(yield = 1, disease = 0.5),
  directions = c(yield = "increase", disease = "decrease"),
  units = NULL,
  method = "smith_hazel"
)
index$scores
#>   id selection_index rank
#> 2 G2        2.944512    1
#> 4 G4        2.890854    2
#> 1 G1        2.790244    3
#> 3 G3        2.431402    4
index$coefficients
#>           trait unit direction coefficient linear_weight expected_response
#> yield     yield <NA>  increase   0.7042683            NA         0.7857940
#> disease disease <NA>  decrease   0.4359756            NA        -0.4428416
#>         realised_response_sd_favourable realised_response_sd_original
#> yield                                NA                            NA
#> disease                              NA                            NA
validate(index)
#>                 check passed                  detail
#> 1     required_fields   TRUE                        
#> 2      schema_version   TRUE                   1.0.0
#> 3              method   TRUE   build_selection_index
#> 4   sample_ids_unique   TRUE  4 sample identifier(s)
#> 5  variant_ids_unique   TRUE 0 variant identifier(s)
#> 6 input_hashes_sha256   TRUE        4 input hash(es)
#> 7       quality_gates   TRUE       4 quality gate(s)
#> 8   validation_status   TRUE                  passed
```

[`build_selection_index()`](https://FAkohoue.github.io/HapBlockR/reference/build_selection_index.md)
also provides deterministic Pesek-Baker, Desired-Gain Selection Index
(DGSI) and Quadratic Genomic Selection Index (QGSI) methods. DGSI and
QGSI call DesiredGainR. DGSI runs independent optimisation replicates
and automatically returns the replicate chosen by its declared holdout
or independent-validation rule; breeders review the compact stability
tables rather than choosing a run manually. QGSI requires an explicit
symmetric quadratic-weight matrix and returns candidate-specific linear,
squared and cross-product contributions.

For every method, `directions` defines favourable orientation.
`economic_weights` and `desired_gains` are non-negative magnitudes in
that oriented space; negative objectives are rejected. DGSI reports both
its model-expected transmitted genetic response in original trait units
and its realised selected-set differential in candidate SD units. DGSI
and QGSI use the selection intensity for the number of candidates
actually selected. QGSI reports model-expected gains using the complete
QGSI variance and retains its linear and quadratic weights separately
because no single global QGSI coefficient vector exists. HapBlockR
converts DesiredGainR’s reported response from its analysis scale back
to the original trait units. Additional QGSI engine controls, such as
explicit `Gamma`, relationship information or trait scaling, can be
supplied through `qgsi_control`. All methods expose a common
coefficient-table schema.

DGSI interprets `desired_gains` in candidate standard deviations even
when `scale_traits = FALSE`; Pesek-Baker uses original trait units. For
example, a DGSI target of 0.5 requests a selected-set shift of half a
candidate standard deviation in the favourable direction. This target is
distinct from the model-expected response transmitted to the next
generation.

Both control lists require exact DesiredGainR argument names. Supply
reference data with the original-unit trait columns and an `id` column.
For QGSI, `Gamma` also uses original trait units, whereas the linear and
quadratic weights refer to the oriented, optionally scaled analysis
space. `genetic_cov` and `phenotypic_cov` remain upstream context for
QGSI; they do not replace `Gamma`.

For a DGSI result, use `coefficients_original_units` when combining
original-unit marker, haplotype or block effects. Candidate scores can
be reconstructed as
`trait_values %*% coefficients_original_units + score_intercept`. The
existing `coefficients$coefficient` field retains DesiredGainR’s
analysis-scale coefficients. Use `engine_result` when passing the fitted
index to DesiredGainR’s comparison tools.

Use
[`fit_multitrait_gblup()`](https://FAkohoue.github.io/HapBlockR/reference/fit_multitrait_gblup.md)
when genetic and residual covariance must be estimated jointly. A
supplied full sampling covariance is the complete record-error
covariance by default; request
`sampling_covariance_mode = "sampling_plus_residual"` only when an
additional residual nugget is intended. Sampling-covariance keys and
diagonal precision are checked before fitting. Use
[`fit_gxe_gblup()`](https://FAkohoue.github.io/HapBlockR/reference/fit_gxe_gblup.md)
with an environment kernel for reaction-norm predictions.
Environment-specific predictions support environment-specific parent
selection. The across-environment output is the internally estimated
genomic main effect; HapBlockR does not combine environment predictions
with breeder-specified environment weights. Both models report REML
likelihoods with residual degrees of freedom; compare those likelihoods
only between models with the same fixed-effect design.

## Operational screening

Parent 1 is the female and parent 2 the male in directed-cross
screening.

``` r
status <- data.frame(
  id = c("G1", "G2", "G3", "G4"),
  female_allowed = c(TRUE, TRUE, TRUE, FALSE),
  male_allowed = c(TRUE, TRUE, FALSE, TRUE),
  fertile = TRUE,
  flowering_start = c(5, 6, 7, 7),
  flowering_end = c(7, 8, 9, 9),
  heterotic_group = c("H1", "H2", "H1", "H2"),
  subpopulation = c("S1", "S1", "S2", "S2"),
  female_capacity = c(30, 20, 20, 0),
  male_capacity = c(20, 30, 0, 20),
  total_capacity = c(30, 30, 20, 20),
  seed_available = c(30, 20, 20, 0),
  pollen_available = c(20, 30, 0, 20)
)

screen_candidate_crosses(
  rbind(c("G1", "G2"), c("G1", "G3"), c("G4", "G1")),
  status,
  require_different_heterotic_groups = TRUE
)
#>   female male feasible flowering_overlap
#> 1     G1   G2     TRUE              TRUE
#> 2     G1   G3    FALSE              TRUE
#> 3     G4   G1    FALSE              TRUE
#>                                              reasons
#> 1                                                   
#> 2 male_role_forbidden;same_heterotic_group_forbidden
#> 3                              female_role_forbidden
```

## Feasibility certificate

``` r
plan <- data.frame(
  female = c("G1", "G3"),
  male = c("G2", "G4"),
  family_size = c(20, 20),
  period = "P1"
)
certificate <- certify_mating_plan(
  plan,
  status,
  min_family_size = 20,
  required_pairs = data.frame(parent1 = "G1", parent2 = "G2"),
  subpopulation_quotas = data.frame(
    subpopulation = c("S1", "S2"),
    min_contribution = c(40, 40)
  ),
  require_different_heterotic_groups = TRUE
)
certificate$certificate
#>   feasible n_crosses total_families n_violations n_binding_constraints
#> 1     TRUE         2             40            0                    12
certificate$binding_constraints
#>          constraint subject observed limit slack
#> 1   female_capacity      G3       20    20     0
#> 2   female_capacity      G4        0     0     0
#> 3     male_capacity      G3        0     0     0
#> 4     male_capacity      G4       20    20     0
#> 5    total_capacity      G3       20    20     0
#> 6    total_capacity      G4       20    20     0
#> 7    seed_available      G3       20    20     0
#> 8    seed_available      G4        0     0     0
#> 9  pollen_available      G3        0     0     0
#> 10 pollen_available      G4       20    20     0
#> 11 min_contribution      S1       40    40     0
#> 12 min_contribution      S2       40    40     0
validate(certificate)
#>                 check passed                  detail
#> 1     required_fields   TRUE                        
#> 2      schema_version   TRUE                   1.0.0
#> 3              method   TRUE     certify_mating_plan
#> 4   sample_ids_unique   TRUE  4 sample identifier(s)
#> 5  variant_ids_unique   TRUE 0 variant identifier(s)
#> 6 input_hashes_sha256   TRUE        5 input hash(es)
#> 7       quality_gates   TRUE       6 quality gate(s)
#> 8   validation_status   TRUE                  passed
```

The certificate must pass before the plan is signed. Keep the violations
and binding constraints even when the plan is feasible.

## Standards-oriented exchange

Remote breeding databases should remain outside the modelling functions.
HapBlockR exports stable tables and mappings that a programme
integration layer can send to BrAPI or align with MIAPPE.

``` r
bundle <- build_breeding_exchange(
  index,
  breeding_program_id = "BP1",
  study_id = "STUDY1",
  trial_id = "TRIAL1",
  environment_id = "ENV1",
  trait_id = "SELECTION_INDEX",
  unit = "index unit",
  ontology_id = "CO_321:0000012"
)
bundle$recommendations
#>   recommendation_id breeding_program_id study_id trial_id environment_id
#> 1         REC000001                 BP1   STUDY1   TRIAL1           ENV1
#> 2         REC000002                 BP1   STUDY1   TRIAL1           ENV1
#> 3         REC000003                 BP1   STUDY1   TRIAL1           ENV1
#> 4         REC000004                 BP1   STUDY1   TRIAL1           ENV1
#>          trait_id    ontology_id       unit                method
#> 1 SELECTION_INDEX CO_321:0000012 index unit build_selection_index
#> 2 SELECTION_INDEX CO_321:0000012 index unit build_selection_index
#> 3 SELECTION_INDEX CO_321:0000012 index unit build_selection_index
#> 4 SELECTION_INDEX CO_321:0000012 index unit build_selection_index
#>   schema_version entity_type entity_id germplasm_id female_germplasm_id
#> 1          1.0.0   germplasm        G2           G2                <NA>
#> 2          1.0.0   germplasm        G4           G4                <NA>
#> 3          1.0.0   germplasm        G1           G1                <NA>
#> 4          1.0.0   germplasm        G3           G3                <NA>
#>   male_germplasm_id rank    score recommendable reason
#> 1              <NA>    1 2.944512          TRUE       
#> 2              <NA>    2 2.890854          TRUE       
#> 3              <NA>    3 2.790244          TRUE       
#> 4              <NA>    4 2.431402          TRUE
bundle$field_mapping
#>        canonical_field                BrAPI_field             MIAPPE_field
#> 1  breeding_program_id        breedingProgramDbId investigation identifier
#> 2             study_id                  studyDbId         study identifier
#> 3             trial_id                  trialDbId         study identifier
#> 4       environment_id        observationUnitDbId      observation unit ID
#> 5             trait_id    observationVariableDbId        observed variable
#> 6         germplasm_id              germplasmDbId   biological material ID
#> 7  female_germplasm_id parent1DbId/additionalInfo        crossing metadata
#> 8    male_germplasm_id parent2DbId/additionalInfo        crossing metadata
#> 9                score                      value           observed value
#> 10                unit                scale/units                     unit
#> 11         ontology_id               ontologyDbId   trait accession number

exchange_path <- tempfile("hapblockr-exchange-")
write_breeding_exchange(bundle, exchange_path)
verified <- read_breeding_exchange(exchange_path)
stopifnot(
  identical(
    verified$recommendations$recommendation_id,
    bundle$recommendations$recommendation_id
  )
)
```

## Sensitivity

Refit the intended threshold, leave-one-environment, or
leave-one-population scenarios explicitly, then compare their result
objects:

``` r
stability <- assess_decision_stability(
  list(
    baseline = baseline_result,
    leave_E1_out = leave_E1_result,
    leave_E2_out = leave_E2_result
  ),
  top_n = 10,
  minimum_selection_frequency = 0.67
)
```

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
