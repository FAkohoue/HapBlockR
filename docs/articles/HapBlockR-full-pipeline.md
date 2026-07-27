# The Full Pipeline: Genotypes to a Certified Mating Plan

## Purpose

Every other vignette in this package is a deep dive into one topic: LD
metrics, phasing, large-scale backends, the parent-selection toolkit,
programme-level operations. This one instead runs the whole thing, start
to finish, as one continuous pipeline on HapBlockR’s own small example
dataset (`ldx_geno`, `ldx_snp_info`, `ldx_blocks`, `ldx_blues`) so the
relationships between stages – what feeds what, what each stage’s output
is actually used for next – are visible in one place rather than left to
be inferred across seven documents.

The stage structure mirrors a real production pipeline built on this
package (a CIAT Peru rice parent-selection script, distributed with the
package source under `scripts/`), organised here into four phases:
evidence, shortlist, mate allocation, sign-off. Each stage is kept short
and points to the vignette or Breeder’s Guide section that covers it in
depth – this is the connective tissue, not a replacement for those.

``` r
data(ldx_geno)
data(ldx_snp_info)
data(ldx_blocks)
data(ldx_blues)

dim(ldx_geno)               # 120 individuals x 230 SNPs
#> [1] 120 230
table(ldx_snp_info$CHR)     # three chromosomes
#> 
#>  1  2  3 
#> 80 80 70
head(ldx_blues)
#>       id     YLD     RES
#> 1 ind001 -0.5175  0.6771
#> 2 ind002  0.7635  1.3764
#> 3 ind003 -1.3093 -0.9946
#> 4 ind004 -1.1162 -1.4089
#> 5 ind005  1.1343 -0.5120
#> 6 ind006  0.9307 -0.4573
```

## Phase 1: Evidence (stages 1-4)

### Stage 1 – Prepare the targets

A real programme’s phenotype rarely arrives as a bare numeric column: it
is the output of a field-trial analysis or genetic evaluation run
outside HapBlockR, and that analysis’s estimand (an adjusted mean? a
BLUP shrunk toward zero? a breeding value?) changes how it should be
used downstream.
[`prepare_breeding_targets()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_breeding_targets.md)
makes that estimand explicit and, for shrunk random-effect predictions,
deregresses them before they re-enter marker or haplotype models. See
“Externally analysed targets” in the *programme operations* vignette and
[`breeding_target_types()`](https://FAkohoue.github.io/HapBlockR/reference/breeding_target_types.md)
for the complete list of accepted input types.

The toy `ldx_blues` already stands in for a completed external analysis
– one row per individual, a single trait (`YLD`) – so here it is used
directly as a named vector, the simplest of the four phenotype formats
[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
accepts:

``` r
blues_vec <- setNames(ldx_blues$YLD, ldx_blues$id)
length(blues_vec)
#> [1] 120
```

### Stage 2 – Genotypes and LD blocks

[`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
is the entry point for real files (VCF, HapMap, GDS, PLINK BED, dosage
tables), with an explicit multiallelic policy and immutable
sample/variant identity preserved throughout. `ldx_geno` is already an
in-memory dosage matrix, so this vignette skips straight to LD block
detection;
[`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md)
(single chromosome) and
[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
(genome-wide) are the tools that would normally produce `ldx_blocks`
from `ldx_geno`/`ldx_snp_info`:

``` r
head(ldx_blocks)
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp n_snps
#> 1     1  25     rs1001   rs1025     1000  25027   1     24028     25
#> 2    31  50     rs1031   rs1050    81064  99022   1     17959     20
#> 3    56  80     rs1056   rs1080   155368 179371   1     24004     25
#> 4    81 110     rs2001   rs2030     1000  30023   2     29024     30
#> 5   116 135     rs2036   rs2055    86236 105290   2     19055     20
#> 6   141 160     rs2061   rs2080   161515 180473   2     18959     20
identical(sort(unique(ldx_blocks$CHR)), sort(unique(ldx_snp_info$CHR)))
#> [1] TRUE
```

See the *intro* vignette for genome-wide detection and block
summarisation, and the *large-scale* vignette for file-backed and
whole-genome-sequencing routes.

### Stage 3 – Fit the prediction model

[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
extracts haplotype alleles per block, builds the shared haplotype-block
GRM, fits marker effects (GBLUP by default), and backsolves each
individual’s per-block local GEBV. The return value carries both the
whole-genome GEBV every shortlist tool below scores candidates on, and
`gebv_uncertainty` – prediction error variance, reliability, and the
`min_reliability` recommendation gate – so “how good is this GEBV,
really” has an answer before any decision is built on top of it.

``` r
pred <- run_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info    = ldx_snp_info,
  blocks      = ldx_blocks,
  blues       = blues_vec,
  min_snps    = 3L,
  seed        = 1,
  verbose     = FALSE
)

head(pred$gebv_uncertainty)
#>       id         gebv PEV reliability min_reliability recommendable
#> 1 ind001 -0.021344120  NA          NA             0.3         FALSE
#> 2 ind002  0.075066051  NA          NA             0.3         FALSE
#> 3 ind003 -0.052453052  NA          NA             0.3         FALSE
#> 4 ind004  0.005175617  NA          NA             0.3         FALSE
#> 5 ind005  0.230417687  NA          NA             0.3         FALSE
#> 6 ind006  0.127726805  NA          NA             0.3         FALSE
dim(pred$local_gebv)
#> [1] 120  39
dim(pred$G)
#> [1] 120 120
```

### Stage 4 – Validate with a deployment-matched design

A predictive-ability estimate is only trustworthy if the
cross-validation design matches how the model will actually be deployed.
`"random"` folds estimate interpolation within an exchangeable
population (used here for speed); `"grouped"` folds (holding whole
families or sites out together) estimate transfer across families or
locations; `"forward"` folds estimate prediction into a not-yet-observed
cycle. See “Leakage-aware validation” in the *workflow* vignette for the
other two designs.

``` r
cv <- cv_haplotype_prediction(
  geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
  blues = blues_vec, k = 3L, validation = "random", seed = 1, verbose = FALSE
)
cv$pa_mean
#>   trait          PA     RMSE PA_sd RMSE_sd
#> 1 trait -0.02368911 1.032983    NA      NA
```

## Phase 2: Shortlist (stages 5-8)

### Stage 5 – Select the blocks that matter

[`select_top_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
reduces `pred$block_importance` to the smallest top-ranked set that
explains most of the variance in local GEBV – here, the smallest set
explaining 90%. The resulting columns of `pred$local_gebv` become
`value_matrix`, the shared input to every founder-selection tool below.

``` r
top_blocks <- select_top_blocks(pred$block_importance, perc_of_total_var = 0.90)
nrow(top_blocks)
#> [1] 6

value_matrix <- pred$local_gebv[, top_blocks$block_id, drop = FALSE]
dim(value_matrix)
#> [1] 120   6
```

### Stage 6 – Truncation baseline

[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
on `pred$gebv` is the reference point every alternative shortlist tool
below should be judged against.

``` r
n_founders <- 15L
ts_sel <- truncation_selection(pred$gebv, n_founders = n_founders)
ts_sel$selected
#>  [1] "ind108" "ind025" "ind088" "ind111" "ind083" "ind104" "ind005" "ind098"
#>  [9] "ind033" "ind032" "ind095" "ind064" "ind062" "ind068" "ind044"
```

### Stage 7 – An alternative shortlist tool

Four alternatives are available depending on which constraint actually
matters for the programme – pure block coverage, coverage plus merit,
guaranteed family or genetic-cluster representation, or diversity itself
as the objective. This vignette runs the GA (pure block-coverage) and
family-quota tools; the *breeding decisions* vignette compares all four
against each other, including a forward-simulation check of realised
genetic gain, in depth.

``` r
ga_sel <- select_parents_ga(
  value_matrix  = value_matrix,
  n_founders    = n_founders,
  strategy      = "no_selfing",
  block_weights = top_blocks$var_scaled,
  popSize = 30L, maxiter = 20L, run = 10L,
  seed = 1, n_reps = 3L, verbose = FALSE
)
ga_sel$selected
#>  [1] "ind029" "ind035" "ind039" "ind041" "ind044" "ind049" "ind073" "ind083"
#>  [9] "ind086" "ind090" "ind091" "ind099" "ind104" "ind108" "ind119"
ga_sel$converged
#> [1] FALSE
```

``` r
# ldx_geno/ldx_blues carry no real pedigree; a synthetic family label (12
# families of 10) stands in here so family_select_mode has something to
# group on. See "Family or genetic-cluster quota selection" in the Breeder's
# Guide for all four family_select_mode variants and the shrinkage-corrected
# family ranking behind this tool.
ids    <- rownames(ldx_geno)
family <- setNames(paste0("Fam", ceiling(seq_along(ids) / 10)), ids)

fam_sel <- select_parents_by_family(
  score = pred$gebv, family = family, n_families = 5L, n_per_family = 3L,
  family_select_mode = "count", variance_method = "anova",
  use_family_relationship = TRUE, G = pred$G, verbose = FALSE
)
fam_sel$selected
#>  [1] "ind005" "ind010" "ind006" "ind015" "ind016" "ind011" "ind025" "ind026"
#>  [9] "ind021" "ind033" "ind032" "ind039" "ind044" "ind049" "ind045"
fam_sel$family_ranking[, c("family", "rank_score", "rank")]
#>    family rank_score rank
#> 1    Fam1  0.1723079    1
#> 2    Fam2  0.1723079    2
#> 3    Fam3  0.1723079    3
#> 4    Fam4  0.1723079    4
#> 5    Fam5  0.1723079    5
#> 6    Fam6  0.1723079    6
#> 7    Fam7  0.1723079    7
#> 8    Fam8  0.1723079    8
#> 9    Fam9  0.1723079    9
#> 10  Fam10  0.1723079   10
#> 11  Fam11  0.1723079   11
#> 12  Fam12  0.1723079   12
```

### Stage 8 – See the realised gain-vs-diversity tradeoff

[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)
sweeps
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
`coancestry_weight` across a grid and Pareto-filters the resulting
(merit, relatedness) points into an empirical frontier, so
`target_degree` in stage 11 can be set from an observed tradeoff rather
than a guess.

``` r
pareto_res <- select_parents_pareto(
  value_matrix = value_matrix, n_founders = n_founders,
  strategy = "no_selfing", block_weights = top_blocks$var_scaled,
  G = pred$G, coancestry_weights = c(0, 1, 2), merit = pred$gebv,
  popSize = 30L, maxiter = 20L, run = 10L, n_reps = 2L,
  seed = 1, verbose = FALSE
)
pareto_res$frontier[, c("coancestry_weight", "mean_merit",
                        "mean_relationship", "pareto_optimal")]
#>   coancestry_weight   mean_merit mean_relationship pareto_optimal
#> 1                 0  0.092417701       -0.00335104           TRUE
#> 2                 2 -0.003483062       -0.01588996           TRUE
#> 3                 1  0.055180348       -0.01020489           TRUE
```

## Phase 3: Mate allocation (stages 9-13)

### Stage 9 – Screen candidate crosses for operational feasibility

A cross that cannot physically be made should never reach a merit-based
ranking.
[`screen_candidate_crosses()`](https://FAkohoue.github.io/HapBlockR/reference/screen_candidate_crosses.md)
checks role, fertility, flowering overlap, heterotic group, quarantine
compatibility, forbidden pairs and reciprocal effects up front.
`ldx_geno` carries no real breeding metadata, so a minimal illustrative
status table stands in here for the GA shortlist; see “Operational
screening” in the *programme operations* vignette for the complete set
of checks.

``` r
ga_status <- data.frame(
  id = ga_sel$selected,
  female_allowed = TRUE,
  male_allowed = TRUE,
  fertile = TRUE,
  flowering_start = 5,
  flowering_end = 9
)
candidate_pairs <- t(utils::combn(ga_sel$selected, 2L))
screened <- screen_candidate_crosses(candidate_pairs, ga_status)
table(screened$feasible)
#> 
#> TRUE 
#>  105

feasible_pairs <- screened[screened$feasible, c("female", "male")]
```

### Stage 10 – Rank feasible crosses by usefulness criterion

[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
scores each feasible pair – only the pairs stage 9 actually cleared,
passed in via `cross_pairs` rather than regenerated from scratch – by
predicted mid-parent value plus predicted segregation variance at the
target blocks. `variance_model = "block_independent"` is used here since
only unphased local GEBV is available; see “Variance-model variants” in
the Breeder’s Guide for the phased, linked and SimpleMating-based
alternatives and what each requires.

``` r
uc_res <- usefulness_criterion(
  cross_pairs = feasible_pairs, gebv = pred$gebv,
  selected_proportion = 0.1, seed = 1,
  variance_model = "block_independent",
  block_importance = top_blocks, local_gebv = pred$local_gebv,
  verbose = FALSE
)
head(uc_res[, c("parent1", "parent2", "mid_parent_gebv", "UC", "rank")])
#>   parent1 parent2 mid_parent_gebv        UC rank
#> 1  ind073  ind083      0.20774221 0.9689435    1
#> 2  ind035  ind073      0.13159565 0.9166926    2
#> 3  ind083  ind099      0.09394776 0.9054839    3
#> 4  ind044  ind086      0.10501938 0.9008342    4
#> 5  ind044  ind083      0.20846077 0.8921764    5
#> 6  ind035  ind044      0.13231421 0.8916564    6
```

### Stage 11 – Solve for actual optimal contributions

[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
solves the real Optimal Contribution Selection problem (Meuwissen 1997):
how much each parent should contribute, and which specific matings
deliver that, under a relatedness ceiling set by `target_degree`
(informed by stage 8’s frontier). It is run here restricted to the stage
5-8 shortlist (`merit`/`G` subset to `ga_sel$selected`), not the full
120-individual population, so its output stays on the same candidate
pool as stage 10’s ranking and stage 13’s exact-validation comparison
below – OCS is solving “how to allocate this already-shortlisted set”,
not re-shortlisting from scratch. This needs the optional `optiSel` or
`SimpleMating` packages; where neither is installed, the usefulness
-criterion ranking from stage 10 stands in as the shortlist for stages
12-13 below.

``` r
shortlist <- ga_sel$selected
ocs_res <- select_parents_ocs(
  merit = pred$gebv[shortlist], G = pred$G[shortlist, shortlist],
  family = family[shortlist],
  engine = "optisel", n_crosses = 10L, max_contrib_per_parent = 3L,
  target_degree = 30, seed = 1, verbose = FALSE
)
ocs_res$mating_plan
#>    parent1 parent2 mean_relationship
#> 1   ind035  ind044        -0.1578575
#> 2   ind041  ind083        -0.1475842
#> 3   ind035  ind086        -0.1327449
#> 4   ind041  ind104        -0.1110568
#> 5   ind044  ind104        -0.1167642
#> 6   ind041  ind108        -0.1692724
#> 7   ind049  ind108        -0.1258960
#> 8   ind104  ind108        -0.1042079
#> 9   ind044  ind119        -0.1578575
#> 10  ind086  ind119        -0.1327449
ocs_res$ok
#> [1] TRUE
```

``` r
message("optiSel not installed -- using the stage 10 usefulness-criterion ",
        "ranking as the mating-plan candidate list instead.")
mating_plan <- uc_res[uc_res$rank <= 10L,
                      c("parent1", "parent2", "UC")]
names(mating_plan)[3] <- "criterion"
mating_plan
```

### Stage 12 – Certify the plan and check its sensitivity

[`certify_mating_plan()`](https://FAkohoue.github.io/HapBlockR/reference/certify_mating_plan.md)
validates the complete plan against family size, parent capacity,
required pairs, period-specific capacity and subpopulation quotas, and
returns every binding constraint alongside the certificate.
[`assess_decision_stability()`](https://FAkohoue.github.io/HapBlockR/reference/assess_decision_stability.md)
(not run here – see “Sensitivity” in the *programme operations*
vignette) then compares the certified plan against threshold, model or
population-scenario variants.

``` r
plan_pairs <- if (have_ocs) ocs_res$mating_plan else mating_plan
plan <- data.frame(
  female = plan_pairs$parent1,
  male = plan_pairs$parent2,
  family_size = 20
)
plan_status <- data.frame(
  id = unique(c(plan$female, plan$male)),
  female_allowed = TRUE, male_allowed = TRUE, fertile = TRUE,
  flowering_start = 5, flowering_end = 9
)
certificate <- certify_mating_plan(plan, plan_status, min_family_size = 20)
certificate$certificate
#>   feasible n_crosses total_families n_violations n_binding_constraints
#> 1     TRUE        10            200            0                     0
```

### Stage 13 – Validate against the true optimum

[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)
solves the same cross-selection problem exactly (binary integer linear
programming) on the culled candidate list and reports the heuristic
plan’s percentage gap below the true optimum – a small-scale sanity
check on stage 11, not a routine replacement for it. This needs the
optional `lpSolve` package.

``` r
exact_res <- validate_crosses_exact(
  data = uc_res, n_cross = min(10L, nrow(uc_res)), max_cross = 3L,
  heuristic_plan = plan_pairs[, c("parent1", "parent2")], verbose = FALSE
)
exact_res$exact_objective
#> [1] 8.874961
exact_res$gap_pct
#> [1] 32.95513
```

## Phase 4: Sign-off (stages 14-16)

### Stage 14 – Cross-check with a diversity-first alternative

[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
on the same relationship matrix, and
[`cluster_selection_groups()`](https://FAkohoue.github.io/HapBlockR/reference/cluster_selection_groups.md)
across every strategy computed above, show whether the merit-driven
selections concentrate on one or two genetic clusters or actually spread
across the population’s real structure.

``` r
core_res <- select_core_collection(
  G = pred$G, n_core = n_founders, type = "relationship",
  strategy = "maximin", merit = pred$gebv, seed = 1, verbose = FALSE
)

selection_groups <- list(
  Truncation  = ts_sel$selected,
  GA          = ga_sel$selected,
  FamilyQuota = fam_sel$selected,
  Core        = core_res$selected
)
clust_res <- cluster_selection_groups(
  G = pred$G, groups = selection_groups,
  variance_threshold = 0.95, method = "hierarchical",
  n_clusters = 3L, verbose = FALSE
)
clust_res$table
#>     Cluster n_total pct_population n_Truncation prop_Truncation n_GA    prop_GA
#> 1 Cluster 1      78       65.00000           13      0.16666667   11 0.14102564
#> 2 Cluster 2      14       11.66667            1      0.07142857    2 0.14285714
#> 3 Cluster 3      28       23.33333            1      0.03571429    2 0.07142857
#>   n_FamilyQuota prop_FamilyQuota n_Core prop_Core
#> 1            11       0.14102564     10 0.1282051
#> 2             2       0.14285714      2 0.1428571
#> 3             2       0.07142857      3 0.1071429
```

### Stage 15 – Validate the result contract

Every breeder-facing result from stage 3 onward carries a
`hapblockr_result` contract;
[`validate()`](https://FAkohoue.github.io/HapBlockR/reference/validate.md)
checks the schema, immutable identifiers, input hashes, quality gates
and decision table before a result is treated as a recommendation. See
“Interpreting quality control and uncertainty” and “Decision sign-off”
in the Breeder’s Guide.

``` r
validate(ga_sel)
#>                 check passed                   detail
#> 1     required_fields   TRUE                         
#> 2      schema_version   TRUE                    1.0.0
#> 3              method   TRUE        select_parents_ga
#> 4   sample_ids_unique   TRUE 120 sample identifier(s)
#> 5  variant_ids_unique   TRUE  6 variant identifier(s)
#> 6 input_hashes_sha256   TRUE         4 input hash(es)
#> 7       quality_gates   TRUE        3 quality gate(s)
#> 8   validation_status   TRUE                   passed
validate(fam_sel)
#>                 check passed                   detail
#> 1     required_fields   TRUE                         
#> 2      schema_version   TRUE                    1.0.0
#> 3              method   TRUE select_parents_by_family
#> 4   sample_ids_unique   TRUE 120 sample identifier(s)
#> 5  variant_ids_unique   TRUE  0 variant identifier(s)
#> 6 input_hashes_sha256   TRUE         2 input hash(es)
#> 7       quality_gates   TRUE        2 quality gate(s)
#> 8   validation_status   TRUE                   passed
validate(core_res)
#>                 check passed                   detail
#> 1     required_fields   TRUE                         
#> 2      schema_version   TRUE                    1.0.0
#> 3              method   TRUE   select_core_collection
#> 4   sample_ids_unique   TRUE 120 sample identifier(s)
#> 5  variant_ids_unique   TRUE  0 variant identifier(s)
#> 6 input_hashes_sha256   TRUE         1 input hash(es)
#> 7       quality_gates   TRUE        2 quality gate(s)
#> 8   validation_status   TRUE                   passed
validate(certificate)
#>                 check passed                  detail
#> 1     required_fields   TRUE                        
#> 2      schema_version   TRUE                   1.0.0
#> 3              method   TRUE     certify_mating_plan
#> 4   sample_ids_unique   TRUE  9 sample identifier(s)
#> 5  variant_ids_unique   TRUE 0 variant identifier(s)
#> 6 input_hashes_sha256   TRUE        5 input hash(es)
#> 7       quality_gates   TRUE       6 quality gate(s)
#> 8   validation_status   TRUE                  passed
validate(exact_res)
#>                 check passed                  detail
#> 1     required_fields   TRUE                        
#> 2      schema_version   TRUE                   1.0.0
#> 3              method   TRUE  validate_crosses_exact
#> 4   sample_ids_unique   TRUE 15 sample identifier(s)
#> 5  variant_ids_unique   TRUE 0 variant identifier(s)
#> 6 input_hashes_sha256   TRUE        1 input hash(es)
#> 7       quality_gates   TRUE       2 quality gate(s)
#> 8   validation_status   TRUE                  passed
```

### Stage 16 – Export for exchange

[`build_breeding_exchange()`](https://FAkohoue.github.io/HapBlockR/reference/build_breeding_exchange.md)
converts a validated result into BrAPI/MIAPPE -oriented tables, and
[`write_breeding_exchange()`](https://FAkohoue.github.io/HapBlockR/reference/write_breeding_exchange.md)
writes them as a checksummed CSV bundle a programme’s own integration
layer can consume. Any `hapblockr_result` qualifies; the diversity-first
core-collection result from stage 14 is exported here because it best
illustrates the cross-clustering check from that stage. See
“Standards-oriented exchange” in the *programme operations* vignette.

``` r
bundle <- build_breeding_exchange(
  core_res,
  breeding_program_id = "DEMO", study_id = "DEMO-CYCLE-1", trial_id = "DEMO-T1"
)
head(bundle$recommendations[, c("entity_id", "rank", "score", "recommendable")])
#>   entity_id rank score recommendable
#> 1    ind013   15     1          TRUE
#> 2    ind003   14     2          TRUE
#> 3    ind093   13     3          TRUE
#> 4    ind015   12     4          TRUE
#> 5    ind025   11     5          TRUE
#> 6    ind029   10     6          TRUE

exchange_path <- tempfile("hapblockr-full-pipeline-")
write_breeding_exchange(bundle, exchange_path)
list.files(exchange_path)
#> [1] "exclusions.csv"      "field_mapping.csv"   "identifiers.csv"    
#> [4] "input_hashes.csv"    "manifest.csv"        "provenance.csv"     
#> [7] "quality_gates.csv"   "recommendations.csv" "uncertainty.csv"
```

## Where each stage is covered in depth

| Stage | Vignette / Guide section |
|:---|:---|
| 1 | *Programme operations* vignette, “Externally analysed targets” |
| 2 | *Intro* vignette; *Large-scale* vignette for file-backed/WGS routes |
| 3-4 | *Workflow* vignette |
| 5-8 | *Breeding decisions* vignette (all four shortlist tools compared) |
| 9, 12 | *Programme operations* vignette, “Operational screening”/“Feasibility certificate” |
| 10-11, 13 | *Breeding decisions* vignette; Breeder’s Guide “The nine decision tools and their variants” |
| 14 | *Breeding decisions* vignette |
| 15 | Breeder’s Guide “Interpreting quality control and uncertainty” |
| 16 | *Programme operations* vignette, “Standards-oriented exchange” |

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
#>  [1] rrBLUP_4.6.3         pspline_1.0-21       xfun_0.57           
#>  [4] bslib_0.11.0         nadiv_2.18.0         htmlwidgets_1.6.4   
#>  [7] lattice_0.22-6       numDeriv_2016.8-1.1  quadprog_1.5-8      
#> [10] vctrs_0.7.3          tools_4.5.0          parallel_4.5.0      
#> [13] rgl_1.3.36           optiSel_2.1.0        pkgconfig_2.0.3     
#> [16] Matrix_1.7-3         data.table_1.18.4    desc_1.4.3          
#> [19] scatterplot3d_0.3-45 alabama_2025.1.0     optiSolve_1.0       
#> [22] lifecycle_1.0.5      compiler_4.5.0       stringr_1.6.0       
#> [25] textshaping_1.0.5    minpack.lm_1.2-4     codetools_0.2-20    
#> [28] kinship2_1.9.6.2     ECOSolveR_0.6.1      htmltools_0.5.9     
#> [31] sass_0.4.10          cccp_0.3-3           yaml_2.3.12         
#> [34] Rttf2pt1_1.3.14      pedigree_1.4.2       pkgdown_2.2.0       
#> [37] nloptr_2.2.1         crayon_1.5.3         jquerylib_0.1.4     
#> [40] extrafontdb_1.1      MASS_7.3-65          fitdistrplus_1.2-6  
#> [43] cachem_1.1.0         iterators_1.0.14     abind_1.4-8         
#> [46] foreach_1.5.2        digest_0.6.39        stringi_1.8.7       
#> [49] reshape2_1.4.5       purrr_1.2.2          magic_1.6-1         
#> [52] splines_4.5.0        extrafont_0.20       fastmap_1.2.0       
#> [55] grid_4.5.0           cli_3.6.6            magrittr_2.0.5      
#> [58] base64enc_0.1-6      survival_3.8-3       rmarkdown_2.31      
#> [61] shapes_1.2.8         igraph_2.3.1         otel_0.2.0          
#> [64] ragg_1.5.2           lpSolve_5.6.23       HaploSim_1.8.4.2    
#> [67] evaluate_1.0.5       GA_3.2.5             knitr_1.51          
#> [70] doParallel_1.0.17    rlang_1.2.0          Rcpp_1.1.1-1.1      
#> [73] glue_1.8.1           rstudioapi_0.18.0    reshape_0.8.10      
#> [76] jsonlite_2.0.0       R6_2.6.1             plyr_1.8.9          
#> [79] systemfonts_1.3.2    fs_2.1.0
```
