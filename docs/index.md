# HapBlockR — Genome-Wide Linkage Disequilibrium (LD) Block Detection, Haplotype Analysis, and Genomic Prediction Features

------------------------------------------------------------------------

## Table of contents

1.  [Motivation](#id_1-motivation)
2.  [Summary](#id_2-summary)
3.  [Relationship to the original Big-LD
    algorithm](#id_3-relationship-to-the-original-big-ld-algorithm)
    - 3.1. [Computational core: R loops replaced by
      C++](#id_31-computational-core-r-loops-replaced-by-c)
    - 3.2. [Memory model:
      never-full-genome](#id_32-memory-model-never-full-genome)
    - 3.3. [Kinship correction: rV²](#id_33-kinship-correction-rv)
    - 3.4. [Singleton SNP handling](#id_34-singleton-snp-handling)
    - 3.5. [Bug fix: zero-row
      assignment](#id_35-bug-fix-zero-row-assignment)
    - 3.6. [Downstream pipeline](#id_36-downstream-pipeline)
    - 3.7. [What is kept exactly](#id_37-what-is-kept-exactly)
4.  [Installation](#id_4-installation)
5.  [Documentation](#id_5-documentation)
6.  [Quick start](#id_6-quick-start)
7.  [Input formats](#id_7-input-formats)
    - 7.1. [Genotype input format](#id_71-genotype-input-format)
    - 7.2. [Phenotype input format](#id_72-phenotype-input-format)
      - 7.2.1. [Format 1 — Named numeric
        vector](#id_721-format-1-named-numeric-vector-single-trait-simplest)
      - 7.2.2. [Format 2 — Data frame, single
        trait](#id_722-format-2-data-frame-single-trait)
      - 7.2.3. [Format 3 — Data frame, multiple
        traits](#id_723-format-3-data-frame-multiple-traits)
      - 7.2.4. [Format 4 — Named
        list](#id_724-format-4-named-list-different-individuals-per-trait)
    - 7.3. [ID matching rules](#id_73-id-matching-rules)
    - 7.4. [Preparing BLUEs from raw phenotype
      data](#id_74-preparing-blues-from-raw-phenotype-data)
8.  [Statistical background](#id_8-statistical-background)
    - 8.1. [MAF filtering](#id_81-maf-filtering)
    - 8.2. [Genotype preparation](#id_82-genotype-preparation)
    - 8.3. [Subsegmentation](#id_83-subsegmentation)
    - 8.4. [Clique detection (CLQD)](#id_84-clique-detection-clqd)
    - 8.5. [Block construction](#id_85-block-construction)
9.  [LD metrics: r² versus rV²](#id_9-ld-metrics-r-versus-rv)
10. [Clique detection mode
    (CLQmode)](#id_10-clique-detection-mode-clqmode)
    - 10.1. [Mode 1 —
      Density](#id_101-mode-1-density-clqmode-density-default)
    - 10.2. [Mode 2 — Maximal](#id_102-mode-2-maximal-clqmode-maximal)
    - 10.3. [Mode 3 — Louvain](#id_103-mode-3-louvain-clqmode-louvain)
    - 10.4. [Mode 4 — Leiden](#id_104-mode-4-leiden-clqmode-leiden)
    - 10.5. [Summary comparison](#id_105-summary-comparison)
    - 10.6. [Recommended
      configurations](#id_106-recommended-configurations)
11. [Haplotype analysis](#id_11-haplotype-analysis)
    - 11.1. [Phase-free haplotype
      extraction](#id_111-phase-free-haplotype-extraction)
    - 11.2. [Haplotype diversity
      metrics](#id_112-haplotype-diversity-metrics)
    - 11.3. [Haplotype feature matrix for genomic
      prediction](#id_113-haplotype-feature-matrix-for-genomic-prediction)
    - 11.4. [Haplotype-based genomic prediction
      pipeline](#id_114-haplotype-based-genomic-prediction-pipeline)
    - 11.5. [Cross-validation and prediction
      accuracy](#id_115-cross-validation-and-prediction-accuracy)
    - 11.6. [Between-population
      comparison](#id_116-between-population-comparison)
    - 11.7. [Haplotype network
      visualisation](#id_117-haplotype-network-visualisation)
    - 11.8. [Multi-environment
      stability](#id_118-multi-environment-stability)
    - 11.9. [Candidate region export](#id_119-candidate-region-export)
    - 11.10. [Per-allele effect
      decomposition](#id_1110-per-allele-effect-decomposition)
    - 11.11. [Sliding-window diversity
      scan](#id_1111-sliding-window-diversity-scan)
    - 11.12. [True diplotype
      inference](#id_1112-true-diplotype-inference)
    - 11.13. [Rare-allele collapsing](#id_1113-rare-allele-collapsing)
    - 11.14. [Cross-panel
      harmonisation](#id_1114-cross-panel-harmonisation)
    - 11.15. [Haplotype association
      testing](#id_1115-haplotype-association-testing)
    - 11.16. [Breeding decision tools](#id_1116-breeding-decision-tools)
    - 11.17. [Cross-population effect concordance
      (haplotype)](#id_1117-cross-population-effect-concordance)
    - 11.18. [Cross-population effect concordance (external
      GWAS)](#id_1118-cross-population-effect-concordance-external-gwas)
    - 11.19. [Within-block and between-block epistasis
      detection](#id_1119-within-block-and-between-block-epistasis-detection)
12. [Parameter auto-tuning](#id_12-parameter-auto-tuning)
13. [Scale strategies and
    backends](#id_13-scale-strategies-and-backends)
    - 13.1. [The HapBlockR_backend
      interface](#id_131-the-hapblockrbackend-interface)
    - 13.2. [Memory requirements by
      configuration](#id_132-memory-requirements-by-configuration)
    - 13.3. [Recommended configurations by dataset
      size](#id_133-recommended-configurations-by-dataset-size)
14. [Full pipeline walkthrough](#id_14-full-pipeline-walkthrough)
15. [Function reference](#id_15-function-reference)
    - 15.1. [Main pipeline](#id_151-main-pipeline)
    - 15.2. [I/O](#id_152-io)
    - 15.3. [LD computation](#id_153-ld-computation)
    - 15.4. [C++ kernels (direct
      access)](#id_154-c-kernels-direct-access)
    - 15.5. [Haplotype analysis](#id_155-haplotype-analysis)
    - 15.6. [Analysis extensions](#id_156-analysis-extensions)
    - 15.7. [True haplotype inference and
      harmonisation](#id_157-true-haplotype-inference-and-harmonisation)
    - 15.8. [Haplotype association
      testing](#id_158-haplotype-association-testing)
    - 15.9. [Breeding decision tools](#id_159-breeding-decision-tools)
    - 15.10. [Utilities](#id_1510-utilities)
16. [Output objects](#id_16-output-objects)
    - 16.1. [`run_Big_LD_all_chr()` — block
      table](#id_161-runbigldallchr-block-table)
    - 16.2. [`run_ldx_pipeline()` — named
      list](#id_162-runldxpipeline-named-list)
    - 16.3. [`tune_LD_params()` — named
      list](#id_163-tuneldparams-named-list)
    - 16.4. [`extract_haplotypes()` — named
      list](#id_164-extracthaplotypes-named-list)
    - 16.5. [`compute_haplotype_diversity()` —
      data.frame](#id_165-computehaplotypediversity-dataframe)
    - 16.6. [`read_geno()` —
      HapBlockR_backend](#id_166-readgeno-hapblockrbackend)
    - 16.7. [`test_block_haplotypes()` —
      HapBlockR_haplotype_assoc](#id_167-test_block_haplotypes--hapblockr_haplotype_assoc)
    - 16.8. [`estimate_diplotype_effects()` —
      HapBlockR_diplotype](#id_168-estimate_diplotype_effects--hapblockr_diplotype)
    - 16.9. [`score_favorable_haplotypes()` — data
      frame](#id_169-score_favorable_haplotypes--data-frame)
    - 16.10. [`summarize_parent_haplotypes()` — data
      frame](#id_1610-summarize_parent_haplotypes--data-frame)
    - 16.11. [`compare_block_effects()` —
      HapBlockR_effect_concordance](#id_1611-compare_block_effects--hapblockr_effect_concordance)
    - 16.12. [`compare_gwas_effects()` —
      HapBlockR_effect_concordance](#id_1612-compare_gwas_effects--hapblockr_effect_concordance)
17. [Memory and performance notes](#id_17-memory-and-performance-notes)
    - 17.1. [C++ core](#id_171-c-core)
    - 17.2. [Never-full-genome memory
      model](#id_172-never-full-genome-memory-model)
    - 17.3. [OpenMP thread count](#id_173-openmp-thread-count)
18. [Citation](#id_18-citation)
19. [Contributing](#id_19-contributing)
20. [License](#id_20-license)
21. [References](#id_21-references)

------------------------------------------------------------------------

## 1. Motivation

Genomic selection tells a breeder *how good* a candidate is. On its own,
it rarely tells them *which specific parents to cross*, *whether the
favourable genetics they carry are actually complementary as a set*, or
*whether that choice will outperform simply taking the top-ranked
individuals*. Answering those questions means going one level below a
single breeding value — down to which haplotype blocks are actually
driving it — and then searching for the founder set that jointly covers
the most favourable blocks, rather than eyeballing a ranked list.

Closing that gap is what HapBlockR’s breeding decision layer is built
for. In short: it decomposes a breeding value you already have (from
your own GBLUP (Genomic Best Linear Unbiased Prediction)/mixed-model
workflow) down to the block level, searches for the founder set that
best covers the important blocks under real crossing-scheme constraints,
benchmarks that founder set against a truncation-selection baseline, and
then carries it forward to a ranked list of crosses or a full
contribution-and-mating plan. Section 2 below breaks down every function
in the layer — from block-level GEBV (Genomic Estimated Breeding Value)
decomposition through to Optimal Contribution Selection — in full. See
the *From Local GEBV to a Crossing Decision* vignette for the worked
example, start to finish: breeding value to crossing block to mating
plan.

Reaching that decision reliably still depends on the block detection
underneath it being correct for the population at hand, which is why LD
block detection remains a foundational step in modern genomic analyses.
Knowing which SNPs (Single Nucleotide Polymorphisms) co-segregate as a
unit determines how GWAS (Genome-Wide Association Study) results are
interpreted, how haplotypes are defined for population genetics, and how
genomic prediction models should be structured. Despite its importance,
most implementations of LD block detection suffer from three problems
that limit their usefulness in practice.

**The kinship problem.** Classical LD estimators (r²) assume
independence between individuals. In livestock, crop, or family-based
human cohorts this assumption is systematically violated. Cryptic
relatedness inflates pairwise correlations, causing LD to appear
stronger than it is and blocks to be drawn too broadly or in the wrong
places. The kinship-adjusted squared correlation rV² (Mangin et
al. 2012, *Heredity* 108:285-291) corrects for this by whitening the
genotype matrix with the inverse square root of the genomic relationship
matrix (GRM):

``` R
rV²ᵢⱼ = Cov(Xᵛᵢ, Xᵛⱼ)² / [Var(Xᵛᵢ) · Var(Xᵛⱼ)]
```

where **X**ᵛ = **V**⁻¹⁄² **G̃** is the kinship-whitened, mean-centred
genotype matrix; **V** = **ZZ**ᵀ / (2 Σⱼ pⱼqⱼ) is the VanRaden (2008)
GRM; and **V**⁻¹⁄² is its inverse square root.

**The scale problem.** The original Big-LD implementation (Kim et
al. 2018) contains no compiled code. For modern whole-genome sequencing
panels with 2-10 million markers the inner loops are prohibitively slow,
and loading the full genotype matrix before detection is impossible on
most workstations. HapBlockR addresses this with a C++/Armadillo
computational core compiled via Rcpp, OpenMP-parallelised LD
computation, a unified multi-format I/O layer, and a strict
never-full-genome memory model.

**The pipeline gap.** The original Big-LD stops at block boundaries.
HapBlockR adds a complete downstream pipeline: statistical phasing via
Beagle 5.x with bigmemory-backed caching, haplotype extraction,
diversity metrics, genomic prediction features, association testing, and
breeding decision tools.

------------------------------------------------------------------------

## 2. Summary

![HapBlockR architecture: inputs and scalable data access; four-stage
core workflow (data preparation, LD block detection, haplotype
reconstruction, feature construction); one haplotype layer feeding seven
analytical pathways (diversity and populations, haplotype-based genomic
prediction, association testing, cross-population concordance, epistasis
and interactions, parent selection, forward-in-time simulation);
computational foundation](reference/figures/HapBlockR_schematic.png)

`HapBlockR` takes a genotyped panel and a breeding value you already
have (from any GBLUP/mixed-model workflow) through to a specific,
defensible set of parents for the next crossing cycle — not just an LD
block table. It extends the Big-LD algorithm of Kim et al. (2018) with a
breeding decision layer built on top of a substantially re-engineered
detection core.

**Breeding decision layer** — see the *From Local GEBV to a Crossing
Decision* vignette for the full worked example:

- **Local GEBV per haplotype block** —
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
  (Tong et al. 2025) decomposes a breeding value into per-block
  contributions, so a breeder can see *where* it comes from, not just
  *how much* of it there is.
- **Block-importance triage** —
  [`select_top_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  and
  [`plot_block_funnel()`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md)
  separate the small set of blocks driving most of the variance from the
  rest.
- **GA-based parent selection** —
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  searches for the founder set that jointly covers the most favourable
  blocks, under five real crossing-scheme constraints
  (no_selfing/selfing/OHS/OPV/ Haploid_OHS), benchmarked against a
  [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
  baseline. Optional `merit_weight` rewards whole-genome merit directly
  inside the search itself (a GA+TS hybrid), guarding against selecting
  on noisy per-block estimates. Don’t want to guess a raw multiplier?
  Set `merit_priority` instead — a plain 0-100 “how much do I care about
  merit vs. coverage” dial, auto-calibrated from your own data (or call
  [`suggest_merit_weight()`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)
  first to see the calibration numbers before committing to one).
  Optional `coancestry_weight`/`G` penalise relatedness directly inside
  the search too; `target_degree` is the easier alternative here — the
  same 0-90 dial
  [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
  uses, converted into a relatedness ceiling from your own data rather
  than a raw penalty multiplier.
- **Family-quota parent selection** —
  [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)
  picks the best-performing families (or, via
  `group_by = "genetic_cluster"`, data-derived genetic clusters) first,
  then the best lines within each, enforcing group balance directly in
  the selection rule. Family ranking is shrinkage-corrected by default
  (small families are pulled toward the pack rather than trusted on a
  raw sample mean alone — see `NEWS.md` for the breaking-change note);
  an optional relatedness ceiling (`within_group_target_degree`)
  balances score against diversity within each group; and an optional
  cross-family haplotype-diversity check (coverage-gain by default, or
  the original dominant-block heuristic) discourages picking
  representatives of different groups that are genomically redundant.
- **Forward-simulation proof of genetic gain** —
  [`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)
  and
  [`plot_ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md)
  simulate recurrent selection over generations for both founder sets,
  so the GA’s answer can be checked against the naive one rather than
  taken on faith.
- **Diversity and family-balance diagnostics** —
  [`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)
  shows where each selection sits in the population’s genetic diversity;
  pairing it with your own pedigree/family table catches an
  over-represented family before the crosses are made.
- **Cross ranking and true Optimal Contribution Selection (OCS)** —
  [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
  ranks candidate crosses by predicted mid-parent value plus segregation
  variance (Schnell & Utz 1975; Bernardo 2003), with four
  variance-prediction modes: unphased block-independent, exact phased
  4-gamete, a native Monte Carlo linkage-aware mode
  (`variance_model = "linked"`, no extra dependency), and a genuine
  multi-locus, linkage-aware mode via
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)
  (Peixoto et al. 2024).
  [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
  solves the classical Meuwissen (1997) contribution- optimisation
  problem via two true-OCS engines (the AlphaMate executable, or
  optiSel’s own `candes()`/`opticont()`/`matings()` solver, Wellmann
  2019), producing an actual contribution-and-mating plan with
  population- wide relatedness under an explicit cap; a third engine
  (SimpleMating’s `planCross()`/`selectCrosses()`) instead does discrete
  greedy cross prediction/selection — a different algorithm class,
  useful when you want fast candidate-pair screening rather than a
  solved contribution optimum.
- **Explicit tradeoffs and validation** —
  [`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)
  sweeps the merit-vs-relatedness tradeoff into a visible Pareto
  frontier instead of one scalar dial;
  [`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)
  checks any heuristic mating plan against a true optimum via exact
  integer linear programming;
  [`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
  covers the reverse case where representing genetic diversity is itself
  the objective (genebank/reference-panel curation), via the classical
  farthest-point maximin heuristic.
- **Haplotype stacking and association tools** —
  [`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md),
  [`summarize_parent_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/summarize_parent_haplotypes.md),
  [`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md),
  and
  [`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
  round out the inference layer feeding into these breeding decisions.

**Detection and computation layer** — 15 core improvements over the
original Big-LD implementation make the layer above possible at
genome-wide scale:

1.  **Dual LD metric** — standard r² (default) and kinship-adjusted rV²
2.  **C++/Armadillo core** — thirteen compiled functions handle all
    expensive operations
3.  **OpenMP parallelism** — outer loop of `compute_r2_cpp()`
    parallelised
4.  **Unified multi-format I/O** — numeric CSV, HapMap, VCF/VCF.gz, GDS,
    PLINK BED, R matrix
5.  **Never-full-genome memory model** — genome never held in RAM at
    once for any format
6.  **MAF filter in C++** — single O(np) pass with NA imputation
7.  **C++ boundary scan** — replaces R inner loop in subsegmentation
8.  **Sparse r² computation** — O(p) cost for large sub-segments
9.  **Automatic parameter tuning** —
    [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md)
    grid search
10. **Haplotype reconstruction** — phase-free and Beagle-phased pathways
11. **Diversity metrics** — richness, He, Shannon entropy, dominant
    haplotype frequency
12. **Prediction feature matrix** — multi-locus dosage columns for
    GBLUP/BayesB/ML
13. **Polynomial community detection** — Louvain and Leiden for WGS
    panels
14. **Sparse LD computation** — `max_bp_distance` restricts pairs to a
    physical window
15. **Memory-mapped genotype store** —
    [`read_geno_bigmemory()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno_bigmemory.md)
    with OS-page access

------------------------------------------------------------------------

## 3. Relationship to the original Big-LD algorithm

### 3.1. Computational core: R loops replaced by C++

| Original R operation | HapBlockR C++ function | Speedup |
|----|----|----|
| `cor(subgeno)` per CLQD call | `compute_r2_cpp()` + OpenMP | ~40× for 1,500-SNP window |
| `apply(Ogeno, 2, ...)` MAF filter | `maf_filter_cpp()` single pass | ~10× for 100k+ SNPs |
| [`cor()`](https://rdrr.io/r/stats/cor.html) inside boundary-scan loop | `boundary_scan_cpp()` compiled | ~20× per chromosome |
| `r2Mat[r2Mat >= CLQcut^2] <- 1` | `build_adj_matrix_cpp()` | eliminates intermediate allocation |
| Single-column correlation | `col_r2_cpp()` | used in boundary scan helper |
| Sparse within-window r² | `compute_r2_sparse_cpp()` | avoids O(p²) for large segments |
| LD-informed overlap resolution | `resolve_overlap_cpp()` | 15,700× per-SNP reduction on chr1 |
| Haplotype string building | `build_hap_strings_cpp()` | ~20-50× per block |
| Block-to-SNP interval lookup | `block_snp_ranges_cpp()` | O(p + n_blocks) single sweep |
| Chromosome haplotype extraction | `extract_chr_haplotypes_cpp()` | strings + freq tabulation in one OpenMP pass |
| Call-rate filter + imputation | `impute_and_filter_cpp()` | single O(n×p) pass |

### 3.2. Memory model: never-full-genome

The full genotype matrix is never held in RAM at once for any format.
Numeric dosage CSV is read in pre-allocated 50,000-row chunks. VCF and
HapMap auto-convert to a streaming GDS cache. GDS and PLINK BED backends
load only the SNP window per CLQD call. `gc(FALSE)` is called after each
chromosome.

### 3.3. Kinship correction: rV²

`method = "rV2"` replaces every pairwise correlation with the
kinship-whitened equivalent. The whitening factor **A** is computed once
per chromosome from the VanRaden (2008) GRM via
[`get_V_inv_sqrt()`](https://FAkohoue.github.io/HapBlockR/reference/get_V_inv_sqrt.md)
(Cholesky or eigendecomposition), then applied to the centred genotype
matrix. In related populations, rV² blocks are typically 10-30% smaller
and more precisely delimited than r² blocks.

### 3.4. Singleton SNP handling

`singleton_as_block = TRUE` collects SNPs that receive `NA` from
[`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md) and
appends them to the block table as single-SNP entries (`start == end`).
Default `FALSE` preserves original behaviour.

### 3.5. Bug fix: zero-row assignment

When a sub-segment contains no valid cliques, `nowLDblocks` has zero
rows. HapBlockR wraps the assignment with `if (nrow(nowLD) > 0L)`,
preventing the backwards-sequence error in the original code.

### 3.6. Downstream pipeline

| Capability | Original Big-LD / gpart | HapBlockR |
|----|----|----|
| Block detection | Yes (core algorithm) | Yes (same algorithm + C++ + rV²) |
| Statistical phasing | No | `run_ldx_pipeline(phase = TRUE)` calls [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md) internally. Phased hap1/hap2/dosage cached as bigmemory backends for fast restart. [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md) reads pre-phased VCF output. |
| Haplotype extraction | No | [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md) — phased and unphased, backend streaming |
| Diversity metrics | No | [`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md) — He, Shannon, richness, f_max |
| Post-GWAS QTL mapping | No | [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md) — pleiotropic block detection |
| Genomic prediction features | No | [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md) |
| Output writers | No | Numeric dosage matrix, nucleotide character matrix, diversity CSV |
| Parameter auto-tuning | No | [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md) — grid search against GWAS marker coverage |
| Multi-format I/O | PLINK, VCF (gpart) | Numeric CSV, HapMap, VCF, GDS, BED, R matrix via unified backend |
| WGS-scale streaming | Partial (gpart GDS) | Full never-full-genome model for all formats |

### 3.7. What is kept exactly

- The interval graph modelling of LD bins (cliques of strong pairwise LD
  SNPs)
- [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md):
  bin vector assignment via maximal clique enumeration and greedy
  density-priority selection
- `constructLDblock()`: maximum-weight independent set via dynamic
  programming
- `appendSGTs()`: rare-SNP appending logic
- `cutsequence.modi()`: boundary-scan logic and forced-split fall-back
- All `CLQmode = "Density"` and `CLQmode = "Maximal"` clique scoring
- The `clstgap` physical distance splitting within cliques
- Block table column format for drop-in compatibility with downstream
  tools

------------------------------------------------------------------------

## 4. Installation

Follow these steps in order. Step 1 matters most on institutional/HPC
servers with restricted or flaky outbound access to GitHub — skipping it
is the most common cause of installation failures reported for this
package.

### Step 1: Install `remotes`

``` r
install.packages("remotes")
```

### Step 2: Install [genomicSimulation](https://github.com/vllrs/genomicSimulation) first

HapBlockR’s forward-simulation function
([`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md))
wraps [genomicSimulation](https://github.com/vllrs/genomicSimulation),
which is not on CRAN and is resolved via the `Remotes:` field in
`DESCRIPTION`. Installing HapBlockR in Step 5 pulls it in automatically,
but letting
[`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html)
fetch it as a *nested* dependency mid-install is where GitHub API
timeouts most often happen — a full HapBlockR install makes several
sequential calls to `api.github.com` in one session (to resolve
HapBlockR itself, then genomicSimulation), and any one of them failing
aborts the whole install. Installing genomicSimulation standalone first,
as a single isolated request, avoids that:

``` r
remotes::install_github("vllrs/genomicSimulation")
```

Confirm it installed successfully before moving on:

``` r
requireNamespace("genomicSimulation", quietly = TRUE)   # should return TRUE
```

### Step 3 (optional, but install now if you’ll use it): Install [optiSel](https://cran.r-project.org/package=optiSel) and/or [SimpleMating](https://github.com/Resende-Lab/SimpleMating)

[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
offers three engines, two of which need one of these packages:
`engine = "optisel"` (true Optimal Contribution Selection via optiSel’s
own `candes()`/`opticont()`/`matings()`/`noffspring()` solver, no
external executable) needs only **optiSel**, which is on CRAN;
`engine = "simplemating"` (discrete greedy cross prediction/selection,
not true OCS) and
`usefulness_criterion(variance_model = "simplemating")` wrap
[SimpleMating](https://github.com/Resende-Lab/SimpleMating), which is
**not on CRAN** and is resolved via the `Remotes:` field in
`DESCRIPTION` (SimpleMating pulls in optiSel itself as a transitive
dependency, but does not use it to implement `engine = "simplemating"`,
which is SimpleMating’s own algorithm). The same nested-dependency
GitHub API timeout risk from Step 2 applies to SimpleMating, for the
same reason: install it standalone first, as its own isolated request,
rather than letting
[`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html)
fetch it mid-install in Step 5.

``` r
install.packages("optiSel")
remotes::install_github("Resende-Lab/SimpleMating")
```

If you only want `engine = "optisel"` (true OCS, no external binary,
recommended default if you don’t have AlphaMate),
`install.packages( "optiSel")` alone is enough — skip the SimpleMating
install entirely.

Confirm optiSel installed successfully:

``` r
requireNamespace("optiSel", quietly = TRUE)                          # should return TRUE
all(c("candes", "opticont", "matings", "noffspring") %in%
      getNamespaceExports("optiSel"))                                # should return TRUE
```

If you also installed SimpleMating, confirm it installed successfully,
**and** that it’s a recent enough version — SimpleMating’s own exported
API has moved between releases, and HapBlockR targets
`planCross()`/`selectCrosses()`/`getUsefA()`, all confirmed present in
`SimpleMating >= 0.2.1`:

``` r
requireNamespace("SimpleMating", quietly = TRUE)                     # should return TRUE
all(c("planCross", "selectCrosses", "getUsefA") %in%
      getNamespaceExports("SimpleMating"))                           # should return TRUE
```

If the second check returns `FALSE`, your installed copy predates the
current API — reinstall with `force = TRUE` (a plain `install_github()`
call is a no-op if a copy is already installed, even an outdated one):

``` r
remotes::install_github("Resende-Lab/SimpleMating", force = TRUE)
```

You can skip this step entirely if you don’t need any true Optimal
Contribution Selection engine or SimpleMating’s cross prediction/
usefulness-criterion mode:
`usefulness_criterion(variance_model = "linked")` is a native,
dependency-free linkage-aware alternative (see the function reference
table), and `select_parents_ocs(engine = "alphamate")` needs the
AlphaMate executable instead, not any R package here.

### Step 4 (only if Step 2 or Step 3 times out): set a GitHub token

A connection timeout to `api.github.com` (rather than an error about the
package itself) means the network path is the problem, not
genomicSimulation or SimpleMating. Setting a GitHub personal access
token raises the GitHub API rate limit from 60/hour to 5000/hour and is
generally routed more reliably:

``` r
Sys.setenv(GITHUB_PAT = "your_token_here")
```

Then retry Step 2 and/or Step 3.

### Step 5: Install HapBlockR

With genomicSimulation (and, if you installed it, SimpleMating) already
present, this step only needs to resolve HapBlockR itself from GitHub —
`remotes` detects the already-satisfied dependencies and skips
re-fetching them:

``` r
remotes::install_github("FAkohoue/HapBlockR",
  build_vignettes = TRUE,
  dependencies    = TRUE
)
```

### Step 6: Required dependencies

These install automatically as part of Step 5; listed here for reference
or in case any need installing manually:

``` r
install.packages(c("Rcpp", "RcppArmadillo", "igraph", "data.table", "dplyr"))
```

### Step 7: Optional dependencies

Install only the ones your workflow needs:

``` r
# GDS backend — required for .gds files; recommended for panels > 2 M SNPs
BiocManager::install("SNPRelate")

# PLINK BED backend
install.packages("BEDMatrix")

# Kinship-adjusted rV² (method = "rV2")
install.packages(c("AGHmatrix", "ASRgenomics"))

# Parallel parameter tuning
install.packages("future.apply")

# LD block visualisation
install.packages("ggplot2")

# Bigmemory-backed genotype store
install.packages("bigmemory")

# GA-based parent selection (select_parents_ga(), select_parents_pareto())
install.packages("GA")

# estimate_marker_effects()'s BayesA/B/C, and
# run_haplotype_prediction(include_dominance = TRUE)'s dual-kernel additive
# + dominance GBLUP
install.packages("BGLR")

# LASSO interaction search for large-block epistasis fine-mapping
install.packages("glmnet")

# select_parents_ocs(engine = "optisel") -- true Optimal Contribution
# Selection, no external binary -- needs only optiSel (on CRAN):
install.packages("optiSel")

# select_parents_ocs(engine = "simplemating") and
# usefulness_criterion(variance_model = "simplemating") need SimpleMating
# -- see Step 3 above (a GitHub-only remote, not a plain install.packages()
# dependency, and worth installing standalone rather than as a nested
# mid-install fetch). variance_model = "linked" needs NEITHER of these --
# it is a native, dependency-free Monte Carlo alternative; see the function
# reference table.

# validate_crosses_exact() — exact ILP cross-selection validation
install.packages("lpSolve")
```

### Providing the AlphaMate executable

`select_parents_ocs(engine = "alphamate")` requires the AlphaMate
executable ([Hickey group / AlphaGenes
suite](https://github.com/AlphaGenes/AlphaMate), Fortran source, MIT
licensed) — a separate, third-party tool, not an R package. **You must
obtain or build it yourself and be entitled to use it under its own
license terms.**

No AlphaMate binary is bundled with or downloaded by HapBlockR:
AlphaGenes does not currently publish pre-built binaries for that
repository (no GitHub Releases), so there is no stable URL to fetch one
from, and CRAN policy prohibits shipping compiled executables in a
source package in any case. Build it yourself (a Fortran compiler such
as `gfortran` is required) or otherwise source a working executable,
then **always pass its full path explicitly** — there is no
auto-detection:

``` r
select_parents_ocs(
  merit         = my_gebv,
  G             = my_grm,
  engine        = "alphamate",
  alphamate_exe = "/full/path/to/AlphaMate",   # or AlphaMate.exe on Windows
  out_dir       = "ocs_results"
)
```

If you don’t have access to AlphaMate, `engine = "optisel"` (true
Optimal Contribution Selection via optiSel’s own solver) requires no
external executable at all — only the `optiSel` R package above — and is
the more portable default, recommended over `engine = "simplemating"`
unless you specifically want SimpleMating’s discrete greedy
cross-selection algorithm rather than a solved contribution optimum.

### A note on `beagle.jar`

`inst/extdata/` also carries a reference copy of a Beagle 5.x release
(`beagle.28Jun21.220.jar`) alongside the package’s existing example
datasets. It is **not** auto-detected by
[`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md)
or `run_ldx_pipeline(phase = TRUE)` — both still expect `beagle.jar` to
be placed in your own `out_dir` or supplied via `beagle_jar` (see
[`vignette("HapBlockR-phasing")`](https://FAkohoue.github.io/HapBlockR/articles/HapBlockR-phasing.md)).
Copy it into your working directory and rename it to `beagle.jar` (or
point `beagle_jar` at it directly) if you want to use that bundled copy.

------------------------------------------------------------------------

## 5. Documentation

Full documentation, function reference, and tutorials:

<https://FAkohoue.github.io/HapBlockR/>

``` r
vignette("HapBlockR-breeding-decisions", package = "HapBlockR")  # start here for parent selection
vignette("HapBlockR-intro",              package = "HapBlockR")
vignette("HapBlockR-workflow",           package = "HapBlockR")
vignette("HapBlockR-phasing",            package = "HapBlockR")
vignette("HapBlockR-ld-metrics",         package = "HapBlockR")
vignette("HapBlockR-large-scale",        package = "HapBlockR")
```

**For breeders and programme managers who don’t run R.** The same eight
parent- and cross-selection tools covered in the *From Local GEBV to a
Crossing Decision* vignette
([`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md),
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md))
are also explained in plain language, with no code, in a standalone
companion document —
[`HapBlockR_Breeder_Guide.pdf`](https://raw.githubusercontent.com/FAkohoue/HapBlockR/master/inst/extdata/HapBlockR_Breeder_Guide.pdf):
what each tool is, when to reach for it, when to be cautious, and how it
works in practice, plus a programme-shape decision guide and glossary.
Distributed as a PDF (not an editable Word document) so it reaches
readers as a fixed reference rather than a document they might
inadvertently edit. Download it directly from the repository (link
above), or, once the package is installed, open it from R:

``` r
HapBlockR::open_breeder_guide()   # opens the .pdf in your default viewer
```

A `.pdf` has no vignette engine, so it isn’t indexed by
[`vignette()`](https://rdrr.io/r/utils/vignette.html) the way the `.Rmd`
vignettes above are — it ships as a static file under `inst/extdata/`
(the same convention this package already uses for `beagle.jar` and its
example datasets), and
[`open_breeder_guide()`](https://FAkohoue.github.io/HapBlockR/reference/open_breeder_guide.md)
is a thin [`system.file()`](https://rdrr.io/r/base/system.file.html)
lookup that locates and opens it.

------------------------------------------------------------------------

## 6. Quick start

``` r
library(HapBlockR)

# ── Option A: end-to-end pipeline ────────────────────────────────────────────
result <- run_ldx_pipeline(
  geno_source    = "mydata.vcf.gz",
  out_dir        = "ldx_results",
  out_blocks     = "ldx_results/blocks.csv",
  out_diversity  = "ldx_results/diversity.csv",
  out_hap_matrix = "ldx_results/hap_matrix.csv",
  CLQcut         = 0.70,
  n_threads      = 8L,
  verbose        = TRUE
)

result$blocks          # data.frame of LD blocks
result$diversity       # per-block diversity table
result$haplotypes      # named list of haplotype strings

# ── Option B: block detection only ───────────────────────────────────────────
be <- read_geno("mydata.vcf.gz")
blocks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.70, n_threads = 8L)
close_backend(be)

# ── Option C: with Beagle phasing ────────────────────────────────────────────
# Place beagle.jar in out_dir first, then:
result <- run_ldx_pipeline(
  geno_source        = "mydata.vcf.gz",
  out_dir            = "ldx_results",
  out_blocks         = "ldx_results/blocks.csv",
  out_hap_matrix     = "ldx_results/hap_matrix.csv",
  phase              = TRUE,
  beagle_jar         = "ldx_results/beagle.jar",
  beagle_threads     = 8L,
  beagle_java_mem_gb = 16L,
  beagle_seed        = 42L,
  CLQcut             = 0.70,
  n_threads          = 8L
)

result$phase_method   # "beagle"
result$phased_vcf     # path to Beagle-phased VCF.gz
```

------------------------------------------------------------------------

## 7. Input formats

[`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
accepts a path to a genotype file (or an in-memory R matrix). The format
is auto-detected from the file extension when `format` is not supplied.

### 7.1. Genotype input format

| Format             | Extension                 | `format =`  |
|--------------------|---------------------------|-------------|
| Numeric dosage     | `.csv`, `.txt`            | `"numeric"` |
| HapMap             | `.hmp.txt`                | `"hapmap"`  |
| VCF / bgzipped VCF | `.vcf`, `.vcf.gz`         | `"vcf"`     |
| SNPRelate GDS      | `.gds`                    | `"gds"`     |
| PLINK binary       | `.bed` (+ `.bim`, `.fam`) | `"bed"`     |
| R matrix           | (in-memory)               | `"matrix"`  |

**Numeric dosage format** — one row per SNP; columns `SNP`, `CHR`,
`POS`, `REF`, `ALT` followed by one column per sample, values in
`{0, 1, 2, NA}`.

**VCF** — standard VCF v4.2. Both phased (`0|1`) and unphased (`0/1`) GT
fields are accepted. Multi-allelic sites use the first ALT allele.
Missing calls (`./.`) become `NA`.

**Chromosome normalisation.** The leading prefix `chr`, `Chr`, or `CHR`
is stripped from chromosome labels at read time. Polyploid sub-genome
labels (`1A`, `2D`) are preserved verbatim.

### 7.2. Phenotype input format

The `blues` argument of
[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
and
[`prepare_gblup_inputs()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md)
accepts four formats:

**Format 1 — Named numeric vector:**

``` r
blues <- c(G001 = 4.21, G002 = 3.87, G003 = 5.14)
```

**Format 2 — Data frame, single trait:**

``` r
blues <- read.csv("blues.csv")   # any column names
res <- run_haplotype_prediction(geno, snp_info, blocks,
                                 blues = blues, id_col = "id", blue_col = "YLD")
```

**Format 3 — Data frame, multiple traits:**

``` r
res <- run_haplotype_prediction(geno, snp_info, blocks,
                                 blues     = blues_df,
                                 id_col    = "id",
                                 blue_cols = c("YLD", "DIS", "PHT"))
```

**Format 4 — Named list (different individuals per trait):**

``` r
blues <- list(
  YLD = c(G001 = 4.21, G002 = 3.87),
  DIS = c(G001 = 0.32, G003 = 0.28)
)
```

### 7.3. ID matching rules

Genotype IDs in phenotype data must match `rownames(geno_matrix)`
exactly (case-sensitive). The function takes the intersection, issues a
message for individuals present in only one source, and errors if no
common individuals are found.

### 7.4. Preparing BLUEs from raw phenotype data

Best Linear Unbiased Estimates (BLUEs) — genotype means adjusted for
design effects such as replicate or block — are typically obtained from
a mixed model fit outside HapBlockR, for example:

``` r
# Example with lme4 (single environment)
library(lme4)
m   <- lmer(YLD ~ (1|id) + (1|rep), data = field_data)
blues_lme4 <- data.frame(id  = rownames(coef(m)$id),
                          YLD = coef(m)$id[, 1])
```

------------------------------------------------------------------------

## 8. Statistical background

### 8.1. MAF filtering

Each marker’s allele frequency (AF) and minor allele frequency (MAF) are
computed as:

``` R
AF_i = (Σⱼ g_ij) / (2 n_i),    MAF_i = min(AF_i, 1 − AF_i)
```

SNPs with MAF_i \< τ_maf (default 0.05) are removed. Both operations run
in a single O(np) C++ pass by `maf_filter_cpp()`.

### 8.2. Genotype preparation

**Standard r² path:** Pearson r² via column standardisation followed by
BLAS-level matrix multiplication in `compute_r2_cpp()`.

**Kinship-adjusted rV² path:** VanRaden (2008) GRM computed via
[`AGHmatrix::Gmatrix()`](https://rdrr.io/pkg/AGHmatrix/man/Gmatrix.html),
bent and conditioned via
[`ASRgenomics::G.tuneup()`](https://rdrr.io/pkg/ASRgenomics/man/G.tuneup.html),
then **A** = **R**⁻¹ (Cholesky) or **Q** Λ⁻¹⁄² **Q**ᵀ (eigen) applied to
mean-centred genotypes before the same `compute_r2_cpp()` kernel.

### 8.3–8.5. Subsegmentation, clique detection, and block construction

See the [full pipeline walkthrough](#id_14-full-pipeline-walkthrough)
and the original Kim et al. (2018) paper for algorithmic details.

------------------------------------------------------------------------

## 9. LD metrics: r² versus rV²

|  | r² | rV² |
|----|----|----|
| **Population type** | Random mating, unrelated | Livestock, inbred lines, family-based cohorts |
| **Computational cost** | O(np) prep + O(p²) per window | O(n²p) GRM + O(n³) Cholesky + O(np) whitening |
| **RAM** | Proportional to one window | n×n GRM + whitening factor held per chromosome |
| **Block accuracy** | Slightly inflated in related populations | Correct for structured populations |
| **External dependencies** | None | AGHmatrix, ASRgenomics |

``` r
blocks_r2  <- run_Big_LD_all_chr(be, method = "r2",  CLQcut = 0.70)
blocks_rv2 <- run_Big_LD_all_chr(be, method = "rV2", CLQcut = 0.70, kin_method = "chol")
```

------------------------------------------------------------------------

## 10. Clique detection mode (CLQmode)

|  | Density | Maximal | Louvain | Leiden |
|----|----|----|----|----|
| **Algorithm** | Bron-Kerbosch + density score | Bron-Kerbosch + size score | Community detection | Community detection |
| **Complexity** | Exponential (worst case) | Exponential (worst case) | O(n log n) | O(n log n) |
| **Connected communities?** | Yes | Yes | Post-processed in HapBlockR | **Yes (formal guarantee)** |
| **WGS feasible?** | Only with `max_bp_distance` | Only with `max_bp_distance` | Yes | **Yes** |
| **Recommended for** | Chip panels (\< 100k SNPs/chr) | Sparse chip panels | Not recommended (use Leiden) | **WGS panels (2M+ SNPs)** |

``` r
# WGS recommended configuration
blocks <- run_Big_LD_all_chr(be, CLQmode = "Leiden", CLQcut = 0.80,
                              max_bp_distance = 500000L,
                              subSegmSize = 500L, leng = 50L,
                              checkLargest = TRUE, n_threads = n_threads)
```

------------------------------------------------------------------------

## 11. Haplotype analysis

### 11.1. Phase-free haplotype extraction

[`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
concatenates each individual’s allele codes (0, 1, or 2) for all SNPs
within a block into a single character string:

    Individual i, block b covering SNPs j1–j4:
      haplotype = paste0(g[i,j1], g[i,j2], g[i,j3], g[i,j4]) = "0120"

> **A note on phasing.** True gametic haplotypes require statistical
> phasing. The strings produced by
> [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
> are *diploid allele strings*, not gametic phases. Within a high-LD
> block these strings are nearly 1:1 with true haplotype classes (Calus
> et al. 2008) and are sufficient for diversity analysis and genomic
> prediction.
>
> If gametic phases are required, set `phase = TRUE` in
> [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md),
> which calls Beagle 5.x after LD block detection and caches phased
> hap1/hap2/dosage as bigmemory backends for fast restart. Place
> `beagle.jar` in `out_dir` before running.

### 11.2. Haplotype diversity metrics

[`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md)
returns four metrics per block:

**Richness (k):** number of unique haplotype strings.

**Expected heterozygosity (Hₑ):**

``` R
Hₑ = n/(n−1) · (1 − Σᵢ pᵢ²)
```

**Shannon entropy (H’):**

``` R
H' = −Σᵢ pᵢ log₂(pᵢ)
```

**Dominant haplotype frequency (f_max):** frequency of the most common
haplotype. Values near 1.0 indicate a selective sweep or strong founder
effect.

### 11.3. Haplotype feature matrix for genomic prediction

[`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md)
converts haplotype strings to a numeric dosage matrix. For each block,
the `top_n` most frequent haplotypes are selected and each individual
receives:

- **Phased data**: 0 (neither gamete), 1 (one gamete — heterozygous), or
  2 (both gametes — homozygous). True allele copy number.
- **Unphased data**: 0 (absent) or 1 (present). The value 2 is never
  produced because the two chromosomes cannot be distinguished — an
  individual homozygous for an allele and one heterozygous for it
  produce different dosage strings and are treated as distinct allele
  classes.

``` r
feat  <- build_haplotype_feature_matrix(haps, scale_features = TRUE)$matrix
G_hap <- tcrossprod(feat) / ncol(feat)   # haplotype GRM
```

### 11.4. Haplotype-based genomic prediction pipeline

[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
produces:

1.  VanRaden GRM from haplotype features
2.  GEBV for all genotyped individuals via REML-based GBLUP
3.  Per-SNP additive effects backsolved from GEBV (Tong et al. 2025)
4.  Local haplotype GEBV per block per individual
5.  Block importance table ranked by scaled Var(local GEBV)

``` r
# Single trait
res <- run_haplotype_prediction(geno, snp_info, blocks,
                                 blues = my_blues, id_col = "id", blue_col = "YLD")

# Multiple traits
res_mt <- run_haplotype_prediction(geno, snp_info, blocks,
                                    blues     = my_blues_df,
                                    id_col    = "id",
                                    blue_cols = c("YLD", "DIS"))

# Integrate GWAS evidence
qtl      <- define_qtl_regions(gwas, blocks, snp_info, p_threshold = 5e-8)
priority <- integrate_gwas_haplotypes(qtl, res, diversity = div)
priority[priority$priority_score == 3, ]
```

### 11.5. Cross-validation and prediction accuracy

[`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md)
estimates prediction accuracy via k-fold cross-validation. Returns a
`HapBlockR_cv` object with `pa_mean` (PA ± SD per trait), `pa_summary`
(per-fold), `k`, and `n_rep`.

### 11.6. Between-population comparison

[`compare_haplotype_populations()`](https://FAkohoue.github.io/HapBlockR/reference/compare_haplotype_populations.md)
computes per-block Weir-Cockerham FST, maximum allele frequency
difference, and a chi-squared test (Monte Carlo, B = 2000). `divergent`
flag: FST \> 0.1 AND p \< 0.05.

### 11.7. Haplotype network visualisation

[`plot_haplotype_network()`](https://FAkohoue.github.io/HapBlockR/reference/plot_haplotype_network.md)
draws a minimum-spanning network using
[`igraph::mst()`](https://r.igraph.org/reference/mst.html). Nodes sized
by frequency; edges weighted by Hamming distance. Optional group
colouring.

### 11.8. Multi-environment stability

[`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
runs Finlay-Wilkinson (1963) regression of per-block local GEBV
contributions against the environmental index. Returns `b` (stability
coefficient), `b_se`, `R²`, `s²d`, and `stable` flag.

### 11.9. Candidate region export

[`export_candidate_regions()`](https://FAkohoue.github.io/HapBlockR/reference/export_candidate_regions.md)
converts
[`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
output to `"bed"` (UCSC-compatible, 0-based), `"csv"`, or `"biomart"`
format.

### 11.10. Per-allele effect decomposition

[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)
aggregates per-SNP additive effects into a per-haplotype-allele effect
table: effect = sum(SNP_effect × allele_dosage).

### 11.11. Sliding-window diversity scan

[`scan_diversity_windows()`](https://FAkohoue.github.io/HapBlockR/reference/scan_diversity_windows.md)
computes He, Shannon entropy, n_eff_alleles, and dominant haplotype
frequency in sliding windows across the genome, independently of LD
block boundaries.

### 11.12. True diplotype inference

[`infer_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/infer_block_haplotypes.md)
converts raw haplotype strings into a structured per-individual,
per-block diplotype table. For phased input (from
`run_ldx_pipeline(phase = TRUE)`), `phase_ambiguous` is always `FALSE`.

### 11.13. Rare-allele collapsing

[`collapse_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/collapse_haplotypes.md)
merges alleles below `min_freq` using `"rare_to_other"`, `"nearest"`
(Hamming-based), or `"tree_based"` (UPGMA) strategies. A `label_map`
attribute records every original → collapsed mapping.

### 11.14. Cross-panel harmonisation

[`harmonize_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/harmonize_haplotypes.md)
anchors allele identity to a reference dictionary, matching by exact
string or nearest Hamming neighbour within `max_hamming`. Novel alleles
receive label `"<novel>"`. A `harmonization_report` attribute reports
match quality per block.

### 11.15. Haplotype association testing

[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
uses a unified Q+K mixed linear model with **simpleM multiple-testing
correction** (Gao et al. 2008, 2010, 2011). The model jointly corrects
for population structure and polygenic kinship:

> y = μ + α · x_hap + Σ β_k · PC_k + g + ε, g ~ MVN(0, σ²_g G)

Key parameters:

- `n_pcs = 0L` (default, EMMAX): GRM random effect absorbs all structure
- `n_pcs = k`: top-k GRM eigenvectors as additional fixed effects (Q+K
  model)
- `n_pcs = NULL`: auto-selected from the GRM scree-plot elbow
- `optimize_pcs = FALSE`: when `TRUE`, auto-selects n_pcs by fitting
  null models for k = 0..`optimize_pcs_max` and minimising the criterion
  set by `optimize_method`. More principled than the elbow heuristic for
  GWAS.
- `optimize_pcs_max = 10L`: maximum PCs evaluated when
  `optimize_pcs = TRUE`
- `optimize_method = c("bic_lambda", "bic", "lambda")`: criterion for PC
  selection (only used when `optimize_pcs = TRUE`):
  - `"bic"` — minimise BIC of the null REML model
  - `"lambda"` — minimise \|λ_GC − 1\| (genomic control calibration)
  - `"bic_lambda"` (**default, recommended**) — hybrid: \|λ−1\| +
    0.01·scaled_BIC. Minimises inflation/deflation while BIC breaks ties
    toward fewer PCs
- `sig_metric`: which p-value drives the `significant` flag:
  - `"p_wald"` — raw Wald p-value (use with a pre-corrected threshold)
  - `"p_fdr"` — Benjamini-Hochberg FDR (recommended for discovery)
  - `"p_simplem"` — simpleM Bonferroni-style: min(p × Meff, 1)
  - `"p_simplem_sidak"` — simpleM Šidák-style: 1−(1−p)^Meff
    (**recommended** for correlated haplotype predictors)
- `meff_scope`: scope for estimating the effective number of independent
  tests (Meff):
  - `"chromosome"` (**recommended**) — separate Meff per chromosome
  - `"global"` — one genome-wide Meff
  - `"block"` — one Meff per LD block (omnibus tests get Meff = 1)
- `meff_percent_cut`: variance threshold for simpleM eigendecomposition
  (default `0.995`)
- `plot = TRUE`: saves **PDF** plots (not PNG). Three files per run:
  `manhattan_<trait>.pdf`, `qq_<trait>.pdf`, `pca_grm.pdf` (GRM PCA
  coloured by phenotype), `grm_scree.pdf` (eigenvalue scree with
  selected PC in red)

All four p-value columns (`p_wald`, `p_fdr`, `p_simplem`,
`p_simplem_sidak`) are **always present** in every output regardless of
`sig_metric`.

[`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
also now supports the full correction set with a new `sig_metric`
parameter. All four p-value columns (`p_omnibus_adj`, `p_omnibus_fdr`,
`p_omnibus_simplem`, `p_omnibus_simplem_sidak`) are always present in
`$omnibus_tests` regardless of the chosen `sig_metric`.

``` r
# FDR-based discovery (EMMAX, default)
assoc_fdr <- test_block_haplotypes(
  haps, blues = blues_vec, blocks = blocks,
  sig_metric = "p_fdr", verbose = FALSE
)

# simpleM Šidák with Q+K correction, chromosome-wise Meff
assoc_sm <- test_block_haplotypes(
  haps, blues = blues_vec, blocks = blocks,
  n_pcs            = 3L,
  sig_metric       = "p_simplem_sidak",
  meff_scope       = "chromosome",
  meff_percent_cut = 0.995,
  verbose          = FALSE
)
assoc_sm$meff$trait$allele$chromosome   # per-chromosome effective test counts
head(assoc_sm$block_tests[order(assoc_sm$block_tests$p_omnibus), ])
```

[`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
decomposes variation at each block into additive (a) and dominance (d)
components. Dominance ratio d/a classifies gene action: 0 = additive, ±1
= complete dominance, \|d/a\| \> 1 = overdominance.

### 11.16. Breeding decision tools

[`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md)
produces a genome-wide stacking index for each individual (sum of
allele_effect × dosage across blocks, normalised to \[0, 1\]).

[`summarize_parent_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/summarize_parent_haplotypes.md)
produces a long-format allele inventory for candidate parents: which
alleles they carry, at what dosage, and with what population frequency.
`is_rare = TRUE` flags alleles with frequency \< 10%.

### 11.17. Cross-population effect concordance

[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)
takes two
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
result objects from independent populations and computes per-block
statistics quantifying how consistently haplotype allele effects
replicate across panels.

**Key output columns (one row per block per trait):**

| Column | What it answers |
|----|----|
| `n_shared_alleles` | Alleles tested in both populations |
| `effect_correlation` | Pearson r of per-allele effects across populations (NA when \< 3 shared alleles) |
| `direction_agreement` | Fraction of shared alleles with the same effect sign |
| `directionally_concordant` | `direction_agreement >= direction_threshold` (default 0.75) |
| `meta_effect` / `meta_SE` / `meta_p` | IVW meta-analytic effect (inverse-variance weighted, same as two-sample MR) |
| `Q_stat` / `Q_p` | Cochran Q: significant Q means effect sizes differ between populations |
| `I2` | I² inconsistency (0–100%); \> 50% = substantial heterogeneity |
| `replicated` | `TRUE` when directionally concordant AND Q_p \> 0.05 AND enough shared alleles |
| `boundary_overlap_ratio` | **Automatically computed output** (not a user-set value): bp(intersection) / bp(union) of the two populations’ block boundaries. Requires `blocks_pop1` / `blocks_pop2` to be supplied; `NA` otherwise |
| `boundary_warning` | `TRUE` when `boundary_overlap_ratio < boundary_overlap_warn` (the **input parameter**, default `0.80`) |
| `match_type` | How the block was matched: `"exact"` (same `block_id` string), `"position"` (matched by genomic overlap), `"pop1_only"` (no Pop2 block overlaps at `overlap_min`), or `NA` when no block tables were supplied |

``` r
# Use Pop A's block boundaries for Pop B (maximises shared alleles)
haps_B_harm <- harmonize_haplotypes(
  extract_haplotypes(geno_B, snp_info, blocks_A),   # same block coords
  reference = haps_A
)
assoc_A <- test_block_haplotypes(haps_A, blues = blues_A, blocks = blocks_A,
                                  sig_metric = "p_simplem_sidak")
assoc_B <- test_block_haplotypes(haps_B_harm, blues = blues_B, blocks = blocks_A,
                                  sig_metric = "p_simplem_sidak")

# block_match = "id" (default): match by block_id string — fast, backward-compatible
# block_match = "position":      match by genomic interval overlap — handles different
#                                 LD block boundaries between populations
conc <- compare_block_effects(
  assoc_A, assoc_B,
  pop1_name             = "PopA",
  pop2_name             = "PopB",
  blocks_pop1           = blocks_A,   # required for boundary_overlap_ratio + position matching
  blocks_pop2           = blocks_B,   # supply Pop B's OWN block table (may differ from A)
  block_match           = "position", # recommended when blocks differ between populations
  overlap_min           = 0.50,       # min IoU for two blocks to be considered the same region
  direction_threshold   = 0.75,
  boundary_overlap_warn = 0.80
)
# $concordance$match_type: "exact" | "position" | "pop1_only"
# "position" rows: same QTL region, different boundary definitions
conc$concordance[conc$concordance$replicated, ]   # replicated blocks
head(conc$shared_alleles)                          # per-allele IVW detail
print(conc)
```

------------------------------------------------------------------------

### 11.18. Cross-population effect concordance (external GWAS)

When GWAS was run externally (GAPIT, TASSEL, FarmCPU, PLINK, or any
other tool),
[`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md)
compares block-level lead-SNP effects between two populations. It
accepts either raw GWAS data frames or the pre-mapped output of
[`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md).

**Key differences from
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md):**

|  | [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md) | [`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md) |
|----|----|----|
| Input | [`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md) results | GWAS tables or [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md) output |
| Unit per block | Multiple haplotype alleles | One lead SNP |
| `effect_correlation` | Pearson r across alleles | Always NA |
| `direction_agreement` | Fraction of alleles with same sign | 0 or 1 |
| `Q_stat` / `I2` | Cochran Q, df = n_alleles − 1 | Always NA (df = 0) |
| `replicated` | dir_concordant AND Q_p \> 0.05 | dir_concordant AND meta_p ≤ 0.05 |

**SE derivation.** Many GWAS tools omit the standard error. When `SE` is
absent,
[`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md)
derives it from the z-score: `SE = |BETA| / |z|`, `z = Φ⁻¹(P/2)`. The
`se_derived_pop1` / `se_derived_pop2` output columns flag which
population required this derivation.

``` r
# Path 1: pre-mapped (recommended) — most auditable
qtl_A <- define_qtl_regions(gwas_A, blocks, snp_info, p_threshold = 5e-8)
qtl_B <- define_qtl_regions(gwas_B, blocks, snp_info, p_threshold = 5e-8)

conc <- compare_gwas_effects(
  qtl_pop1      = qtl_A,
  qtl_pop2      = qtl_B,
  blocks_pop1   = blocks_A,       # Pop A's block table
  blocks_pop2   = blocks_B,       # Pop B's block table (may differ from A)
  block_match   = "position",     # match by genomic overlap — handles different boundaries
  overlap_min   = 0.50,           # min IoU for a valid match
  pop1_name     = "PopA",
  pop2_name     = "PopB"
)

# Path 2: raw GWAS + blocks (convenience — calls define_qtl_regions internally)
conc2 <- compare_gwas_effects(
  gwas_pop1     = gwas_A,
  gwas_pop2     = gwas_B,
  blocks_pop1   = blocks_A,
  blocks_pop2   = blocks_B,
  snp_info_pop1 = snp_info,
  # snp_info_pop2 = NULL: reuses snp_info_pop1 when marker panels are shared
  pop1_name     = "PopA",
  pop2_name     = "PopB",
  block_match   = "position",  # match blocks by genomic interval overlap
  overlap_min   = 0.50,        # blocks with IoU < 0.5 become "pop1_only"
  p_threshold   = 5e-8,
  beta_col      = "BETA",      # column name — change to match your GWAS tool
  se_col        = "SE",        # NULL or absent: derived from z-score automatically
  p_col         = "P"
)

# Output uses the same HapBlockR_effect_concordance class
# GWAS-specific extra columns in $concordance:
#   lead_snp_pop1, lead_snp_pop2  — which SNP tagged the block in each pop
#   lead_p_pop1, lead_p_pop2      — lead SNP p-values
#   se_derived_pop1/pop2          — TRUE when SE was derived from z-score
#   both_pleiotropic              — TRUE when block is pleiotropic in both pops
conc$concordance[conc$concordance$replicated, ]
print(conc)
```

> **Shared marker panel.** When both populations were genotyped on the
> same array or sequenced with the same reference, set
> `snp_info_pop2 = NULL` (default) to reuse `snp_info_pop1`. When
> populations have different marker sets, supply `snp_info_pop2`
> explicitly.

### 11.19. Within-block and between-block epistasis detection

HapBlockR provides three functions for epistasis detection that extend
the haplotype association framework. All operate on GRM-corrected REML
residuals from the same null model as
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md),
ensuring population-structure-corrected tests throughout.

**[`scan_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_epistasis.md)**
tests all C(p,2) SNP pairs within each significant block for pairwise
interaction on GRM-corrected REML residuals. The model for each pair is:

> y = μ + aᵢxᵢ + aⱼxⱼ + aaᵢⱼ(xᵢ × xⱼ) + ε

Restricting to significant blocks avoids the genome-wide explosion: for
15 significant blocks with ~200 SNPs each, the total number of tests is
~300,000 rather than ~4.4 billion. Multiple testing is corrected by
Bonferroni and simpleM Sidak within each block, where Meff is estimated
from the eigenspectrum of the pairwise interaction column matrix.

**[`scan_block_by_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_by_block_epistasis.md)**
is a trans-haplotype epistasis scan that tests each significant
haplotype allele against every allele at all other blocks. This is
conceptually analogous to a trans-eQTL scan but for phenotypic haplotype
interactions — it detects genetic background dependence where a
resistance haplotype at one locus only functions in the presence of a
specific background at another locus. For 25 significant alleles ×
17,943 total alleles, the scan involves ~450,000 tests corrected by
Bonferroni.

**[`fine_map_epistasis_block()`](https://FAkohoue.github.io/HapBlockR/reference/fine_map_epistasis_block.md)**
fine-maps a single block by identifying the specific interacting SNP
pairs. For blocks with p ≤ 200 SNPs it runs an exhaustive pairwise scan.
For larger blocks it uses LASSO with pairwise interaction terms
([`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
`lambda.1se`), which avoids the multiple-testing burden while
identifying the most influential pairs.

``` r
# Within-block epistasis scan (significant blocks only)
epi_within <- scan_block_epistasis(
  assoc              = assoc,          # test_block_haplotypes() result
  geno_matrix        = res$geno_matrix,
  snp_info           = snp_info,
  blocks             = blocks,
  blues              = blues_list,
  haplotypes         = haps,
  trait              = "BL",
  sig_blocks         = NULL,           # NULL = use significant_omnibus blocks
  max_snps_per_block = 300L,
  sig_metric         = "p_simplem_sidak",
  sig_threshold      = 0.05
)
print(epi_within)
epi_within$results[epi_within$results$significant, ]
epi_within$scan_summary

# Between-block trans-haplotype epistasis scan
epi_between <- scan_block_by_block_epistasis(
  assoc         = assoc,
  haplotypes    = haps,
  blues         = blues_list,
  blocks        = blocks,
  trait         = "BL",
  sig_alleles   = NULL,    # NULL = significant alleles from assoc
  sig_threshold = 0.05
)
print(epi_between)
epi_between$results[epi_between$results$significant, ]

# Fine-map a single block (auto-dispatches: pairwise or LASSO)
fine <- fine_map_epistasis_block(
  block_id   = "block_12_1054210_1086071",
  geno_matrix = res$geno_matrix,
  snp_info    = snp_info,
  blocks      = blocks,
  y_resid     = my_reml_residuals,   # from .fit_null_reml() or null model
  method      = "auto",              # pairwise <= 200 SNPs, lasso otherwise
  sig_threshold = 0.05
)
head(fine)
```

------------------------------------------------------------------------

## 12. Parameter auto-tuning

[`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md)
selects `CLQcut` (and optionally other parameters) minimising in
priority order: unassigned GWAS markers, forced assignments, number of
blocks, deviation from target median block size.

``` r
result <- tune_LD_params(
  geno_matrix    = my_geno,
  snp_info       = my_snp_info,
  gwas_df        = my_gwas,
  prefer_perfect = TRUE,
  target_bp_band = c(5e4, 5e5),
  parallel       = FALSE,
  seed           = 42L
)
result$best_params
result$gwas_assigned
```

------------------------------------------------------------------------

## 13. Scale strategies and backends

### 13.1. The HapBlockR_backend interface

``` r
be    <- read_geno("mydata.bed")       # opens PLINK BED
chunk <- read_chunk(be, col_idx)       # n_samples × length(col_idx) matrix
close_backend(be)                      # release file handle
```

### 13.2. Memory requirements by configuration

| Configuration                                | Peak RAM                   |
|----------------------------------------------|----------------------------|
| Plain matrix, r², all in RAM                 | ~40 GB (500 ind, 10M SNPs) |
| GDS or BED backend, r², `subSegmSize = 1500` | ~300 MB                    |
| Any backend, rV², `method = "rV2"`           | ~4 GB                      |
| bigmemory, r², `subSegmSize = 500`           | ~0.8 MB/window             |

### 13.3. Recommended configurations by dataset size

| Markers | Individuals | Recommended configuration |
|----|----|----|
| \< 100 k | any | `method = "r2"`, `format = "matrix"` |
| 100 k – 500 k | \< 5,000 | `method = "rV2"` for structured populations |
| 100 k – 2 M | any | `method = "r2"`, `format = "vcf"` or `"bed"` |
| 2 M – 10 M | any | `CLQmode = "Leiden"`, `max_bp_distance = 500000L`, `format = "gds"` |
| \> 10 M | any | `CLQmode = "Leiden"`, `max_bp_distance = 500000L`, [`read_geno_bigmemory()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno_bigmemory.md) |

------------------------------------------------------------------------

## 14. Full pipeline walkthrough

| Step | Action | Key parameter(s) |
|----|----|----|
| 1 | Accept genotype backend or wrap plain matrix | [`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md), `format =` |
| 2 | Extract per-chromosome genotype slice | `read_chunk(backend, chr_idx)` |
| 3 | MAF filter + monomorphic removal in C++ | `MAFcut`, `maf_filter_cpp()` |
| 4 | Centre (r²) or centre + whiten (rV²) | `method`, `kin_method`, [`prepare_geno()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_geno.md) |
| 4b | **Optional: Beagle phasing** (`phase = TRUE`). Calls [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md) on the original VCF after imputation/filtering. Reads the phased VCF via [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md), aligns samples and SNPs via CHR+POS+REF+ALT+SNP composite key, and caches hap1/hap2/dosage as bigmemory backends (`hapblockr_bm_phased_hap1`, `_hap2`, `_dos`) in `bigmemory_path`. On restart the VCF is **not re-read** — backends reattach from disk via fingerprint match (VCF path + mtime + n_snps). Haplotype extraction then uses gametic strings (`g1\|g2`). | `phase`, `beagle_jar`, `beagle_threads`, `beagle_java_mem_gb`, `beagle_seed`, `beagle_ref_panel`, `beagle_map_file` |
| 5 | C++ boundary scan — find weak-LD cut points | `leng`, `subSegmSize`, `boundary_scan_cpp()` |
| 6 | Divide chromosome into sub-segments | `subSegmSize` |
| 7 | Per sub-segment: compute r² or rV² matrix in C++ | `CLQcut`, `compute_r2_cpp()` |
| 8 | Build binary adjacency matrix in C++ | `CLQcut`, `build_adj_matrix_cpp()` |
| 9 | Find communities/cliques | `CLQmode`, `checkLargest`, `max_bp_distance` |
| 10 | Greedy clique assignment → bin vector | `split`, `clstgap`, [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md) |
| 11 | MWIS block construction | internal |
| 12 | Re-merge across forced cut-points | automatic |
| 13 | Merge overlapping blocks | automatic |
| 14 | Map indices to bp position and rsID | `SNPinfo` |
| 15 | Re-index over full SNP set including monomorphics | automatic |
| 16 | Optionally append rare SNPs | `appendrare` |

**Detailed example with all parameters:**

``` r
blocks <- run_Big_LD_all_chr(
  # ── Genotype input ──────────────────────────────────────────────────────────
  geno_matrix  = be,

  # ── LD metric ───────────────────────────────────────────────────────────────
  method       = "r2",          # "r2" (default) or "rV2"
  kin_method   = "chol",        # "chol" (default) or "eigen" for rV2

  # ── Clique detection ────────────────────────────────────────────────────────
  CLQcut       = 0.70,
  CLQmode      = "Density",     # "Density"/"Maximal"/"Louvain"/"Leiden"
  clstgap      = 40000L,
  split        = FALSE,

  # ── Subsegmentation ─────────────────────────────────────────────────────────
  leng         = 200L,
  subSegmSize  = 1500L,

  # ── Filtering ───────────────────────────────────────────────────────────────
  MAFcut       = 0.05,
  appendrare   = FALSE,

  # ── Large window heuristics ─────────────────────────────────────────────────
  checkLargest = FALSE,

  # ── Parallelism ─────────────────────────────────────────────────────────────
  n_threads    = 8L,
  digits       = -1L,

  # ── Chromosome minimum ──────────────────────────────────────────────────────
  min_snps_chr = 10L,

  # ── Reproducibility ─────────────────────────────────────────────────────────
  seed         = 42L,
  verbose      = TRUE
)
```

**Full
[`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)
example with phasing:**

``` r
result <- run_ldx_pipeline(
  # ── Genotype input ──────────────────────────────────────────────────────────
  geno_source    = "mydata.vcf.gz",   # VCF required for phase = TRUE
  out_dir        = "ldx_results",

  # ── Output paths ────────────────────────────────────────────────────────────
  out_blocks     = "ldx_results/blocks.csv",
  out_diversity  = "ldx_results/diversity.csv",
  out_hap_matrix = "ldx_results/hap_matrix.csv",
  hap_format     = "numeric",

  # ── Phasing (optional — requires VCF input and beagle.jar in out_dir) ───────
  phase              = TRUE,
  beagle_jar         = "ldx_results/beagle.jar",
  beagle_threads     = 8L,
  beagle_java_mem_gb = 16L,          # -Xmx16g; NULL = JVM default
  beagle_seed        = 42L,          # integer seed for reproducibility
  beagle_ref_panel   = NULL,         # phased reference VCF (optional)
  beagle_map_file    = NULL,         # genetic map for improved accuracy

  # ── Filtering & imputation ───────────────────────────────────────────────────
  maf_cut        = 0.05,
  impute         = "mean_rounded",

  # ── LD block detection ───────────────────────────────────────────────────────
  CLQcut         = 0.70,
  method         = "r2",
  CLQmode        = "Leiden",
  subSegmSize    = 500L,
  leng           = 50L,
  max_bp_distance = 500000L,
  n_threads      = 8L,

  # ── Haplotype extraction ─────────────────────────────────────────────────────
  min_snps_block = 3L,
  top_n          = 5L,

  # ── Bigmemory caching (restart-safe) ─────────────────────────────────────────
  use_bigmemory  = TRUE,
  bigmemory_path = "ldx_results/bm_cache",
  bigmemory_type = "char",

  # ── General ──────────────────────────────────────────────────────────────────
  verbose = TRUE
)
```

------------------------------------------------------------------------

## 15. Function reference

### 15.1. Main pipeline

| Function | Description |
|----|----|
| [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md) | **Recommended end-to-end pipeline.** One call from genotype source to haplotype matrix, diversity table, and block CSV. Accepts `phase = TRUE` for Beagle phasing (place `beagle.jar` in `out_dir`). Phased data cached as bigmemory backends for fast restart. |
| [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md) | Chromosome-wise LD block detection. Accepts both plain matrices and `HapBlockR_backend` objects. |
| [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md) | Core per-chromosome segmentation. Called internally by [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md). Not exported. |
| [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md) | Grid-search auto-tuner minimising unassigned GWAS marker placements. |

### 15.2. I/O

| Function | Description |
|----|----|
| [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md) | See Main pipeline. Also handles all I/O for the end-to-end run. |
| [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md) | Statistical phasing via Beagle 5.x. Called internally by `run_ldx_pipeline(phase = TRUE)`. Exposes all Beagle parameters: `java_mem_gb`, `map_file`, `chrom`, `seed`, `burnin`, `iterations`, `window`, `overlap`. Log written to `out_prefix.log`. Requires `beagle.jar`. |
| [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md) | Read a phased VCF (pipe-separated `0\|1` GT fields) into `list(hap1, hap2, dosage, snp_info, sample_ids, phased=TRUE)`. Called internally after Beagle phasing. |
| [`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md) | Auto-dispatch genotype reader. Returns an `HapBlockR_backend` object. |
| [`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md) | Extract a genotype slice (n_samples × width) from any backend type. |
| [`close_backend()`](https://FAkohoue.github.io/HapBlockR/reference/close_backend.md) | Release file handles. No-op for in-memory backends. |
| [`read_geno_bigmemory()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno_bigmemory.md) | Create a file-backed memory-mapped store from any source. |

### 15.3. LD computation

| Function | Description |
|----|----|
| [`compute_r2()`](https://FAkohoue.github.io/HapBlockR/reference/compute_r2.md) | Standard r² matrix via C++ Armadillo + optional OpenMP. |
| [`compute_rV2()`](https://FAkohoue.github.io/HapBlockR/reference/compute_rV2.md) | rV² on a pre-whitened matrix. |
| [`prepare_geno()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_geno.md) | Centre (r²) or centre + whiten (rV²). |
| [`get_V_inv_sqrt()`](https://FAkohoue.github.io/HapBlockR/reference/get_V_inv_sqrt.md) | Whitening factor **A** such that AVA’ = I. |

### 15.4. C++ kernels (direct access)

| Function | Description |
|----|----|
| `compute_r2_cpp()` | Full r² matrix. OpenMP outer loop. |
| `maf_filter_cpp()` | MAF + monomorphic filter in one O(np) C++ pass. |
| `build_adj_matrix_cpp()` | LD threshold → 0/1 integer adjacency matrix. |
| `col_r2_cpp()` | r² of one query column against all others. |
| `compute_r2_sparse_cpp()` | Sparse r² for pairs within a bp distance window. |
| `boundary_scan_cpp()` | Cross-boundary LD scan. Returns 0/1 vector of valid cut positions. |
| `resolve_overlap_cpp()` | BLAS DGEMM scoring of overlapping block boundaries. |

### 15.5. Haplotype analysis

| Function | Description |
|----|----|
| [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md) | Phase-free diploid allele strings per block × individual. For phased input (from `run_ldx_pipeline(phase = TRUE)`), produces gametic strings `g1\|g2`. |
| [`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md) | Per-block richness, He, n_eff_alleles, Shannon, sweep_flag. |
| [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md) | Haplotype dosage matrix for genomic prediction. Phased: 0/1/2. Unphased: 0/1. |
| [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md) | VanRaden GRM from haplotype feature matrix. |
| [`compute_dominance_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md) | Vitezica et al. (2013) dominance relationship matrix from a 0/1/2 genotype matrix. Does not fit a dual-kernel G_A+G_D model itself — pass the result to `sommer`/`ASReml-R` alongside [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)’s additive GRM. |
| [`decode_haplotype_strings()`](https://FAkohoue.github.io/HapBlockR/reference/decode_haplotype_strings.md) | Decode dosage strings to nucleotide sequences. |
| [`write_haplotype_numeric()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_numeric.md) | Write haplotype dosage matrix to file. |
| [`write_haplotype_character()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_character.md) | Write nucleotide character matrix to file. |
| [`write_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_diversity.md) | Write diversity table to CSV. |
| [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md) | Map GWAS hits to LD blocks; detect pleiotropic blocks. |
| [`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md) | Derive per-SNP effects from GEBV (Tong et al. 2025). |
| [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md) | Local haplotype GEBV per block per individual. |
| [`prepare_gblup_inputs()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md) | Align phenotype data and haplotype GRM for external GBLUP solvers. |
| [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md) | Single or multi-trait Tong et al. (2025) haplotype stacking pipeline. Optional `include_dominance = TRUE` (requires `marker_effect_method = "gblup"` and `BGLR`) fits a dual-kernel additive+dominance GBLUP (Vitezica et al. 2013) and adds `dominance_deviation`/`total_genetic_value`/`G_dominance` to the output; default `FALSE` leaves additive-only behaviour unchanged. |
| [`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/integrate_gwas_haplotypes.md) | Combine GWAS, variance, and diversity evidence per block. |
| [`rank_haplotype_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/rank_haplotype_blocks.md) | Unified block ranking across 3 use cases. |

### 15.6. Analysis extensions

| Function | Description |
|----|----|
| [`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md) | K-fold cross-validation for the haplotype GBLUP model. |
| [`compare_haplotype_populations()`](https://FAkohoue.github.io/HapBlockR/reference/compare_haplotype_populations.md) | Per-block Weir-Cockerham FST between two sample groups. |
| [`plot_haplotype_network()`](https://FAkohoue.github.io/HapBlockR/reference/plot_haplotype_network.md) | Minimum-spanning network of haplotype alleles for one block. |
| [`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md) | Finlay-Wilkinson stability regression across environments. |
| [`export_candidate_regions()`](https://FAkohoue.github.io/HapBlockR/reference/export_candidate_regions.md) | Convert QTL regions to BED, CSV, or biomaRt format. |
| [`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md) | Per-haplotype-allele effect table from per-SNP additive effects. |
| [`scan_diversity_windows()`](https://FAkohoue.github.io/HapBlockR/reference/scan_diversity_windows.md) | Sliding-window diversity scan across the genome. |

### 15.7. True haplotype inference and harmonisation

| Function | Description |
|----|----|
| [`infer_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/infer_block_haplotypes.md) | Structured diplotype table from raw haplotype strings. Phased input: `phase_ambiguous = FALSE`. |
| [`collapse_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/collapse_haplotypes.md) | Merge rare alleles via `"rare_to_other"`, `"nearest"`, or `"tree_based"`. |
| [`harmonize_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/harmonize_haplotypes.md) | Cross-panel allele label harmonisation using reference dictionary. |

### 15.8. Haplotype association testing

| Function | Description |
|----|----|
| [`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md) | Block-level association tests via Q+K mixed model with simpleM multiple-testing correction (Gao et al. 2008, 2010, 2011). Parameters `sig_metric`, `meff_scope`, `meff_percent_cut`, `meff_max_cols` control correction method and Meff estimation scope. All four p-value flavours always present in output. Returns `HapBlockR_haplotype_assoc`. |
| [`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md) | Additive (a) and dominance (d) effects per block. Returns `HapBlockR_diplotype`. |
| [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md) | Cross-population haplotype effect concordance. Computes IVW meta-analytic effects, Cochran Q heterogeneity, I² inconsistency, direction agreement, and replication flag per block. New `block_match` parameter: `"id"` (default, backward-compatible) matches by `block_id` string; `"position"` matches by genomic interval overlap (IoU ≥ `overlap_min`, default 0.50) — handles different LD block boundaries between populations. `boundary_overlap_ratio` is automatically computed (not user-set); `boundary_overlap_warn` (default 0.80) controls `boundary_warning`. Output `match_type` column: `"exact"` / `"position"` / `"pop1_only"`. Returns `HapBlockR_effect_concordance`. |
| [`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md) | Cross-population concordance from **external GWAS results** (GAPIT, TASSEL, FarmCPU, PLINK, etc.). Accepts raw GWAS data frames or pre-mapped [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md) output. Derives SE from z-score when absent. `block_match = "position"` handles different block boundaries between populations (same as [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)). One lead SNP per block; `effect_correlation` and Cochran Q are NA. Same output class as [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md). |
| [`scan_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_epistasis.md) | Within-block pairwise SNP epistasis scan. Tests all C(p,2) SNP pairs within significant blocks on GRM-corrected REML residuals. Model: y = mu + a_i*x_i + a_j*x_j + aa_ij*(x_i*x_j) + e. Corrected via Bonferroni and simpleM Sidak within each block (Meff from interaction column eigenspectrum). `sig_metric` controls which correction drives the `significant` flag. Returns `HapBlockR_epistasis`. |
| [`scan_block_by_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_by_block_epistasis.md) | Trans-haplotype between-block epistasis scan. Tests significant haplotype alleles against all other block alleles: O(n_sig x n_total_alleles) tests corrected by Bonferroni. Identifies genetic background dependence that single-block analyses cannot detect. Returns `HapBlockR_block_epistasis`. |
| [`fine_map_epistasis_block()`](https://FAkohoue.github.io/HapBlockR/reference/fine_map_epistasis_block.md) | Single-block epistasis fine-mapping. Dispatches to exhaustive pairwise scan (p \<= 200 SNPs) or LASSO interaction search via `glmnet` (p \> 200 SNPs). Requires pre-computed REML residuals (`y_resid`). |

### 15.9. Breeding decision tools

See the *From Local GEBV to a Crossing Decision* vignette for how these
fit together end to end: local GEBV → block triage → parent selection →
forward simulation → crossing shortlist.

| Function | Description |
|----|----|
| [`select_top_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md) | Smallest top-ranked set of blocks by fixed count, fixed percentage of blocks, or cumulative share of variance explained. |
| [`plot_block_funnel()`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md) | Funnel plot of per-block local GEBV vs. scaled variance, highlighting important blocks. |
| [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md) | GA search for the founder set (`GA::ga(type = "binary")`) jointly maximising coverage of favourable per-block values, under five crossing-scheme strategies (`no_selfing`/`selfing`/`OHS`/`OPV`/`Haploid_OHS`). Optional `coancestry_weight` soft-penalises mean relatedness of the chosen set (requires `G`); optional `merit_weight` (GA+TS hybrid) adds whole-genome merit directly into the fitness function, alongside the pre-search `min_sel_value` floor, so merit keeps pulling on the search throughout rather than only at its entry gate — guards against selecting on noisy per-block estimates. `merit_priority` (0-100) is the easy alternative to a raw `merit_weight` — auto-calibrated from your data (see [`suggest_merit_weight()`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)). `target_degree` (0-90) is the easy alternative to a raw `coancestry_weight` — the same dial [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md) uses, converting your preference into a relatedness ceiling interpolated between a greedy high-coverage reference group and `select_core_collection(strategy = "maximin")`’s diversity-maximising one, enforced by a self-scaled penalty (no separate tuning constant needed). |
| [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md) | Top-n-by-score baseline, used as the comparison arm for [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md). |
| [`suggest_merit_weight()`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md) | Suggests a starting `merit_weight` for [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md) from your own data — estimates each term’s realistic best-vs-worst achievable spread (coverage’s reference point built via a greedy, provably near-optimal reachable group, not an inflated unreachable ceiling; both references trimmed against single-outlier distortion) and returns the ratio, optionally converted into a literal `merit_weight` for a given `merit_priority` (0-100). Same calculation [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s own `merit_priority` argument uses internally — call this first to inspect the numbers before committing. |
| [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md) | Family-quota parent selection: the `n_families` best-performing families (ranked by each family’s own top-`n_per_family` mean), then the `n_per_family` best lines within each — family balance enforced directly by the selection rule rather than checked afterward. Optional `ensure_haplotype_diversity` adjusts within-family picks so representatives of different chosen families avoid sharing the same dominant target block, with an allele-level fallback auto-derived from `haplotypes` (the direct [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md) output — no manual matrix construction needed), so crossing across families is more likely to combine different favourable haplotypes rather than duplicate one. |
| [`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md) | PCA of a genomic/haplotype relationship matrix, coloured by GA-selected/TS-selected/both/neither, for a diversity-space sanity check of either selection. |
| [`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md) | Forward-in-time recurrent-selection simulation (wraps `genomicSimulation`) comparing founder sets/mating schemes over generations. Requires phased haplotypes. Legacy `ga_selected`/`ts_selected` args still give the original 2-scheme (GA vs. TS) `"truncation"`-mating comparison; a new `schemes` argument generalises to any named list of schemes, each with its own `mating_scheme` — `"truncation"` (unchanged), `"ocs"` (per-generation [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)-chosen mating plan), or `"uc"` (per-generation [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)-ranked top crosses; requires `blocks`). |
| [`plot_ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md) | Realised mean/max breeding-population GEBV per generation, for every scheme in the simulation (2 or more). |
| [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md) | Usefulness Criterion (UC, Schnell & Utz 1975; Bernardo 2003) ranking candidate **crosses** by predicted mid-parent value plus selection-intensity-scaled segregation variance. Four `variance_model` modes: `"block_independent"` (unphased, local GEBV), `"phased"` (exact 4-gamete enumeration, target blocks summed independently), `"linked"` (native, no extra dependency: Monte Carlo progeny simulation via Haldane-mapped block-to-block recombination, closing `"phased"`’s independent-blocks gap; requires a `genetic_map`), `"simplemating"` (genuine multi-locus linkage-aware variance via [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html), requires 0/2-coded DH/RIL genotypes). Optional finite-population selection-intensity correction (`n_progeny`, Monte Carlo). |
| [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md) | Optimal Contribution Selection / mate allocation under an explicit relatedness cap (`target_degree`). Three engines: `"alphamate"` (AlphaMate executable, true OCS via full evolutionary-algorithm allocation), `"optisel"` (true OCS via optiSel’s own `candes()`/`opticont()`/`matings()` solver, Meuwissen 1997, Wellmann 2019, no external binary), and `"simplemating"` (not true OCS — discrete greedy cross prediction/selection via SimpleMating’s `planCross()`/`selectCrosses()`, no external binary; honors `max_contrib_per_parent` directly and guarantees no repeated matings by construction). |
| [`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md) | Sweeps [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s `coancestry_weight` across a grid and Pareto-filters the results into an explicit merit-vs-relatedness tradeoff frontier, via [`pareto_front()`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md). |
| [`pareto_front()`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md) | General non-dominated-sort utility (with NSGA-II crowding distance) for any data frame of objective columns — not specific to parent selection. |
| [`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md) | Exact binary integer linear programming (`lpSolve`) solution to the cross-selection problem, for sanity-checking a heuristic mating plan’s optimality gap on a small-enough candidate set (`max_vars`). |
| [`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md) | Core-collection / diversity-maximising subset selection (Schoen & Brown 1993) via farthest-point maximin greedy traversal (Gonzalez 1985) or a mean-distance variant — diversity itself as the objective, not a constraint, with an optional merit floor. |
| [`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md) | Genome-wide haplotype stacking index per individual. |
| [`summarize_parent_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/summarize_parent_haplotypes.md) | Long-format allele inventory for candidate parents. |

**Choosing between
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
and
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md).**
These answer genuinely different questions, and neither is a strict
improvement on the other — run both and compare rather than picking one
on faith.

|  | [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md) | [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md) |
|----|----|----|
| Optimises for | A single whole-genome score, top-n | Joint coverage of favourable value across a set of target blocks |
| Ignores | Which haplotype blocks each individual carries — its top-n list can concentrate on the same favourable blocks and leave others uncovered | Whole-genome merit, unless `min_sel_value` (a hard pre-search floor) and/or `merit_weight` (a soft, continuous term rewarded throughout the search) are set (without either, a poor overall performer who uniquely covers one target block will still be selected) |
| Output | Flat top-n list | Flat n-founder set (no differential contributions, unless `coancestry_weight`/`target_degree` is used to soft-penalise relatedness) |
| Use when | You want the simplest defensible shortlist, or a baseline to beat | You want to check whether “best by value” is actually “best coverage of what matters,” and to engineer complementarity |

**A third shortlist strategy: family (or genetic-cluster) quotas.**
[`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)
sits alongside both of the above rather than replacing either — it
answers neither “best by value” nor “best joint coverage,” but “our best
few groups, and our best few lines out of each,” a shape many real
programme shortlists already take. A “group” is either your own
pedigree/family labels (`group_by = "family"`, the default) or a
data-derived genetic cluster built directly from a relationship matrix
(`group_by = "genetic_cluster"`, hierarchical Ward’s-linkage clustering
— useful when pedigree labels aren’t a fully faithful proxy for actual
relatedness). It enforces group balance directly in the selection rule
(a fixed, or per-group uneven, quota) instead of checking it as a
diagnostic afterward; groups are ranked by a shrinkage-corrected
function of each group’s own top members by default (small groups pulled
toward the pack rather than trusted on a raw sample mean — an
empirical-Bayes/BLUP correction standard in quantitative genetics family
evaluation); an optional `within_group_target_degree` (the same 0-90
convention as
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)/[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md))
balances score against a relatedness ceiling within each group; and its
own optional `ensure_haplotype_diversity` mode brings a lighter-weight
version of
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
complementarity idea to a group-quota shortlist (two methods:
coverage-gain, the default, or the original dominant-block heuristic),
without requiring a GA search. See
[`?select_parents_by_family`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)
for the full decision guide on when to reach for it instead of the two
functions above.

**From a founder shortlist to a full crossing decision.**
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)/
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)/[`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)
all return a flat parent set with no notion of pairing, differential
contributions, or a population-wide relatedness cap.
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
and
[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
close that gap — ranking candidate crosses, solving the classical
Meuwissen (1997) contribution- optimisation problem via two true-OCS
engines (AlphaMate; optiSel’s own `candes()`/`opticont()`/`matings()`
solver) or, alternatively, discrete greedy cross prediction/selection
via a third engine (SimpleMating’s `planCross()`/`selectCrosses()`),
making the gain-vs-diversity tradeoff and a heuristic plan’s optimality
gap explicit, and covering the diversity-first case respectively. See
the *From Local GEBV to a Crossing Decision* vignette (Sections 8-13)
for the full worked example end to end, from a founder shortlist through
to an actual mating plan, including a real GA-vs-
truncation-vs-OCS-vs-Pareto comparison and a decision table mapping each
tool to the breeding question it answers.

**How these five functions relate, precisely.**
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
and
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
(and
[`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md))
sit at the *same* stage and are alternatives to each other, not a
sequence — each answers “who are my candidate parents” by different
logic, and the normal workflow is to run more than one and compare, not
to feed one’s output into another.
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
and
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
sit one stage downstream of *whichever* shortlist you settled on, and
neither builds its own candidate pool: both are written to consume any
flat ID list plus its merit values (and, for
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
`G`) — the vignette’s `ga_sel$selected` in its worked examples is
illustrative, not a requirement; `ts_sel$selected`, `fam_sel$selected`,
or a hand-picked list work identically.
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
and
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
are themselves siblings, not a sequence, either — cross ranking scores
candidate pairs independently with no notion of an overall plan, while
OCS solves a different, population-level problem (contributions and
matings under an explicit relatedness cap) directly from merit and `G`,
without needing cross ranking’s table as an input.

**A qualification on “downstream” for
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md).**
That function can never select a parent absent from the `merit`/`G` list
you supply — in that sense it’s a strict downstream consumer. But none
of its three engines is obliged to give every supplied candidate a
nonzero contribution, and `n_parents_max` makes that narrowing explicit
and engine-dependent: `engine = "alphamate"` enforces it *natively
inside* its own evolutionary search (jointly deciding which subset of
your candidates is used, how much each contributes, and who mates with
whom, in one optimisation — a real selection step happening inside the
mating-plan tool itself); `engine = "optisel"` approximates it by
solving the unconstrained contribution optimum over your whole list
first, then zeroing all but the `n_parents_max` largest contributions as
a post-hoc, non-re-optimised patch; `engine = "simplemating"` does not
support it at all (ignored, with a message). Your shortlist is always a
ceiling on who *can* appear in the final plan, not a guarantee everyone
on it *will* — use `engine = "alphamate"` if you want that further
narrowing done well, as one integrated optimisation rather than
shortlist-then-truncate.

### 15.10. Utilities

| Function | Description |
|----|----|
| [`summarise_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/summarise_blocks.md) | Per-chromosome and genome-wide block size summary statistics. |
| [`plot_ld_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ld_blocks.md) | ggplot2 block diagram coloured by block size or chromosome. |

------------------------------------------------------------------------

## 16. Output objects

### 16.1. `run_Big_LD_all_chr()` — block table

| Column       | Type      | Description                                         |
|--------------|-----------|-----------------------------------------------------|
| `start`      | integer   | Index of the first SNP in the block (full SNP set). |
| `end`        | integer   | Index of the last SNP.                              |
| `start.rsID` | character | SNP identifier at the block start.                  |
| `end.rsID`   | character | SNP identifier at the block end.                    |
| `start.bp`   | numeric   | Base-pair position of the block start.              |
| `end.bp`     | numeric   | Base-pair position of the block end.                |
| `CHR`        | character | Chromosome label (normalised, no `chr` prefix).     |
| `length_bp`  | integer   | `end.bp - start.bp + 1`.                            |

### 16.2. `run_ldx_pipeline()` — named list

| Element | Type | Description |
|----|----|----|
| `blocks` | data.frame | LD block table (same columns as Section 16.1). |
| `diversity` | data.frame | Per-block diversity metrics (same columns as Section 16.4). |
| `hap_matrix` | matrix | Haplotype feature matrix (individuals × haplotype columns). |
| `hap_matrix_info` | data.frame | Column metadata for `hap_matrix`. |
| `haplotypes` | named list | Raw haplotype strings per block (same structure as Section 16.3). |
| `geno_matrix` | matrix | Cleaned, imputed genotype matrix (individuals × SNPs). |
| `snp_info_filtered` | data.frame | SNP metadata after MAF filtering. |
| `phased_vcf` | character | Path to Beagle-phased VCF.gz in `out_dir`. `NULL` when `phase = FALSE`. |
| `phased_backend_desc` | character | Path to `hapblockr_bm_phased_dos.desc` for reattachment. `NULL` when `use_bigmemory = FALSE` or `phase = FALSE`. |
| `phase_method` | character | `"beagle"` when Beagle was run; `"unphased"` otherwise. |
| `n_blocks` | integer | Total number of LD blocks detected. |
| `n_hap_columns` | integer | Number of columns in `hap_matrix`. |

**Bigmemory phased backends** (written to `bigmemory_path` when
`use_bigmemory = TRUE` and `phase = TRUE`):

| File | Content | Type |
|----|----|----|
| `hapblockr_bm_phased_hap1.bin/.desc` | Gamete 1 alleles (individuals × SNPs) | `char`, 0/1 |
| `hapblockr_bm_phased_hap2.bin/.desc` | Gamete 2 alleles (individuals × SNPs) | `char`, 0/1 |
| `hapblockr_bm_phased_dos.bin/.desc` | Dosage = hap1 + hap2 (individuals × SNPs) | `char`, 0/1/2 |
| `hapblockr_bm_phased_snpinfo.rds` | Aligned SNP metadata | — |
| `hapblockr_bm_phased_sampleids.rds` | Sample IDs in aligned order | — |
| `hapblockr_bm_phased_params.rds` | Fingerprint for cache validation | — |

On subsequent runs with the same `bigmemory_path` and unchanged phased
VCF, backends are reattached from disk without re-reading the VCF.

### 16.3. `tune_LD_params()` — named list

| Element | Type | Description |
|----|----|----|
| `best_params` | named list | Selected parameter values. |
| `score_table` | data.frame | All grid combinations and scores. |
| `perfect_table` | data.frame or NULL | Combinations with n_unassigned = 0 and n_forced = 0. |
| `final_blocks` | data.table | Block table from `best_params`, all chromosomes. |
| `gwas_assigned` | data.frame | Input GWAS data with `LD_block` column added. Entries ending in `*` denote forced (nearest-block) assignments. |

### 16.4. `extract_haplotypes()` — named list

| Element | Type | Description |
|----|----|----|
| `block_<start>_<end>` | character vector | One haplotype string per individual. |
| `attr(., "block_info")` | data.frame | block_id, CHR, start_bp, end_bp, n_snps. |

### 16.5. `compute_haplotype_diversity()` — data.frame

| Column | Description |
|----|----|
| `block_id` | Block name matching `names(haplotypes)`. |
| `CHR` | Chromosome (normalised). |
| `start_bp`, `end_bp` | Block coordinates in base pairs. |
| `n_snps` | Number of SNPs in the block. |
| `n_ind` | Individuals with non-missing haplotypes. |
| `n_haplotypes` | Richness: number of unique haplotype strings. |
| `He` | Expected heterozygosity (Nei 1973), sample-size corrected. |
| `Shannon` | Shannon entropy in bits. |
| `n_eff_alleles` | Effective number of alleles = 1/Σpᵢ². |
| `freq_dominant` | Frequency of the most common haplotype. |
| `sweep_flag` | TRUE when freq_dominant ≥ 0.90. |
| `phased` | Logical: was phased input used? |

### 16.6. `read_geno()` — HapBlockR_backend

| Element | Type | Description |
|----|----|----|
| `type` | character | `"numeric"`, `"hapmap"`, `"vcf"`, `"gds"`, `"bed"`, or `"matrix"`. |
| `n_samples` | integer | Number of individuals. |
| `n_snps` | integer | Number of SNPs. |
| `sample_ids` | character | Individual identifiers. |
| `snp_info` | data.frame | SNP metadata: SNP, CHR, POS, REF, ALT. |

### 16.7. `test_block_haplotypes()` — HapBlockR_haplotype_assoc

**`$allele_tests`** — one row per allele per block per trait:

| Column | Type | Description |
|----|----|----|
| `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` | — | Identifiers |
| `allele` | character | Haplotype allele string identifier |
| `allele_freq_tested` | numeric \[0,1\] | Allele frequency among tested individuals |
| `effect` | numeric | Additive effect on de-regressed scale |
| `SE` | numeric | Standard error |
| `t_stat` | numeric | t-statistic |
| `p_wald` | numeric | Two-sided Wald p-value (raw) |
| `p_fdr` | numeric | Benjamini-Hochberg FDR-adjusted p-value |
| `Meff` | numeric | Effective number of independent tests (simpleM, scope = `meff_scope`) |
| `alpha_simplem` | numeric | simpleM Bonferroni-style significance threshold: α / Meff |
| `alpha_simplem_sidak` | numeric | simpleM Šidák-style threshold: 1 − (1−α)^(1/Meff) |
| `p_simplem` | numeric | simpleM Bonferroni-style adjusted p-value: min(p × Meff, 1) |
| `p_simplem_sidak` | numeric | simpleM Šidák-style adjusted p-value: 1 − (1−p)^Meff |
| `significant` | logical | TRUE when the p-value chosen by `sig_metric` ≤ `sig_threshold` |

All p-value columns are always present regardless of `sig_metric`. The
`significant` flag is the only column that changes based on
`sig_metric`.

**`$block_tests`** — one row per block per trait:

| Column | Type | Description |
|----|----|----|
| `n_alleles_tested` | integer | Alleles above `min_freq` tested jointly |
| `F_stat` | numeric | Omnibus F-statistic |
| `df_LRT` | integer | Numerator df = n_alleles_tested |
| `p_omnibus` | numeric | Raw omnibus p-value |
| `p_omnibus_fdr` | numeric | BH FDR-adjusted omnibus p-value |
| `p_omnibus_adj` | numeric | Bonferroni-style per-trait adjustment (backward-compat) |
| `var_explained` | numeric \[0,1\] | Proportion of de-regressed variance explained |
| `Meff` | numeric | Block-level effective test count (from block-summary PC1 matrix) |
| `alpha_simplem` | numeric | Block-level simpleM Bonferroni threshold |
| `alpha_simplem_sidak` | numeric | Block-level simpleM Šidák threshold |
| `p_omnibus_simplem` | numeric | simpleM Bonferroni-style adjusted omnibus p-value |
| `p_omnibus_simplem_sidak` | numeric | simpleM Šidák-style adjusted omnibus p-value |
| `significant_omnibus` | logical | TRUE when the p-value chosen by `sig_metric` ≤ `sig_threshold` |

**Additional return list elements:**

| Element | Description |
|----|----|
| `sig_metric` | Which p-value was used for `significant` / `significant_omnibus` |
| `meff_scope` | Scope used for Meff estimation (`"chromosome"`, `"global"`, or `"block"`) |
| `meff_percent_cut` | Variance threshold used for simpleM eigendecomposition |
| `meff` | Named list of Meff summaries per trait: `$allele$global`, `$allele$chromosome`, `$allele$block`, `$block$global`, `$block$chromosome` |
| `pc_model_selection` | `data.frame` with one row per k tested when `optimize_pcs = TRUE`: `n_pcs`, `BIC`, `lambda_gc`, `score`, `selected`. `NULL` when `optimize_pcs = FALSE` |

**Plots produced when `plot = TRUE`** (all saved as PDF):

| File | Description |
|----|----|
| `manhattan_<trait>.pdf` | Manhattan plot — one per trait |
| `qq_<trait>.pdf` | QQ plot — one per trait |
| `pca_grm.pdf` | PCA of GRM eigenvectors in PC1 × PC2, coloured by first trait phenotype |
| `grm_scree.pdf` | Scree plot of GRM eigenvalues (up to PC30); red bar and vertical line at selected k |

### 16.8. `estimate_diplotype_effects()` — HapBlockR_diplotype

**`$omnibus_tests`** — one row per LD block per trait (updated columns):

| Column | Type | Description |
|----|----|----|
| `block_id`, `trait` | — | Identifiers |
| `n_diplotypes` | integer | Diplotype classes with ≥ `min_n_diplotype` individuals |
| `F_stat` | numeric | F-statistic (one-way ANOVA on GRM-corrected residuals) |
| `df1`, `df2` | integer | Numerator and denominator degrees of freedom |
| `p_omnibus` | numeric | Raw p-value |
| `p_omnibus_adj` | numeric | Plain Bonferroni × n_blocks_per_trait (backward-compat) |
| `p_omnibus_fdr` | numeric | BH-FDR per trait |
| `p_omnibus_simplem` | numeric | simpleM Bonferroni-style: min(p × Meff, 1) |
| `p_omnibus_simplem_sidak` | numeric | simpleM Šidák-style: 1 − (1−p)^Meff (**recommended**) |
| `Meff` | numeric | Effective number of blocks (block-summary PC1 eigenspectrum) |
| `significant` | logical | TRUE when the p-value chosen by `sig_metric` \< `sig_threshold` |

All five p-value columns are always present regardless of `sig_metric`.

**`$dominance_table`** — one row per allele pair per block per trait:

| Column | Type | Description |
|----|----|----|
| `allele_A`, `allele_B` | character | The two alleles (A ≤ B alphabetically) |
| `mean_AA`, `mean_AB`, `mean_BB` | numeric | Diplotype class means |
| `a` | numeric | Additive effect: (mean_BB − mean_AA) / 2 |
| `d` | numeric | Dominance deviation: mean_AB − midpoint |
| `d_over_a` | numeric or NA | Dominance ratio. 0 = additive; ±1 = complete dominance; \|d/a\| \> 1 = overdominance |
| `overdominance` | logical | TRUE when \|d/a\| \> 1 |

### 16.9. `score_favorable_haplotypes()` — data frame

| Column | Type | Description |
|----|----|----|
| `id` | character | Individual identifier |
| `stacking_index` | numeric \[0,1\] | Genome-wide normalised score |
| `n_blocks_scored` | integer | Blocks with at least one matched allele effect |
| `mean_block_score` | numeric | Raw mean per-block score |
| `rank` | integer | Rank by stacking_index (1 = highest) |
| `score_<block_id>` | numeric | Per-block score columns |

### 16.10. `summarize_parent_haplotypes()` — data frame

| Column | Type | Description |
|----|----|----|
| `id`, `block_id`, `CHR`, `start_bp`, `end_bp` | — | Identifiers |
| `allele` | character | Haplotype allele string |
| `dosage` | integer {0,1,2} | Copies carried. Phased: 0/1/2. Unphased: 0/1. |
| `allele_freq` | numeric \[0,1\] | Population frequency in full panel |
| `allele_effect` | numeric or NA | Effect from allele_effects argument |
| `is_rare` | logical | TRUE when allele_freq \< 0.10 |

### 16.11. `compare_block_effects()` — HapBlockR_effect_concordance

**`$concordance`** — one row per block per trait:

| Column | Type | Description |
|----|----|----|
| `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` | — | Identifiers |
| `n_alleles_pop1`, `n_alleles_pop2` | integer | Alleles tested in each population before intersection |
| `n_shared_alleles` | integer | Alleles present (by string match) in both populations |
| `enough_shared` | logical | `n_shared_alleles >= min_shared_alleles` |
| `effect_correlation` | numeric | Pearson r of per-allele effects across populations (NA when n_shared \< 3) |
| `direction_agreement` | numeric \[0,1\] | Fraction of shared alleles with the same effect sign |
| `directionally_concordant` | logical | `direction_agreement >= direction_threshold` |
| `meta_effect` | numeric | IVW meta-analytic effect (weighted mean of per-allele IVW estimates) |
| `meta_SE` | numeric | SE of IVW estimate |
| `meta_z` | numeric | meta-analytic z-score |
| `meta_p` | numeric | Two-sided p-value of meta-analytic effect |
| `Q_stat` | numeric | Cochran Q heterogeneity statistic |
| `Q_df` | integer | Degrees of freedom (n_shared_alleles − 1) |
| `Q_p` | numeric | p-value of Q under chi-squared; significant = heterogeneity between populations |
| `I2` | numeric \[0,100\] | I² inconsistency: max(0, (Q − df) / Q × 100). \> 50% = substantial heterogeneity |
| `replicated` | logical | `enough_shared AND directionally_concordant AND Q_p > 0.05` |
| `boundary_overlap_ratio` | numeric \[0,1\] | **Automatically computed output** from `blocks_pop1` / `blocks_pop2`: bp(intersection) / bp(union). `NA` when block tables not supplied |
| `boundary_warning` | logical | `TRUE` when `boundary_overlap_ratio < boundary_overlap_warn` (the **input parameter**, default 0.80); flags blocks where differing LD structure may make haplotype strings non-comparable |

**`$shared_alleles`** — one row per shared allele per block per trait:

| Column | Type | Description |
|----|----|----|
| `allele` | character | Haplotype allele string |
| `effect_pop1`, `SE_pop1`, `p_wald_pop1` | numeric | Effect, SE, and p-value from population 1 |
| `effect_pop2`, `SE_pop2`, `p_wald_pop2` | numeric | Effect, SE, and p-value from population 2 |
| `direction_agree` | logical | Same sign in both populations? |
| `ivw_effect` | numeric | Per-allele IVW combined effect (weighted mean of pop1 and pop2 by 1/SE²) |
| `ivw_SE` | numeric | SE of IVW combined effect |
| `match_type` | character | How this block was matched: `"exact"` (same block_id), `"position"` (matched by genomic IoU ≥ `overlap_min`), `"pop1_only"` (no Pop2 block overlapped), or `NA` when no block tables supplied. Set when `block_match = "position"`. |

------------------------------------------------------------------------

### 16.12. `compare_gwas_effects()` — HapBlockR_effect_concordance

Same output class as
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)
(Section 16.11). Additional columns in `$concordance` specific to
external GWAS input:

| Column | Type | Description |
|----|----|----|
| `lead_snp_pop1`, `lead_snp_pop2` | character | Lead SNP ID from each population (same SNP = same LD tag; different = different proxies of the same QTL region) |
| `lead_p_pop1`, `lead_p_pop2` | numeric | Lead SNP p-values |
| `se_derived_pop1`, `se_derived_pop2` | logical | `TRUE` when SE was derived from `BETA` and `P` via z-score rather than read directly |
| `both_pleiotropic` | logical | `TRUE` when the block is pleiotropic in both populations (requires `pleiotropic` column from [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)) |

Columns that differ in meaning from
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md):

| Column | In [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md) | In [`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md) |
|----|----|----|
| `n_shared_alleles` | Haplotype alleles in both pops | Always 1 (one lead SNP per block) |
| `effect_correlation` | Pearson r across alleles | Always `NA` (needs ≥ 3 alleles) |
| `direction_agreement` | Fraction of alleles with same sign | 0 or 1 only |
| `Q_stat`, `Q_p`, `I2` | Cochran Q heterogeneity test | Always `NA` (df = 0) |
| `replicated` | dir_concordant AND Q_p \> 0.05 | dir_concordant AND meta_p ≤ 0.05 |
| `match_type` | `"exact"` / `"position"` / `"pop1_only"` / `NA` | Same semantics as in [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md). `block_match = "position"` applies to both functions. |

The `$shared_alleles` data frame contains one row per block per trait
with `lead_snp_pop1`, `lead_snp_pop2` instead of `allele`.

## 17. Memory and performance notes

### 17.1. C++ core

The thirteen compiled functions in `src/ld_core.cpp` replace the most
expensive R operations: `compute_r2_cpp`, `compute_rV2_cpp`,
`maf_filter_cpp`, `build_adj_matrix_cpp`, `col_r2_cpp`,
`compute_r2_sparse_cpp`, `boundary_scan_cpp`, `build_hap_strings_cpp`,
`resolve_overlap_cpp`, `block_snp_ranges_cpp`,
`extract_chr_haplotypes_cpp`, `extract_chr_haplotypes_phased_cpp`, and
`impute_and_filter_cpp`. Key speedups:

- **`compute_r2_cpp()`**: ~40× over pure R for a 1,500 × 1,500 window
  with 500 individuals.
- **`boundary_scan_cpp()`**: eliminates ~150,000 small R-level matrix
  operations for a 50,000-SNP chromosome.
- **`maf_filter_cpp()`**: ~10× faster for panels \> 100,000 SNPs.
- **`resolve_overlap_cpp()`**: 15,700× per-SNP reduction on chr1 via
  BLAS DGEMM scoring with a lazy column cache.

### 17.2. Never-full-genome memory model

- **Numeric dosage CSV**: two-pass chunked reading. Peak RAM = one
  50,000-row chunk.
- **VCF and HapMap**: auto-converted to GDS cache; all access streaming
  via
  [`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md).
- **GDS and PLINK BED**:
  [`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md)
  called once per sub-segment per chromosome.
- **Chromosome loop**: [`rm()`](https://rdrr.io/r/base/rm.html) and
  `gc(FALSE)` called after each chromosome.

### 17.3. OpenMP thread count

``` r
n_thr  <- parallel::detectCores(logical = FALSE)
blocks <- run_Big_LD_all_chr(be, n_threads = n_thr)
```

Efficient scaling up to approximately 8-16 threads for
`subSegmSize = 1500`.

------------------------------------------------------------------------

## 18. Citation

    HapBlockR Development Team (2025).
    HapBlockR: Genome-Wide LD Block Detection, Haplotype Analysis, and Genomic
    Prediction Features with Kinship-Adjusted Correlations.
    R package version 0.3.1.
    https://github.com/FAkohoue/HapBlockR

Please also cite:

    Kim S-A et al. (2018). Bioinformatics 34(4):588-596.
    VanRaden PM (2008). Journal of Dairy Science 91(11):4414-4423.
    Calus MPL et al. (2008). Genetics 178(1):553-561.
    Nei M (1973). PNAS 70(12):3321-3323.
    Blondel VD et al. (2008). J. Stat. Mech. P10008.
    Traag VA et al. (2019). Scientific Reports 9:5233.

------------------------------------------------------------------------

## 19. Contributing

<https://github.com/FAkohoue/HapBlockR/issues>

Before opening a pull request: run
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
with zero errors/warnings, add tests in `tests/testthat/test-core.R`,
rebuild docs with
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html).

------------------------------------------------------------------------

## 20. License

MIT + file LICENSE © Félicien Akohoue

------------------------------------------------------------------------

## 21. References

Gao X, Starmer J, Martin ER (2008). A multiple testing correction method
for genetic association studies using correlated single nucleotide
polymorphisms. *Genetic Epidemiology* **32**:361-369.
<doi:10.1002/gepi.20310>

Gao X, Becker LC, Becker DM, Starmer JD, Province MA (2010). Avoiding
the high Bonferroni penalty in genome-wide association studies. *Genetic
Epidemiology* **34**:100-105. <doi:10.1002/gepi.20430>

Gao X (2011). Multiple testing corrections for imputed SNPs. *Genetic
Epidemiology* **35**:154-158. <doi:10.1002/gepi.20563>

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009). *Introduction
to Meta-Analysis*. Wiley.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine* **21**(11):1539-1558.
<doi:10.1002/sim.1186>

Kim S-A, Cho C-S, Kim S-R, Bull SB, Yoo Y-J (2018). A new haplotype
block detection method for dense genome sequencing data based on
interval graph modeling and dynamic programming. *Bioinformatics*
**34**(4):588-596.

VanRaden PM (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science* **91**(11):4414-4423.

Mangin B et al. (2012). Novel measures of linkage disequilibrium that
correct the bias due to population structure and relatedness. *Heredity*
**108**(3):285-291.

Calus MPL et al. (2008). Accuracy of genomic selection using different
methods to define haplotypes. *Genetics* **178**(1):553-561.

de Roos APW, Hayes BJ, Goddard ME (2009). Reliability of genomic
predictions across multiple populations. *Genetics*
**183**(4):1545-1553.

Nei M (1973). Analysis of gene diversity in subdivided populations.
*PNAS* **70**(12):3321-3323.

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection. *Theor. Appl. Genet.* **137**:274.

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theor. Appl. Genet.* **138**:267.

Blondel VD et al. (2008). Fast unfolding of communities in large
networks. *J. Stat. Mech.* P10008.

Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden:
guaranteeing well-connected communities. *Scientific Reports*
**9**:5233.
