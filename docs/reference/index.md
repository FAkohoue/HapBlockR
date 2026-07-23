# Package index

## Package overview

Main documentation entry point for HapBlockR, and the standalone
Breeder’s Guide covering all parent-selection and mate-allocation
strategies in narrative form.

- [`HapBlockR-package`](https://FAkohoue.github.io/HapBlockR/reference/HapBlockR-package.md)
  [`HapBlockR`](https://FAkohoue.github.io/HapBlockR/reference/HapBlockR-package.md)
  : HapBlockR: Genome-Wide LD Block Detection, Haplotype Analysis,
  Genomic Prediction, and Breeding Decision Support
- [`open_breeder_guide()`](https://FAkohoue.github.io/HapBlockR/reference/open_breeder_guide.md)
  : Locate or Open the HapBlockR Breeder's Guide

## Main workflows

High-level functions covering the full HapBlockR pipeline, from genotype
input to LD block detection, haplotype construction, and genomic
prediction.

- [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md)
  : LD Block Segmentation (r^2 or rV^2, C++ accelerated)
- [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)
  : End-to-End Haplotype Block Pipeline
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
  : Genome-Wide LD Block Detection by Chromosome
- [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md)
  : Auto-Tune LD Block Detection Parameters

## Genotype input and backend

Import genotype data from multiple formats and manage the HapBlockR
memory-efficient backend for large genomic datasets.

- [`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
  : Read Genotype Data into an HapBlockR Backend
- [`read_geno_bigmemory()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno_bigmemory.md)
  : Open a bigmemory-backed Genotype Store
- [`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md)
  : Extract a Genotype Slice from an HapBlockR Backend
- [`prepare_geno()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_geno.md)
  : Prepare Genotype Matrix for LD Computation
- [`close_backend()`](https://FAkohoue.github.io/HapBlockR/reference/close_backend.md)
  : Close an HapBlockR Backend and Release File Handles
- [`print(`*`<HapBlockR_backend>`*`)`](https://FAkohoue.github.io/HapBlockR/reference/print.HapBlockR_backend.md)
  : Print Method for HapBlockR Backend
- [`summary(`*`<HapBlockR_backend>`*`)`](https://FAkohoue.github.io/HapBlockR/reference/summary.HapBlockR_backend.md)
  : Summary Method for HapBlockR Backend

## LD computation and block detection

Core algorithms for LD estimation, LD decay analysis, and genome-wide LD
block segmentation using clique-based and community-detection
approaches.

- [`compute_r2()`](https://FAkohoue.github.io/HapBlockR/reference/compute_r2.md)
  : Compute Standard r^2 LD Matrix
- [`compute_rV2()`](https://FAkohoue.github.io/HapBlockR/reference/compute_rV2.md)
  : Compute Kinship-Adjusted rV^2 LD Matrix
- [`get_V_inv_sqrt()`](https://FAkohoue.github.io/HapBlockR/reference/get_V_inv_sqrt.md)
  : Compute the Inverse Square Root (Whitening Factor) of a Kinship
  Matrix
- [`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  : Compute LD Decay and Chromosome-Specific Decay Distances
- [`plot_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ld_decay.md)
  : Plot LD Decay Curve
- [`summarise_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/summarise_blocks.md)
  : Summarise LD Block Characteristics
- [`plot_ld_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ld_blocks.md)
  : Plot LD Block Structure Across Chromosomes

## Phasing tools

Functions for obtaining gametic phase from WGS data (Beagle) or
pre-phased VCF files.

- [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md)
  : Read Pre-Phased VCF
- [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md)
  : Statistical Phasing via Beagle 5.x

## Haplotype construction and diversity

Extract haplotypes within LD blocks, decode nucleotide sequences,
compute haplotype diversity statistics, and construct genomic prediction
features.

- [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  : Extract Haplotype Dosage Strings from LD Blocks
- [`decode_haplotype_strings()`](https://FAkohoue.github.io/HapBlockR/reference/decode_haplotype_strings.md)
  : Decode Haplotype Strings to Nucleotide Sequences
- [`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md)
  : Compute Haplotype Diversity Per Block
- [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md)
  : Build Haplotype Dosage Matrix for Genomic Prediction
- [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
  : Compute Haplotype-Based Genomic Relationship Matrix
- [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
  : Map GWAS Hits to LD Blocks (Post-GWAS QTL Region Definition)

## Haplotype harmonisation

Harmonise haplotypes across populations, infer diplotypes, and collapse
rare alleles into biologically meaningful groups.

- [`infer_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/infer_block_haplotypes.md)
  : Infer Structured Block-Level Diplotypes Per Individual
- [`collapse_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/collapse_haplotypes.md)
  : Collapse Rare Haplotype Alleles Into Biologically Meaningful Groups
- [`harmonize_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/harmonize_haplotypes.md)
  : Harmonize Haplotype Allele Labels Across Panels or Analysis Runs

## Haplotype export

Export haplotype dosage matrices and diversity summaries for downstream
statistical analyses.

- [`write_haplotype_numeric()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_numeric.md)
  : Write Haplotype Feature Matrix as Numeric Dosage Table
- [`write_haplotype_character()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_character.md)
  : Write Haplotype Character (Nucleotide) Matrix
- [`write_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_diversity.md)
  : Write Haplotype Diversity Table

## Association and QTL interpretation

Block-level association tests (Q+K mixed model with simpleM
multiple-testing correction), diplotype effect decomposition,
cross-population effect concordance from haplotype association results
(compare_block_effects) or external GWAS tools (compare_gwas_effects),
GWAS integration, and epistasis detection.

- [`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
  : Block-Level Haplotype Association Testing (Q+K Mixed Linear Model
  with simpleM Multiple-Testing Correction)
- [`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
  : Estimate Diplotype Effects and Dominance Deviations Per LD Block
- [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)
  : Cross-Population Haplotype Effect Concordance
- [`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md)
  : Cross-Population GWAS Effect Concordance from External Results
- [`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/integrate_gwas_haplotypes.md)
  : Integrate GWAS QTL Regions with Haplotype Prediction Results

## Genomic prediction

Haplotype-based genomic prediction following Tong et al. (2024-2025),
including GBLUP integration, local GEBV estimation, and additive +
dominance GBLUP via the Vitezica et al. (2013) dominance relationship
matrix (compute_dominance_grm, paired with run_haplotype_prediction’s
include_dominance = TRUE).

- [`prepare_gblup_inputs()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md)
  : Prepare Genomic Prediction Inputs for External GBLUP Software
- [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
  : Haplotype Prediction and Block Importance from Pre-Adjusted
  Phenotypes
- [`estimate_marker_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
  : Estimate Per-SNP Marker Effects via GBLUP, rrBLUP (SNP-BLUP), or
  Bayesian Regression
- [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
  : Compute Local Haplotype GEBV per Block (Tong et al. 2025)
- [`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  : Backsolve SNP Effects from GEBV (Tong et al. 2025)
- [`compute_dominance_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md)
  : Genomic Dominance Relationship Matrix (Vitezica et al. 2013)

## Epistasis detection

Within-block pairwise SNP interaction scan (scan_block_epistasis),
trans-haplotype between-block epistasis scan
(scan_block_by_block_epistasis), and single-block fine-mapping with
pairwise or LASSO dispatch (fine_map_epistasis_block). All functions
operate on GRM-corrected REML residuals from test_block_haplotypes().

- [`scan_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_epistasis.md)
  : Within-Block Pairwise SNP Epistasis Scan
- [`scan_block_by_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_by_block_epistasis.md)
  : Between-Block Haplotype Allele Epistasis Scan
- [`fine_map_epistasis_block()`](https://FAkohoue.github.io/HapBlockR/reference/fine_map_epistasis_block.md)
  : Fine-Map Epistatic SNP Pairs Within a Single Block

## Breeding decision support

Identify favourable haplotypes and rank LD blocks for haplotype stacking
strategies.

- [`rank_haplotype_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/rank_haplotype_blocks.md)
  : Rank Haplotype Blocks by Evidence Strength
- [`select_top_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  : Select Top-Ranked Blocks by Count, Percentage, or Cumulative
  Variance
- [`plot_block_funnel()`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md)
  : Funnel Plot of Block Importance
- [`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md)
  : Score Individual Haplotype Portfolios Against Known Allele Effects
- [`summarize_parent_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/summarize_parent_haplotypes.md)
  : Summarise Haplotype Allele Inventory Per Candidate Parent

## Parent selection and mate allocation

Genetic-algorithm founder-parent optimiser (select_parents_ga, via the
GA package) searching for the combination of individuals maximising
coverage of favourable local GEBV/haplotype values across target blocks,
under crossing-scheme constraints (no_selfing/selfing/OHS/
OPV/Haploid_OHS), with a truncation-selection baseline, a calibration
helper for its optional merit-weighted fitness (suggest_merit_weight),
and a PCA visualisation for comparing parent sets. Family- (or
genetic-cluster-) quota selection (select_parents_by_family): best-
performing groups first, then the best individual lines within each, via
four mutually exclusive rules for how many lines a group takes
(family_select_mode = count/percentage/sd_threshold/check_relative),
with an optional GBLUP-corrected group ranking and a within-group
haplotype-diversity tiebreak. Beyond founder-set selection: an explicit
merit-vs-diversity Pareto tradeoff frontier (select_parents_pareto,
pareto_front); true Optimal Contribution Selection producing an actual
contribution-and-mating plan (select_parents_ocs, via the AlphaMate
executable, optiSel, or SimpleMating’s planCross()/selectCrosses());
cross ranking by predicted mid-parent value plus segregation variance
across four variance-prediction modes (usefulness_criterion); exact
binary integer linear programming validation of a heuristic mating
plan’s optimality gap (validate_crosses_exact); and diversity-maximising
core-collection selection for germplasm-bank/founder-diversity use cases
(select_core_collection).

- [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  : Genetic-Algorithm Founder-Parent Selection
- [`suggest_merit_weight()`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)
  : Suggest a Starting merit_weight for select_parents_ga()
- [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
  : Truncation Selection (Top-n by a Single Score)
- [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)
  : Family- or Genetic-Cluster-Quota Parent Selection (Best Groups, Then
  Best Lines Within Them)
- [`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)
  : PCA of Parent-Selection Groups
- [`cluster_selection_groups()`](https://FAkohoue.github.io/HapBlockR/reference/cluster_selection_groups.md)
  : Cluster Individuals on Retained Principal Components and
  Cross-Tabulate Selection Coverage
- [`plot_selection_clusters()`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md)
  : Plot Individuals Coloured by Genetic Cluster
- [`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)
  : Pareto Frontier of Parent Sets: Merit vs. Relatedness
- [`pareto_front()`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md)
  : Pareto (Non-Dominated) Front of a Set of Candidates
- [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
  : True Optimal Contribution Selection (OCS) and Mate Allocation
- [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
  : Usefulness Criterion (UC) / Genomic Mating: Rank Candidate Crosses
- [`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)
  : Exact (ILP) Cross Selection: A Validation Check for Heuristic Mating
  Plans
- [`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
  : Core-Collection / Diversity-Maximizing Subset Selection

## Forward-in-time simulation

Recurrent-selection simulator (ga_vs_ts_simulation), built as a wrapper
around the genomicSimulation package, comparing GA-selected against
truncation-selected founders over generations of realistic meiosis and
recombination, tracking realised mean/max breeding- population GEBV per
generation.

- [`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)
  : Simulate GA-Selected vs. Truncation-Selected Founders Over
  Generations
- [`plot_ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md)
  : Plot Realised Genetic Gain: GA vs. Truncation Selection (and Beyond)

## Population and stability analysis

Compare haplotype distributions across populations and evaluate
stability across environments.

- [`compare_haplotype_populations()`](https://FAkohoue.github.io/HapBlockR/reference/compare_haplotype_populations.md)
  : Compare Haplotype Allele Frequencies Between Two Population Groups
- [`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
  : Finlay-Wilkinson Stability Analysis of Haplotype Effects Across
  Environments
- [`scan_diversity_windows()`](https://FAkohoue.github.io/HapBlockR/reference/scan_diversity_windows.md)
  : Sliding-Window Genome-Wide Diversity Scan
- [`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)
  : Decompose Per-SNP Effects into Per-Haplotype-Allele Effect Table
- [`export_candidate_regions()`](https://FAkohoue.github.io/HapBlockR/reference/export_candidate_regions.md)
  : Export Candidate Gene Regions to BED, CSV, or biomaRt Format

## Cross-validation

Evaluate predictive ability of haplotype-based genomic models.

- [`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md)
  : K-Fold Cross-Validation for Haplotype-Based Genomic Prediction

## Visualisation

Plot LD structure and haplotype relationships.

- [`plot_haplotype_network()`](https://FAkohoue.github.io/HapBlockR/reference/plot_haplotype_network.md)
  : Plot a Minimum-Spanning Haplotype Network for One LD Block

## Example datasets

Simulated datasets for tutorials and reproducible examples.

- [`ldx_geno`](https://FAkohoue.github.io/HapBlockR/reference/ldx_geno.md)
  : Example Genotype Matrix
- [`ldx_snp_info`](https://FAkohoue.github.io/HapBlockR/reference/ldx_snp_info.md)
  : Example SNP Information Table
- [`ldx_blocks`](https://FAkohoue.github.io/HapBlockR/reference/ldx_blocks.md)
  : Example LD Block Table
- [`ldx_gwas`](https://FAkohoue.github.io/HapBlockR/reference/ldx_gwas.md)
  : Example GWAS Marker Table
- [`ldx_blues`](https://FAkohoue.github.io/HapBlockR/reference/ldx_blues.md)
  : Pre-Adjusted Phenotype Means (BLUEs) for Genomic Prediction Examples
- [`ldx_blues_list`](https://FAkohoue.github.io/HapBlockR/reference/ldx_blues_list.md)
  : Per-Environment BLUEs for Stability Analysis Examples
