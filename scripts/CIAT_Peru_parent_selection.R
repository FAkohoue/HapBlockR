# ==============================================================================
# CIAT Peru Rice Breeding Program -- Parent Selection for Next Cycle
#
# Uses HapBlockR to decompose an already-estimated Selection Index (GCA-based,
# ASReml GBLUP, missing testcrosses imputed) onto haplotype blocks, cross-
# validates that decomposition, then uses BOTH of HapBlockR's parent-selection
# layers to go from block importance to an actual crossing decision:
#   (a) founder-SET selection (GA vs. truncation) -- "which individuals";
#   (b) mate ALLOCATION (Optimal Contribution Selection + cross ranking +
#       exact validation) -- "how much does each contribute, and which
#       specific matings" -- the true OCS layer this program already runs
#       manually in AlphaMate, now reproducible end-to-end in HapBlockR.
#
# Every argument of every HapBlockR function called below is passed
# explicitly (including ones left at their package default) so nothing here
# depends on a default you can't see. Defaults were verified against the
# current source in R/pipeline.R, R/haplotypes.R, R/haplotype_analysis.R,
# R/parent_selection.R, R/pareto_selection.R, R/ocs.R, R/genomic_mating.R,
# R/exact_validation.R, R/core_collection.R.
#
# Two things this script deliberately does NOT do, with reasons (not silent
# omissions):
#   - Forward-in-time simulation (ga_vs_ts_simulation()) is NOT included.
#     It requires phased hap1/hap2 haplotype matrices; this script's data
#     source is unphased numeric dosage (phase = FALSE below), so there is
#     nothing to feed it. If you later have phased data (a VCF phased with
#     Beagle, phase = TRUE, or a pre-phased VCF via read_phased_vcf()), the
#     forward-simulation vignette section shows how to compare this script's
#     candidate mating schemes (truncation / OCS / UC-informed) over several
#     generations before committing to one for real.
#   - Additive + dominance GBLUP (compute_dominance_grm(),
#     run_haplotype_prediction(include_dominance = TRUE)) is NOT used. Your
#     Selection Index is explicitly GCA-based (General Combining Ability),
#     which is by definition an additive-effects composite derived from
#     testcross means -- layering a dominance relationship matrix on top of
#     an already-additive index does not have a clean interpretation the way
#     it would if you were predicting raw phenotypes or SCA/hybrid vigor
#     directly. Reconsider this if this program starts predicting hybrid
#     performance from raw (not GCA-decomposed) phenotypes.
#
# Pipeline stages used (this list is kept in sync with the numbered section
# comments below -- if you add/remove a stage, update both):
#   1.  (Selection Index)            read Index_Peru.csv
#   2.  run_ldx_pipeline()           genotypes -> LD blocks + geno_matrix
#   3.  (ID check)                   Selection Index IDs vs. genotyped IDs
#   4.  run_haplotype_prediction()   Selection Index -> per-block local GEBV
#   5.  cv_haplotype_prediction()    k-fold predictive-ability check on (4)
#   6.  plot_block_funnel()          visual check of block importance
#   7.  select_top_blocks()          keep the most informative blocks
#   7b. truncation_selection()       truncation-selection baseline
#   7c. (must_include)               force in top whole-genome performers
#   8.  select_parents_ga()          founder-SET selection: GA parent search
#   9.  (family_summary())           family/pedigree balance check
#   10. plot_parent_selection_pca()  diversity-space check -- both genome-
#                                    wide (G) and target-block feature-space
#                                    (feature_matrix) views
#   11. select_parents_pareto()      explicit merit-vs-diversity frontier,
#                                    to inform the OCS target_degree below
#   12. select_parents_ocs()         true Optimal Contribution Selection:
#                                    contributions + an actual mating plan
#                                    (the AlphaMate-equivalent step, now
#                                    reproducible without leaving R)
#   13. usefulness_criterion()       rank the OCS-selected parents' candidate
#                                    crosses by predicted mid-parent value
#                                    plus segregation variance
#                                    (variance_model = "block_independent")
#   13b. usefulness_criterion()      cross-check: same crosses, ranked by
#        (simplemating)              variance_model = "simplemating" -- a
#                                    genuinely linkage-aware, genome-wide
#                                    model (SimpleMating::getUsefA()), IF
#                                    this program's genotypes turn out to be
#                                    strictly homozygous (0/2/NA); errors
#                                    clearly and is skipped, not silently
#                                    ignored, otherwise
#   13c. usefulness_criterion()      cross-check: variance_model = "linked"
#        (linked, not run)           -- NOT run; this program's genotypes are
#                                    unphased, and phasing is the one
#                                    requirement the new ld_matrix fallback
#                                    (added this cycle) does not remove. See
#                                    the commented-out ready-to-use call at
#                                    this stage for exactly what to run once
#                                    phased haplotypes are available.
#   14. validate_crosses_exact()     exact ILP check of the OCS mating plan's
#                                    optimality gap on the ranked cross list
#   15. select_core_collection()     diversity-first cross-check: how does a
#                                    pure-diversity selection of the same
#                                    size compare to the merit-aware picks?
#   16. cluster_selection_groups() / genetic-cluster representation check,
#       plot_selection_clusters()    across ALL FOUR selection strategies
#                                    (GA/TS/OCS/core-collection) at once --
#                                    the strategic population-improvement-
#                                    vs.-fast-release decision support layer
#   17. (final outputs)              write the founder shortlist + mating
#                                    plan CSVs for the next crossing cycle
# ==============================================================================

remove.packages("HapBlockR")

unlink(file.path(.libPaths()[1], "00LOCK-HapBlockR"), recursive = TRUE)

# 2. Delete the compiled DLL from the source tree (src/*.so / src/*.dll).
#    Do this BEFORE the restart so there is nothing to unload conflicts with.
#devtools::clean_dll()

# 3. Restart R.  All loaded DLLs are released, file locks are cleared.
.rs.restartR()
# -- after restart, continue in SESSION B --------------------------------------


#remotes::install_github("vllrs/genomicSimulation")

#install.packages("optiSel")
# remotes::install_github("Resende-Lab/SimpleMating")
#
#
# remotes::install_github("FAkohoue/HapBlockR",
#                         build_vignettes = TRUE,
#                         dependencies    = TRUE
# )

library(HapBlockR)

# ------------------------------------------------------------------------
# 0. USER CONFIG -- edit these before running
# ------------------------------------------------------------------------

geno_file  <- "/opt/Mega-GWAS/1.Data/Libraries Tests/HapBlockR/SNP_data_NUM.csv"   # numeric dosage file -- see note below
index_file <- "/opt/Mega-GWAS/1.Data/Libraries Tests/HapBlockR/Index_Peru.csv"

id_col     <- "GenoID"          # column in Index_Peru.csv with genotype IDs
# (must match the sample column names used
# as headers in SNP_data_NUM, e.g. "Gen490")
index_col  <- "SelectionIndex"    # desired-gain index column (higher = better)
family_col <- "Family"            # family/pedigree grouping column


n_founders <- 30                  # number of parents wanted for the next
# crossing cycle's SHORTLIST (stage 6 --
# GA/truncation) -- set to your program's
# actual crossing-block capacity
strategy   <- "no_selfing"        # crossing-scheme assumption for the GA
# search: "no_selfing" requires each
# favourable block to be covered by TWO
# distinct parents (standard biparental
# crossing block). See ?select_parents_ga
# for "selfing" / "OHS" / "OPV" /
# "Haploid_OHS" alternatives.

min_sel_value <- 0.5               # whole-genome merit floor, a native
# select_parents_ga()/truncation_selection()/
# select_core_collection() argument.
# select_parents_ga() optimises block
# coverage only -- it has no term that
# penalises poor genome-wide merit, so a
# low-performing individual who happens to
# uniquely carry one favourable block will
# still be selected unless excluded here.
min_sel_mode  <- "percentile"      # "value" (min_sel_value is an absolute
# SelectionIndex cutoff -- requires knowing
# its scale up front), "percentile"
# (min_sel_value in (0,1] = fraction of
# candidates kept from the top, e.g. 0.5 =
# top half -- self-scaling, used here),
# or "sd_below_mean" (min_sel_value = SDs
# below the mean). Set min_sel_value <-
# NULL for no floor at all.

n_reps <- 5                        # GA replicates run by select_parents_ga()
# (package default). A single fixed-seed
# GA run gives no evidence its answer is a
# robust optimum rather than one of several
# near-equally-good solutions -- n_reps > 1
# runs independent replicates and reports
# per-individual selection frequency in
# $stability. See ?select_parents_ga's "GA
# rigour" section.

n_must_include <- 2               # always force these many of the very top
# whole-genome-value candidates into the
# final GA set, regardless of block
# coverage. select_parents_ga() has no
# "must include" argument -- its fitness
# function only rewards being top-1/top-2
# AT A TARGET BLOCK, so it can and does
# leave out your single best overall line
# if that line's value is spread broadly
# across many blocks rather than
# concentrated in the highest-variance
# ones. That is mathematically consistent
# with the objective but not what most
# breeders want. Set to 0 to disable and
# let the GA search completely freely.

# -- new: mate-allocation (OCS/UC/Pareto/core-collection) config ------------

n_crosses_ocs <- 60L               # number of actual matings select_parents_
# ocs() should produce for the next
# crossing block -- this is a CAPACITY
# number (how many crosses your program can
# physically make), independent of
# n_founders above (a shortlist size).
# Edit to your program's real crossing
# capacity.
target_degree_ocs <- 30            # select_parents_ocs()'s diversity-vs-gain
# lever, [0, 90], matching AlphaMate's own
# TargetDegree convention: 0 = prioritise
# maximising merit (accept more relatedness),
# 90 = prioritise minimising relatedness.
# 30 (the package default) is a
# starting point ONLY -- stage 8's Pareto
# frontier below is what this program
# should actually use to tune this value
# deliberately, not guess it. There is no
# exact numeric equivalence between this
# and select_parents_ga()'s
# coancestry_weight (stage 8) -- they are
# different mechanisms on different scales
# (see ?select_parents_ocs's "Engine
# differences" section) -- use the Pareto
# frontier to build intuition for how much
# diversity costs THIS population in merit
# terms, then adjust target_degree_ocs by
# trial from there.
max_contrib_per_parent_ocs <- 4L   # (package default) maximum number of
# matings any single parent can
# participate in -- honoured directly by
# all three select_parents_ocs() engines.
ocs_engine <- "simplemating"       # Pinned explicitly (not "auto") to
# preserve this script's existing, already-
# tuned behaviour: as of the HapBlockR
# engine rename/addition, "optisel" now
# means a DIFFERENT true-OCS engine (real
# optiSel::candes()/opticont() solver, new)
# and "auto" falls back to THAT rather than
# to this SimpleMating-based cross-selection
# engine -- so this line must stay pinned to
# "simplemating" for this script's tuning
# (n_crosses/target_degree/etc below) to
# keep meaning what it always has. If you
# want to try the new true-OCS optiSel
# engine here instead, set this to "optisel"
# deliberately and re-check your
# target_degree tuning (see ?select_parents_ocs
# "Engine differences" -- the two engines'
# target_degree map onto different
# mechanisms, not numerically identical).
alphamate_exe <- NULL #"/opt/Mega-GWAS/1.Data/Libraries Tests/HapBlockR/AlphaMate.exe"               # Path to your own AlphaMate executable
# (Hickey group / AlphaGenes suite) -- a
# separate, third-party binary you must
# obtain and be entitled to use yourself;
# see README's "Providing the AlphaMate
# executable" section. Only used if
# ocs_engine above is set to "alphamate" (or
# "auto" and this resolves to a real file).
# In a real deployment, set this to your
# own installed AlphaMate path.

selected_proportion_uc <- 0.1      # (package default) usefulness_criterion()
# -- proportion of each cross's progeny
# you intend to keep, used to weight
# predicted variance against the
# mid-parent mean. Smaller = more weight
# on variance (favour higher-upside,
# riskier crosses).

max_vars_exact <- 2000L            # (package default) validate_crosses_
# exact()'s safety cap on candidate
# crosses considered -- C(length(ocs
# parents), 2) from stage 9 will be well
# under this for a program this size.

n_core <- n_founders               # select_core_collection() core-collection
# size -- matched to n_founders for a
# direct comparison against the GA/
# truncation shortlists (stage 6).

# -- new: strategic decision support (stage 16) ------------------------------

variance_threshold_diversity <- 0.95  # cluster_selection_groups()'s
# cumulative-variance threshold: how
# many leading PCs to retain before
# clustering (unlike stage 10's PCA
# plot, which only ever looks at
# PC1/PC2). Higher = more PCs
# retained = finer-grained clusters,
# at the cost of more noise. 0.95
# (package default) is a reasonable
# starting point.
n_clusters_diversity <- 3L            # number of genetically distinct
# sub-groups to cluster the merit-
# eligible candidate pool into. No
# universal default -- this package
# does not silently guess it for you
# (see ?cluster_selection_groups's
# "Choosing n_clusters" section). If
# this cycle's parents come from a
# known number of source families/
# breeding lines, start there;
# otherwise try a small range (2-6)
# and compare dendrograms
# (clust_res$cluster_fit below).

out_dir <- "/opt/Mega-GWAS/1.Data/Libraries Tests/HapBlockR/ldx_results"
dir.create(out_dir, showWarnings = FALSE)

# IMPORTANT -- file format auto-detection:
# run_ldx_pipeline() detects genotype file FORMAT from its file extension
# only (no explicit format= argument). Your numeric dosage file must be
# saved with a ".csv" extension exactly as in your snippet (SNP, CHR, POS,
# REF, ALT, then one column per sample, values 0/1/2/NA) to be read as
# format = "numeric". If you want to use SNP_data_HAP instead, it must be
# named ending in ".hmp.txt" (not just ".txt") to auto-detect as HapMap,
# OR be read separately via read_geno(SNP_data_HAP, format = "hapmap") and
# the resulting object passed in as geno_source.


# ------------------------------------------------------------------------
# 1. Selection Index (your externally-estimated GCA composite, from ASReml)
# ------------------------------------------------------------------------
index_df <- read.csv(index_file, stringsAsFactors = FALSE)

stopifnot(all(c(id_col, index_col, family_col) %in% names(index_df)))
index_df[[id_col]] <- as.character(index_df[[id_col]])

cat("Selection Index file:", nrow(index_df), "parents,",
    length(unique(index_df[[family_col]])), "families\n")


# ------------------------------------------------------------------------
# 2. LD block detection + genotype preparation
#    method = "rV2" (kinship-adjusted LD, Mangin et al. 2012): your 180
#    parents come from multiple families, so LD should not be estimated
#    assuming an unrelated panel -- rV2 whitens for relatedness before
#    computing LD, which is the metric HapBlockR recommends for
#    family-structured breeding populations (see README, "when to use rV2").
#    Requires the AGHmatrix and ASRgenomics packages.
#
#    Every argument of run_ldx_pipeline() is listed explicitly below, grouped
#    exactly as in its own source (Phasing / Filtering & imputation / LD
#    block detection / Haplotype extraction / General). Values that are just
#    the package default are commented "(default)"; everything else is a
#    deliberate choice for this analysis.
# ------------------------------------------------------------------------
result <- run_ldx_pipeline(
  # -- required ------------------------------------------------------------
  geno_source    = geno_file,
  out_dir        = out_dir,
  out_blocks     = file.path(out_dir, "blocks.csv"),
  out_diversity  = file.path(out_dir, "diversity.csv"),
  out_hap_matrix = file.path(out_dir, "hap_matrix.csv"),
  hap_format     = "numeric",        # (default) additive 0/1/2/NA hap matrix

  # -- Phasing (unused here: geno_file is unphased numeric dosage, not VCF) -
  phase              = FALSE,        # (default) no statistical phasing --
                                      # this is WHY forward simulation
                                      # (ga_vs_ts_simulation(), which needs
                                      # hap1/hap2) is out of scope, see the
                                      # file header
  beagle_jar         = NULL,         # (default) unused when phase = FALSE
  beagle_threads     = 1L,           # (default)
  java_path          = "java",       # (default)
  beagle_java_mem_gb = NULL,         # (default)
  beagle_args        = "",           # (default)
  beagle_ref_panel   = NULL,         # (default)
  beagle_map_file    = NULL,         # (default)
  beagle_chrom       = NULL,         # (default)
  beagle_seed        = NULL,         # (default)
  beagle_burnin      = NULL,         # (default)
  beagle_iterations  = NULL,         # (default)
  beagle_window      = NULL,         # (default)
  beagle_overlap     = NULL,         # (default)

  # -- Filtering & imputation ------------------------------------------------
  maf_cut        = 0.05,             # (default) drop SNPs with MAF < 5%
  impute         = "mean_rounded",   # (default) mean-impute then round to
                                      # nearest dosage class
  min_callrate   = 0.0,              # (default) no per-SNP call-rate filter

  # -- LD block detection -----------------------------------------------------
  CLQcut             = 0.70,         # raised from the 0.5 default -- tighter
                                      # LD threshold for cleaner haploblocks,
                                      # per the README quick-start example
  method             = "rV2",        # kinship-adjusted LD (see note above)
  kin_method         = "chol",       # (default) Cholesky-based kinship
                                      # correction for rV2
  CLQmode            = "Density",    # (default) clique-detection algorithm
  leng               = 200L,         # (default)
  subSegmSize        = 1500L,        # (default)
  clstgap            = 40000L,       # (default) max bp gap within a clique
  split              = FALSE,        # (default) don't split cliques at gaps
  appendrare         = FALSE,        # (default) don't append rare SNPs to
                                      # nearest block
  singleton_as_block = FALSE,        # (default) isolated SNPs are NOT listed
                                      # as their own rows in `blocks` -- but
                                      # they are NOT dropped from geno_matrix
                                      # either (filtering below is MAF/call-
                                      # rate only), so run_haplotype_
                                      # prediction()'s complete_decomposition
                                      # recovers their GEBV afterward (step 4)
  close_gaps_with_snps = TRUE,       # (default) expand adjacent blocks to
                                      # absorb unassigned gap SNPs
  gap_k_rep          = 2L,           # (default)
  checkLargest       = FALSE,        # (default)
  digits             = -1L,          # (default) no rounding of LD values
  n_threads          = 4L,           # raised from the 1L default -- adjust
                                      # to the cores available on your machine
  min_snps_chr       = 10L,          # (default) skip chromosomes with fewer
                                      # than 10 SNPs
  max_bp_distance    = 0L,           # (default) 0 = no sparse distance cap
                                      # (all pairs on a chromosome considered)

  # -- Haplotype extraction ---------------------------------------------------
  min_snps_block   = 3L,             # (default) drop blocks with < 3 SNPs
                                      # from the `blocks` table (recovered as
                                      # singleton pseudo-blocks in step 4)
  top_n            = NULL,           # (default) keep all haplotype alleles
                                      # (no top-n truncation at this stage)
  min_freq         = 0.01,           # (default) drop haplotype alleles with
                                      # frequency < 1%
  scale_hap_matrix = FALSE,          # (default) hap matrix left on its raw
                                      # dosage/count scale

  # -- General ------------------------------------------------------------
  chr             = NULL,            # (default) process all chromosomes
  clean_malformed = FALSE,           # (default) don't attempt to repair
                                      # malformed input rows
  use_bigmemory   = FALSE,           # (default) in-memory matrix (180
                                      # individuals is small; no need for a
                                      # file-backed store)
  bigmemory_path  = tempdir(),       # (default) unused when use_bigmemory=FALSE
  bigmemory_type  = "char",          # (default) unused when use_bigmemory=FALSE
  verbose         = TRUE             # print progress
)

geno_matrix <- result$geno_matrix        # individuals x SNPs, imputed
snp_info    <- result$snp_info_filtered
blocks      <- result$blocks

cat(nrow(blocks), "LD blocks detected across",
    length(unique(snp_info$CHR)), "chromosomes\n")


# ------------------------------------------------------------------------
# 3. ID check -- Selection Index IDs must match geno_matrix row names exactly
# ------------------------------------------------------------------------
missing_geno  <- setdiff(index_df[[id_col]], rownames(geno_matrix))
missing_index <- setdiff(rownames(geno_matrix), index_df[[id_col]])

if (length(missing_geno))
  message(length(missing_geno), " Index_Peru parent(s) have no genotype: ",
          paste(head(missing_geno, 10), collapse = ", "))
if (length(missing_index))
  message(length(missing_index), " genotyped individual(s) have no Selection ",
          "Index value and will be predicted from relationship information ",
          "only: ", paste(head(missing_index, 10), collapse = ", "))


# ------------------------------------------------------------------------
# 4. Haplotype-based decomposition of the Selection Index
#    marker_effect_method = "gblup" (default): fits rrBLUP::kin.blup() on
#    the haplotype-block GRM using your Selection Index as the response,
#    then backsolves per-SNP effects and sums them into per-block local
#    GEBV (Tong et al. 2025). This is a second-stage decomposition of your
#    already-estimated GCA values onto haplotype blocks, not a
#    re-estimation of GCA from raw testcross phenotypes -- appropriate
#    since ASReml already handled the testcross design and missing-cross
#    imputation.
#
#    complete_decomposition = TRUE (default): isolated SNPs / sub-min_snps
#    blocks dropped from `blocks` above are still present in geno_matrix
#    (MAF/call-rate filtering only, not block-membership), so this recovers
#    each one as its own one-SNP "singleton_CHR_POS" pseudo-block, scaled on
#    the same footing as multi-SNP blocks in block_importance -- their GEBV
#    contribution is therefore available to select_top_blocks() /
#    select_parents_ga() below, not lost.
#
#    include_dominance = FALSE (default, left explicit): see the file
#    header's reasoning -- your Selection Index is a GCA (additive-effects)
#    composite, so a dominance relationship matrix has no clean role here.
#
#    Every argument is listed explicitly; all but blues/id_col/blue_col/
#    include_dominance are package defaults for this single-trait GBLUP use
#    case.
# ------------------------------------------------------------------------
pred <- run_haplotype_prediction(
  # -- required --------------------------------------------------------------
  geno_matrix = geno_matrix,
  snp_info    = snp_info,
  blocks      = blocks,
  blues       = index_df,

  # -- phenotype parsing -------------------------------------------------
  id_col      = id_col,               # Index_Peru.csv genotype-ID column
  blue_col    = index_col,            # Index_Peru.csv SelectionIndex column
  blue_cols   = NULL,                 # (default) unused -- single trait via
                                       # blue_col above, not multi-trait

  # -- block / allele filtering -------------------------------------------
  importance_rule = "any",            # (default) a block is "important" if
                                       # ANY trait's scaled variance clears
                                       # importance_threshold (moot here with
                                       # one trait, but stated for clarity)
  top_n           = NULL,             # (default) no top-n truncation of
                                       # haplotype alleles at this stage
  min_freq        = 0.01,             # (default) drop haplotype alleles
                                       # with frequency < 1%
  min_snps        = 3L,               # (default) blocks with < 3 SNPs are
                                       # excluded from `blocks`-derived
                                       # importance rows (recovered as
                                       # singletons below)
  bend            = TRUE,             # (default) bend the haplotype GRM to
                                       # be positive-definite before kin.blup()

  # -- marker-effect model -------------------------------------------------
  marker_effect_method   = "gblup",   # (default) rrBLUP::kin.blup() on the
                                       # haplotype GRM, then backsolve to SNP
                                       # effects -- appropriate for a single
                                       # pre-estimated Selection Index
  complete_decomposition = TRUE,      # (default) keep singleton SNPs' GEBV,
                                       # see note above
  importance_threshold   = 0.9,       # (default) scaled-variance cutoff used
                                       # to flag block_importance$important
                                       # (informational only; step 5 below
                                       # uses select_top_blocks() instead)
  include_dominance      = FALSE,     # (default, kept explicit) -- see file
                                       # header: your Index is GCA/additive
                                       # by construction

  # -- Bayesian-method-only arguments (unused under marker_effect_method =
  #    "gblup"; listed for completeness since they are real arguments) ------
  n_iter          = 6000L,            # (default) BGLR iterations (BayesA/B/C)
  burn_in         = 1000L,            # (default) BGLR burn-in (BayesA/B/C)

  # -- general ---------------------------------------------------------------
  seed            = NULL,             # (default) no fixed seed
  verbose         = TRUE,             # print progress
  ploidy          = 2L                # rice is diploid (default)
)

cat(pred$n_train, "parents used to train the model,",
    pred$n_predict, "predicted from relationship information alone\n")
cat(sum(pred$block_importance$important, na.rm = TRUE),
    "of", pred$n_blocks, "blocks flagged important",
    "(scaled variance >= 0.9)\n")


# ------------------------------------------------------------------------
# 5. K-fold cross-validation of the prediction model itself
#    Before ranking blocks/parents on top of `pred`, this asks how much to
#    trust it: k-fold predictive ability of the SAME haplotype-GRM GBLUP
#    model used in step 4, holding out 1/k of parents' Selection Index
#    values each fold and correlating predicted vs. observed. Low PA here
#    is a warning sign for every downstream decision, not just this step.
#    Every argument listed explicitly; only geno_matrix/snp_info/blocks/
#    blues/id_col/blue_col are specific to this analysis.
# ------------------------------------------------------------------------
cv <- cv_haplotype_prediction(
  # -- required --------------------------------------------------------------
  geno_matrix = geno_matrix,
  snp_info    = snp_info,
  blocks      = blocks,
  blues       = index_df,

  # -- CV design ---------------------------------------------------------
  k           = 5L,                   # (default) 5-fold CV
  n_rep       = 1L,                   # (default) single fold assignment --
                                       # raise for a more stable PA estimate
                                       # at the cost of runtime (k x n_rep
                                       # total model fits)

  # -- haplotype-allele filtering (independent of step 4's own min_freq --
  #    this function's own default is 0.05, not 0.01; listed explicitly so
  #    that difference isn't silently hidden) ---------------------------
  top_n       = NULL,                 # (default) all alleles above min_freq
  min_freq    = 0.05,                 # (default -- note: NOT the same
                                       # default as run_haplotype_prediction()
                                       # above, which uses 0.01)
  min_snps    = 3L,                   # (default)

  # -- phenotype parsing -------------------------------------------------
  id_col      = id_col,               # same Index_Peru.csv ID column as step 4
  blue_col    = index_col,            # same Index_Peru.csv index column
  blue_cols   = NULL,                 # (default) unused -- single trait

  # -- general ---------------------------------------------------------------
  seed        = 42L,                  # (default)
  verbose     = TRUE
)

cat("\nCross-validated predictive ability (Pearson r, out-of-fold):\n")
print(cv$pa_mean)
if (any(cv$pa_mean$PA < 0.3, na.rm = TRUE))
  warning("Predictive ability below 0.3 for at least one trait/fold summary ",
          "-- downstream block/parent rankings should be treated with ",
          "extra caution.", call. = FALSE)


# ------------------------------------------------------------------------
# 6. Visual check of block importance before choosing a cutoff
# ------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  funnel_plot <- plot_block_funnel(
    local_gebv          = pred$local_gebv,
    block_importance    = pred$block_importance,
    highlight_threshold = 0.9,     # (default) matches importance_threshold
                                    # used above, for a consistent cutoff line
    max_blocks          = 2000L    # (default) plot at most 2000 blocks
  )
  print(funnel_plot)
  ggplot2::ggsave(file.path(out_dir, "block_funnel.png"), funnel_plot,
                   width = 7, height = 5, dpi = 300)
}


# ------------------------------------------------------------------------
# 7. Keep the smallest top-ranked set of blocks explaining 90% of the
#    variance in local GEBV. Data-driven; cross-check against the funnel
#    plot above and swap for a fixed count/percentage if you prefer, e.g.
#    select_top_blocks(pred$block_importance, n = 50, perc_total = NULL,
#    perc_of_total_var = NULL, var_col = "var_scaled").
# ------------------------------------------------------------------------
top_blocks <- select_top_blocks(
  block_importance   = pred$block_importance,
  n                  = NULL,           # not using fixed-count mode
  perc_total         = NULL,           # not using fixed-percentage-of-blocks
                                        # mode
  perc_of_total_var  = 0.90,           # keep smallest top-ranked set
                                        # explaining >= 90% cumulative variance
  var_col            = "var_scaled"    # explicit -- this is what var_col =
                                        # NULL resolves to internally, since
                                        # block_importance carries a
                                        # var_scaled column
)

cat(nrow(top_blocks), "of", pred$n_blocks,
    "blocks retained (top", round(100 * max(top_blocks$cum_var_share), 1),
    "% cumulative variance)\n")

value_matrix <- pred$local_gebv[, top_blocks$block_id, drop = FALSE]


# ------------------------------------------------------------------------
# 7b. Truncation-selection baseline, with the merit floor now applied
#     natively by truncation_selection() itself (min_sel_value/min_sel_mode
#     are package arguments -- see ?truncation_selection).
# ------------------------------------------------------------------------
ts_sel <- truncation_selection(
  score         = pred$gebv,
  n_founders    = n_founders,
  min_sel_value = min_sel_value,
  min_sel_mode  = min_sel_mode
)

cat("Merit-floor cutoff applied:", round(ts_sel$cutoff, 4),
    "(mode =", min_sel_mode, ")\n")


# ------------------------------------------------------------------------
# 7c. Force in the top n_must_include performers by whole-genome merit.
#
#     select_parents_ga()'s fitness function only rewards being top-1/top-2
#     AT A TARGET BLOCK -- it has no term for "generally excellent
#     everywhere," even with the merit floor applied (the floor excludes
#     poor performers, it does not force in the best ones). A candidate
#     whose high value comes from being solidly good across many blocks,
#     without topping any single high-variance target block, contributes
#     ZERO marginal fitness and can be left out entirely -- even the single
#     best individual in the whole panel. This is not a package feature
#     (select_parents_ga() has no "must include" argument), so it stays
#     script-level: force the top performers in, then let the GA search only
#     the remaining pool for the complementary founders.
#
#     `merit_ok` (the full merit-floor-passing candidate pool) is also the
#     candidate set fed to the mate-allocation stages (8-12) below -- those
#     are independent of the specific GA shortlist, since OCS's job is
#     precisely to decide who contributes and how much from a broad pool,
#     not to refine an already-shrunk shortlist.
# ------------------------------------------------------------------------
merit_ok <- names(pred$gebv)[is.finite(pred$gebv) & pred$gebv >= ts_sel$cutoff]

if (n_must_include > 0) {
  must_include <- names(sort(pred$gebv[merit_ok], decreasing = TRUE))[seq_len(n_must_include)]
  cat("\nForced into the final set by whole-genome merit (n_must_include =",
      n_must_include, "):", paste(must_include, collapse = ", "), "\n")
} else {
  must_include <- character(0)
}
n_founders_ga <- n_founders - length(must_include)
if (n_founders_ga < 1)
  stop("n_must_include (", n_must_include, ") leaves no founders for the GA ",
       "to search for. Lower n_must_include or raise n_founders.")

value_matrix_remaining <- value_matrix[setdiff(rownames(value_matrix), must_include), ,
                                       drop = FALSE]

cat(length(merit_ok), "of", length(pred$gebv),
    "individuals clear the merit floor and form the candidate pool for",
    "stages 8-12 (mate allocation) below\n")


# ------------------------------------------------------------------------
# 8. GA-based founder-SET selection. Every select_parents_ga() argument is
#    listed explicitly. merit_score/min_sel_value/min_sel_mode apply the
#    same floor as step 7b (natively, package-level); n_reps runs multiple
#    replicates and reports stability -- see ?select_parents_ga's "GA
#    rigour" section. must_include from step 7c is prepended afterward.
# ------------------------------------------------------------------------
ga_sel_search <- select_parents_ga(
  value_matrix   = value_matrix_remaining,
  n_founders     = n_founders_ga,
  strategy       = strategy,           # "no_selfing" -- see USER CONFIG note
  block_weights  = top_blocks$var_scaled, # weight each retained block by its
                                        # own scaled variance rather than
                                        # equally (package default: NULL ->
                                        # equal weight 1 for every block)
  top_candidates = NULL,               # (default) search the full remaining
                                        # candidate pool, no further prefilter
  popSize        = 100L,               # (default) GA::ga() population size
  maxiter        = 200L,               # (default) GA::ga() max generations
  run            = 50L,                # (default) stop if best fitness hasn't
                                        # improved for this many generations
  pmutation      = 0.1,                # (default) GA::ga() mutation rate
  pcrossover     = 0.8,                # (default) GA::ga() crossover rate
  penalty_weight = NULL,               # (default) auto-scaled from
                                        # block_weights (10 * sum(block_weights)
                                        # / n_blocks)
  seed           = 1L,                 # fixed for reproducibility (default:
                                        # NULL, no fixed seed) -- also seeds
                                        # replicates 2..n_reps deterministically
  verbose        = FALSE,              # (default) suppress GA::ga() iteration
                                        # monitor
  merit_score    = pred$gebv,          # whole-genome value for the floor below
  min_sel_value  = min_sel_value,      # same floor as truncation_selection()
                                        # in step 7b, for an apples-to-apples
                                        # comparison
  min_sel_mode   = min_sel_mode,
  n_reps         = n_reps              # GA replicates; see USER CONFIG note
)

cat("\nGA search converged:", ga_sel_search$converged,
    "| fitness range across", n_reps, "replicates:",
    paste(round(ga_sel_search$stability$fitness_range, 3), collapse = " to "), "\n")
cat("Candidates selected in every replicate (selection_freq == 1):\n")
print(names(ga_sel_search$stability$selection_freq)[
  ga_sel_search$stability$selection_freq == 1])

# Combine: must_include parents are added by merit, not by GA search, so
# ga_sel$fitness/per_block/ga_fit/stability below describe the complementary
# search over value_matrix_remaining only, not the full final set.
ga_sel <- ga_sel_search
ga_sel$selected <- c(must_include, ga_sel_search$selected)

cat("\nGA-selected parents (", length(ga_sel$selected), "):\n", sep = "")
print(ga_sel$selected)
cat("\nTruncation-selected parents (", length(ts_sel$selected), "):\n", sep = "")
print(ts_sel$selected)
cat("\nOverlap between the two sets:",
    length(intersect(ga_sel$selected, ts_sel$selected)), "parents\n")


# ------------------------------------------------------------------------
# 9. Family balance check
#    select_parents_ga() has NO built-in family/pedigree constraint -- it
#    optimises purely on haplotype-block coverage. Since your 180 parents
#    span multiple families, review the resulting set here before
#    finalising crosses; re-run with a different n_founders or manually
#    swap a parent if one family dominates. Reused below (stage 12b) for
#    the OCS- and core-collection-derived sets too.
# ------------------------------------------------------------------------
family_lookup <- stats::setNames(index_df[[family_col]], index_df[[id_col]])

family_summary <- function(selected_ids, label) {
  tab <- sort(table(family_lookup[selected_ids]), decreasing = TRUE)
  cat("\n", label, " -- family representation:\n", sep = "")
  print(tab)
  top_share <- max(tab) / length(selected_ids)
  if (top_share > 0.5)
    warning(label, ": a single family ('", names(tab)[1], "') makes up ",
            round(100 * top_share), "% of the selected parents.",
            call. = FALSE)
  invisible(tab)
}

family_summary(ga_sel$selected, "GA selection")
family_summary(ts_sel$selected, "Truncation selection")


# ------------------------------------------------------------------------
# 10. Diversity-space visualisation of the GA and truncation selections.
#     plot_parent_selection_pca() is GA-vs-TS-specific by design (two fixed
#     group labels), so the OCS/core-collection sets are checked separately
#     in stage 16 below (which is not limited to 2 groups), not forced into
#     this same two-group plot.
#
#     Two views, both listing every argument explicitly:
#       (a) G = pred$G          -- genome-wide diversity space (as before).
#       (b) feature_matrix = value_matrix -- the LITERAL target-block local-
#           GEBV space select_parents_ga() searched (stage 8). (a) is a
#           useful cross-check ("do these picks also look diversity-
#           covering genome-wide?") but is not the space GA was actually
#           asked to spread out in; (b) answers that directly. See
#           ?plot_parent_selection_pca's "feature_matrix" argument.
# ------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  pca_plot <- plot_parent_selection_pca(
    G              = pred$G,
    ga_selected    = ga_sel$selected,
    ts_selected    = ts_sel$selected,
    feature_matrix = NULL              # (default) -- using G above, not
                                        # feature_matrix, for this call
  )
  print(pca_plot)
  ggplot2::ggsave(file.path(out_dir, "parent_selection_pca.png"), pca_plot,
                   width = 7, height = 5.5, dpi = 300)

  pca_plot_blocks <- plot_parent_selection_pca(
    G              = NULL,             # exactly one of G/feature_matrix --
                                        # using feature_matrix here instead
    ga_selected    = ga_sel$selected,
    ts_selected    = ts_sel$selected,
    feature_matrix = value_matrix      # the actual space GA searched --
                                        # see stage 7/8 above
  )
  print(pca_plot_blocks)
  ggplot2::ggsave(file.path(out_dir, "parent_selection_pca_blockspace.png"),
                   pca_plot_blocks, width = 7, height = 5.5, dpi = 300)
}


# ==========================================================================
# MATE ALLOCATION -- stages 11-15
#
# Everything above (6-10) answers "which individuals should be parents"
# (founder-SET selection). Everything below answers "how much should each
# parent contribute, and which specific matings should be made" -- the
# actual Optimal Contribution Selection problem this program already solves
# manually in AlphaMate today. All four stages operate on `merit_ok` (the
# full merit-floor-passing candidate pool from step 7c), independently of
# the specific 30-parent GA/truncation shortlists above.
# ==========================================================================

# ------------------------------------------------------------------------
# 11. Merit-vs-diversity Pareto frontier, to inform target_degree_ocs below
#     BEFORE committing to a single value. Sweeps select_parents_ga()'s
#     coancestry_weight across a grid and Pareto-filters the results --
#     shows the actual gain-vs-diversity tradeoff curve for this candidate
#     pool, rather than guessing target_degree_ocs and hoping it was right.
#     Every argument listed explicitly.
# ------------------------------------------------------------------------
pareto_res <- select_parents_pareto(
  value_matrix       = value_matrix[merit_ok, , drop = FALSE],
  n_founders         = n_founders,
  strategy           = strategy,                    # same crossing-scheme
                                                      # assumption as stage 8
  block_weights      = top_blocks$var_scaled,        # same weighting as stage 8
  top_candidates     = NULL,                         # (default)
  G                  = pred$G[merit_ok, merit_ok, drop = FALSE],
  coancestry_weights = c(0, 0.25, 0.5, 1, 2, 4),     # (default) broad first pass
  merit              = pred$gebv[merit_ok],          # whole-genome merit axis
                                                      # for the frontier (not
                                                      # the GA's own block-
                                                      # coverage fitness)
  popSize            = 100L,                         # (default)
  maxiter            = 200L,                         # (default)
  run                = 50L,                          # (default)
  pmutation          = 0.1,                          # (default)
  pcrossover         = 0.8,                          # (default)
  n_reps             = 3L,                           # (default; lower than
                                                      # stage 8's 5, to keep
                                                      # the 6-grid-point sweep's
                                                      # total runtime reasonable)
  seed               = 1L,
  verbose            = TRUE
)

cat("\nMerit-vs-diversity Pareto frontier ",
    "(review before finalising target_degree_ocs):\n", sep = "")
print(pareto_res$frontier[, c("coancestry_weight", "mean_merit",
                              "mean_relationship", "n_selected",
                              "pareto_optimal")])


# ------------------------------------------------------------------------
# 12. True Optimal Contribution Selection (Meuwissen 1997): the actual
#     production mate-allocation step -- contributions AND a specific
#     mating plan under an explicit relatedness cap, replacing manual
#     AlphaMate use with a reproducible in-R call (or still calling
#     AlphaMate itself, if alphamate_exe is set -- see USER CONFIG).
#     `family` is reporting/labelling only (neither engine enforces family
#     representation as a hard constraint -- see stage 15 below, where
#     family_summary(ocs_parents, ...) runs the post-hoc check). Every
#     argument listed explicitly.
# ------------------------------------------------------------------------
ocs_res <- select_parents_ocs(
  merit                  = pred$gebv[merit_ok],
  G                      = pred$G[merit_ok, merit_ok, drop = FALSE],
  family                 = family_lookup[merit_ok],
  engine                 = ocs_engine,                # "simplemating" -- see USER CONFIG
  n_crosses              = n_crosses_ocs,
  n_parents_max          = NULL,                      # (default) not capped --
                                                        # engine = "simplemating"
                                                        # cannot enforce this
                                                        # directly anyway (see
                                                        # ?select_parents_ocs)
  max_contrib_per_parent = max_contrib_per_parent_ocs,
  allow_selfing          = FALSE,                     # (default) matches
                                                        # strategy = "no_selfing"
                                                        # above
  allow_repeated_matings = FALSE,                     # (default)
  target_degree          = target_degree_ocs,         # see stage 11's frontier
  rescale_nrm            = TRUE,                      # (default) rescale G to
                                                        # a numerator-relationship-
                                                        # matrix-like mean before
                                                        # either engine uses it
  alphamate_exe          = alphamate_exe,             # NULL = auto-detect; see
                                                        # USER CONFIG
  out_dir                = file.path(out_dir, "ocs"),
  seed                   = 1L,
  verbose                = TRUE
)

cat("\nOCS engine used:", ocs_res$engine_used,
    "| constraints satisfied:", ocs_res$ok, "\n")
cat(nrow(ocs_res$mating_plan), "of", n_crosses_ocs, "requested crosses produced\n")
print(head(ocs_res$mating_plan, 10))

ocs_parents <- if (!is.null(ocs_res$contributors) && nrow(ocs_res$contributors))
  ocs_res$contributors$id else
  unique(c(ocs_res$mating_plan$parent1, ocs_res$mating_plan$parent2))

if (!length(ocs_parents))
  stop("select_parents_ocs() produced no usable mating plan -- check the ",
       "engine's console output above (out_dir = ",
       file.path(out_dir, "ocs"), ") before continuing to stages 13-14.")

cat(length(ocs_parents), "distinct parents contribute to the OCS mating plan\n")


# ------------------------------------------------------------------------
# 13. Cross ranking among the OCS-selected parents by predicted mid-parent
#     value plus segregation variance (Schnell & Utz 1975; Bernardo 2003).
#     variance_model = "block_independent": the only mode this pipeline's
#     data supports -- no phased hap1/hap2 (see file header), and the
#     Selection Index is a composite score rather than strictly 0/2/NA
#     homozygous dosage, so "simplemating" does not apply either. Every
#     argument listed explicitly, including ones unused under
#     "block_independent" (kept visible rather than silently omitted, per
#     this script's own convention).
# ------------------------------------------------------------------------
uc_res <- usefulness_criterion(
  # -- candidate crosses ---------------------------------------------------
  parent_ids           = ocs_parents,     # all pairwise crosses among the
                                           # OCS-selected parents, generated
                                           # via combn()
  cross_pairs           = NULL,           # (default) -- overridden by
                                           # parent_ids above
  gebv                  = pred$gebv,      # whole-genome merit, mid-parent term

  # -- selection-intensity parameters --------------------------------------
  selected_proportion   = selected_proportion_uc,
  n_progeny             = NULL,           # (default) classical infinite-
                                           # population asymptotic i_sel;
                                           # supply a realistic per-cross
                                           # progeny count here for the
                                           # finite-population Monte Carlo
                                           # correction instead
  n_sim                 = 20000L,         # (default) unused since
                                           # n_progeny = NULL
  seed                  = 1L,

  # -- variance model --------------------------------------------------------
  variance_model         = "block_independent",
  block_importance       = top_blocks,    # same top-ranked blocks as stages
                                           # 7-11
  block_ids              = NULL,          # (default) use every block in
                                           # top_blocks as-is
  local_gebv              = pred$local_gebv,
  segregation_factor      = 0.5,          # (default)

  # -- phased/linked-only arguments (unused under "block_independent") -----
  haplotypes              = NULL,         # (default)
  snp_info                = NULL,         # (default)
  snp_effects              = NULL,        # (default)

  # -- simplemating-only arguments (unused here) ---------------------------
  geno_matrix              = NULL,        # (default)
  G                         = NULL,       # (default)
  genetic_map               = NULL,       # (default)
  ld_matrix                  = NULL,      # (default)
  type                        = "RIL",    # (default)
  generation                  = 1L,       # (default)
  n_threads                   = 1L,       # (default)

  # -- linked-only argument (unused here) ------------------------------------
  n_sim_linked                 = 2000L,   # (default)

  verbose                       = TRUE
)

cat("\nTop 10 candidate crosses among OCS-selected parents by UC:\n")
print(head(uc_res, 10))


# ------------------------------------------------------------------------
# 13b. Cross-check: same OCS-selected parents' candidate crosses, ranked
#      instead by variance_model = "simplemating" (Peixoto et al. 2024,
#      SimpleMating::getUsefA()) -- a genuinely linkage-aware, genome-wide
#      multi-locus Mendelian-sampling model, rather than stage 13's
#      independent-per-target-block sum.
#
#      Two real data gaps, handled rather than glossed over:
#        - snp_effects: not available from this program's phenotype-only
#          Selection Index (no marker-effect model was fit directly) --
#          backsolved instead from the whole-genome GEBV already computed
#          in stage 4, via backsolve_snp_effects() (Tong et al. 2025). Same
#          scale as pred$gebv, no re-fitting needed.
#        - genetic_map: this program has no real genetic map, only physical
#          bp positions (snp_info$POS) -- so ld_matrix (an LD-based
#          recombination-fraction PROXY, not a validated genetic distance;
#          see ?usefulness_criterion) is supplied instead of genetic_map.
#      SNP set restricted to the SAME target-block SNPs as the rest of this
#      script (top_blocks, stage 7) purely for computational tractability --
#      a genome-wide SNP x SNP LD matrix would be far too large to build or
#      hold in memory. This restriction is a practical necessity for THIS
#      script, not a requirement of simplemating itself, which is otherwise
#      a genome-wide model unrelated to HapBlockR's own target blocks.
#
#      Real constraint, handled automatically rather than worked around:
#      SimpleMating::getUsefA() requires geno_matrix coded strictly 0/2/NA
#      (fully homozygous DH/RIL-style calls). This program's genotypes are
#      012-coded and may still carry residual heterozygosity (dosage = 1)
#      -- usefulness_criterion()'s default het_to_na = TRUE (left at its
#      default below) treats any such calls as missing (NA) rather than
#      erroring or fabricating a homozygous call, and reports how many via
#      message() when verbose = TRUE. The tryCatch() below is kept as a
#      safety net for any OTHER failure mode (e.g. SimpleMating not
#      installed, an unrecognised dosage value, or a version mismatch) --
#      caught and reported rather than aborting the rest of this script
#      (stage 13 above and stage 14 below do not depend on this cross-check
#      succeeding).
# ------------------------------------------------------------------------
target_snp_ids <- unique(unlist(lapply(seq_len(nrow(top_blocks)), function(i) {
  snp_info$SNP[snp_info$CHR == top_blocks$CHR[i] &
               snp_info$POS >= top_blocks$start_bp[i] &
               snp_info$POS <= top_blocks$end_bp[i]]
})))
cat("\n", length(target_snp_ids), " SNPs across the ", nrow(top_blocks),
    " target blocks used for the simplemating cross-check below\n", sep = "")

snp_effects_bs <- backsolve_snp_effects(
  geno_matrix = geno_matrix,   # individuals x SNPs, from stage 2 -- ALL
                                # SNPs; backsolve_snp_effects() itself
                                # doesn't need restricting to target blocks
  gebv        = pred$gebv,     # whole-genome GEBV -- Tong et al. (2025)
                                # backsolves per-SNP effects from these
                                # directly, without refitting a marker model
  G           = pred$G,        # same GRM already used throughout this
                                # script (not recomputed)
  ploidy      = 2L             # (default) diploid rice
)

ld_matrix_targets <- compute_r2(
  geno_matrix[, target_snp_ids, drop = FALSE],  # individuals x target-block
                                                 # SNPs only (see note above)
  digits    = 6L,
  n_threads = 1L                                # raise if you have many
                                                 # target-block SNPs and
                                                 # multiple cores available
)
dimnames(ld_matrix_targets) <- list(target_snp_ids, target_snp_ids)

uc_res_simplemating <- tryCatch(
  usefulness_criterion(
    # -- candidate crosses ----------------------------------------------------
    parent_ids           = ocs_parents,     # SAME parent set and combn()
                                             # order as stage 13, so the two
                                             # results' parent1/parent2
                                             # columns line up directly for
                                             # the rank comparison below
    cross_pairs           = NULL,           # (default)
    gebv                  = pred$gebv,

    # -- selection-intensity parameters ----------------------------------------
    selected_proportion   = selected_proportion_uc,
    n_progeny             = NULL,           # (default)
    n_sim                 = 20000L,         # (default) unused since
                                             # n_progeny = NULL
    seed                  = 1L,

    # -- variance model ----------------------------------------------------------
    variance_model         = "simplemating",
    block_importance       = NULL,          # (default) not used by
                                             # simplemating
    block_ids              = NULL,          # (default) not used by
                                             # simplemating
    local_gebv              = NULL,         # (default) not used by
                                             # simplemating
    segregation_factor      = 0.5,          # (default) not used by
                                             # simplemating

    # -- phased-only argument (unused here) -----------------------------------
    haplotypes              = NULL,         # (default)
    snp_info                = NULL,         # (default) simplemating uses
                                             # genetic_map/ld_matrix instead,
                                             # not snp_info -- see
                                             # ?usefulness_criterion
    snp_effects              = snp_effects_bs, # backsolved above -- REQUIRED
                                             # by simplemating (this argument
                                             # is phased/linked-only in
                                             # stage 13 above, but is also
                                             # required here)

    # -- simplemating arguments --------------------------------------------------
    geno_matrix              = geno_matrix[, target_snp_ids, drop = FALSE],
    het_to_na                 = TRUE,     # (default) auto-convert heterozygous
                                           # (dosage = 1) calls to NA rather
                                           # than erroring -- see this stage's
                                           # header comment; set FALSE to
                                           # error instead
    G                         = pred$G,
    genetic_map               = NULL,       # no real genetic map for this
                                             # program -- see note above
    ld_matrix                  = ld_matrix_targets,
    type                        = "RIL",    # rice, advanced by selfing --
                                             # matches this program's actual
                                             # population structure
    generation                  = 1L,       # (default)
    n_threads                   = 1L,       # (default)

    # -- linked-only argument (unused here) --------------------------------------
    n_sim_linked                 = 2000L,   # (default)

    verbose                       = TRUE
  ),
  error = function(e) {
    message("[stage 13b] variance_model = 'simplemating' could not be run: ",
            conditionMessage(e),
            "\nHeterozygous calls are handled automatically (het_to_na = ",
            "TRUE above), so this is NOT the usual cause -- more likely ",
            "SimpleMating is not installed/too old, an unrecognised dosage ",
            "value exists in geno_matrix, or a SNP-ID mismatch between ",
            "geno_matrix/snp_effects/ld_matrix. See the message above for ",
            "the specific cause. Skipping this cross-check; stage 13's ",
            "block_independent ranking above and stage 14's exact ",
            "validation below are unaffected.")
    NULL
  }
)

if (!is.null(uc_res_simplemating)) {
  cat("\nTop 10 candidate crosses among OCS-selected parents by UC ",
      "(simplemating cross-check):\n", sep = "")
  print(head(uc_res_simplemating, 10))

  rank_agree <- merge(
    uc_res[, c("parent1", "parent2", "rank")],
    uc_res_simplemating[, c("parent1", "parent2", "rank")],
    by = c("parent1", "parent2"), suffixes = c("_indep", "_simplemating")
  )
  if (nrow(rank_agree) >= 3L) {
    rho <- suppressWarnings(stats::cor(
      rank_agree$rank_indep, rank_agree$rank_simplemating,
      method = "spearman", use = "complete.obs"
    ))
    cat("\nRank agreement (Spearman's rho) between block_independent and ",
        "simplemating cross rankings: ", round(rho, 3), " (over ",
        nrow(rank_agree), " crosses in common; closer to 1 means the two ",
        "models agree on which crosses look best; closer to 0 means they ",
        "disagree and the choice of variance model matters for this ",
        "program's decision)\n", sep = "")
  }
}


# ------------------------------------------------------------------------
# 13c. Cross-check: variance_model = "linked" -- NOT run here, by data
#      necessity rather than choice.
#
#      "linked" now accepts ld_matrix as an alternative to genetic_map
#      (added this cycle -- see NEWS.md), which removes the "this program
#      has no real genetic map" barrier that used to block it entirely.
#      It does NOT remove the OTHER requirement: phased hap1/hap2
#      haplotypes (same as variance_model = "phased"), which this script's
#      data source still does not provide (unphased numeric dosage,
#      phase = FALSE throughout -- see this file's header). Population-
#      level LD (what ld_matrix carries) cannot substitute for phase
#      (which specific alleles a given individual's two chromosomes carry)
#      -- there is no proxy for that second requirement the way there now
#      is for the genetic-map one.
#
#      If this program later phases its genotypes (e.g. via
#      phase_with_beagle() / run_ldx_pipeline(phase = TRUE) -- see this
#      file's header's "forward-in-time simulation" note and README's
#      phasing section), re-run stages 2-13b above on the phased data, then
#      the call below becomes usable as-is (every argument already listed
#      explicitly; ld_matrix_targets and snp_effects_bs reused from stage
#      13b above, no need to recompute them):
#
#        uc_res_linked <- usefulness_criterion(
#          parent_ids          = ocs_parents,
#          gebv                = pred$gebv,
#          selected_proportion = selected_proportion_uc,
#          seed                = 1L,
#          variance_model      = "linked",
#          block_importance    = top_blocks,
#          haplotypes          = <phased extract_haplotypes() output>,
#          snp_info            = snp_info,
#          snp_effects         = snp_effects_bs,    # from stage 13b
#          ld_matrix           = ld_matrix_targets, # from stage 13b -- or a
#                                                    # real genetic_map
#                                                    # instead, if you
#                                                    # obtain one (preferred
#                                                    # when available -- see
#                                                    # ?usefulness_criterion)
#          n_sim_linked        = 2000L,
#          verbose             = TRUE
#        )
# ------------------------------------------------------------------------
cat("\n[stage 13c] variance_model = 'linked' skipped: this program's ",
    "genotypes are unphased (phase = FALSE), and phasing is the one ",
    "requirement the new ld_matrix fallback does NOT remove. See the ",
    "commented-out call in this stage's header comment for the exact code ",
    "to run once phased haplotypes are available.\n", sep = "")


# ------------------------------------------------------------------------
# 14. Exact ILP validation of the OCS mating plan's optimality gap.
#     Solves the SAME cross-selection problem (n_crosses_ocs crosses,
#     max_contrib_per_parent_ocs per parent) exactly via binary integer
#     linear programming, over the UC-ranked candidate list from stage 13 --
#     a small enough set (C(length(ocs_parents), 2) candidates) for exact
#     solving to be practical. Every argument listed explicitly.
# ------------------------------------------------------------------------
exact_res <- validate_crosses_exact(
  data                = uc_res,
  n_cross              = n_crosses_ocs,       # same plan size as stage 12,
                                               # for an apples-to-apples
                                               # comparison
  max_cross            = max_contrib_per_parent_ocs,
  culling_pairwise_k   = NULL,                # (default) no additional
                                               # culling here -- OCS's
                                               # target_degree already shaped
                                               # the candidate parent pool
                                               # upstream (stage 12)
  parent1_col          = "parent1",           # (default) -- usefulness_
                                               # criterion()'s own column name
  parent2_col          = "parent2",           # (default)
  criterion_col        = "UC",                # (default) -- maximise UC,
                                               # not just mid-parent value
  relatedness_col      = NULL,                # (default)
  G                    = NULL,                # (default) unused since
                                               # culling_pairwise_k = NULL
  heuristic_plan       = ocs_res$mating_plan, # compare OCS's heuristic plan
                                               # against the true optimum
  max_vars             = max_vars_exact,
  verbose              = TRUE
)

cat("\nExact optimum total UC:", round(exact_res$exact_objective, 3), "\n")
if (!is.null(exact_res$heuristic_objective))
  cat("OCS heuristic plan's total UC:", round(exact_res$heuristic_objective, 3),
      "(", round(exact_res$gap_pct, 1), "% below the exact optimum)\n")


# ------------------------------------------------------------------------
# 15. Diversity-first cross-check: select_core_collection() picks n_core
#     individuals maximising genetic diversity itself (not merit-subject-
#     to-a-diversity-cap the way GA/OCS above do), from the SAME merit-
#     floor-passing candidate pool. Useful as an independent sanity check
#     on whether the merit-aware selections (GA/OCS) are eroding the
#     genetic base more than a pure-diversity baseline would. Every
#     argument listed explicitly.
# ------------------------------------------------------------------------
core_res <- select_core_collection(
  G             = pred$G[merit_ok, merit_ok, drop = FALSE],
  n_core        = n_core,
  type          = "relationship",       # (default) G above is a relationship
                                         # matrix, not a distance matrix
  strategy      = "maximin",            # (default) farthest-point greedy
                                         # traversal (Gonzalez 1985) --
                                         # guards against near-duplicate
                                         # individuals; see ?select_core_
                                         # collection for "mean_distance"
  merit         = pred$gebv[merit_ok],  # tie-break/pre-filter only -- this
                                         # function optimises diversity, not
                                         # merit, among eligible candidates
  min_sel_value = NULL,                 # (default) no additional floor here
                                         # -- merit_ok already applied one
                                         # upstream (stage 7c)
  min_sel_mode  = "value",              # (default) unused since
                                         # min_sel_value = NULL
  seed          = 1L,
  verbose       = TRUE
)

cat("\nCore-collection selection (", length(core_res$selected),
    "):\n", sep = "")
print(core_res$selected)
family_summary(core_res$selected, "Core-collection selection")
family_summary(ocs_parents, "OCS mating-plan parents")

cat("\nOverlap: GA vs OCS parents:",
    length(intersect(ga_sel$selected, ocs_parents)), "| ",
    "GA vs core-collection:",
    length(intersect(ga_sel$selected, core_res$selected)), "| ",
    "OCS vs core-collection:",
    length(intersect(ocs_parents, core_res$selected)), "\n")


# ==========================================================================
# STRATEGIC DECISION SUPPORT -- stage 16
#
# Everything above answers "who did each method pick." This stage answers
# a different, program-level question: does EVERY method represent every
# real genetic sub-group of the population, or does it concentrate on just
# one or two? That is precisely the axis this program's two recurring
# scenarios sit on:
#   - Population improvement (broad crossing block, every family
#     represented into the next cycle, OCS in AlphaMate balancing gain vs.
#     diversity) wants HIGH prop_* coverage spread across every cluster
#     below.
#   - Fast-track release (3-4 best families, one best line each, pushed
#     hard toward homozygosity for rapid recurrent genomic selection)
#     deliberately wants selection CONCENTRATED in only the best 1-2
#     clusters.
# This table is how you check which one a given cycle's selections
# actually look like, rather than assuming. Neither pattern is "wrong" --
# they answer different program goals.
#
# plot_parent_selection_pca() (stage 10) only ever visualises PC1/PC2, a
# potentially small fraction of total variance in a genetically complex,
# multi-family population. cluster_selection_groups() instead retains
# enough PCs for a real cumulative-variance threshold, clusters the WHOLE
# merit-eligible candidate pool, and cross-tabulates ALL FOUR selection
# strategies used above (GA, TS, OCS, core-collection) against that
# structure at once -- something the 2-group plot_parent_selection_pca()
# cannot do. Both available clustering algorithms are run and reported
# side by side, not just one: hierarchical (Ward, deterministic) and
# k-means (needs a seed) can disagree on borderline individuals, and
# cross-checking both is cheap insurance against over-trusting a single
# algorithm's idiosyncrasies.
# ==========================================================================

# ------------------------------------------------------------------------
# 16. Genetic-cluster representation across GA/TS/OCS/core-collection.
#     Clustered in genome-wide diversity space (G), the same space stage
#     10's family-representation question lives in -- not target-block
#     feature space, since the question here is population-level genetic
#     structure, not GA's specific optimisation objective.
#     Every cluster_selection_groups()/plot_selection_clusters() argument
#     listed explicitly. Individuals selected by more than one strategy
#     (e.g. both GA and OCS) are counted under each independently -- see
#     ?cluster_selection_groups.
# ------------------------------------------------------------------------
selection_groups <- list(
  GA   = ga_sel$selected,
  TS   = ts_sel$selected,
  OCS  = ocs_parents,
  Core = core_res$selected
)
# Defensive: every group above should already be a subset of merit_ok (each
# was itself computed from a merit_ok-restricted G/value_matrix upstream),
# but intersecting explicitly guards against a stale/mismatched variable
# rather than erroring deep inside cluster_selection_groups().
selection_groups <- lapply(selection_groups, intersect, merit_ok)

cluster_results <- list()
for (clust_method in c("hierarchical", "kmeans")) {
  cat("\n--- Genetic-cluster representation (method = '", clust_method,
      "', n_clusters = ", n_clusters_diversity, ") ---\n", sep = "")
  clust_res <- cluster_selection_groups(
    G                  = pred$G[merit_ok, merit_ok, drop = FALSE],
    feature_matrix     = NULL,                    # (default) -- using G,
                                                    # see note above
    groups             = selection_groups,
    variance_threshold = variance_threshold_diversity,
    method             = clust_method,
    n_clusters         = n_clusters_diversity,
    seed               = 1L,                       # only used by "kmeans";
                                                    # ignored (with a
                                                    # message) for
                                                    # "hierarchical"
    verbose            = TRUE
  )
  cluster_results[[clust_method]] <- clust_res
  print(clust_res$table)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_selection_clusters(
      cluster_res = clust_res,
      save_path   = file.path(out_dir,
                              paste0("selection_clusters_", clust_method, ".pdf")),
      width       = 7,
      height      = 5.5
    )
  }
}

# "Cluster 1"/"Cluster 2"/... labels are arbitrary per method (hierarchical's
# "Cluster 1" need not correspond to k-means's "Cluster 1" even when both
# recover the SAME grouping), so compare the induced PARTITIONS, not the raw
# labels: agreement means the cross-tab collapses to a one-to-one mapping
# (every hierarchical cluster maps to exactly one k-means cluster and vice
# versa), i.e. every row and every column of the cross-tab has exactly one
# non-zero cell.
agree_tab <- table(hierarchical = cluster_results$hierarchical$cluster,
                   kmeans       = cluster_results$kmeans$cluster)
same_partition <- all(rowSums(agree_tab > 0) == 1) && all(colSums(agree_tab > 0) == 1)
if (!same_partition)
  message("[stage 16] Hierarchical and k-means disagree on at least one ",
          "individual's cluster grouping -- review both tables above ",
          "rather than trusting either algorithm alone for borderline ",
          "cases.")


# ------------------------------------------------------------------------
# 17. Final outputs for the next crossing cycle
#     Two deliverables: (a) the founder shortlist (GA-selected, stage 8) for
#     continuity with prior cycles' reporting, and (b) the actual, actionable
#     mating plan (OCS-selected parents x specific pairs, UC-ranked, stage
#     12-13) -- the latter is what stage-11's Pareto frontier and stage-14's
#     exact validation were done in service of.
# ------------------------------------------------------------------------
final_parents <- data.frame(
  Genotype        = ga_sel$selected,
  SelectionIndex  = index_df[[index_col]][match(ga_sel$selected, index_df[[id_col]])],
  Family          = index_df[[family_col]][match(ga_sel$selected, index_df[[id_col]])],
  WholeGenomeGEBV = pred$gebv[ga_sel$selected],
  stringsAsFactors = FALSE
)
final_parents <- final_parents[order(-final_parents$SelectionIndex), ]

write.csv(final_parents, file.path(out_dir, "selected_parents_next_cycle.csv"),
          row.names = FALSE)
print(final_parents)

# NOTE: usefulness_criterion() built uc_res's cross_pairs via combn(unique(
# ocs_parents), 2), so a given pair's parent1/parent2 there reflects each ID's
# position in ocs_parents -- NOT necessarily the parent1/parent2 direction
# ocs_res$mating_plan (from SimpleMating's own selectCrosses(), engine =
# "simplemating" per USER CONFIG above) assigned
# to the same unordered pair. A direct merge(by = c("parent1", "parent2")) is
# therefore direction-sensitive and silently drops any pair where the two
# sources disagree on which parent is listed first (in practice this can be
# most rows). Match on an order-independent pair key instead, keeping
# ocs_res$mating_plan's own parent1/parent2 labels for the output.
.pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "___")

ocs_plan_keyed <- ocs_res$mating_plan
ocs_plan_keyed$.key <- .pair_key(ocs_plan_keyed$parent1, ocs_plan_keyed$parent2)

uc_keyed <- uc_res[, c("parent1", "parent2", "mid_parent_gebv",
                       "predicted_variance", "UC", "rank")]
uc_keyed$.key <- .pair_key(uc_keyed$parent1, uc_keyed$parent2)
uc_keyed$parent1 <- NULL
uc_keyed$parent2 <- NULL

final_mating_plan <- merge(ocs_plan_keyed, uc_keyed, by = ".key", all.x = TRUE)
final_mating_plan$.key <- NULL
final_mating_plan <- final_mating_plan[order(final_mating_plan$rank), ]

write.csv(final_mating_plan, file.path(out_dir, "mating_plan_next_cycle.csv"),
          row.names = FALSE)
print(head(final_mating_plan, 20))

cat("\nDone. Outputs written to:", normalizePath(out_dir), "\n",
    "  - selected_parents_next_cycle.csv        (GA founder shortlist, stage 8)\n",
    "  - mating_plan_next_cycle.csv             (OCS mating plan + UC ranking,",
    "stages 12-13)\n",
    "  - parent_selection_pca.png               (genome-wide diversity view, stage 10)\n",
    "  - parent_selection_pca_blockspace.png    (target-block feature-space view, stage 10)\n",
    "  - selection_clusters_hierarchical.pdf    (GA/TS/OCS/Core cluster representation, stage 16)\n",
    "  - selection_clusters_kmeans.pdf          (same, k-means cross-check, stage 16)\n",
    sep = "")
