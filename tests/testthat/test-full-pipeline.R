## tests/testthat/test-full-pipeline.R
## -----------------------------------------------------------------------------
## End-to-end integration test: chains every breeder-facing decision tool this
## package offers into ONE pipeline, mirroring the structure of
## scripts/CIAT_Peru_parent_selection.R (a real production pipeline) but run
## on the package's own small toy dataset (ldx_geno/ldx_snp_info/ldx_blocks/
## ldx_blues) so it stays fast enough for CI.
##
## Every other test-*.R file in this suite unit-tests one function at a time,
## in isolation, against hand-built fixtures. This file instead checks that
## real output from one stage stays USABLE as input to the next, the way an
## actual breeding decision flows in practice:
##
##   prediction (+ CV)  ->  block selection  ->  founder-SET selection
##   (truncation / GA / family-quota, all four family_select_mode variants,
##   both grouping modes)  ->  mate ALLOCATION (Pareto frontier -> OCS ->
##   usefulness criterion -> exact ILP validation)  ->  diversity cross-checks
##   (core collection, multi-strategy cluster representation)
##
## A bug that only shows up when two functions' outputs are chained together
## (an ID subset that silently drops names, a column-name mismatch, a matrix
## that loses its dimnames) would not necessarily be caught by any single
## function's own isolated unit tests -- that is precisely the gap this file
## closes.
##
## Every optional-dependency stage is individually guarded with
## skip_if_not_installed()/skip_if_simplemating_too_old()/
## skip_if_optisel_missing_fn() (the latter two from helper.R), consistent
## with the rest of the suite's convention, so this file degrades gracefully
## rather than failing on a minimal install. rrBLUP is a hard Import (always
## available), so stage 1-4 (prediction through truncation) always run.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

test_that("full pipeline: prediction through mate allocation and diversity checks runs end-to-end and stays internally consistent", {
  data(ldx_geno,     package = "HapBlockR")
  data(ldx_snp_info, package = "HapBlockR")
  data(ldx_blocks,   package = "HapBlockR")
  data(ldx_blues,    package = "HapBlockR")

  ids <- rownames(ldx_geno)

  # Synthetic family labels: 12 families of 10 individuals each, in genotype
  # order. ldx_geno/ldx_blues carry no real pedigree, so this grouping exists
  # purely to exercise family_select_mode with multiple members per family.
  family <- stats::setNames(paste0("Fam", ceiling(seq_along(ids) / 10)), ids)

  # ==========================================================================
  # STAGE 1-2: haplotype-based decomposition of a real BLUE (YLD), + its
  # own k-fold cross-validation
  # ==========================================================================
  pred <- run_haplotype_prediction(
    geno_matrix = ldx_geno,
    snp_info    = ldx_snp_info,
    blocks      = ldx_blocks,
    blues       = ldx_blues,
    id_col      = "id",
    blue_col    = "YLD",
    marker_effect_method = "gblup",
    min_snps    = 3L,
    seed        = 1L,
    verbose     = FALSE
  )
  expect_true(all(c("id", "gebv", "PEV", "reliability",
                    "recommendable") %in% names(pred$gebv_uncertainty)))
  expect_equal(pred$gebv_uncertainty$id, names(pred$gebv))
  expect_type(pred$gebv, "double")
  expect_false(is.null(names(pred$gebv)))
  expect_true(all(ids %in% names(pred$gebv)))
  expect_true(is.matrix(pred$G))
  expect_equal(dim(pred$G), c(length(ids), length(ids)))
  expect_true(is.matrix(pred$local_gebv))
  expect_type(pred$haplotypes, "list")

  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = ldx_blues, id_col = "id", blue_col = "YLD",
    k = 3L, n_rep = 1L, min_snps = 3L, seed = 1L, verbose = FALSE
  )
  expect_true("PA" %in% names(cv$pa_mean))

  # ==========================================================================
  # STAGE 3: keep the smallest top-ranked set of blocks explaining most of
  # the variance in local GEBV -- output feeds every founder-selection stage
  # below via value_matrix
  # ==========================================================================
  top_blocks <- select_top_blocks(pred$block_importance, perc_of_total_var = 0.90)
  expect_true(nrow(top_blocks) >= 1L)
  expect_true(nrow(top_blocks) <= pred$n_blocks)

  value_matrix <- pred$local_gebv[, top_blocks$block_id, drop = FALSE]
  expect_equal(nrow(value_matrix), length(ids))
  expect_equal(rownames(value_matrix), ids)

  n_founders <- 15L

  # ==========================================================================
  # STAGE 4: truncation-selection baseline (no optional dependency)
  # ==========================================================================
  ts_sel <- truncation_selection(pred$gebv, n_founders = n_founders)
  expect_equal(length(ts_sel$selected), n_founders)
  expect_true(all(ts_sel$selected %in% ids))

  # ==========================================================================
  # STAGE 5: GA-based founder-SET selection -- which set jointly covers the
  # most favourable blocks (skipped cleanly if GA is not installed)
  # ==========================================================================
  has_ga <- requireNamespace("GA", quietly = TRUE)
  if (has_ga) {
    ga_sel <- select_parents_ga(
      value_matrix   = value_matrix,
      n_founders     = n_founders,
      strategy       = "no_selfing",
      block_weights  = top_blocks$var_scaled,
      popSize        = 20L,
      maxiter        = 15L,
      run            = 8L,
      seed           = 1L,
      verbose        = FALSE,
      n_reps         = 1L
    )
    expect_equal(length(ga_sel$selected), n_founders)
    expect_true(all(ga_sel$selected %in% ids))
  }

  # ==========================================================================
  # STAGE 6: family- (and genetic-cluster-) quota selection -- ALL FOUR
  # family_select_mode variants, plus both grouping modes and the GBLUP-
  # corrected / haplotype-diversity-aware options
  # ==========================================================================

  # "count": fixed 2 lines from each of the top 4 families
  fam_count <- select_parents_by_family(
    score = pred$gebv, family = family, n_families = 4L, n_per_family = 2L,
    family_select_mode = "count", variance_method = "anova", verbose = FALSE
  )
  expect_equal(length(fam_count$selected), 8L)
  expect_equal(length(unique(family[fam_count$selected])), 4L)
  expect_true(all(fam_count$selected %in% ids))

  # "percentage": top 30% of each of the top 3 families -- self-scaling to
  # each family's own size, unlike "count" above
  fam_pct <- select_parents_by_family(
    score = pred$gebv, family = family, n_families = 3L,
    family_select_mode = "percentage", pct_per_family = 0.3,
    variance_method = "anova", verbose = FALSE
  )
  expect_true(length(unique(family[fam_pct$selected])) <= 3L)
  expect_true(all(fam_pct$selected %in% ids))

  # "sd_threshold": members of retained families scoring at least one SD
  # above the eligible population mean -- a threshold rule, not a fixed
  # count/percentage
  fam_sd <- select_parents_by_family(
    score = pred$gebv, family = family, n_families = 3L,
    family_select_mode = "sd_threshold", sd_threshold = 1,
    variance_method = "anova", verbose = FALSE
  )
  expect_true(all(fam_sd$selected %in% ids))

  # "check_relative": every eligible member scoring at least 5% above a
  # named check individual's own score (a mid-ranked check, not the very
  # best, so the margin is meaningful rather than trivially near-empty)
  check_id <- names(sort(pred$gebv, decreasing = TRUE))[round(length(pred$gebv) / 2)]
  fam_check <- select_parents_by_family(
    score = pred$gebv, family = family, n_families = 4L,
    family_select_mode = "check_relative",
    check_id = check_id, check_margin_pct = 5,
    variance_method = "anova", verbose = FALSE
  )
  expect_true(all(fam_check$selected %in% ids))

  # Family-quota with the GBLUP-corrected family ranking (use_family_
  # relationship = TRUE, needs G) AND a haplotype-diversity tiebreak within
  # each family (needs value_matrix/haplotypes/block_weights) -- the "more
  # robust approach to correct family mean" combination documented in the
  # Breeder's Guide
  fam_robust <- select_parents_by_family(
    score = pred$gebv, family = family, n_families = 3L, n_per_family = 2L,
    family_select_mode = "count", variance_method = "anova",
    use_family_relationship = TRUE, G = pred$G,
    ensure_haplotype_diversity = TRUE, value_matrix = value_matrix,
    haplotypes = pred$haplotypes, block_weights = top_blocks$var_scaled,
    verbose = FALSE
  )
  expect_equal(length(fam_robust$selected), 6L)
  expect_true(all(fam_robust$selected %in% ids))

  # group_by = "genetic_cluster": groups derived directly from G rather than
  # the pedigree/family label above -- the other half of "family- (or
  # genetic-cluster-) quota selection"
  fam_cluster <- select_parents_by_family(
    score = pred$gebv, family = NULL, n_families = 3L, n_per_family = 2L,
    family_select_mode = "count", group_by = "genetic_cluster",
    G = pred$G, n_clusters = 4L, cluster_method = "ward.D2",
    variance_method = "anova", verbose = FALSE
  )
  expect_equal(length(fam_cluster$selected), 6L)
  expect_true(all(fam_cluster$selected %in% ids))

  # ==========================================================================
  # STAGE 7: merit-vs-diversity Pareto frontier -- informs OCS's
  # target_degree below rather than guessing it (skipped cleanly if GA is
  # not installed, same dependency as select_parents_ga())
  # ==========================================================================
  if (has_ga) {
    pareto_res <- select_parents_pareto(
      value_matrix = value_matrix, n_founders = n_founders,
      strategy = "no_selfing", block_weights = top_blocks$var_scaled,
      G = pred$G, coancestry_weights = c(0, 1), merit = pred$gebv,
      popSize = 20L, maxiter = 15L, run = 8L, n_reps = 1L,
      seed = 1L, verbose = FALSE
    )
    expect_true(is.data.frame(pareto_res$frontier))
    expect_true(nrow(pareto_res$frontier) >= 1L)
    expect_true(all(c("mean_merit", "mean_relationship", "pareto_optimal")
                   %in% names(pareto_res$frontier)))
  }

  # ==========================================================================
  # STAGE 8-10: true Optimal Contribution Selection -> cross ranking by
  # usefulness criterion -> exact ILP validation of the OCS heuristic plan.
  # engine = "simplemating" needs BOTH SimpleMating and optiSel installed
  # (see helper.R and test-ocs.R); exact validation additionally needs
  # lpSolve. All skipped cleanly, independently, if unavailable.
  # ==========================================================================
  has_ocs <- requireNamespace("SimpleMating", quietly = TRUE) &&
    requireNamespace("optiSel", quietly = TRUE)
  if (has_ocs) skip_if_simplemating_too_old("selectCrosses")

  ocs_parents <- NULL
  if (has_ocs) {
    ocs_res <- select_parents_ocs(
      merit = pred$gebv, G = pred$G, family = family,
      engine = "simplemating", n_crosses = 10L, max_contrib_per_parent = 3L,
      allow_selfing = FALSE, target_degree = 30, seed = 1L, verbose = FALSE
    )
    expect_true(all(c("mating_plan", "contributors", "engine_used", "ok")
                   %in% names(ocs_res)))
    expect_equal(ocs_res$engine_used, "simplemating")

    if (!is.null(ocs_res$mating_plan) && nrow(ocs_res$mating_plan)) {
      expect_true(all(c("parent1", "parent2") %in% names(ocs_res$mating_plan)))
      ocs_parents <- unique(c(ocs_res$mating_plan$parent1,
                              ocs_res$mating_plan$parent2))
      expect_true(all(ocs_parents %in% ids))

      if (length(ocs_parents) >= 2L) {
        uc_res <- usefulness_criterion(
          parent_ids = ocs_parents, gebv = pred$gebv,
          selected_proportion = 0.1, seed = 1L,
          variance_model = "block_independent",
          block_importance = top_blocks, local_gebv = pred$local_gebv,
          verbose = FALSE
        )
        expect_true(all(c("parent1", "parent2", "UC", "rank") %in% names(uc_res)))
        expect_true(nrow(uc_res) >= 1L)
        expect_equal(uc_res$rank, seq_len(nrow(uc_res)))

        if (requireNamespace("lpSolve", quietly = TRUE) && nrow(uc_res) >= 1L) {
          exact_res <- validate_crosses_exact(
            data = uc_res, n_cross = min(10L, nrow(uc_res)), max_cross = 3L,
            heuristic_plan = ocs_res$mating_plan, verbose = FALSE
          )
          expect_true(is.finite(exact_res$exact_objective))
        }
      }
    }
  }

  # ==========================================================================
  # STAGE 11: diversity-first cross-check -- select_core_collection() picks
  # n_founders individuals maximising genetic diversity itself, from the
  # same candidate pool, as an independent sanity check on the merit-aware
  # selections above (no optional dependency)
  # ==========================================================================
  core_res <- select_core_collection(
    G = pred$G, n_core = n_founders, type = "relationship",
    strategy = "maximin", merit = pred$gebv, seed = 1L, verbose = FALSE
  )
  expect_equal(length(core_res$selected), n_founders)
  expect_true(all(core_res$selected %in% ids))

  # ==========================================================================
  # STAGE 12: genetic-cluster representation across EVERY strategy computed
  # above -- does each method spread its picks across the population's real
  # genetic structure, or concentrate on just one or two clusters?
  # cluster_selection_groups() is not limited to two groups (unlike
  # plot_parent_selection_pca()), so every available strategy is included.
  # ==========================================================================
  selection_groups <- list(
    TS          = ts_sel$selected,
    FamilyQuota = fam_count$selected,
    Core        = core_res$selected
  )
  if (has_ga)                      selection_groups$GA  <- ga_sel$selected
  if (!is.null(ocs_parents))       selection_groups$OCS <- ocs_parents

  clust_res <- cluster_selection_groups(
    G = pred$G, groups = selection_groups,
    variance_threshold = 0.95, method = "hierarchical",
    n_clusters = 3L, verbose = FALSE
  )
  expect_true(is.data.frame(clust_res$table))
  expect_true("Cluster" %in% names(clust_res$table))
  expect_true(all(paste0("n_", names(selection_groups)) %in% names(clust_res$table)))
  expect_equal(sum(clust_res$table$n_total), length(ids))
})
