## tests/testthat/test-parent-selection.R
## -----------------------------------------------------------------------------
## Tests for R/parent_selection.R (HapSelect gap-analysis, part 2):
##   truncation_selection(), select_parents_ga(), select_parents_ga_ts(),
##   plot_parent_selection_pca()
##
## GA parent-selection/plot_parent_selection_pca() tests are individually
## skipped when GA/ggplot2 are unavailable; truncation_selection() has no
## dependency and always runs.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# ==============================================================================
# truncation_selection()
# ==============================================================================

test_that("truncation_selection(): selects the top n_founders by descending score", {
  score <- setNames(c(5, 1, 9, 3, 7, 2, 8, 4, 6, 0), paste0("ind", 1:10))
  res <- truncation_selection(score, n_founders = 3)
  expect_equal(res$selected, c("ind3", "ind7", "ind5"))  # scores 9, 8, 7
  expect_equal(res$score, c(9, 8, 7))
})

test_that("truncation_selection(): unnamed score errors", {
  expect_error(truncation_selection(1:10, n_founders = 3), "named")
})

test_that("truncation_selection(): n_founders < 1 errors", {
  score <- setNames(1:5, paste0("ind", 1:5))
  expect_error(truncation_selection(score, n_founders = 0), "n_founders")
})

test_that("truncation_selection(): NA/non-finite scores are excluded before ranking", {
  score <- setNames(c(5, NA, 9, Inf, 7), paste0("ind", 1:5))
  res <- truncation_selection(score, n_founders = 2)
  expect_equal(res$selected, c("ind3", "ind5"))  # 9, 7 (NA and Inf dropped)
})

test_that("truncation_selection(): warns and returns all available when fewer finite scores than requested", {
  score <- setNames(c(5, NA, 9), paste0("ind", 1:3))
  expect_warning(res <- truncation_selection(score, n_founders = 5), "Only")
  expect_equal(length(res$selected), 2L)
})

# -- min_sel_value / min_sel_mode --------------------------------------------
# score = 1..10 for ind1..ind10, so each mode's cutoff is hand-computable.

test_that("truncation_selection(): min_sel_mode = 'value' keeps scores >= the literal value", {
  score <- setNames(1:10, paste0("ind", 1:10))
  # n_founders=10 exceeds the 5 individuals that clear the floor, so this
  # (correctly) warns "Only 5 individuals clear..."; suppressed here since
  # that warning behaviour is already covered by its own dedicated test above.
  res <- suppressWarnings(truncation_selection(
    score, n_founders = 10L, min_sel_value = 6, min_sel_mode = "value"
  ))
  # expect_setequal(), not sort() + expect_equal(): R's default string sort
  # is lexicographic ("ind10" < "ind6"), so it would not match a
  # numerically-ordered expected vector; selection order doesn't matter here.
  expect_setequal(res$selected, paste0("ind", 6:10))  # scores 6..10
})

test_that("truncation_selection(): min_sel_mode = 'percentile' keeps the top fraction by score", {
  score <- setNames(1:10, paste0("ind", 1:10))
  res <- suppressWarnings(truncation_selection(
    score, n_founders = 10L, min_sel_value = 0.3, min_sel_mode = "percentile"
  ))
  expect_equal(length(res$selected), 3L)  # top 30% of 10 = 3
})

test_that("truncation_selection(): relaxed_pool matches the broad-pool cutoff", {
  score <- setNames(1:10, paste0("ind", 1:10))  # mean = 5.5, sd = ~3.0277
  res <- suppressWarnings(truncation_selection(
    score, n_founders = 10L, min_sel_value = 1, min_sel_mode = "relaxed_pool"
  ))
  cutoff <- mean(score) - 1 * sd(score)
  # expect_setequal(), not sort() + expect_equal() -- see the lexicographic-
  # sort note in the min_sel_mode = 'value' test above.
  expect_setequal(res$selected, paste0("ind", which(score >= cutoff)))
})

test_that("truncation_selection(): sd_above_mean retains superior directional scores", {
  score <- setNames(1:10, paste0("ind", 1:10))
  res <- suppressWarnings(truncation_selection(
    score, n_founders = 10L, min_sel_value = 1,
    min_sel_mode = "sd_above_mean"
  ))
  cutoff <- mean(score) + sd(score)
  expect_equal(res$cutoff, cutoff)
  expect_setequal(res$selected, names(score)[score >= cutoff])

  lower_is_better <- setNames(10:1, paste0("ind", 1:10))
  directional_merit <- -lower_is_better
  res_directional <- suppressWarnings(truncation_selection(
    directional_merit, n_founders = 10L, min_sel_value = 0,
    min_sel_mode = "sd_above_mean"
  ))
  expect_true(all(
    lower_is_better[res_directional$selected] <= mean(lower_is_better)
  ))
})

test_that("truncation_selection(): sd_above_mean requires a non-negative SD distance", {
  score <- setNames(1:10, paste0("ind", 1:10))
  expect_error(
    truncation_selection(
      score, n_founders = 5L, min_sel_value = -1,
      min_sel_mode = "sd_above_mean"
    ),
    "non-negative"
  )
})

test_that("truncation_selection(): min_sel_value with no candidates clearing the floor errors", {
  score <- setNames(1:10, paste0("ind", 1:10))
  expect_error(
    truncation_selection(score, n_founders = 5L,
                         min_sel_value = 100, min_sel_mode = "value"),
    "min_sel_value"
  )
})

# ==============================================================================
# select_parents_ga()
# ==============================================================================

make_vmat <- function(n = 12, p = 5, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), nrow = n, ncol = p,
             dimnames = list(paste0("ind", seq_len(n)), paste0("blk", seq_len(p))))
  m
}

test_that("GA and GA+TS public APIs expose distinct objective arguments", {
  base_args <- names(formals(select_parents_ga))
  hybrid_args <- names(formals(select_parents_ga_ts))
  merit_args <- c(
    "merit_score", "merit_priority", "merit_weight",
    "min_sel_value", "min_sel_mode"
  )

  expect_false(any(merit_args %in% base_args))
  expect_true(all(merit_args %in% hybrid_args))
  expect_lt(length(base_args), length(hybrid_args))
})

test_that("select_parents_ga(): requires GA and errors with a clear message when absent", {
  skip_if(requireNamespace("GA", quietly = TRUE),
         "GA is installed; this test targets the absent-dependency error path only")
  expect_error(select_parents_ga(make_vmat(), n_founders = 4), "GA")
})

test_that("select_parents_ga(): runs and returns the documented structure", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  res <- select_parents_ga(vmat, n_founders = 4, popSize = 20, maxiter = 15,
                           run = 10, seed = 1)
  expect_type(res, "list")
  expect_true(all(c("selected", "fitness", "per_block", "strategy", "ga_fit") %in% names(res)))
  expect_equal(length(res$selected), 4L)
  expect_true(all(res$selected %in% rownames(vmat)))
  expect_equal(nrow(res$per_block), ncol(vmat))
  expect_true(all(c("block_id", "best_value", "contributor_1", "contributor_2") %in%
                 names(res$per_block)))
  expect_equal(res$strategy, "no_selfing")
  expect_identical(res$selection_method, "select_parents_ga")
  expect_identical(res$result_contract$method, "select_parents_ga")
  expect_true(is.na(res$mean_merit))
  expect_equal(res$merit_weight, 0)
})

test_that("select_parents_ga(): public API does not accept hybrid merit arguments", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  merit <- setNames(rnorm(nrow(vmat)), rownames(vmat))
  expect_error(
    select_parents_ga(vmat, n_founders = 4, merit_score = merit),
    "unused argument"
  )
})

test_that("select_parents_ga(): all five strategies run without error", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  for (s in c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS")) {
    expect_no_error(
      select_parents_ga(vmat, n_founders = 4, strategy = s,
                        popSize = 20, maxiter = 15, run = 10, seed = 1)
    )
  }
})

test_that("select_parents_ga(): 'selfing'/'OPV'/'Haploid_OHS' give identical results (documented equivalence)", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  r_self <- select_parents_ga(vmat, n_founders = 4, strategy = "selfing",
                              popSize = 20, maxiter = 15, run = 10, seed = 1)
  r_opv  <- select_parents_ga(vmat, n_founders = 4, strategy = "OPV",
                              popSize = 20, maxiter = 15, run = 10, seed = 1)
  r_hohs <- select_parents_ga(vmat, n_founders = 4, strategy = "Haploid_OHS",
                              popSize = 20, maxiter = 15, run = 10, seed = 1)

  expect_equal(sort(r_self$selected), sort(r_opv$selected))
  expect_equal(sort(r_self$selected), sort(r_hohs$selected))
  expect_equal(r_self$fitness, r_opv$fitness, tolerance = 1e-8)
  expect_equal(r_self$fitness, r_hohs$fitness, tolerance = 1e-8)
})

test_that("select_parents_ga(): 'no_selfing' and 'OHS' give identical results (same distinct-pair rule)", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  r_ns  <- select_parents_ga(vmat, n_founders = 4, strategy = "no_selfing",
                             popSize = 20, maxiter = 15, run = 10, seed = 1)
  r_ohs <- select_parents_ga(vmat, n_founders = 4, strategy = "OHS",
                             popSize = 20, maxiter = 15, run = 10, seed = 1)
  expect_equal(sort(r_ns$selected), sort(r_ohs$selected))
  expect_equal(r_ns$fitness, r_ohs$fitness, tolerance = 1e-8)
})

test_that("select_parents_ga(): value_matrix without row names errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  dimnames(vmat) <- NULL
  expect_error(select_parents_ga(vmat, n_founders = 4), "row names")
})

test_that("select_parents_ga(): n_founders exceeding candidate pool errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 5)
  expect_error(select_parents_ga(vmat, n_founders = 10), "exceeds")
})

test_that("select_parents_ga(): block_weights of wrong length errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(select_parents_ga(vmat, n_founders = 4, block_weights = c(1, 2)),
              "block_weights")
})

test_that("select_parents_ga(): top_candidates smaller than n_founders errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(select_parents_ga(vmat, n_founders = 4, top_candidates = 2),
              "top_candidates")
})

test_that("select_parents_ga(): top_candidates prefilters the candidate pool", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 20)
  res <- select_parents_ga(vmat, n_founders = 4, top_candidates = 8,
                           popSize = 20, maxiter = 15, run = 10, seed = 1)
  row_best <- apply(vmat, 1L, max)
  top8 <- names(sort(row_best, decreasing = TRUE))[1:8]
  expect_true(all(res$selected %in% top8))
})

# ==============================================================================
# select_parents_ga_ts()
# ==============================================================================

# -- min_sel_value / min_sel_mode (merit_score floor) -------------------------

test_that("select_parents_ga_ts(): merit_score is required", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4),
    "merit_score"
  )
})

test_that("select_parents_ga_ts(): min_sel_value restricts selection to candidates clearing the floor", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 16)
  merit <- setNames(seq_len(16), rownames(vmat))  # ind1 lowest merit .. ind16 highest
  res <- select_parents_ga_ts(
    vmat, n_founders = 4, merit_score = merit,
    min_sel_value = 0.5, min_sel_mode = "percentile",
    popSize = 20, maxiter = 15, run = 10, seed = 1
  )
  top8 <- names(sort(merit, decreasing = TRUE))[1:8]  # top 50% by merit
  expect_true(all(res$selected %in% top8))
})

test_that("select_parents_ga_ts(): min_sel_value with no candidates clearing the floor errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  merit <- setNames(rep(1, nrow(vmat)), rownames(vmat))
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit,
                         min_sel_value = 100, min_sel_mode = "value"),
    "min_sel_value"
  )
})

# -- coancestry_weight / G ----------------------------------------------------

test_that("select_parents_ga(): coancestry_weight > 0 without G errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(
    select_parents_ga(vmat, n_founders = 4, coancestry_weight = 1),
    "coancestry_weight"
  )
})

test_that("select_parents_ga(): G with candidates missing from its dimnames errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 6)
  Gmat <- diag(4)
  dimnames(Gmat) <- list(rownames(vmat)[1:4], rownames(vmat)[1:4])  # missing ind5/ind6
  expect_error(
    select_parents_ga(vmat, n_founders = 4, G = Gmat, coancestry_weight = 1),
    "not present in G"
  )
})

test_that("select_parents_ga(): mean_relationship is reported (non-NA) whenever G is supplied", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 10)
  set.seed(6L)
  X <- matrix(rnorm(10 * 5), 10, 5)
  Gmat <- tcrossprod(scale(X, scale = FALSE)) / 5
  diag(Gmat) <- diag(Gmat) + 0.5
  dimnames(Gmat) <- list(rownames(vmat), rownames(vmat))

  res_no_weight <- select_parents_ga(vmat, n_founders = 4, G = Gmat,
                                     coancestry_weight = 0,
                                     popSize = 20, maxiter = 15, run = 10, seed = 1)
  expect_false(is.na(res_no_weight$mean_relationship))

  res_no_G <- select_parents_ga(vmat, n_founders = 4,
                                popSize = 20, maxiter = 15, run = 10, seed = 1)
  expect_true(is.na(res_no_G$mean_relationship))
})

test_that("select_parents_ga(): a strong coancestry_weight overrides a flat coverage landscape to minimise relatedness", {
  skip_if_not_installed("GA")
  # 12 candidates in two tight clusters (ind1-6 highly related to each other,
  # ind7-12 highly related to each other, the two clusters nearly unrelated
  # to one another), with NEAR-IDENTICAL block values across every candidate
  # (so raw merit gives the GA no reason to prefer one cluster over the
  # other). With coancestry_weight = 0 the fitness landscape is flat and the
  # selection could land anywhere; with a large coancestry_weight, the
  # penalty term becomes the ONLY thing that meaningfully differentiates
  # candidate sets, so a working GA search should reliably find a set spread
  # across both clusters (mean pairwise relationship well below the
  # within-cluster level of ~0.95).
  set.seed(8L)
  ids <- paste0("ind", 1:12)
  vmat <- matrix(rnorm(12 * 3, mean = 5, sd = 0.05), nrow = 12, ncol = 3,
                dimnames = list(ids, paste0("blk", 1:3)))

  Gmat <- matrix(0.02, 12, 12, dimnames = list(ids, ids))
  Gmat[1:6, 1:6]   <- 0.95
  Gmat[7:12, 7:12] <- 0.95
  diag(Gmat) <- 1

  res <- select_parents_ga(vmat, n_founders = 4, G = Gmat,
                           coancestry_weight = 50, popSize = 60, maxiter = 60,
                           run = 30, seed = 3, n_reps = 1)
  expect_lt(res$mean_relationship, 0.5)
})

# -- merit_weight (GA+TS hybrid) ----------------------------------------------

test_that("select_parents_ga_ts(): merit_weight > 0 without merit_score errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_weight = 1),
    "merit_score"
  )
})

test_that("select_parents_ga_ts(): merit_score missing candidates errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 6)
  merit <- setNames(1:4, rownames(vmat)[1:4])  # missing ind5/ind6
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit,
                         merit_weight = 1),
    "merit_score"
  )
})

test_that("select_parents_ga_ts(): mean_merit and method are reported", {
  skip_if_not_installed("GA")
  vmat  <- make_vmat(n = 10)
  merit <- setNames(rnorm(10), rownames(vmat))

  res <- select_parents_ga_ts(
    vmat, n_founders = 4, merit_score = merit, merit_priority = 50,
    popSize = 20, maxiter = 15, run = 10, seed = 1
  )
  expect_false(is.na(res$mean_merit))
  expect_gt(res$merit_weight, 0)
  expect_identical(res$selection_method, "select_parents_ga_ts")
  expect_identical(
    res$result_contract$method,
    "select_parents_ga_ts"
  )
})

test_that("select_parents_ga_ts(): a strong merit_weight prefers high merit on a flat coverage landscape", {
  skip_if_not_installed("GA")
  # 12 candidates with NEAR-IDENTICAL block values (flat coverage landscape,
  # so raw block coverage gives the GA no reason to prefer one candidate over
  # another) but a large, two-group merit split: ind1-6 low merit, ind7-12
  # high merit. With merit_weight = 0 the fitness landscape is flat and
  # selection could land anywhere; with a large merit_weight, the merit term
  # becomes the dominant signal, so a working GA search should reliably
  # select from the high-merit group (mirrors the coancestry_weight
  # override test above, same logic applied to the merit term instead).
  set.seed(9L)
  ids  <- paste0("ind", 1:12)
  vmat <- matrix(rnorm(12 * 3, mean = 5, sd = 0.05), nrow = 12, ncol = 3,
                dimnames = list(ids, paste0("blk", 1:3)))
  merit <- setNames(c(rep(0, 6), rep(100, 6)), ids)

  res <- select_parents_ga_ts(
    vmat, n_founders = 4, merit_score = merit,
    merit_weight = 50, popSize = 60, maxiter = 60,
    run = 30, seed = 3, n_reps = 1
  )
  expect_true(all(res$selected %in% ids[7:12]))
  expect_gt(res$mean_merit, 50)
})

# -- merit_priority / suggest_merit_weight() calibration ----------------------
# Hand-verifiable fixture: 6 candidates, 2 blocks, n_founders = 2,
# strategy = "OPV" (single-best-value, no pairing, simplest to hand-compute).
#   A: blk1=10 blk2=0   merit=100
#   B: blk1=0  blk2=10  merit=90
#   C: blk1=5  blk2=5   merit=50
#   D: blk1=1  blk2=1   merit=10
#   E: blk1=2  blk2=2   merit=20
#   F: blk1=3  blk2=3   merit=30
# Greedy coverage build: step 1 ties at sum=10 among A/B/C, first in row
# order (A) wins; step 2, adding B to {A} gives blk1=max(10,0)=10,
# blk2=max(0,10)=10 -> sum=20, beating every other pairing -- so greedy (and
# the TRUE optimum here) is {A,B}, coverage_ceiling = 20.
# Coverage floor (n_cand=6 < 20 -> trim=0): row_best ascending is
# D=1,E=2,F=3,C=5,A=10,B=10 -- bottom 2 = {D,E}, actual coverage
# blk1=max(1,2)=2, blk2=max(1,2)=2 -> coverage_floor = 4. coverage_span = 16.
# Merit ceiling: top 2 by merit = {A,B}, mean = 95. Merit floor (trim=0):
# bottom 2 by merit = {D,E} (10,20), mean = 15. merit_span = 80.
# scale_factor = coverage_span / merit_span = 16 / 80 = 0.2.

.calib_ids   <- c("A", "B", "C", "D", "E", "F")
.calib_vmat  <- matrix(c(10, 0, 0, 10, 5, 5, 1, 1, 2, 2, 3, 3),
                      nrow = 6, ncol = 2, byrow = TRUE,
                      dimnames = list(.calib_ids, c("blk1", "blk2")))
.calib_merit <- setNames(c(100, 90, 50, 10, 20, 30), .calib_ids)

test_that(".greedy_coverage_group(): hand-verified greedy pick matches the true optimum here", {
  chosen <- HapBlockR:::.greedy_coverage_group(.calib_vmat, block_weights = c(1, 1),
                                               strategy = "OPV", n_founders = 2)
  expect_setequal(chosen, c("A", "B"))
})

test_that(".calibrate_merit_scale(): spans and scale_factor match the hand-computed fixture", {
  cal <- HapBlockR:::.calibrate_merit_scale(.calib_vmat, .calib_merit,
                                            block_weights = c(1, 1),
                                            strategy = "OPV", n_founders = 2)
  expect_equal(cal$coverage_ceiling, 20)
  expect_equal(cal$coverage_floor, 4)
  expect_equal(cal$coverage_span, 16)
  expect_equal(cal$merit_ceiling, 95)
  expect_equal(cal$merit_floor, 15)
  expect_equal(cal$merit_span, 80)
  expect_equal(cal$scale_factor, 0.2, tolerance = 1e-8)
  expect_true(cal$ok)
  expect_equal(cal$trim_n, 0L)  # n_cand = 6 < 20
})

test_that(".calibrate_merit_scale(): errors on negative block_weights", {
  expect_error(
    HapBlockR:::.calibrate_merit_scale(.calib_vmat, .calib_merit,
                                       block_weights = c(-1, 1),
                                       strategy = "OPV", n_founders = 2),
    "block_weights"
  )
})

test_that(".calibrate_merit_scale(): ok = FALSE and NA scale_factor when merit has no spread", {
  flat_merit <- setNames(rep(5, 6), .calib_ids)
  cal <- HapBlockR:::.calibrate_merit_scale(.calib_vmat, flat_merit,
                                            block_weights = c(1, 1),
                                            strategy = "OPV", n_founders = 2)
  expect_false(cal$ok)
  expect_true(is.na(cal$scale_factor))
})

test_that("suggest_merit_weight(): returns the hand-computed diagnostic numbers and suggested weight", {
  res <- suggest_merit_weight(.calib_vmat, .calib_merit, n_founders = 2,
                              strategy = "OPV", merit_priority = 100)
  expect_equal(res$coverage_span, 16)
  expect_equal(res$merit_span, 80)
  expect_equal(res$scale_factor, 0.2, tolerance = 1e-8)
  expect_equal(res$suggested_merit_weight, 0.2, tolerance = 1e-8)  # priority 100 -> full scale_factor

  res0 <- suggest_merit_weight(.calib_vmat, .calib_merit, n_founders = 2,
                               strategy = "OPV", merit_priority = 0)
  expect_equal(res0$suggested_merit_weight, 0)

  res_half <- suggest_merit_weight(.calib_vmat, .calib_merit, n_founders = 2,
                                   strategy = "OPV", merit_priority = 50)
  expect_equal(res_half$suggested_merit_weight, 0.1, tolerance = 1e-8)  # half of 0.2
})

test_that("suggest_merit_weight(): merit_priority = NULL returns diagnostics only, no suggestion", {
  res <- suggest_merit_weight(.calib_vmat, .calib_merit, n_founders = 2, strategy = "OPV")
  expect_null(res$merit_priority)
  expect_null(res$suggested_merit_weight)
  expect_equal(res$scale_factor, 0.2, tolerance = 1e-8)  # diagnostics still computed
})

test_that("suggest_merit_weight(): merit_priority outside [0, 100] errors", {
  expect_error(
    suggest_merit_weight(.calib_vmat, .calib_merit, n_founders = 2,
                         strategy = "OPV", merit_priority = 150),
    "merit_priority"
  )
})

test_that("suggest_merit_weight(): missing merit_score entries error clearly", {
  merit_sub <- .calib_merit[1:4]  # missing E, F
  expect_error(
    suggest_merit_weight(.calib_vmat, merit_sub, n_founders = 2, strategy = "OPV"),
    "merit_score"
  )
})

test_that("select_parents_ga_ts(): merit_priority and merit_weight together errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  merit <- setNames(rnorm(nrow(vmat)), rownames(vmat))
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit,
                         merit_weight = 1, merit_priority = 50),
    "exactly one"
  )
})

test_that("select_parents_ga_ts(): merit_priority outside (0, 100] errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  merit <- setNames(rnorm(nrow(vmat)), rownames(vmat))
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit,
                         merit_priority = 200),
    "merit_priority"
  )
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit,
                         merit_priority = 0),
    "coverage-only"
  )
})

test_that("select_parents_ga_ts(): non-finite merit_score errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  merit <- setNames(rnorm(nrow(vmat)), rownames(vmat))
  merit[1] <- NA_real_
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit),
    "finite"
  )
})

test_that("select_parents_ga_ts(): merit_score must cover candidates before calibration", {
  skip_if_not_installed("GA")
  vmat <- make_vmat(n = 6)
  merit <- setNames(1:4, rownames(vmat)[1:4])  # missing ind5/ind6
  expect_error(
    select_parents_ga_ts(vmat, n_founders = 4, merit_score = merit,
                         merit_priority = 50),
    "merit_score"
  )
})

test_that("select_parents_ga_ts(): merit_priority matches suggest_merit_weight() calibration", {
  skip_if_not_installed("GA")
  res <- select_parents_ga_ts(
    .calib_vmat, n_founders = 2, strategy = "OPV",
    merit_score = .calib_merit, merit_priority = 50,
    popSize = 20, maxiter = 15, run = 10, seed = 1
  )
  expect_equal(res$merit_priority, 50)
  # The calibration itself is deterministic (independent of the GA's own
  # stochastic search), so the merit_weight actually used must match
  # suggest_merit_weight()'s independently-computed number exactly.
  expect_equal(res$merit_weight, 0.1, tolerance = 1e-8)
})

test_that("select_parents_ga_ts(): flat merit cannot silently become coverage-only GA", {
  skip_if_not_installed("GA")
  flat_merit <- setNames(rep(1, nrow(.calib_vmat)), rownames(.calib_vmat))
  expect_error(
    select_parents_ga_ts(
      .calib_vmat, n_founders = 2, strategy = "OPV",
      merit_score = flat_merit,
      popSize = 20, maxiter = 15, run = 10, seed = 1
    ),
    "cannot be calibrated"
  )
})

# -- target_degree (easy alternative to coancestry_weight) --------------------
#
# .calib_G: relationship matrix over the same A-F fixture, diag = 1
# throughout, off-diagonals = 0.1 everywhere except G_AB = 0.5 (the greedy
# coverage pick {A,B} is deliberately made the MOST related pair) and
# G_DF = -0.3 (deliberately the LEAST related, i.e. most distant, pair).
# Since D_ij = G_ii + G_jj - 2*G_ij = 2 - 2*G_ij with diag = 1 throughout,
# distance is a strictly decreasing function of G_ij here, so:
#   - gain end (merit_weight = 0): .greedy_coverage_group() picks {A,B}
#     (already hand-verified above) -> gain_end = G_AB = 0.5.
#   - diversity end: select_core_collection(strategy = "maximin") with
#     n_core = 2 picks the single globally-most-distant pair directly (no
#     growth loop needed for n_core = 2) = the pair with the LOWEST G, i.e.
#     {D,F} -> diversity_end = G_DF = -0.3.
#   - span = 0.5 - (-0.3) = 0.8; ceiling(0) = 0.5, ceiling(90) = -0.3,
#     ceiling(45) = 0.5 - 0.5*0.8 = 0.1.
.calib_G <- matrix(0.1, nrow = 6, ncol = 6, dimnames = list(.calib_ids, .calib_ids))
diag(.calib_G) <- 1
.calib_G["A", "B"] <- .calib_G["B", "A"] <- 0.5
.calib_G["D", "F"] <- .calib_G["F", "D"] <- -0.3

test_that(".coverage_span_reference(): refactor still matches the hand-computed fixture", {
  cov <- HapBlockR:::.coverage_span_reference(.calib_vmat, block_weights = c(1, 1),
                                              strategy = "OPV", n_founders = 2, trim = 0L)
  expect_equal(cov$ceiling, 20)
  expect_equal(cov$floor, 4)
  expect_equal(cov$span, 16)
})

test_that(".greedy_coverage_group(): merit_weight = 0 default is unchanged (backward compatible)", {
  chosen <- HapBlockR:::.greedy_coverage_group(.calib_vmat, block_weights = c(1, 1),
                                               strategy = "OPV", n_founders = 2,
                                               merit_weight = 0, merit_score = .calib_merit)
  expect_setequal(chosen, c("A", "B"))
})

test_that(".greedy_coverage_group(): merit_weight > 0 can override the coverage-only pick", {
  # P has higher coverage (10 vs 9) but Q has far higher merit (100 vs 0);
  # with merit_weight = 1, objective(P) = 10 + 0 = 10, objective(Q) =
  # 9 + 100 = 109, so the merit-aware greedy must pick Q over P.
  vmat_pq  <- matrix(c(10, 9), nrow = 2, dimnames = list(c("P", "Q"), "blk1"))
  merit_pq <- setNames(c(0, 100), c("P", "Q"))
  chosen0 <- HapBlockR:::.greedy_coverage_group(vmat_pq, block_weights = 1,
                                                strategy = "OPV", n_founders = 1,
                                                merit_weight = 0, merit_score = merit_pq)
  chosen1 <- HapBlockR:::.greedy_coverage_group(vmat_pq, block_weights = 1,
                                                strategy = "OPV", n_founders = 1,
                                                merit_weight = 1, merit_score = merit_pq)
  expect_equal(chosen0, "P")
  expect_equal(chosen1, "Q")
})

test_that(".calibrate_relatedness_ceiling(): endpoints and interpolated ceiling match the hand-computed fixture", {
  cal0 <- HapBlockR:::.calibrate_relatedness_ceiling(
    .calib_vmat, .calib_G, block_weights = c(1, 1), strategy = "OPV",
    n_founders = 2, target_degree = 0
  )
  expect_setequal(cal0$gain_ids, c("A", "B"))
  expect_setequal(cal0$diversity_ids, c("D", "F"))
  expect_equal(cal0$gain_end, 0.5)
  expect_equal(cal0$diversity_end, -0.3)
  expect_equal(cal0$span, 0.8, tolerance = 1e-8)
  expect_equal(cal0$ceiling, 0.5, tolerance = 1e-8)
  expect_equal(cal0$coverage_span, 16)

  cal90 <- HapBlockR:::.calibrate_relatedness_ceiling(
    .calib_vmat, .calib_G, block_weights = c(1, 1), strategy = "OPV",
    n_founders = 2, target_degree = 90
  )
  expect_equal(cal90$ceiling, -0.3, tolerance = 1e-8)

  cal45 <- HapBlockR:::.calibrate_relatedness_ceiling(
    .calib_vmat, .calib_G, block_weights = c(1, 1), strategy = "OPV",
    n_founders = 2, target_degree = 45
  )
  expect_equal(cal45$ceiling, 0.1, tolerance = 1e-8)
})

test_that("select_parents_ga(): target_degree and explicit coancestry_weight together errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  G <- diag(nrow(vmat)); dimnames(G) <- list(rownames(vmat), rownames(vmat))
  expect_error(
    select_parents_ga(vmat, n_founders = 4, G = G,
                      coancestry_weight = 1, target_degree = 30),
    "coancestry_weight or target_degree"
  )
})

test_that("select_parents_ga(): target_degree outside [0, 90] errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  G <- diag(nrow(vmat)); dimnames(G) <- list(rownames(vmat), rownames(vmat))
  expect_error(
    select_parents_ga(vmat, n_founders = 4, G = G, target_degree = 120),
    "target_degree"
  )
})

test_that("select_parents_ga(): target_degree without G errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(
    select_parents_ga(vmat, n_founders = 4, target_degree = 30),
    "target_degree"
  )
})

test_that("select_parents_ga(): target_degree end-to-end matches the hand-computed calibration", {
  skip_if_not_installed("GA")
  res0 <- select_parents_ga(.calib_vmat, n_founders = 2, strategy = "OPV",
                            G = .calib_G, target_degree = 0,
                            popSize = 20, maxiter = 15, run = 10, seed = 1)
  expect_equal(res0$target_degree, 0)
  expect_equal(res0$relatedness_ceiling, 0.5, tolerance = 1e-8)

  res90 <- select_parents_ga(.calib_vmat, n_founders = 2, strategy = "OPV",
                             G = .calib_G, target_degree = 90,
                             popSize = 20, maxiter = 15, run = 10, seed = 1)
  expect_equal(res90$relatedness_ceiling, -0.3, tolerance = 1e-8)

  res45 <- select_parents_ga(.calib_vmat, n_founders = 2, strategy = "OPV",
                             G = .calib_G, target_degree = 45,
                             popSize = 20, maxiter = 15, run = 10, seed = 1)
  expect_equal(res45$relatedness_ceiling, 0.1, tolerance = 1e-8)
})

test_that("select_parents_ga(): founder count and target_degree ceiling are hard constraints", {
  skip_if_not_installed("GA")
  res <- select_parents_ga(
    .calib_vmat, n_founders = 2, strategy = "OPV",
    G = .calib_G, target_degree = 90,
    popSize = 10, maxiter = 3, run = 2, seed = 11, n_reps = 3
  )

  expect_length(res$selected, 2L)
  expect_true(res$feasible)
  expect_equal(unname(res$constraint_violations[["founder_count"]]), 0)
  expect_lte(res$mean_relationship, res$relatedness_ceiling + 1e-10)
  expect_lte(unname(res$constraint_violations[["relatedness_excess"]]), 1e-10)
  expect_true(all(res$stability$feasible))
})

test_that("select_parents_ga_ts(): repeated runs are ranked by the complete objective", {
  skip_if_not_installed("GA")
  merit <- setNames(c(10, 9, 1, 0, -1, -2), rownames(.calib_vmat))
  res <- select_parents_ga_ts(
    .calib_vmat, n_founders = 2, strategy = "OPV",
    G = .calib_G, coancestry_weight = 0.75,
    merit_score = merit, merit_weight = 1.5,
    popSize = 12, maxiter = 5, run = 3, seed = 21, n_reps = 4
  )

  expect_equal(res$fitness, max(res$stability$fitness_values))
  expect_equal(unname(res$objective_components[["total"]]), res$fitness)
  expect_equal(
    unname(res$objective_components[["total"]]),
    unname(res$objective_components[["coverage"]]) +
      unname(res$objective_components[["merit_bonus"]]) -
      unname(res$objective_components[["coancestry_penalty"]]),
    tolerance = 1e-12
  )
  expect_equal(res$run_id, res$stability$best_rep)
})

# -- n_reps / $stability / $converged -----------------------------------------

test_that("select_parents_ga(): n_reps controls the number of independent replicates in $stability", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  res <- select_parents_ga(vmat, n_founders = 4, popSize = 20, maxiter = 15,
                           run = 10, seed = 1, n_reps = 3)
  expect_equal(res$stability$n_reps, 3L)
  expect_equal(length(res$stability$fitness_values), 3L)
  expect_equal(length(res$stability$converged), 3L)
  expect_equal(length(res$stability$mean_relationship), 3L)
  expect_equal(length(res$stability$mean_merit), 3L)
  expect_true(all(res$stability$selection_freq >= 0 & res$stability$selection_freq <= 1))
  expect_true(is.logical(res$converged))
})

test_that("select_parents_ga(): n_reps < 1 errors", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  expect_error(
    select_parents_ga(vmat, n_founders = 4, n_reps = 0),
    "n_reps"
  )
})

test_that("select_parents_ga(): every candidate in $stability$selection_freq was selected in at least one replicate", {
  skip_if_not_installed("GA")
  vmat <- make_vmat()
  res <- select_parents_ga(vmat, n_founders = 4, popSize = 20, maxiter = 15,
                           run = 10, seed = 1, n_reps = 3)
  expect_true(all(names(res$stability$selection_freq) %in% rownames(vmat)))
  expect_true(res$selected[1] %in% names(res$stability$selection_freq))
})

# ==============================================================================
# plot_parent_selection_pca()
# ==============================================================================

test_that("plot_parent_selection_pca(): runs and returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  n <- 15
  X <- matrix(rnorm(n * 6), n, 6)
  G <- tcrossprod(scale(X, scale = FALSE)) / 6
  diag(G) <- diag(G) + 1e-6
  ids <- paste0("ind", seq_len(n))
  dimnames(G) <- list(ids, ids)

  p <- plot_parent_selection_pca(G, ga_selected = ids[1:3], ts_selected = ids[3:5])
  expect_s3_class(p, "ggplot")
})

test_that("plot_parent_selection_pca(): errors on unnamed G", {
  skip_if_not_installed("ggplot2")
  G <- matrix(rnorm(9), 3, 3)
  expect_error(plot_parent_selection_pca(G, character(0), character(0)), "row/column names")
})
