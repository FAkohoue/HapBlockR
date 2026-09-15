## tests/testthat/test-exact-validation.R
## -----------------------------------------------------------------------------
## Tests for R/exact_validation.R:
##   validate_crosses_exact() -- exact binary ILP cross selection (lpSolve)
##
## All tests use a small, fully hand-verified candidate-cross table (5 crosses
## among 5 parents) where the true optimum under each constraint combination
## was worked out by exhaustive enumeration by hand (documented inline), so
## exact_objective/exact_plan/gap_pct can be checked against known values
## rather than only structural properties. Guarded by skip_if_not_installed
## ("lpSolve") throughout, since lpSolve is a Suggests dependency.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# -- Shared fixture ---------------------------------------------------------
# Candidate crosses (parent1, parent2, UC):
#   A: P1-P2, UC=10   B: P1-P3, UC=9   C: P2-P3, UC=8
#   D: P1-P4, UC=1    E: P4-P5, UC=7
#
# Hand-verified optima (exhaustive enumeration over all C(5,2)=10 pairs):
#   n_cross=2, no max_cross      -> {A,B}, objective = 19 (simple top-2-by-UC)
#   n_cross=2, max_cross=1       -> {A,E}, objective = 17 (the only 2-cross
#                                    combinations respecting "no parent in
#                                    more than 1 selected cross" are
#                                    {A,E}=17, {B,E}=16, {C,D}=9, {C,E}=15;
#                                    every combination involving two of
#                                    {A,B,C,D} shares a parent (P1 in A/B/D,
#                                    P2 in A/C, P3 in B/C) and so is
#                                    infeasible)
.evc_data <- data.frame(
  parent1 = c("P1", "P1", "P2", "P1", "P4"),
  parent2 = c("P2", "P3", "P3", "P4", "P5"),
  UC      = c(10,   9,    8,    1,    7),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 1. Exact optimum: hand-verified values
# ==============================================================================

test_that("validate_crosses_exact: unconstrained n_cross=2 picks the top 2 by UC", {
  skip_if_not_installed("lpSolve")
  res <- validate_crosses_exact(data = .evc_data, n_cross = 2L, verbose = FALSE)
  expect_equal(res$exact_objective, 19)
  expect_equal(sort(res$exact_plan$UC), c(9, 10))
  expect_equal(res$status, 0)
  expect_equal(res$n_candidates, 5L)
})

test_that("validate_crosses_exact: max_cross=1 finds the hand-verified constrained optimum", {
  skip_if_not_installed("lpSolve")
  res <- validate_crosses_exact(data = .evc_data, n_cross = 2L, max_cross = 1L,
                                verbose = FALSE)
  expect_equal(res$exact_objective, 17)
  expect_equal(sort(res$exact_plan$UC), c(7, 10))  # crosses A (10) and E (7)
  # No parent appears in more than one selected cross:
  parents_used <- c(res$exact_plan$parent1, res$exact_plan$parent2)
  expect_equal(anyDuplicated(parents_used), 0L)
})

test_that("validate_crosses_exact: exact_plan always contains exactly n_cross rows", {
  skip_if_not_installed("lpSolve")
  res <- validate_crosses_exact(data = .evc_data, n_cross = 3L, verbose = FALSE)
  expect_equal(nrow(res$exact_plan), 3L)
})

# ==============================================================================
# 2. Heuristic-plan comparison (gap_pct)
# ==============================================================================

test_that("validate_crosses_exact: gap_pct is negative when heuristic_plan violates max_cross", {
  skip_if_not_installed("lpSolve")
  # Naive "top 2 by UC" heuristic picks {A, B} = 19, which is INFEASIBLE
  # under max_cross=1 (both crosses use P1) -- the exact optimum {A,E}=17
  # is lower in raw UC sum but respects the constraint. gap_pct should be
  # negative here (heuristic_objective > exact_objective) precisely because
  # the "heuristic" plan is comparing an infeasible total against the
  # feasible optimum: 100 * (17 - 19) / 17.
  heuristic <- .evc_data[.evc_data$UC %in% c(10, 9), ]  # A and B
  res <- validate_crosses_exact(data = .evc_data, n_cross = 2L, max_cross = 1L,
                                heuristic_plan = heuristic, verbose = FALSE)
  expect_equal(res$heuristic_objective, 19)
  expect_equal(res$gap_pct, 100 * (17 - 19) / 17, tolerance = 1e-8)
  expect_lt(res$gap_pct, 0)
})

test_that("validate_crosses_exact: gap_pct is exactly zero when heuristic_plan equals the exact optimum", {
  skip_if_not_installed("lpSolve")
  optimal_plan <- .evc_data[.evc_data$UC %in% c(10, 7), ]  # A and E
  res <- validate_crosses_exact(data = .evc_data, n_cross = 2L, max_cross = 1L,
                                heuristic_plan = optimal_plan, verbose = FALSE)
  expect_equal(res$gap_pct, 0, tolerance = 1e-8)
})

test_that("validate_crosses_exact: unmatched heuristic_plan crosses are excluded with a message", {
  skip_if_not_installed("lpSolve")
  bogus <- data.frame(parent1 = "Px", parent2 = "Py", UC = 99)
  expect_message(
    res <- validate_crosses_exact(data = .evc_data, n_cross = 2L,
                                  heuristic_plan = bogus, verbose = TRUE),
    "not found"
  )
  expect_equal(res$heuristic_objective, 0)
})

# ==============================================================================
# 3. culling_pairwise_k
# ==============================================================================

test_that("validate_crosses_exact: culling_pairwise_k excludes crosses above the threshold", {
  skip_if_not_installed("lpSolve")
  dat <- .evc_data
  dat$related <- c(0.5, 0.1, 0.1, 0.1, 0.1)  # cross A (highest UC) is highly related
  res <- validate_crosses_exact(data = dat, n_cross = 2L,
                                culling_pairwise_k = 0.2,
                                relatedness_col = "related", verbose = FALSE)
  expect_equal(res$n_candidates, 4L)
  expect_false(10 %in% res$exact_plan$UC)  # cross A was culled out
})

test_that("validate_crosses_exact: culling_pairwise_k without relatedness info errors", {
  skip_if_not_installed("lpSolve")
  expect_error(
    validate_crosses_exact(data = .evc_data, n_cross = 2L,
                           culling_pairwise_k = 0.2, verbose = FALSE),
    "culling_pairwise_k"
  )
})

# ==============================================================================
# 4. Error paths
# ==============================================================================

test_that("validate_crosses_exact: errors on missing required columns", {
  skip_if_not_installed("lpSolve")
  bad <- .evc_data[, c("parent1", "parent2")]  # no UC column
  expect_error(validate_crosses_exact(data = bad, n_cross = 2L, verbose = FALSE),
              "missing column")
})

test_that("validate_crosses_exact: errors when n_cross exceeds available candidates", {
  skip_if_not_installed("lpSolve")
  expect_error(
    validate_crosses_exact(data = .evc_data, n_cross = 10L, verbose = FALSE),
    "exceeds"
  )
})

test_that("validate_crosses_exact: errors when candidate count exceeds max_vars", {
  skip_if_not_installed("lpSolve")
  expect_error(
    validate_crosses_exact(data = .evc_data, n_cross = 2L, max_vars = 2L,
                           verbose = FALSE),
    "max_vars"
  )
})

test_that("validate_crosses_exact: NA criterion_col rows are dropped before solving", {
  skip_if_not_installed("lpSolve")
  dat <- .evc_data
  dat$UC[1] <- NA
  res <- validate_crosses_exact(data = dat, n_cross = 2L, verbose = FALSE)
  expect_equal(res$n_candidates, 4L)
})

test_that("validate_crosses_exact returns a valid common result contract", {
  skip_if_not_installed("lpSolve")
  result <- validate_crosses_exact(
    data = .evc_data, n_cross = 2L, max_cross = 1L,
    verbose = FALSE
  )

  expect_s3_class(result, "HapBlockR_exact_cross_validation")
  expect_s3_class(result, "hapblockr_result")
  expect_identical(
    result$result_contract$method,
    "validate_crosses_exact"
  )
  expect_equal(
    nrow(result$result_contract$decision_table),
    nrow(result$exact_plan)
  )
  expect_no_error(validate(result))
})
