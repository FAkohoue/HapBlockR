## tests/testthat/test-pareto.R
## -----------------------------------------------------------------------------
## Tests for R/pareto_selection.R:
##   pareto_front()          -- general non-dominated-sort utility, tested with
##                               hand-verified dominance examples (both
##                               directions = "max" and a mixed "max"/"min" case)
##   select_parents_pareto() -- sweeps select_parents_ga() across
##                               coancestry_weight; smoke-tested structurally
##                               (requires GA; skipped otherwise)
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# Local value-matrix fixture (mirrors test-parent-selection.R's make_vmat();
# not shared across test files, so redefined locally here).
.make_vmat_pareto <- function(n = 12, p = 5, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(n * p), nrow = n, ncol = p,
        dimnames = list(paste0("ind", seq_len(n)), paste0("blk", seq_len(p))))
}

# ==============================================================================
# 1. pareto_front(): hand-verified dominance, directions = "max" (default)
# ==============================================================================

test_that("pareto_front: correctly flags a hand-verified dominance example (both max)", {
  # A=(1,1) dominated by B=(2,2); D=(0,0) dominated by A/B/C.
  # B=(2,2) and C=(3,1) are mutually non-dominated (different tradeoff points).
  df <- data.frame(
    id   = c("A", "B", "C", "D"),
    obj1 = c(1, 2, 3, 0),
    obj2 = c(1, 2, 1, 0)
  )
  res <- pareto_front(df, objectives = c("obj1", "obj2"), directions = "max")
  flags <- setNames(res$pareto_optimal, res$id)
  expect_false(flags[["A"]])
  expect_true(flags[["B"]])
  expect_true(flags[["C"]])
  expect_false(flags[["D"]])
})

test_that("pareto_front: with exactly 2 non-dominated points, both get Inf crowding distance", {
  df <- data.frame(id = c("A", "B", "C", "D"),
                   obj1 = c(1, 2, 3, 0), obj2 = c(1, 2, 1, 0))
  res <- pareto_front(df, objectives = c("obj1", "obj2"), directions = "max")
  front <- res[res$pareto_optimal, ]
  expect_equal(nrow(front), 2L)
  expect_true(all(is.infinite(front$crowding_distance)))
})

test_that("pareto_front: correctly flags a hand-verified dominance example (max + min)", {
  # E=(merit=5, related=0.5) is dominated by F=(merit=5, related=0.3):
  # same merit, strictly lower relatedness -> F dominates E.
  # F and G=(merit=3, related=0.1) are mutually non-dominated tradeoff points.
  df <- data.frame(
    id      = c("E", "F", "G"),
    merit   = c(5, 5, 3),
    related = c(0.5, 0.3, 0.1)
  )
  res <- pareto_front(df, objectives = c("merit", "related"),
                      directions = c("max", "min"))
  flags <- setNames(res$pareto_optimal, res$id)
  expect_false(flags[["E"]])
  expect_true(flags[["F"]])
  expect_true(flags[["G"]])
})

test_that("pareto_front: pareto-optimal rows are sorted first", {
  df <- data.frame(id = c("A", "B", "C", "D"),
                   obj1 = c(1, 2, 3, 0), obj2 = c(1, 2, 1, 0))
  res <- pareto_front(df, objectives = c("obj1", "obj2"), directions = "max")
  n_opt <- sum(res$pareto_optimal)
  expect_true(all(res$pareto_optimal[seq_len(n_opt)]))
  if (n_opt < nrow(res))
    expect_true(all(!res$pareto_optimal[(n_opt + 1L):nrow(res)]))
})

test_that("pareto_front: dominated rows have NA crowding_distance", {
  df <- data.frame(id = c("A", "B", "C", "D"),
                   obj1 = c(1, 2, 3, 0), obj2 = c(1, 2, 1, 0))
  res <- pareto_front(df, objectives = c("obj1", "obj2"), directions = "max")
  expect_true(all(is.na(res$crowding_distance[!res$pareto_optimal])))
})

test_that("pareto_front: errors on unknown objective column", {
  df <- data.frame(id = "A", obj1 = 1)
  expect_error(pareto_front(df, objectives = c("obj1", "not_a_column")),
              "not_a_column")
})

test_that("pareto_front: errors on directions of the wrong length", {
  df <- data.frame(obj1 = c(1, 2), obj2 = c(1, 2), obj3 = c(1, 2))
  expect_error(
    pareto_front(df, objectives = c("obj1", "obj2", "obj3"),
                directions = c("max", "min")),
    "directions"
  )
})

test_that("pareto_front: errors on an invalid direction value", {
  df <- data.frame(obj1 = c(1, 2))
  expect_error(pareto_front(df, objectives = "obj1", directions = "up"),
              "'max' or 'min'")
})

test_that("pareto_front: a single dominant point is the only Pareto-optimal one", {
  df <- data.frame(id = c("A", "B", "C"),
                   obj1 = c(1, 2, 5), obj2 = c(1, 2, 5))
  res <- pareto_front(df, objectives = c("obj1", "obj2"), directions = "max")
  expect_equal(res$id[res$pareto_optimal], "C")
})

# ==============================================================================
# 2. select_parents_pareto()
# ==============================================================================

test_that("select_parents_pareto: requires G", {
  skip_if_not_installed("GA")
  vmat <- .make_vmat_pareto()
  expect_error(
    select_parents_pareto(value_matrix = vmat, n_founders = 4L, verbose = FALSE),
    "G is required"
  )
})

test_that("select_parents_pareto: runs and returns the documented structure", {
  skip_if_not_installed("GA")
  vmat <- .make_vmat_pareto(n = 15)
  set.seed(4L)
  X <- matrix(rnorm(15 * 6), 15, 6)
  Gmat <- tcrossprod(scale(X, scale = FALSE)) / 6
  diag(Gmat) <- diag(Gmat) + 0.5
  dimnames(Gmat) <- list(rownames(vmat), rownames(vmat))

  res <- select_parents_pareto(
    value_matrix = vmat, n_founders = 4L, G = Gmat,
    coancestry_weights = c(0, 1), merit = NULL,
    popSize = 20L, maxiter = 15L, run = 10L, n_reps = 1L,
    seed = 1L, verbose = FALSE
  )
  expect_true(all(c("frontier", "runs") %in% names(res)))
  expect_equal(nrow(res$frontier), 2L)
  expect_equal(length(res$runs), 2L)
  req <- c("coancestry_weight", "fitness", "mean_relationship", "mean_merit",
           "n_selected", "converged", "run_index", "pareto_optimal",
           "crowding_distance")
  expect_true(all(req %in% names(res$frontier)))
})

test_that("select_parents_pareto: run_index correctly points back into $runs", {
  skip_if_not_installed("GA")
  vmat <- .make_vmat_pareto(n = 15)
  set.seed(4L)
  X <- matrix(rnorm(15 * 6), 15, 6)
  Gmat <- tcrossprod(scale(X, scale = FALSE)) / 6
  diag(Gmat) <- diag(Gmat) + 0.5
  dimnames(Gmat) <- list(rownames(vmat), rownames(vmat))

  res <- select_parents_pareto(
    value_matrix = vmat, n_founders = 4L, G = Gmat,
    coancestry_weights = c(0, 2), popSize = 20L, maxiter = 15L, run = 10L,
    n_reps = 1L, seed = 1L, verbose = FALSE
  )
  for (i in seq_len(nrow(res$frontier))) {
    run_i <- res$runs[[res$frontier$run_index[i]]]
    expect_equal(run_i$fitness, res$frontier$fitness[i])
  }
})

test_that("select_parents_pareto: errors when coancestry_weights is empty", {
  skip_if_not_installed("GA")
  vmat <- .make_vmat_pareto()
  Gmat <- diag(nrow(vmat)); dimnames(Gmat) <- list(rownames(vmat), rownames(vmat))
  expect_error(
    select_parents_pareto(value_matrix = vmat, n_founders = 4L, G = Gmat,
                          coancestry_weights = numeric(0), verbose = FALSE),
    "coancestry_weights"
  )
})
