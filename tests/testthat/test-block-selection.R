## tests/testthat/test-block-selection.R
## -----------------------------------------------------------------------------
## Tests for select_top_blocks() and plot_block_funnel() (HapSelect
## gap-analysis, part 1 -- block-selection ergonomics).
##
## Coverage:
##   1. select_top_blocks(): exactly-one-mode enforcement
##   2. select_top_blocks(): n / perc_total / perc_of_total_var modes,
##      checked against a hand-computed reference (not just "does it run")
##   3. select_top_blocks(): rank and cum_var_share columns
##   4. select_top_blocks(): var_col auto-detection and override, and errors
##   5. plot_block_funnel(): runs and returns a ggplot object (skipped if
##      ggplot2 unavailable)
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# -- Shared fixture: a block-importance table with known, hand-checkable
# variance values (descending, so ranking is unambiguous) -------------------
make_bi <- function() {
  data.frame(
    block_id   = paste0("b", 1:10),
    var_scaled = c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
    stringsAsFactors = FALSE
  )
}

# -- 1. Exactly-one-mode enforcement ------------------------------------------

test_that("select_top_blocks(): errors when zero modes supplied", {
  expect_error(select_top_blocks(make_bi()), "Exactly one")
})

test_that("select_top_blocks(): errors when more than one mode supplied", {
  expect_error(select_top_blocks(make_bi(), n = 3, perc_total = 0.5),
              "Exactly one")
})

# -- 2. n mode -----------------------------------------------------------------

test_that("select_top_blocks(): n mode returns exactly n rows, sorted descending", {
  out <- select_top_blocks(make_bi(), n = 4)
  expect_equal(nrow(out), 4L)
  expect_equal(out$block_id, c("b1", "b2", "b3", "b4"))
  expect_equal(out$rank, 1:4)
})

test_that("select_top_blocks(): n larger than nrow(block_importance) is capped", {
  out <- select_top_blocks(make_bi(), n = 100)
  expect_equal(nrow(out), 10L)
})

# -- 3. perc_total mode ----------------------------------------------------

test_that("select_top_blocks(): perc_total mode keeps ceiling(perc_total * n) blocks", {
  out <- select_top_blocks(make_bi(), perc_total = 0.3)
  expect_equal(nrow(out), 3L)   # ceiling(0.3 * 10) = 3
  expect_equal(out$block_id, c("b1", "b2", "b3"))
})

test_that("select_top_blocks(): perc_total out of (0,1] errors", {
  expect_error(select_top_blocks(make_bi(), perc_total = 0), "perc_total")
  expect_error(select_top_blocks(make_bi(), perc_total = 1.5), "perc_total")
})

# -- 4. perc_of_total_var mode -------------------------------------------------

test_that("select_top_blocks(): perc_of_total_var picks the smallest set reaching the target", {
  # total_var = sum(1.0..0.1) = 5.5
  # cumsum/5.5 = 0.1818 0.3455 0.4909 0.6182 0.7273 0.8182 0.8909 0.9455 0.9818 1.0
  # first index >= 0.9 is index 8 (0.9455)
  out <- select_top_blocks(make_bi(), perc_of_total_var = 0.9)
  expect_equal(nrow(out), 8L)
  expect_true(out$cum_var_share[nrow(out)] >= 0.9)
  expect_true(out$cum_var_share[nrow(out) - 1L] < 0.9)
})

test_that("select_top_blocks(): perc_of_total_var = 1 keeps all blocks", {
  out <- select_top_blocks(make_bi(), perc_of_total_var = 1)
  expect_equal(nrow(out), 10L)
  expect_equal(out$cum_var_share[nrow(out)], 1, tolerance = 1e-8)
})

test_that("select_top_blocks(): perc_of_total_var out of (0,1] errors", {
  expect_error(select_top_blocks(make_bi(), perc_of_total_var = 0), "perc_of_total_var")
  expect_error(select_top_blocks(make_bi(), perc_of_total_var = 1.1), "perc_of_total_var")
})

test_that("select_top_blocks(): cum_var_share is computed over the FULL input, not the subset", {
  # With n = 3 (only 3 rows kept), the cum_var_share on those rows should
  # still reflect the full 10-block total (5.5), not a total recomputed from
  # just the 3 kept rows.
  out <- select_top_blocks(make_bi(), n = 3)
  expect_equal(out$cum_var_share, cumsum(c(1.0, 0.9, 0.8)) / 5.5, tolerance = 1e-8)
})

# -- 5. var_col auto-detection and override ------------------------------------

test_that("select_top_blocks(): prefers var_scaled over var_local_gebv when both present", {
  bi <- make_bi()
  bi$var_local_gebv <- rev(bi$var_scaled)  # deliberately different ordering
  out <- select_top_blocks(bi, n = 3)
  expect_equal(out$block_id, c("b1", "b2", "b3"))  # ranked by var_scaled
})

test_that("select_top_blocks(): falls back to var_local_gebv when var_scaled absent", {
  bi <- make_bi()
  bi$var_local_gebv <- bi$var_scaled
  bi$var_scaled <- NULL
  out <- select_top_blocks(bi, n = 3)
  expect_equal(out$block_id, c("b1", "b2", "b3"))
})

test_that("select_top_blocks(): var_col override selects a different ranking column", {
  bi <- make_bi()
  bi$alt_score <- rev(bi$var_scaled)  # b10 has the highest alt_score
  out <- select_top_blocks(bi, n = 1, var_col = "alt_score")
  expect_equal(out$block_id, "b10")
})

test_that("select_top_blocks(): missing var_col errors", {
  expect_error(select_top_blocks(make_bi(), n = 3, var_col = "nope"), "nope")
})

test_that("select_top_blocks(): missing var_scaled/var_local_gebv errors", {
  bi <- data.frame(block_id = paste0("b", 1:3))
  expect_error(select_top_blocks(bi, n = 1), "var_scaled")
})

test_that("select_top_blocks(): non-data.frame input errors", {
  expect_error(select_top_blocks(as.list(make_bi()), n = 1), "data frame")
})

# -- 6. plot_block_funnel() ----------------------------------------------------

test_that("plot_block_funnel(): runs and returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  bi <- make_bi()
  lg <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10,
              dimnames = list(paste0("ind", 1:20), bi$block_id))
  p <- plot_block_funnel(lg, bi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_block_funnel(): highlight_threshold marks the correct blocks", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  bi <- make_bi()
  lg <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10,
              dimnames = list(paste0("ind", 1:20), bi$block_id))
  p <- plot_block_funnel(lg, bi, highlight_threshold = 0.65)
  # blocks b1-b4 have var_scaled >= 0.65 (1.0, 0.9, 0.8, 0.7); b5..b10 do not
  expect_equal(sum(p$data$important), 4L * 20L)  # 4 blocks x 20 individuals each
})

test_that("plot_block_funnel(): errors on non-matrix local_gebv", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_block_funnel(as.data.frame(matrix(1:4, 2, 2)), make_bi()),
              "matrix")
})

test_that("plot_block_funnel(): errors when block_importance lacks block_id", {
  skip_if_not_installed("ggplot2")
  lg <- matrix(rnorm(4), 2, 2, dimnames = list(c("i1","i2"), c("b1","b2")))
  expect_error(plot_block_funnel(lg, data.frame(x = 1:2)), "block_id")
})

test_that("plot_block_funnel(): errors when no block_id overlap with local_gebv columns", {
  skip_if_not_installed("ggplot2")
  lg <- matrix(rnorm(4), 2, 2, dimnames = list(c("i1","i2"), c("bA","bB")))
  bi <- data.frame(block_id = c("bX", "bY"), var_scaled = c(0.5, 0.6))
  expect_error(plot_block_funnel(lg, bi), "common")
})
