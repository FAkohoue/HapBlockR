## tests/testthat/test-core-collection.R
## -----------------------------------------------------------------------------
## Tests for R/core_collection.R:
##   select_core_collection() -- maximin/mean_distance farthest-point greedy
##                                selection, hand-verified against a small
##                                1-D "coordinate" fixture where every pairwise
##                                distance (and hence every greedy step) is
##                                known exactly by construction.
##
## Fixture: 5 points on a number line, coordinates A=0, B=1, C=2, D=8, E=9.
## Pairwise |distance|:
##   A-B=1  A-C=2  A-D=8  A-E=9
##   B-C=1  B-D=7  B-E=8
##   C-D=6  C-E=7
##   D-E=1
## The unique farthest pair is (A,E) [distance 9] -- select_core_collection()
## always starts there. Given sel={A,E}, each remaining point's distance to
## the selected set (min of its distances to A and E) is: B=min(1,8)=1,
## C=min(2,7)=2, D=min(8,1)=1 -- so C (distance 2) is the unique 3rd pick
## under the "maximin" strategy.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

.cc_ids <- c("A", "B", "C", "D", "E")
.cc_x   <- c(A = 0, B = 1, C = 2, D = 8, E = 9)

# Distance matrix directly (type = "distance")
.cc_D <- as.matrix(dist(.cc_x, method = "euclidean"))
dimnames(.cc_D) <- list(.cc_ids, .cc_ids)

# Relationship-style Gram matrix (type = "relationship"): G_ij = x_i * x_j,
# so D_ij = G_ii + G_jj - 2*G_ij = x_i^2 + x_j^2 - 2*x_i*x_j = (x_i - x_j)^2
# (squared distance -- same relative ordering as .cc_D, so the same greedy
# picks result, with squared distance VALUES).
.cc_G <- outer(.cc_x, .cc_x)
dimnames(.cc_G) <- list(.cc_ids, .cc_ids)

# ==============================================================================
# 1. maximin strategy, type = "distance" (hand-verified)
# ==============================================================================

test_that("select_core_collection (maximin, distance): n_core=2 picks the globally farthest pair", {
  res <- select_core_collection(G = .cc_D, n_core = 2L, type = "distance",
                                strategy = "maximin", verbose = FALSE)
  expect_equal(sort(res$selected), c("A", "E"))
  expect_equal(res$mean_distance, 9)
  expect_equal(res$min_distance, 9)
})

test_that("select_core_collection (maximin, distance): n_core=3 adds the hand-verified 3rd pick (C)", {
  res <- select_core_collection(G = .cc_D, n_core = 3L, type = "distance",
                                strategy = "maximin", verbose = FALSE)
  expect_equal(sort(res$selected), c("A", "C", "E"))
  expect_equal(res$mean_distance, mean(c(2, 9, 7)))  # A-C=2, A-E=9, C-E=7
  expect_equal(res$min_distance, 2)
})

test_that("select_core_collection (maximin, distance): trace has one row per selection step, in order", {
  res <- select_core_collection(G = .cc_D, n_core = 3L, type = "distance",
                                strategy = "maximin", verbose = FALSE)
  expect_equal(nrow(res$trace), 3L)
  expect_equal(res$trace$step, 1:3)
  expect_true(is.na(res$trace$criterion[1]))
  expect_equal(res$trace$criterion[3], 2)  # the 3rd pick's maximin distance
})

test_that("select_core_collection (maximin, distance): n_core=5 selects every individual", {
  res <- select_core_collection(G = .cc_D, n_core = 5L, type = "distance",
                                strategy = "maximin", verbose = FALSE)
  expect_equal(sort(res$selected), sort(.cc_ids))
})

# ==============================================================================
# 2. type = "relationship" (D_ij = G_ii + G_jj - 2*G_ij identity)
# ==============================================================================

test_that("select_core_collection (relationship): applies the exact distance identity", {
  res <- select_core_collection(G = .cc_G, n_core = 3L, type = "relationship",
                                strategy = "maximin", verbose = FALSE)
  expect_equal(sort(res$selected), c("A", "C", "E"))  # same picks (squaring preserves order here)
  expect_equal(res$mean_distance, mean(c(4, 81, 49)), tolerance = 1e-8)  # squared distances
  expect_equal(res$min_distance, 4)
})

# ==============================================================================
# 3. mean_distance strategy (smoke-tested; full-selection edge case verified exactly)
# ==============================================================================

test_that("select_core_collection (mean_distance): n_core=5 selects every individual", {
  res <- select_core_collection(G = .cc_D, n_core = 5L, type = "distance",
                                strategy = "mean_distance", verbose = FALSE)
  expect_equal(sort(res$selected), sort(.cc_ids))
})

test_that("select_core_collection (mean_distance): starts from the same farthest pair as maximin", {
  res_md <- select_core_collection(G = .cc_D, n_core = 2L, type = "distance",
                                   strategy = "mean_distance", verbose = FALSE)
  expect_equal(sort(res_md$selected), c("A", "E"))
})

# ==============================================================================
# 4. n_core = 1 edge case
# ==============================================================================

test_that("select_core_collection: n_core=1 with merit picks the single highest-merit individual", {
  merit <- setNames(c(1, 2, 10, 3, 4), .cc_ids)
  res <- select_core_collection(G = .cc_D, n_core = 1L, type = "distance",
                                merit = merit, verbose = FALSE)
  expect_equal(res$selected, "C")
  expect_true(is.na(res$mean_distance))
  expect_true(is.na(res$min_distance))
})

# ==============================================================================
# 5. Merit floor
# ==============================================================================

test_that("select_core_collection: min_sel_value/min_sel_mode restricts the eligible pool", {
  merit <- setNames(c(1, 2, 3, 4, 5), .cc_ids)  # A lowest, E highest
  expect_message(
    res <- select_core_collection(G = .cc_D, n_core = 3L, type = "distance",
                                  merit = merit, min_sel_value = 2.5,
                                  min_sel_mode = "value", verbose = TRUE),
    "eligible"
  )
  expect_true(all(res$selected %in% c("C", "D", "E")))  # merit > 2.5
})

# ==============================================================================
# 6. Errors
# ==============================================================================

test_that("select_core_collection: errors on undimnamed G", {
  expect_error(select_core_collection(G = unname(.cc_D), n_core = 2L, verbose = FALSE),
              "dimnamed matrix")
})

test_that("select_core_collection: errors when G rownames != colnames", {
  G_bad <- .cc_D
  colnames(G_bad) <- rev(colnames(G_bad))
  expect_error(select_core_collection(G = G_bad, n_core = 2L, verbose = FALSE),
              "row and column names")
})

test_that("select_core_collection: errors when n_core < 1", {
  expect_error(select_core_collection(G = .cc_D, n_core = 0L, verbose = FALSE),
              "n_core")
})

test_that("select_core_collection: errors when n_core exceeds the number of eligible individuals", {
  expect_error(select_core_collection(G = .cc_D, n_core = 10L, verbose = FALSE),
              "exceeds")
})

test_that("select_core_collection: errors when all pairwise distances are zero", {
  ids3 <- c("X", "Y", "Z")
  G_flat <- matrix(1, 3, 3, dimnames = list(ids3, ids3))  # identical -> all distances 0
  expect_error(select_core_collection(G = G_flat, n_core = 2L, verbose = FALSE),
              "identical")
})

test_that("select_core_collection: errors with fewer than 2 individuals after a merit floor", {
  merit <- setNames(c(1, 1, 1, 1, 100), .cc_ids)
  expect_error(
    select_core_collection(G = .cc_D, n_core = 2L, merit = merit,
                           min_sel_value = 50, min_sel_mode = "value", verbose = FALSE),
    "Fewer than 2"
  )
})
