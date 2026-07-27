## tests/testthat/test-compute-ld.R
## -----------------------------------------------------------------------------
## Unit tests for the R-level exported API in R/compute_ld.R:
##   - compute_r2()      -- thin wrapper around compute_r2_cpp()
##   - compute_rV2()     -- thin wrapper around compute_rV2_cpp()
##   - get_V_inv_sqrt()  -- whitening factor (chol / eigen)
## These wrap the already-tested C++ kernels in test-cpp.R; the point here is
## to exercise the R-level argument handling (matrix coercion, defaults,
## rounding) that test-cpp.R -- which calls the _cpp functions directly --
## never touches.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# ==============================================================================
# compute_r2()
# ==============================================================================

test_that("compute_r2(): matches compute_r2_cpp() directly for a plain matrix", {
  set.seed(55)
  G <- matrix(sample(0:2, 40 * 15, replace = TRUE), 40, 15)
  expect_equal(compute_r2(G), compute_r2_cpp(G, digits = -1L, n_threads = 1L))
})

test_that("compute_r2(): coerces a data.frame input via as.matrix()", {
  set.seed(56)
  G  <- matrix(sample(0:2, 30 * 10, replace = TRUE), 30, 10)
  df <- as.data.frame(G)
  expect_equal(compute_r2(df), compute_r2(G))
})

test_that("compute_r2(): digits argument rounds the output", {
  set.seed(57)
  G        <- matrix(sample(0:2, 30 * 10, replace = TRUE), 30, 10)
  r2_full  <- compute_r2(G, digits = -1L)
  r2_round <- compute_r2(G, digits = 2L)
  expect_equal(r2_round, round(r2_full, 2))
})

test_that("compute_r2(): output is symmetric, diagonal 0, values in [0,1]", {
  set.seed(58)
  G  <- matrix(sample(0:2, 50 * 20, replace = TRUE), 50, 20)
  r2 <- compute_r2(G)
  expect_true(isSymmetric(r2, tol = 1e-10))
  expect_equal(diag(r2), rep(0, 20))
  expect_true(all(r2 >= 0 & r2 <= 1 + 1e-8))
})

test_that("compute_r2(): n_threads argument does not change the result", {
  set.seed(62)
  G <- matrix(sample(0:2, 60 * 25, replace = TRUE), 60, 25)
  expect_equal(compute_r2(G, n_threads = 1L), compute_r2(G, n_threads = 2L))
})

# ==============================================================================
# compute_rV2()
# ==============================================================================

test_that("compute_rV2(): matches compute_rV2_cpp() directly for a plain matrix", {
  set.seed(59)
  G <- matrix(sample(0:2, 40 * 15, replace = TRUE), 40, 15)
  expect_equal(compute_rV2(G), compute_rV2_cpp(G, digits = -1L, n_threads = 1L))
})

test_that("compute_rV2(): coerces a data.frame input via as.matrix()", {
  set.seed(60)
  G  <- matrix(sample(0:2, 30 * 10, replace = TRUE), 30, 10)
  df <- as.data.frame(G)
  expect_equal(compute_rV2(df), compute_rV2(G))
})

test_that("compute_rV2(): digits argument rounds the output", {
  set.seed(63)
  G         <- matrix(sample(0:2, 30 * 10, replace = TRUE), 30, 10)
  rv2_full  <- compute_rV2(G, digits = -1L)
  rv2_round <- compute_rV2(G, digits = 2L)
  expect_equal(rv2_round, round(rv2_full, 2))
})

test_that("compute_rV2(): applied to the same matrix equals compute_r2() (documented equivalence)", {
  # Per ?compute_rV2: "Mathematically identical to compute_r2() applied to X
  # ... The C++ kernel is the same -- the distinction is purely in the
  # preparation step." This confirms both wrappers dispatch to numerically
  # identical kernels; the whitening step itself (prepare_geno()) is tested
  # separately and needs AGHmatrix/ASRgenomics, which this test deliberately
  # avoids depending on.
  set.seed(61)
  G  <- matrix(sample(0:2, 40 * 12, replace = TRUE), 40, 12)
  Xc <- scale(G, center = TRUE, scale = FALSE)
  expect_equal(compute_r2(Xc), compute_rV2(Xc))
})

# ==============================================================================
# get_V_inv_sqrt()
# ==============================================================================

test_that("get_V_inv_sqrt(): chol method gives the exact closed-form result for a diagonal matrix", {
  V <- diag(c(4, 9, 16))
  A <- get_V_inv_sqrt(V, method = "chol")
  expect_equal(A, diag(c(1/2, 1/3, 1/4)), tolerance = 1e-10)
})

test_that("get_V_inv_sqrt(): both methods satisfy A V A^T = I", {
  V <- diag(c(4, 9, 16))
  for (m in c("chol", "eigen")) {
    A <- get_V_inv_sqrt(V, method = m)
    expect_equal(A %*% V %*% t(A), diag(3), tolerance = 1e-8, label = m)
  }
})

test_that("get_V_inv_sqrt(): both methods satisfy A V A^T = I for a NON-diagonal matrix", {
  # A diagonal V has a diagonal (hence trivially symmetric) Cholesky factor,
  # so R^-1 and t(R^-1) coincide and a transpose bug in the "chol" branch is
  # invisible to every test above. A genuine relationship/kinship matrix is
  # never diagonal -- that is the entire point of rV2's kinship whitening --
  # so this is the case that actually exercises the whitening property.
  set.seed(321)
  n <- 5L
  M <- matrix(stats::rnorm(n * n), n, n)
  V <- M %*% t(M) + diag(n)   # symmetric positive-definite, non-diagonal
  for (m in c("chol", "eigen")) {
    A <- get_V_inv_sqrt(V, method = m)
    expect_equal(A %*% V %*% t(A), diag(n), tolerance = 1e-8, label = m)
  }
})

test_that("get_V_inv_sqrt(): chol and eigen agree on a well-conditioned diagonal matrix", {
  V <- diag(c(4, 9, 16))
  A_chol  <- get_V_inv_sqrt(V, method = "chol")
  A_eigen <- get_V_inv_sqrt(V, method = "eigen")
  expect_equal(A_chol, A_eigen, tolerance = 1e-8)
})

test_that("get_V_inv_sqrt(): default method is 'chol'", {
  V <- diag(c(4, 9, 16))
  expect_equal(get_V_inv_sqrt(V), get_V_inv_sqrt(V, method = "chol"))
})

test_that("get_V_inv_sqrt(): eigenvalue floor engages for method='eigen' but not for method='chol' on a near-singular matrix", {
  # chol() has no floor: the tiny eigenvalue produces a very large entry.
  # eigen() clamps eigenvalues to a 1e-6 floor before inverting, so its
  # corresponding entry is capped far lower -- this is the documented
  # "more robust for near-singular matrices" behaviour.
  V <- diag(c(1e-10, 1, 2))
  A_chol  <- get_V_inv_sqrt(V, method = "chol")
  A_eigen <- get_V_inv_sqrt(V, method = "eigen")
  expect_equal(max(diag(A_chol)),  1 / sqrt(1e-10), tolerance = 1e-4)
  expect_equal(max(diag(A_eigen)), 1 / sqrt(1e-6),  tolerance = 1e-4)
  expect_true(max(diag(A_eigen)) < max(diag(A_chol)))
})

test_that("get_V_inv_sqrt(): errors on an unrecognised method", {
  V <- diag(c(4, 9, 16))
  expect_error(get_V_inv_sqrt(V, method = "bogus"))
})
