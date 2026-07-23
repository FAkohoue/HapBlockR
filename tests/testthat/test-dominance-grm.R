## tests/testthat/test-dominance-grm.R
## -----------------------------------------------------------------------------
## Tests for R/haplotypes.R: compute_dominance_grm()
## (Vitezica et al. 2013 dominance relationship matrix).
##
## Hand-verified fixture: 4 individuals x 2 SNPs, both SNPs at p=q=0.5
## (dosages chosen so mean = 1 exactly), so the dominance coefficients are
## the same simple numbers for every column:
##   homozygous (dose 0 or 2): S = -2*q^2 = -2*p^2 = -0.5   (equal since p=q)
##   heterozygous (dose 1):    S =  2*p*q =  0.5
##
##   SNP1 dosage: ind1=0 ind2=1 ind3=1 ind4=2
##   SNP2 dosage: ind1=0 ind2=0 ind3=2 ind4=2
##
## S matrix:
##   ind1 = (-0.5, -0.5)   ind2 = (0.5, -0.5)
##   ind3 = (0.5, -0.5)    ind4 = (-0.5, -0.5)
##
## denom = sum((2pq)^2) = 0.5^2 + 0.5^2 = 0.5
## D = tcrossprod(S) / denom gives the fully hand-verified matrix below
## (worked out by hand, not by running the package's own code).
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

.dgrm_geno <- matrix(
  c(0, 1, 1, 2,   # SNP1
    0, 0, 2, 2),  # SNP2
  nrow = 4, ncol = 2,
  dimnames = list(paste0("ind", 1:4), c("rs1", "rs2"))
)

.dgrm_expected <- matrix(
  c(1, 0, 0, 1,
    0, 1, 1, 0,
    0, 1, 1, 0,
    1, 0, 0, 1),
  nrow = 4, byrow = TRUE,
  dimnames = list(paste0("ind", 1:4), paste0("ind", 1:4))
)

# ==============================================================================
# 1. Hand-verified values
# ==============================================================================

test_that("compute_dominance_grm: matches the hand-derived matrix exactly (bend = FALSE)", {
  D <- compute_dominance_grm(.dgrm_geno, bend = FALSE)
  expect_equal(unname(D), unname(.dgrm_expected), tolerance = 1e-8)
})

test_that("compute_dominance_grm: bend = TRUE adds exactly 0.001 to the diagonal only", {
  D_raw  <- compute_dominance_grm(.dgrm_geno, bend = FALSE)
  D_bent <- compute_dominance_grm(.dgrm_geno, bend = TRUE)
  expect_equal(diag(D_bent), diag(D_raw) + 0.001, tolerance = 1e-10)
  offdiag_idx <- upper.tri(D_raw)
  expect_equal(D_bent[offdiag_idx], D_raw[offdiag_idx], tolerance = 1e-10)
})

test_that("compute_dominance_grm: individuals with the same genotype-class pattern have D=1", {
  D <- compute_dominance_grm(.dgrm_geno, bend = FALSE)
  expect_equal(D["ind1", "ind4"], 1, tolerance = 1e-8)  # both homozygous at both SNPs
  expect_equal(D["ind2", "ind3"], 1, tolerance = 1e-8)  # both het at SNP1, homo-ref at SNP2
})

# ==============================================================================
# 2. Structure
# ==============================================================================

test_that("compute_dominance_grm: result is symmetric and dimnamed by individual ID", {
  D <- compute_dominance_grm(.dgrm_geno)
  expect_equal(D, t(D), tolerance = 1e-10)
  expect_equal(rownames(D), rownames(.dgrm_geno))
  expect_equal(colnames(D), rownames(.dgrm_geno))
})

test_that("compute_dominance_grm: works on a larger random 0/1/2 matrix without error", {
  set.seed(21L)
  G <- matrix(sample(0:2, 40 * 15, replace = TRUE), 40, 15,
             dimnames = list(paste0("ind", 1:40), paste0("rs", 1:15)))
  D <- compute_dominance_grm(G)
  expect_equal(dim(D), c(40L, 40L))
  expect_equal(D, t(D), tolerance = 1e-8)
})

test_that("compute_dominance_grm: NA values are mean-imputed rather than erroring", {
  G <- .dgrm_geno
  G[1, 1] <- NA
  expect_no_error(D <- compute_dominance_grm(G))
  expect_false(anyNA(D))
})

test_that("compute_dominance_grm: monomorphic columns are silently dropped", {
  G <- cbind(.dgrm_geno, rs3 = c(1, 1, 1, 1))  # monomorphic
  D <- compute_dominance_grm(G, bend = FALSE)
  expect_equal(unname(D), unname(.dgrm_expected), tolerance = 1e-8)
})

# ==============================================================================
# 3. Errors
# ==============================================================================

test_that("compute_dominance_grm: errors on ploidy != 2", {
  expect_error(compute_dominance_grm(.dgrm_geno, ploidy = 4L), "ploidy = 2")
})

test_that("compute_dominance_grm: errors when geno_matrix has no row names", {
  G <- .dgrm_geno; rownames(G) <- NULL
  expect_error(compute_dominance_grm(G), "row names")
})

test_that("compute_dominance_grm: errors on dosage values outside [0, 2]", {
  G <- .dgrm_geno; G[1, 1] <- 3
  expect_error(compute_dominance_grm(G), "0, 2")
})

test_that("compute_dominance_grm: errors when every column is monomorphic", {
  G <- matrix(1, 4, 2, dimnames = list(paste0("ind", 1:4), c("rs1", "rs2")))
  expect_error(compute_dominance_grm(G), "polymorphic")
})
