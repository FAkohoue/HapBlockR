## tests/testthat/test-marker-effects.R
## -----------------------------------------------------------------------------
## Tests for this session's marker-effect / local-GEBV additions:
##   - estimate_marker_effects()  (gblup / rrblup / bayesb / bayesc / bayesr)
##   - backsolve_snp_effects()    ploidy generalisation
##   - compute_local_gebv()       complete_decomposition, importance_threshold,
##                                ploidy, and the pre-filtered-geno_matrix guard
##   - compute_haplotype_grm()    ploidy generalisation
##
## Where possible, tests check an exact hand-computed reference rather than
## just "does it run" -- see comments before each block for the identity
## being verified.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# -- Local helper: a simple VanRaden-style GRM for method = "gblup" tests -----
# (Deliberately re-derived here rather than calling compute_haplotype_grm(),
# so this test file does not depend on that function's correctness.)
toy_grm <- function(G, ploidy = 2L) {
  p <- pmax(pmin(colMeans(G) / ploidy, 1 - 1e-8), 1e-8)
  M <- sweep(G, 2, ploidy * p, "-")
  denom <- ploidy * sum(p * (1 - p))
  Gr <- tcrossprod(M) / denom
  diag(Gr) <- diag(Gr) + 1e-6
  dimnames(Gr) <- list(rownames(G), rownames(G))
  Gr
}

# ==============================================================================
# estimate_marker_effects()
# ==============================================================================

test_that("estimate_marker_effects(): method='rrblup' runs and returns full structure", {
  G <- make_geno(n = 40, p = 15, seed = 2L)
  y <- setNames(rnorm(nrow(G)), rownames(G))
  res <- estimate_marker_effects(G, y, method = "rrblup")
  expect_type(res, "list")
  expect_true(all(c("snp_effects", "gebv", "method", "fit") %in% names(res)))
  expect_equal(res$method, "rrblup")
  expect_equal(length(res$snp_effects), ncol(G))
  expect_equal(names(res$snp_effects), colnames(G))
  expect_equal(length(res$gebv), nrow(G))
  expect_equal(names(res$gebv), rownames(G))
})

test_that("estimate_marker_effects(): method='gblup' runs with a supplied G", {
  G <- make_geno(n = 40, p = 15, seed = 3L)
  y <- setNames(rnorm(nrow(G)), rownames(G))
  Gr <- toy_grm(G)
  res <- estimate_marker_effects(G, y, method = "gblup", G = Gr)
  expect_equal(res$method, "gblup")
  expect_equal(length(res$snp_effects), ncol(G))
})

test_that("estimate_marker_effects(): method='gblup' without G errors", {
  G <- make_geno(n = 40, p = 15, seed = 3L)
  y <- setNames(rnorm(nrow(G)), rownames(G))
  expect_error(estimate_marker_effects(G, y, method = "gblup"), "G")
})

test_that("estimate_marker_effects(): prepared targets apply precision weights", {
  G <- make_geno(n = 40, p = 15, seed = 112L)
  ids <- rownames(G)
  set.seed(113)
  target_data <- data.frame(
    id = ids,
    trait = "yield",
    value = rnorm(length(ids)),
    SE = seq(0.15, 0.60, length.out = length(ids))
  )
  targets <- prepare_breeding_targets(
    target_data,
    input_type = "BLUE",
    se_col = "SE"
  )
  result <- estimate_marker_effects(G, targets, method = "rrblup")
  expect_equal(mean(result$precision_weights), 1)
  expect_gt(stats::sd(result$precision_weights), 0)
  expect_equal(names(result$precision_weights), ids)
})

test_that("estimate_marker_effects(): fewer than 10 phenotyped individuals errors", {
  G <- make_geno(n = 40, p = 15, seed = 4L)
  y <- setNames(rep(NA_real_, nrow(G)), rownames(G))
  y[1:5] <- rnorm(5)
  expect_error(estimate_marker_effects(G, y, method = "rrblup"), "10")
})

test_that("estimate_marker_effects(): invalid ploidy errors", {
  G <- make_geno(n = 40, p = 15, seed = 5L)
  y <- setNames(rnorm(nrow(G)), rownames(G))
  expect_error(estimate_marker_effects(G, y, method = "rrblup", ploidy = 1L),
              "ploidy")
  expect_error(estimate_marker_effects(G, y, method = "rrblup", ploidy = c(2L, 3L)),
              "ploidy")
})

test_that("estimate_marker_effects(): method='rrblup' -- ploidy leaves the centred GEBV formula's *value* unchanged (mathematical identity, not a bug)", {
  # CORRECTED from an earlier, wrong version of this test: I originally
  # expected ploidy to change the reported `gebv`. It does not, and that
  # is CORRECT behaviour, not a wiring bug -- caught by actually running
  # this test rather than by inspection.
  #
  # The centring term is ploidy * p_hat, where p_hat = colMeans(geno)/ploidy.
  # Algebraically, ploidy * (colMeans(geno) / ploidy) == colMeans(geno) for
  # ANY ploidy (the ploidy cancels), as long as the pmax/pmin clamp to
  # [1e-8, 1-1e-8] does not engage (it does not, for realistic dosage
  # ranges). So the centred matrix M = geno - colMeans(geno) is ploidy
  # invariant, and since rrBLUP::mixed.solve() fits snp_effects on the RAW
  # (uncentred) genotype matrix -- also ploidy independent -- the reported
  # `gebv = M %*% snp_effects` ends up IDENTICAL across ploidy values for
  # this method. Ploidy's real, discriminating effect is on `denom` in
  # backsolve_snp_effects()/compute_haplotype_grm() (see those tests below),
  # not on this centring step.
  G <- make_geno(n = 40, p = 12, seed = 6L)
  y <- setNames(rnorm(nrow(G)), rownames(G))

  me2 <- estimate_marker_effects(G, y, method = "rrblup", ploidy = 2L)
  me4 <- estimate_marker_effects(G, y, method = "rrblup", ploidy = 4L)

  expect_equal(me2$snp_effects, me4$snp_effects, tolerance = 1e-8)
  expect_equal(unname(me2$gebv), unname(me4$gebv), tolerance = 1e-6)

  expected <- as.numeric(sweep(G, 2, colMeans(G), "-") %*% me2$snp_effects)
  expect_equal(unname(me2$gebv), expected, tolerance = 1e-6)
})

test_that("estimate_marker_effects(): Bayesian methods run when BGLR is available", {
  skip_if_not_installed("BGLR")
  # Note: "bayesa" is used here, not "bayesr" -- BGLR does not implement a
  # "BayesR" model (an earlier version of this test used "bayesr" and it
  # failed at runtime with BGLR's own "model BayesR not implemented" error;
  # see NEWS.md's 0.3.9.9000 entry).
  G <- make_geno(n = 40, p = 12, seed = 7L)
  y <- setNames(rnorm(nrow(G)), rownames(G))
  for (m in c("bayesb", "bayesc", "bayesa")) {
    res <- estimate_marker_effects(G, y, method = m, n_iter = 200L,
                                   burn_in = 50L, seed = 1L, verbose = FALSE)
    expect_equal(res$method, m)
    expect_equal(length(res$snp_effects), ncol(G))
    expect_equal(length(res$gebv), nrow(G))
  }
})

# ==============================================================================
# backsolve_snp_effects() -- ploidy generalisation
# ==============================================================================

test_that("backsolve_snp_effects(): recovered alpha exactly reconstructs the input GEBV (n <= p, ploidy=4)", {
  # Mathematical basis: alpha_hat = M' (M M' / denom)^{-1} gebv / denom * denom
  #   = M' (M M')^{-1} gebv is the minimum-norm solution to M %*% alpha = gebv.
  # For M with full row rank (virtually certain here, n=12 < p=20, random
  # dosage), M %*% alpha_hat reproduces gebv EXACTLY (up to the tiny 1e-6
  # ridge added to G's diagonal for invertibility) -- this holds regardless
  # of whether alpha_hat equals the *true* simulated alpha (which it need
  # not, when p > n). This checks the ploidy-centering is wired correctly
  # end-to-end: if `ploidy` were silently ignored internally, the internal M
  # would differ from the M computed here and this identity would break.
  set.seed(21)
  n <- 12; p <- 20
  G <- matrix(sample(0:4, n * p, replace = TRUE), n, p)
  rownames(G) <- paste0("ind", seq_len(n)); colnames(G) <- paste0("rs", seq_len(p))

  p_hat <- pmax(pmin(colMeans(G) / 4, 1 - 1e-8), 1e-8)
  M     <- sweep(G, 2, 4 * p_hat, "-")

  true_alpha <- setNames(rnorm(p), colnames(G))
  gebv       <- setNames(as.numeric(M %*% true_alpha), rownames(G))

  alpha_hat <- backsolve_snp_effects(G, gebv, ploidy = 4L)
  recon     <- as.numeric(M %*% alpha_hat[colnames(G)])

  expect_equal(recon, unname(gebv), tolerance = 1e-3)
})

test_that("backsolve_snp_effects(): ploidy=2 (default) matches explicit ploidy=2L", {
  G <- make_geno(n = 30, p = 10, seed = 8L)
  gebv <- setNames(rnorm(nrow(G)), rownames(G))
  a_default <- backsolve_snp_effects(G, gebv)
  a_explicit <- backsolve_snp_effects(G, gebv, ploidy = 2L)
  expect_equal(a_default, a_explicit)
})

test_that("backsolve_snp_effects(): invalid ploidy errors", {
  G <- make_geno(n = 30, p = 10, seed = 8L)
  gebv <- setNames(rnorm(nrow(G)), rownames(G))
  expect_error(backsolve_snp_effects(G, gebv, ploidy = 1L), "ploidy")
})

test_that("backsolve_snp_effects(): no matching individuals errors", {
  G <- make_geno(n = 10, p = 5, seed = 9L)
  gebv <- setNames(rnorm(3), c("nobody1", "nobody2", "nobody3"))
  expect_error(backsolve_snp_effects(G, gebv), "matching")
})

# ==============================================================================
# compute_haplotype_grm() -- ploidy generalisation
# ==============================================================================

test_that("compute_haplotype_grm(): ploidy=2 default matches hand-computed VanRaden G", {
  set.seed(9)
  n <- 20; p <- 8
  H <- matrix(sample(0:2, n * p, replace = TRUE), n, p)
  rownames(H) <- paste0("ind", seq_len(n)); colnames(H) <- paste0("h", seq_len(p))

  G1 <- compute_haplotype_grm(H, phased = TRUE)
  G2 <- compute_haplotype_grm(H, phased = TRUE, ploidy = 2L)
  expect_equal(G1, G2)

  p_hat <- pmax(pmin(colMeans(H) / 2, 1 - 1e-6), 1e-6)
  Z     <- sweep(H, 2, 2 * p_hat, "-")
  denom <- 2 * sum(p_hat * (1 - p_hat))
  expected <- tcrossprod(Z) / denom
  dimnames(expected) <- dimnames(G1)
  expect_equal(G1, expected, tolerance = 1e-8)
})

test_that("compute_haplotype_grm(): ploidy=4 rescales centering/denominator correctly", {
  set.seed(9)
  n <- 20; p <- 8
  H <- matrix(sample(0:4, n * p, replace = TRUE), n, p)
  rownames(H) <- paste0("ind", seq_len(n)); colnames(H) <- paste0("h", seq_len(p))

  G4 <- compute_haplotype_grm(H, phased = TRUE, ploidy = 4L)

  p_hat <- pmax(pmin(colMeans(H) / 4, 1 - 1e-6), 1e-6)
  Z     <- sweep(H, 2, 4 * p_hat, "-")
  denom <- 4 * sum(p_hat * (1 - p_hat))
  expected <- tcrossprod(Z) / denom
  dimnames(expected) <- dimnames(G4)
  expect_equal(G4, expected, tolerance = 1e-8)
})

test_that("compute_haplotype_grm(): invalid ploidy errors", {
  H <- make_geno(n = 20, p = 8, seed = 10L)
  expect_error(compute_haplotype_grm(H, ploidy = 1L), "ploidy")
})

# ==============================================================================
# compute_local_gebv() -- complete_decomposition, importance_threshold, ploidy
# ==============================================================================

test_that("compute_local_gebv(): ploidy=4 local GEBVs sum to the hand-computed genome-wide GEBV", {
  set.seed(5)
  n <- 24; p <- 12
  G <- matrix(sample(0:4, n * p, replace = TRUE), n, p)
  rownames(G) <- paste0("ind", seq_len(n)); colnames(G) <- paste0("rs", seq_len(p))
  snp_info <- data.frame(SNP = colnames(G), CHR = "1",
                         POS = seq(1000L, by = 2000L, length.out = p),
                         stringsAsFactors = FALSE)
  blocks <- make_blocks(snp_info, n_blocks = 3L)
  alpha  <- setNames(rnorm(p), colnames(G))

  res <- compute_local_gebv(G, snp_info, blocks, alpha, ploidy = 4L,
                            complete_decomposition = TRUE)

  p_hat <- colMeans(G) / 4
  M     <- sweep(G, 2, 4 * p_hat, "-")
  expected_gebv <- setNames(as.numeric(M %*% alpha), rownames(G))

  recon <- rowSums(res$local_gebv)
  expect_equal(unname(recon[rownames(G)]), unname(expected_gebv), tolerance = 1e-8)
})

test_that("compute_local_gebv(): propagates marker-effect standard errors", {
  G <- matrix(c(
    0, 0,
    1, 2,
    2, 1
  ), nrow = 3L, byrow = TRUE,
  dimnames = list(paste0("i", 1:3), c("s1", "s2")))
  snp_info <- data.frame(
    SNP = c("s1", "s2"), CHR = "1", POS = c(100L, 200L)
  )
  blocks <- data.frame(
    CHR = "1", start.bp = 50L, end.bp = 250L
  )
  alpha <- c(s1 = 0.4, s2 = -0.2)
  alpha_se <- c(s1 = 0.10, s2 = 0.20)
  result <- compute_local_gebv(
    G, snp_info, blocks, alpha, snp_effect_se = alpha_se
  )
  centred <- sweep(G, 2L, colMeans(G), "-")
  expected_se <- sqrt(rowSums(
    sweep(centred^2, 2L, alpha_se^2, "*")
  ))
  expect_equal(as.numeric(result$local_gebv_se[, 1L]), unname(expected_se))
  expect_identical(
    result$uncertainty_assumption,
    "independent_marker_effect_errors"
  )
  expect_true(all(is.finite(result$block_importance$mean_local_gebv_se)))
})

test_that("compute_local_gebv(): complete_decomposition=TRUE creates singleton blocks for SNPs outside every block window", {
  n <- 20; p <- 12
  G <- make_geno(n = n, p = p, seed = 11L)
  snp_info <- make_snpinfo(p = p)
  # Blocks only cover the first 3 SNPs (by position) -- SNPs 4..12 are
  # deliberately left uncovered by any block window.
  blocks <- data.frame(
    start = 1L, end = 3L,
    start.rsID = snp_info$SNP[1L], end.rsID = snp_info$SNP[3L],
    start.bp = snp_info$POS[1L], end.bp = snp_info$POS[3L],
    CHR = "1", length_bp = snp_info$POS[3L] - snp_info$POS[1L] + 1L,
    stringsAsFactors = FALSE
  )
  alpha <- setNames(rnorm(p), colnames(G))

  res_complete <- compute_local_gebv(G, snp_info, blocks, alpha,
                                     complete_decomposition = TRUE)
  res_partial  <- compute_local_gebv(G, snp_info, blocks, alpha,
                                     complete_decomposition = FALSE)

  # 1 real block + 9 singleton blocks (SNPs 4..12) = 10 columns
  expect_equal(ncol(res_complete$local_gebv), 10L)
  expect_true(sum(res_complete$block_importance$singleton) == 9L)

  # Without complete_decomposition, only the 1 real block survives
  expect_equal(ncol(res_partial$local_gebv), 1L)

  # rowSums under complete_decomposition reconstructs the full GEBV
  p_hat <- pmax(pmin(colMeans(G) / 2, 1 - 1e-8), 1e-8)
  M     <- sweep(G, 2, 2 * p_hat, "-")
  expected_gebv <- as.numeric(M %*% alpha)
  expect_equal(unname(rowSums(res_complete$local_gebv)), expected_gebv,
              tolerance = 1e-8)
})

test_that("compute_local_gebv(): importance_threshold controls the 'important' flag", {
  n <- 20; p <- 12
  G <- make_geno(n = n, p = p, seed = 12L)
  snp_info <- make_snpinfo(p = p)
  blocks <- make_blocks(snp_info, n_blocks = 3L)
  alpha <- setNames(rnorm(p), colnames(G))

  res_90 <- compute_local_gebv(G, snp_info, blocks, alpha, importance_threshold = 0.9)
  res_50 <- compute_local_gebv(G, snp_info, blocks, alpha, importance_threshold = 0.5)

  expect_equal(res_90$block_importance$important,
              res_90$block_importance$var_scaled >= 0.9)
  expect_equal(res_50$block_importance$important,
              res_50$block_importance$var_scaled >= 0.5)
  # A lower threshold can never flag FEWER blocks as important.
  expect_true(sum(res_50$block_importance$important) >=
             sum(res_90$block_importance$important))
})

test_that("compute_local_gebv(): warns when geno_matrix looks pre-filtered to block-member SNPs only", {
  n <- 20; p <- 12
  G <- make_geno(n = n, p = p, seed = 13L)
  snp_info <- make_snpinfo(p = p)
  # Block covers only the first 3 SNPs.
  blocks <- data.frame(
    start = 1L, end = 3L,
    start.rsID = snp_info$SNP[1L], end.rsID = snp_info$SNP[3L],
    start.bp = snp_info$POS[1L], end.bp = snp_info$POS[3L],
    CHR = "1", length_bp = snp_info$POS[3L] - snp_info$POS[1L] + 1L,
    stringsAsFactors = FALSE
  )
  alpha <- setNames(rnorm(p), colnames(G))

  # geno_matrix pre-filtered to ONLY the 3 block-member SNPs -- the other 9
  # SNPs in snp_info have no corresponding column here at all.
  G_filtered <- G[, snp_info$SNP[1:3], drop = FALSE]

  expect_warning(
    compute_local_gebv(G_filtered, snp_info, blocks, alpha,
                      complete_decomposition = TRUE),
    "pre-filtered|already been filtered"
  )
})

test_that("compute_local_gebv(): invalid ploidy errors", {
  n <- 20; p <- 12
  G <- make_geno(n = n, p = p, seed = 14L)
  snp_info <- make_snpinfo(p = p)
  blocks <- make_blocks(snp_info, n_blocks = 3L)
  alpha <- setNames(rnorm(p), colnames(G))
  expect_error(compute_local_gebv(G, snp_info, blocks, alpha, ploidy = 1L),
              "ploidy")
})

# ==============================================================================
# prepare_gblup_inputs()
# ==============================================================================
# hap_matrix just needs to be a numeric matrix with individual-ID rownames --
# make_geno()'s 0/1/2 dosage matrix is a convenient, already-trusted stand-in
# for a real haplotype feature matrix (e.g. from build_haplotype_feature_matrix()).

test_that("prepare_gblup_inputs(): full overlap, single trait -- dimensions, names, n_train/n_predict", {
  hap <- make_geno(n = 20, p = 6, seed = 101L)
  colnames(hap) <- paste0("hap", 1:6)
  pheno <- data.frame(id = rownames(hap), YLD = rnorm(20), stringsAsFactors = FALSE)

  # NOTE: expect_warning()/expect_message() return the captured CONDITION
  # object, not the wrapped expression's value -- the assignment must
  # happen INSIDE the wrapped expression, not around the whole
  # expect_message() call, or `out` ends up being the message itself.
  expect_message(
    out <- prepare_gblup_inputs(hap, pheno, id_col = "id", trait_col = "YLD"),
    "training individuals"
  )
  expect_type(out, "list")
  expect_true(all(c("G", "pheno_df", "y_vec", "n_train", "n_predict") %in% names(out)))
  expect_equal(dim(out$G), c(20L, 20L))
  expect_equal(rownames(out$G), rownames(hap))
  expect_equal(nrow(out$pheno_df), 20L)
  expect_equal(length(out$y_vec), 20L)
  expect_equal(names(out$y_vec), rownames(hap))
  expect_equal(out$n_train, 20L)
  expect_equal(out$n_predict, 0L)
})

test_that("prepare_gblup_inputs(): trait_col = NULL leaves y_vec NULL and n_train/n_predict NA", {
  hap <- make_geno(n = 15, p = 5, seed = 102L)
  pheno <- data.frame(id = rownames(hap), stringsAsFactors = FALSE)
  out <- prepare_gblup_inputs(hap, pheno, id_col = "id")
  expect_null(out$y_vec)
  expect_true(is.na(out$n_train))
  expect_true(is.na(out$n_predict))
  expect_equal(nrow(out$G), 15L)
})

test_that("prepare_gblup_inputs(): NA values in trait_col are kept and counted as prediction candidates", {
  hap <- make_geno(n = 12, p = 5, seed = 108L)
  y <- rnorm(12); y[c(2, 5, 9)] <- NA
  pheno <- data.frame(id = rownames(hap), YLD = y, stringsAsFactors = FALSE)
  out <- suppressMessages(
    prepare_gblup_inputs(hap, pheno, id_col = "id", trait_col = "YLD")
  )
  expect_equal(out$n_train, 9L)
  expect_equal(out$n_predict, 3L)
  expect_equal(sum(is.na(out$y_vec)), 3L)
  expect_equal(length(out$y_vec), 12L)
})

test_that("prepare_gblup_inputs(): errors when id_col is not in pheno_df", {
  hap <- make_geno(n = 10, p = 4, seed = 103L)
  pheno <- data.frame(sample_id = rownames(hap), stringsAsFactors = FALSE)
  expect_error(prepare_gblup_inputs(hap, pheno, id_col = "id"), "id_col")
})

test_that("prepare_gblup_inputs(): errors when trait_col is not in pheno_df", {
  hap <- make_geno(n = 10, p = 4, seed = 104L)
  pheno <- data.frame(id = rownames(hap), stringsAsFactors = FALSE)
  expect_error(
    prepare_gblup_inputs(hap, pheno, id_col = "id", trait_col = "nope"),
    "trait_col"
  )
})

test_that("prepare_gblup_inputs(): errors when there is no ID overlap at all", {
  hap <- make_geno(n = 10, p = 4, seed = 105L)
  pheno <- data.frame(id = paste0("other", 1:10), YLD = rnorm(10),
                      stringsAsFactors = FALSE)
  expect_error(
    prepare_gblup_inputs(hap, pheno, id_col = "id", trait_col = "YLD"),
    "No matching individual IDs"
  )
})

test_that("prepare_gblup_inputs(): one-sided partial overlap excludes unmatched genotyped individuals and reports it", {
  hap <- make_geno(n = 20, p = 6, seed = 106L)   # ind1..ind20
  pheno <- data.frame(id = rownames(hap)[11:20], YLD = rnorm(10),
                      stringsAsFactors = FALSE)  # a subset of hap's IDs

  msgs <- testthat::capture_messages(
    out <- prepare_gblup_inputs(hap, pheno, id_col = "id", trait_col = "YLD")
  )
  expect_true(any(grepl("excluded from G", msgs)))
  expect_true(any(grepl("training individuals", msgs)))
  expect_false(any(grepl("excluded from output", msgs)))  # pheno_df has no unmatched IDs here

  expect_equal(nrow(out$G), 10L)
  expect_equal(sort(rownames(out$G)), sort(rownames(hap)[11:20]))
  expect_equal(nrow(out$pheno_df), 10L)
  expect_equal(out$n_train, 10L)
  expect_equal(out$n_predict, 0L)
})

test_that("prepare_gblup_inputs(): two-sided partial overlap excludes from both sides and reports both", {
  hap <- make_geno(n = 20, p = 6, seed = 109L)                    # ind1..ind20
  pheno <- data.frame(id = paste0("ind", 15:34),                  # ind15..ind34
                      YLD = rnorm(20), stringsAsFactors = FALSE)
  # Overlap is ind15..ind20 (6 individuals); hap loses ind1..14, pheno loses ind21..34.
  msgs <- testthat::capture_messages(
    out <- prepare_gblup_inputs(hap, pheno, id_col = "id", trait_col = "YLD")
  )
  expect_true(any(grepl("excluded from G", msgs)))
  expect_true(any(grepl("excluded from output", msgs)))
  expect_equal(nrow(out$G), 6L)
  expect_equal(sort(rownames(out$G)), paste0("ind", 15:20))
  expect_equal(nrow(out$pheno_df), 6L)
})

test_that("prepare_gblup_inputs(): bend adds exactly 0.001 to G's diagonal and leaves off-diagonal entries unchanged", {
  hap <- make_geno(n = 12, p = 5, seed = 107L)
  pheno <- data.frame(id = rownames(hap), stringsAsFactors = FALSE)
  out_bend   <- prepare_gblup_inputs(hap, pheno, id_col = "id", bend = TRUE)
  out_nobend <- prepare_gblup_inputs(hap, pheno, id_col = "id", bend = FALSE)
  expect_equal(unname(diag(out_bend$G) - diag(out_nobend$G)), rep(0.001, 12),
              tolerance = 1e-10)
  expect_equal(out_bend$G[upper.tri(out_bend$G)],
              out_nobend$G[upper.tri(out_nobend$G)], tolerance = 1e-10)
})
