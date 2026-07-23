## tests/testthat/test-genomic-mating.R
## -----------------------------------------------------------------------------
## Tests for R/genomic_mating.R:
##   usefulness_criterion()          -- all four variance_model modes
##   .selection_intensity()          -- internal, asymptotic + finite-population
##   .block_contrib_independent()    -- internal, exact formula
##   .block_contrib_phased()         -- internal, exact formula
##   .phased_allele_effect()         -- internal, exact formula
##   .block_snp_order()              -- internal
##   .haldane_r()                    -- internal, exact formula (linked mode)
##   .adjacent_r()                   -- internal, exact formula (linked mode)
##   .block_genetic_positions()      -- internal (linked mode, genetic_map)
##   .block_physical_positions()     -- internal (linked mode, ld_matrix
##                                       fallback -- see section 6c)
##   .adjacent_r_from_ld()           -- internal, exact formula (linked mode,
##                                       ld_matrix fallback -- see section 6c)
##   .block_contrib_linked_mc()      -- internal, Monte Carlo (linked mode),
##                                       hand-verified against 2 exact
##                                       closed-form boundary cases (r -> 0.5
##                                       and r = 0) -- see section 6 below
##
## Internal (dot-prefixed) helpers are accessed via HapBlockR::: for direct,
## exact-value unit testing of the underlying arithmetic, in addition to the
## end-to-end usefulness_criterion() tests below. This mirrors the package's
## own file-header comments explaining these are simple, hand-verifiable
## formulas (as opposed to the external-tool-wrapping code, which is tested
## structurally instead since its correctness depends on the external
## package's own behaviour).
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# -- Shared fixtures ------------------------------------------------------------

.si20  <- make_snpinfo(p = 20, chr = "1")
.blk20 <- make_blocks(.si20, n_blocks = 4L)
.phased <- make_phased(n = 40, p = 20, seed = 11L)  # hap1/hap2: 20 SNPs x 40 ind

.haps_ph <- extract_haplotypes(.phased, .si20, .blk20, min_snps = 3L)
.bi_ph   <- attr(.haps_ph, "block_info")

.gebv40 <- setNames(rnorm(40, 0, 1), paste0("ind", 1:40))
.snpfx20 <- setNames(rnorm(20, 0, 0.5), .si20$SNP)

# Small unphased fixture for block_independent mode
.G_u   <- make_geno(n = 30, p = 24, seed = 5L)
.si_u  <- make_snpinfo(p = 24, chr = "1")
.blk_u <- make_blocks(.si_u, n_blocks = 4L)
.pred_u <- run_haplotype_prediction(
  geno_matrix = .G_u, snp_info = .si_u, blocks = .blk_u,
  blues = make_blues(.G_u, seed = 5L, format = "data.frame"),
  id_col = "id", blue_col = "blue",
  marker_effect_method = "gblup", verbose = FALSE
)

# ==============================================================================
# 1. .selection_intensity() (internal)
# ==============================================================================

test_that(".selection_intensity: asymptotic formula matches the classical p=0.1 value (~1.755)", {
  i_sel <- HapBlockR:::.selection_intensity(0.1)
  expect_equal(i_sel, 1.755, tolerance = 0.01)
})

test_that(".selection_intensity: p=1 gives zero intensity (keeping everyone)", {
  expect_equal(HapBlockR:::.selection_intensity(1), 0)
})

test_that(".selection_intensity: errors on p outside (0, 1]", {
  expect_error(HapBlockR:::.selection_intensity(0), "selected_proportion")
  expect_error(HapBlockR:::.selection_intensity(1.5), "selected_proportion")
})

test_that(".selection_intensity: finite n_progeny gives a SMALLER intensity than the asymptotic formula", {
  i_inf    <- HapBlockR:::.selection_intensity(0.1)
  i_finite <- HapBlockR:::.selection_intensity(0.1, n_progeny = 30L, n_sim = 2000L, seed = 1L)
  expect_lt(i_finite, i_inf)
})

test_that(".selection_intensity: finite-population estimate is reproducible with the same seed", {
  a <- HapBlockR:::.selection_intensity(0.1, n_progeny = 40L, n_sim = 1000L, seed = 42L)
  b <- HapBlockR:::.selection_intensity(0.1, n_progeny = 40L, n_sim = 1000L, seed = 42L)
  expect_equal(a, b)
})

test_that(".selection_intensity: does not disturb the caller's global RNG state", {
  set.seed(99L)
  before <- runif(1)
  set.seed(99L)
  invisible(HapBlockR:::.selection_intensity(0.1, n_progeny = 20L, n_sim = 500L, seed = 7L))
  after <- runif(1)
  expect_equal(before, after)
})

test_that(".selection_intensity: k >= n_progeny (keeping ~everyone) returns zero", {
  # k <- max(1, round(p * n_progeny)); the function returns 0 early when
  # k >= n_progeny. p = 0.9, n_progeny = 5 does NOT hit this branch: R's
  # round() uses round-half-to-even ("banker's rounding"), so
  # round(0.9 * 5) = round(4.5) = 4 (rounds to the nearest EVEN integer,
  # not always up), giving k = 4 < n_progeny = 5 -- the Monte Carlo path
  # runs instead (selecting the top 4 of 5), which is a real, nonzero
  # intensity, not a bug. p = 0.99 / n_progeny = 10 avoids the .5-tie
  # rounding case entirely (round(9.9) = 10 unambiguously) and reliably
  # exercises the k >= n_progeny early-return branch this test targets.
  expect_equal(HapBlockR:::.selection_intensity(0.99, n_progeny = 10L, n_sim = 500L), 0)
})

# ==============================================================================
# 2. .block_contrib_independent() / .block_contrib_phased() (internal, exact)
# ==============================================================================

test_that(".block_contrib_independent: matches the hand-derived formula exactly", {
  res <- HapBlockR:::.block_contrib_independent(vi = 2, vj = 0, segregation_factor = 0.5)
  expect_equal(unname(res["mean"]), 1)
  expect_equal(unname(res["var"]),  0.5)
})

test_that(".block_contrib_independent: identical parent values give zero variance", {
  res <- HapBlockR:::.block_contrib_independent(vi = 3, vj = 3, segregation_factor = 0.5)
  expect_equal(unname(res["var"]), 0)
  expect_equal(unname(res["mean"]), 3)
})

test_that(".block_contrib_phased: matches the hand-derived 4-gamete formula exactly", {
  res <- HapBlockR:::.block_contrib_phased(eff_i1 = 0, eff_i2 = 2, eff_j1 = 0, eff_j2 = 2)
  expect_equal(unname(res["mean"]), 2)
  expect_equal(unname(res["var"]),  2)
})

test_that(".block_contrib_phased: identical gametes on both sides give zero variance", {
  res <- HapBlockR:::.block_contrib_phased(eff_i1 = 1, eff_i2 = 1, eff_j1 = 1, eff_j2 = 1)
  expect_equal(unname(res["var"]), 0)
})

# ==============================================================================
# 3. .phased_allele_effect() / .block_snp_order() (internal, exact)
# ==============================================================================

test_that(".phased_allele_effect: sums dose x effect correctly", {
  snp_ids <- c("s1", "s2", "s3")
  fx <- c(s1 = 1, s2 = 2, s3 = 3)
  res <- HapBlockR:::.phased_allele_effect("101", snp_ids, fx)
  expect_equal(res, 1 * 1 + 0 * 2 + 1 * 3)  # = 4
})

test_that(".phased_allele_effect: errors on length mismatch", {
  snp_ids <- c("s1", "s2", "s3")
  fx <- c(s1 = 1, s2 = 2, s3 = 3)
  expect_error(HapBlockR:::.phased_allele_effect("10", snp_ids, fx), "length")
})

test_that(".phased_allele_effect: errors on non-numeric characters", {
  snp_ids <- c("s1", "s2")
  fx <- c(s1 = 1, s2 = 2)
  expect_error(HapBlockR:::.phased_allele_effect(".1", snp_ids, fx), "Non-numeric")
})

test_that(".block_snp_order: returns SNP IDs sorted by position within range", {
  res <- HapBlockR:::.block_snp_order(.si20, chr = "1",
                                      start_bp = .si20$POS[1], end_bp = .si20$POS[5])
  expect_equal(res, .si20$SNP[order(.si20$POS)][1:5])
})

# ==============================================================================
# 4. usefulness_criterion(): variance_model = "block_independent"
# ==============================================================================

test_that("usefulness_criterion (block_independent): returns the documented columns", {
  parents <- rownames(.pred_u$local_gebv)[1:6]
  uc <- usefulness_criterion(
    parent_ids          = parents,
    gebv                 = .pred_u$gebv,
    variance_model        = "block_independent",
    block_importance      = .pred_u$block_importance,
    local_gebv            = .pred_u$local_gebv,
    verbose               = FALSE
  )
  req <- c("parent1", "parent2", "mid_parent_gebv", "predicted_variance",
           "selection_intensity", "UC", "rank")
  expect_true(all(req %in% names(uc)))
  expect_equal(nrow(uc), choose(length(parents), 2L))
})

test_that("usefulness_criterion (block_independent): UC = mid_parent_gebv + i_sel * sqrt(predicted_variance)", {
  parents <- rownames(.pred_u$local_gebv)[1:6]
  uc <- usefulness_criterion(
    parent_ids          = parents,
    gebv                 = .pred_u$gebv,
    variance_model        = "block_independent",
    block_importance      = .pred_u$block_importance,
    local_gebv            = .pred_u$local_gebv,
    verbose               = FALSE
  )
  expected <- uc$mid_parent_gebv + uc$selection_intensity * sqrt(pmax(uc$predicted_variance, 0))
  expect_equal(uc$UC, expected, tolerance = 1e-8)
})

test_that("usefulness_criterion (block_independent): output is ranked by descending UC", {
  parents <- rownames(.pred_u$local_gebv)[1:8]
  uc <- usefulness_criterion(
    parent_ids          = parents,
    gebv                 = .pred_u$gebv,
    variance_model        = "block_independent",
    block_importance      = .pred_u$block_importance,
    local_gebv            = .pred_u$local_gebv,
    verbose               = FALSE
  )
  expect_equal(uc$rank, seq_len(nrow(uc)))
  expect_true(all(diff(uc$UC) <= 1e-8))  # non-increasing
})

test_that("usefulness_criterion: cross_pairs argument scores a specific, non-exhaustive list", {
  parents <- rownames(.pred_u$local_gebv)[1:6]
  pairs <- data.frame(parent1 = parents[1:3], parent2 = parents[4:6],
                      stringsAsFactors = FALSE)
  uc <- usefulness_criterion(
    cross_pairs          = pairs,
    gebv                 = .pred_u$gebv,
    variance_model        = "block_independent",
    block_importance      = .pred_u$block_importance,
    local_gebv            = .pred_u$local_gebv,
    verbose               = FALSE
  )
  expect_equal(nrow(uc), 3L)
})

test_that("usefulness_criterion: errors on unnamed gebv", {
  parents <- rownames(.pred_u$local_gebv)[1:4]
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = unname(.pred_u$gebv),
                         block_importance = .pred_u$block_importance,
                         local_gebv = .pred_u$local_gebv, verbose = FALSE),
    "gebv"
  )
})

test_that("usefulness_criterion: errors with fewer than 2 unique parent_ids", {
  expect_error(
    usefulness_criterion(parent_ids = "ind1", gebv = .pred_u$gebv,
                         block_importance = .pred_u$block_importance,
                         local_gebv = .pred_u$local_gebv, verbose = FALSE),
    "parent_ids"
  )
})

test_that("usefulness_criterion (block_independent): errors when block_importance is missing required columns", {
  parents <- rownames(.pred_u$local_gebv)[1:4]
  bad_bi <- data.frame(block_id = "x")  # missing CHR
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .pred_u$gebv,
                         block_importance = bad_bi,
                         local_gebv = .pred_u$local_gebv, verbose = FALSE),
    "block_importance"
  )
})

test_that("usefulness_criterion (block_independent): selected_proportion closer to 1 shrinks UC toward mid_parent_gebv", {
  parents <- rownames(.pred_u$local_gebv)[1:6]
  uc_tight <- usefulness_criterion(
    parent_ids = parents, gebv = .pred_u$gebv, selected_proportion = 0.9,
    block_importance = .pred_u$block_importance, local_gebv = .pred_u$local_gebv,
    verbose = FALSE
  )
  uc_wide <- usefulness_criterion(
    parent_ids = parents, gebv = .pred_u$gebv, selected_proportion = 0.05,
    block_importance = .pred_u$block_importance, local_gebv = .pred_u$local_gebv,
    verbose = FALSE
  )
  expect_lt(mean(uc_tight$selection_intensity), mean(uc_wide$selection_intensity))
})

# ==============================================================================
# 5. usefulness_criterion(): variance_model = "phased"
# ==============================================================================

test_that("usefulness_criterion (phased): runs on phased haplotypes and returns documented columns", {
  parents <- paste0("ind", 1:6)
  uc <- usefulness_criterion(
    parent_ids       = parents,
    gebv             = .gebv40,
    variance_model   = "phased",
    block_importance = .bi_ph,
    haplotypes       = .haps_ph,
    snp_info         = .si20,
    snp_effects      = .snpfx20,
    verbose          = FALSE
  )
  req <- c("parent1", "parent2", "mid_parent_gebv", "predicted_variance",
           "selection_intensity", "UC", "rank")
  expect_true(all(req %in% names(uc)))
  expect_equal(nrow(uc), choose(length(parents), 2L))
})

test_that("usefulness_criterion (phased): errors when haplotypes/snp_info/snp_effects are missing", {
  parents <- paste0("ind", 1:4)
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "phased", block_importance = .bi_ph,
                         verbose = FALSE),
    "phased"
  )
})

# ==============================================================================
# 6. .block_contrib_linked_mc() (internal): Monte Carlo linked-variance model,
#    hand-verified against two EXACT closed-form boundary cases
# ==============================================================================

# Fixture effect vectors: 3 blocks, arbitrary but fixed additive effects for
# each of the two haplotype copies of each of the two parents (i, j). Used by
# every test in this section. Hand-computed reference values (see comments
# below) were derived directly from these numbers.
.eff_i1_lk <- c(1.0, 0.5, -0.3)
.eff_i2_lk <- c(-0.4, 0.2, 0.6)
.eff_j1_lk <- c(0.3, -0.6, 0.1)
.eff_j2_lk <- c(0.7, 0.4, -0.2)

test_that(".block_contrib_linked_mc: r -> 0.5 (fully independent) converges to the sum of per-block .block_contrib_phased() variances", {
  # When every adjacent pair of blocks has recombination fraction 0.5, blocks
  # assort completely independently -- exactly the assumption variance_model
  # = "phased" already makes. So the linked MC model's total mean/variance
  # must converge (as n_sim grows) to the independent SUM of each block's own
  # exact 4-combo .block_contrib_phased() mean/variance.
  per_block <- lapply(1:3, function(k) {
    HapBlockR:::.block_contrib_phased(.eff_i1_lk[k], .eff_i2_lk[k],
                                      .eff_j1_lk[k], .eff_j2_lk[k])
  })
  expected_mean <- sum(vapply(per_block, `[`, numeric(1), "mean"))
  expected_var  <- sum(vapply(per_block, `[`, numeric(1), "var"))
  # Hand check: block1 mean=0.8/var=0.53, block2 mean=0.25/var=0.2725,
  # block3 mean=0.1/var=0.225 -> totals mean=1.15, var=1.0275.
  expect_equal(expected_mean, 1.15, tolerance = 1e-6)
  expect_equal(expected_var, 1.0275, tolerance = 1e-6)

  contrib <- HapBlockR:::.block_contrib_linked_mc(
    .eff_i1_lk, .eff_i2_lk, .eff_j1_lk, .eff_j2_lk,
    r_adjacent = c(0.5, 0.5), n_sim = 200000L, seed = 1L
  )

  expect_equal(unname(contrib["mean"]), expected_mean, tolerance = 0.02)
  expect_equal(unname(contrib["var"]), expected_var, tolerance = 0.05)
})

test_that(".block_contrib_linked_mc: r = 0 (fully linked) converges to an exact 4-combo enumeration of the fused super-locus", {
  # When every adjacent recombination fraction is exactly 0, each parent
  # transmits ALL of hap1 or ALL of hap2 across every block (no crossover
  # ever occurs) -- so the 3 blocks behave as one fused super-locus with a
  # single summed effect per haplotype copy. The exact distribution is then
  # the same 4-combo enumeration .block_contrib_phased() uses for a single
  # block, applied to the block-summed haplotype effects.
  expected <- HapBlockR:::.block_contrib_phased(
    sum(.eff_i1_lk), sum(.eff_i2_lk), sum(.eff_j1_lk), sum(.eff_j2_lk)
  )
  # Hand check: summed effects (1.2, 0.4, -0.2, 0.9) -> combos
  # (1.0, 2.1, 0.2, 1.3), mean=1.15, var=0.4625.
  expect_equal(unname(expected["mean"]), 1.15, tolerance = 1e-6)
  expect_equal(unname(expected["var"]), 0.4625, tolerance = 1e-6)

  contrib <- HapBlockR:::.block_contrib_linked_mc(
    .eff_i1_lk, .eff_i2_lk, .eff_j1_lk, .eff_j2_lk,
    r_adjacent = c(0, 0), n_sim = 200000L, seed = 1L
  )

  expect_equal(unname(contrib["mean"]), unname(expected["mean"]), tolerance = 0.02)
  expect_equal(unname(contrib["var"]), unname(expected["var"]), tolerance = 0.05)
})

test_that(".block_contrib_linked_mc: intermediate r gives variance between the r=0 and r=0.5 boundary cases", {
  # Sanity check that intermediate linkage produces an answer bounded by the
  # two exact boundary cases -- catches gross implementation errors (e.g.
  # r_adjacent applied backwards) that the two boundary cases alone, both
  # using symmetric r_adjacent vectors, might not surface.
  var_r0   <- unname(HapBlockR:::.block_contrib_linked_mc(
    .eff_i1_lk, .eff_i2_lk, .eff_j1_lk, .eff_j2_lk,
    r_adjacent = c(0, 0), n_sim = 100000L, seed = 2L)["var"])
  var_r05  <- unname(HapBlockR:::.block_contrib_linked_mc(
    .eff_i1_lk, .eff_i2_lk, .eff_j1_lk, .eff_j2_lk,
    r_adjacent = c(0.5, 0.5), n_sim = 100000L, seed = 2L)["var"])
  var_rmid <- unname(HapBlockR:::.block_contrib_linked_mc(
    .eff_i1_lk, .eff_i2_lk, .eff_j1_lk, .eff_j2_lk,
    r_adjacent = c(0.1, 0.1), n_sim = 100000L, seed = 2L)["var"])
  expect_gte(var_rmid, min(var_r0, var_r05) - 0.05)
  expect_lte(var_rmid, max(var_r0, var_r05) + 0.05)
})

test_that(".block_contrib_linked_mc: does not disturb the caller's global RNG state when seed is supplied", {
  set.seed(99L)
  before <- runif(1)
  set.seed(99L)
  invisible(HapBlockR:::.block_contrib_linked_mc(
    .eff_i1_lk, .eff_i2_lk, .eff_j1_lk, .eff_j2_lk,
    r_adjacent = c(0.2, 0.2), n_sim = 5000L, seed = 7L))
  after <- runif(1)
  expect_equal(before, after)
})

# ==============================================================================
# 6b. .haldane_r() / .adjacent_r() / .block_genetic_positions() (internal, exact)
# ==============================================================================

test_that(".haldane_r: r = 0 at zero genetic distance, converges to 0.5 far away", {
  expect_equal(HapBlockR:::.haldane_r(0), 0)
  expect_equal(HapBlockR:::.haldane_r(1000), 0.5, tolerance = 1e-6)
})

test_that(".haldane_r: matches the closed-form value at a known distance (50 cM)", {
  # Haldane 1919: r = 0.5 * (1 - exp(-2d/100)); at d = 50 cM, r = 0.5*(1-exp(-1))
  expect_equal(HapBlockR:::.haldane_r(50), 0.5 * (1 - exp(-1)), tolerance = 1e-8)
})

test_that(".adjacent_r: uses r = 0.5 across a chromosome boundary regardless of cM", {
  bp <- data.frame(block_id = c("b1", "b2"), CHR = c("1", "2"), cM = c(0, 0))
  expect_equal(HapBlockR:::.adjacent_r(bp), 0.5)
})

test_that(".adjacent_r: returns length(block_order) - 1, matching Haldane's formula within a chromosome", {
  bp <- data.frame(block_id = c("b1", "b2", "b3"), CHR = "1", cM = c(0, 10, 10))
  r <- HapBlockR:::.adjacent_r(bp)
  expect_length(r, 2L)
  expect_equal(r[1], HapBlockR:::.haldane_r(10))
  expect_equal(r[2], 0)  # zero genetic distance between b2 and b3
})

test_that(".block_genetic_positions: drops blocks with no SNP in common with genetic_map and sorts by (CHR, cM)", {
  block_coords <- data.frame(block_id = c("bA", "bB", "bC"), CHR = c("1", "1", "1"),
                             start_bp = c(1, 100, 200), end_bp = c(50, 150, 250))
  snp_info <- data.frame(SNP = c("s1", "s2", "s3"), CHR = "1", POS = c(10, 110, 999))
  genetic_map <- data.frame(SNP = c("s1", "s2"), CHR = "1", cM = c(5, 1))
  pos <- HapBlockR:::.block_genetic_positions(block_coords, snp_info, genetic_map)
  expect_equal(pos$block_id, c("bB", "bA"))  # sorted by cM ascending: bB(1) before bA(5)
  expect_false("bC" %in% pos$block_id)       # s3 not in genetic_map -> bC dropped
})

# ==============================================================================
# 6c. .block_physical_positions() / .adjacent_r_from_ld() (internal, exact) --
#     the ld_matrix fallback for variance_model = "linked" when no
#     genetic_map is supplied
# ==============================================================================

test_that(".block_physical_positions: sorts by (CHR, start_bp) and drops nothing (no genetic_map needed)", {
  block_coords <- data.frame(block_id = c("bA", "bB", "bC"), CHR = c("2", "1", "1"),
                             start_bp = c(5, 200, 10), end_bp = c(50, 250, 60))
  pos <- HapBlockR:::.block_physical_positions(block_coords)
  expect_equal(pos$block_id, c("bC", "bB", "bA"))  # CHR "1" before "2"; within CHR1, 10 before 200
  expect_equal(nrow(pos), 3L)  # physical position always known -- nothing dropped
})

test_that(".adjacent_r_from_ld: uses r = 0.5 across a chromosome boundary regardless of LD", {
  block_pos <- data.frame(block_id = c("b1", "b2"), CHR = c("1", "2"), start_bp = c(1, 1))
  snp_ids_by_block <- list(b1 = "s1", b2 = "s2")
  ld <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("s1", "s2"), c("s1", "s2")))
  expect_equal(HapBlockR:::.adjacent_r_from_ld(block_pos, snp_ids_by_block, ld), 0.5)
})

test_that(".adjacent_r_from_ld: r = 0 at LD = 1 (perfect association), r = 0.5 at LD = 0 (no association)", {
  block_pos <- data.frame(block_id = c("b1", "b2"), CHR = c("1", "1"), start_bp = c(1, 100))
  snp_ids_by_block <- list(b1 = "s1", b2 = "s2")

  ld_full <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("s1", "s2"), c("s1", "s2")))
  expect_equal(HapBlockR:::.adjacent_r_from_ld(block_pos, snp_ids_by_block, ld_full), 0)

  ld_zero <- matrix(0, 2, 2, dimnames = list(c("s1", "s2"), c("s1", "s2")))
  expect_equal(HapBlockR:::.adjacent_r_from_ld(block_pos, snp_ids_by_block, ld_zero), 0.5)
})

test_that(".adjacent_r_from_ld: matches the closed-form 0.5 * (1 - mean(LD)) at a known intermediate value", {
  block_pos <- data.frame(block_id = c("b1", "b2"), CHR = c("1", "1"), start_bp = c(1, 100))
  # 2 SNPs in b1, 1 SNP in b2 -> mean LD across the 2 cross-block pairs.
  snp_ids_by_block <- list(b1 = c("s1", "s2"), b2 = "s3")
  ld <- matrix(0, 3, 3, dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3")))
  ld["s1", "s3"] <- ld["s3", "s1"] <- 0.8
  ld["s2", "s3"] <- ld["s3", "s2"] <- 0.4
  # mean(0.8, 0.4) = 0.6 -> r = 0.5 * (1 - 0.6) = 0.2
  expect_equal(HapBlockR:::.adjacent_r_from_ld(block_pos, snp_ids_by_block, ld), 0.2, tolerance = 1e-8)
})

test_that(".adjacent_r_from_ld: falls back to r = 0.5 when neither adjacent block has any SNP covered by ld_matrix", {
  block_pos <- data.frame(block_id = c("b1", "b2"), CHR = c("1", "1"), start_bp = c(1, 100))
  snp_ids_by_block <- list(b1 = "sX", b2 = "sY")  # neither present in ld_matrix below
  ld <- matrix(1, 2, 2, dimnames = list(c("s1", "s2"), c("s1", "s2")))
  expect_equal(HapBlockR:::.adjacent_r_from_ld(block_pos, snp_ids_by_block, ld), 0.5)
})

# ==============================================================================
# 7. usefulness_criterion(): variance_model = "linked"
# ==============================================================================

.gmap20 <- data.frame(SNP = .si20$SNP, CHR = .si20$CHR,
                      cM = (.si20$POS - min(.si20$POS)) / 1000)

test_that("usefulness_criterion (linked): runs on phased haplotypes + genetic_map and returns documented columns", {
  parents <- paste0("ind", 1:6)
  uc <- usefulness_criterion(
    parent_ids       = parents,
    gebv             = .gebv40,
    variance_model   = "linked",
    block_importance = .bi_ph,
    haplotypes       = .haps_ph,
    snp_info         = .si20,
    snp_effects      = .snpfx20,
    genetic_map      = .gmap20,
    n_sim_linked     = 500L,
    seed             = 1L,
    verbose          = FALSE
  )
  req <- c("parent1", "parent2", "mid_parent_gebv", "predicted_variance",
           "selection_intensity", "UC", "rank")
  expect_true(all(req %in% names(uc)))
  expect_equal(nrow(uc), choose(length(parents), 2L))
  expect_true(all(is.na(uc$predicted_variance) | uc$predicted_variance >= 0))
})

test_that("usefulness_criterion (linked): errors when genetic_map is missing", {
  parents <- paste0("ind", 1:4)
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "linked", block_importance = .bi_ph,
                         haplotypes = .haps_ph, snp_info = .si20,
                         snp_effects = .snpfx20, verbose = FALSE),
    "genetic_map"
  )
})

test_that("usefulness_criterion (linked): errors when haplotypes/snp_info/snp_effects are missing", {
  parents <- paste0("ind", 1:4)
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "linked", block_importance = .bi_ph,
                         genetic_map = .gmap20, verbose = FALSE),
    "linked"
  )
})

test_that("usefulness_criterion (linked): errors on n_sim_linked below the minimum", {
  parents <- paste0("ind", 1:4)
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "linked", block_importance = .bi_ph,
                         haplotypes = .haps_ph, snp_info = .si20,
                         snp_effects = .snpfx20, genetic_map = .gmap20,
                         n_sim_linked = 10L, verbose = FALSE),
    "n_sim_linked"
  )
})

test_that("usefulness_criterion (linked): is reproducible with the same seed", {
  parents <- paste0("ind", 1:6)
  args <- list(parent_ids = parents, gebv = .gebv40, variance_model = "linked",
              block_importance = .bi_ph, haplotypes = .haps_ph, snp_info = .si20,
              snp_effects = .snpfx20, genetic_map = .gmap20, n_sim_linked = 500L,
              seed = 3L, verbose = FALSE)
  uc1 <- do.call(usefulness_criterion, args)
  uc2 <- do.call(usefulness_criterion, args)
  expect_equal(uc1$predicted_variance, uc2$predicted_variance)
})

# -- 7b. variance_model = "linked" via ld_matrix instead of genetic_map --------
# (the fallback added so "linked" is usable without a real genetic map --
# still requires phased haplotypes either way, see .adjacent_r_from_ld() in
# R/genomic_mating.R for the recombination-fraction proxy formula, exact-
# value-verified in section 6c above)

set.seed(21L)
.ld20 <- matrix(stats::runif(20 * 20, 0, 1), 20, 20)
.ld20 <- (.ld20 + t(.ld20)) / 2
diag(.ld20) <- 1
dimnames(.ld20) <- list(.si20$SNP, .si20$SNP)

test_that("usefulness_criterion (linked): runs on phased haplotypes + ld_matrix (no genetic_map) and returns documented columns", {
  parents <- paste0("ind", 1:6)
  uc <- usefulness_criterion(
    parent_ids       = parents,
    gebv             = .gebv40,
    variance_model   = "linked",
    block_importance = .bi_ph,
    haplotypes       = .haps_ph,
    snp_info         = .si20,
    snp_effects      = .snpfx20,
    ld_matrix        = .ld20,
    n_sim_linked     = 500L,
    seed             = 1L,
    verbose          = FALSE
  )
  req <- c("parent1", "parent2", "mid_parent_gebv", "predicted_variance",
           "selection_intensity", "UC", "rank")
  expect_true(all(req %in% names(uc)))
  expect_equal(nrow(uc), choose(length(parents), 2L))
  expect_true(all(is.na(uc$predicted_variance) | uc$predicted_variance >= 0))
})

test_that("usefulness_criterion (linked): errors when neither genetic_map nor ld_matrix is supplied", {
  parents <- paste0("ind", 1:4)
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "linked", block_importance = .bi_ph,
                         haplotypes = .haps_ph, snp_info = .si20,
                         snp_effects = .snpfx20, verbose = FALSE),
    "genetic_map"
  )
})

test_that("usefulness_criterion (linked): errors on a malformed ld_matrix (missing dimnames)", {
  parents <- paste0("ind", 1:4)
  bad_ld <- matrix(1, 20, 20)  # no dimnames at all
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "linked", block_importance = .bi_ph,
                         haplotypes = .haps_ph, snp_info = .si20,
                         snp_effects = .snpfx20, ld_matrix = bad_ld,
                         verbose = FALSE),
    "ld_matrix"
  )
})

test_that("usefulness_criterion (linked): still requires phased haplotypes even when ld_matrix is supplied", {
  parents <- paste0("ind", 1:4)
  expect_error(
    usefulness_criterion(parent_ids = parents, gebv = .gebv40,
                         variance_model = "linked", block_importance = .bi_ph,
                         ld_matrix = .ld20, verbose = FALSE),
    "haplotypes"
  )
})

test_that("usefulness_criterion (linked): ld_matrix path is reproducible with the same seed", {
  parents <- paste0("ind", 1:6)
  args <- list(parent_ids = parents, gebv = .gebv40, variance_model = "linked",
              block_importance = .bi_ph, haplotypes = .haps_ph, snp_info = .si20,
              snp_effects = .snpfx20, ld_matrix = .ld20, n_sim_linked = 500L,
              seed = 3L, verbose = FALSE)
  uc1 <- do.call(usefulness_criterion, args)
  uc2 <- do.call(usefulness_criterion, args)
  expect_equal(uc1$predicted_variance, uc2$predicted_variance)
})

# ==============================================================================
# 8. usefulness_criterion(): variance_model = "simplemating"
# ==============================================================================

test_that("usefulness_criterion (simplemating): requires SimpleMating and errors clearly when absent", {
  skip_if(requireNamespace("SimpleMating", quietly = TRUE),
         "SimpleMating is installed; this test targets the absent-dependency error path only")
  expect_error(
    usefulness_criterion(parent_ids = paste0("ind", 1:4), gebv = .gebv40,
                         variance_model = "simplemating", verbose = FALSE),
    "SimpleMating"
  )
})

test_that("usefulness_criterion (simplemating): by default (het_to_na = TRUE), heterozygous (dose=1) calls are converted to NA and it runs successfully", {
  skip_if_not_installed("SimpleMating")
  skip_if_simplemating_too_old("getUsefA")
  geno02 <- matrix(0, nrow = 6, ncol = 5,
                   dimnames = list(paste0("ind", 1:6), paste0("rs", 1:5)))
  geno02[1, 1] <- 1  # heterozygous call -- auto-converted to NA by default
  Gmat <- diag(6); dimnames(Gmat) <- list(rownames(geno02), rownames(geno02))
  gm <- data.frame(SNP = paste0("rs", 1:5), CHR = "1", cM = seq(0, 4, by = 1))
  fx <- setNames(rnorm(5), paste0("rs", 1:5))
  expect_message(
    uc <- usefulness_criterion(
      parent_ids = rownames(geno02)[1:4], gebv = setNames(rnorm(6), rownames(geno02)),
      variance_model = "simplemating", geno_matrix = geno02, G = Gmat,
      genetic_map = gm, snp_effects = fx, verbose = TRUE
      # het_to_na left at its default (TRUE)
    ),
    "heterozygous"
  )
  req <- c("parent1", "parent2", "mid_parent_gebv", "predicted_variance",
           "selection_intensity", "UC", "rank")
  expect_true(all(req %in% names(uc)))
})

test_that("usefulness_criterion (simplemating): het_to_na = FALSE restores the strict error-on-heterozygous behaviour", {
  skip_if_not_installed("SimpleMating")
  skip_if_simplemating_too_old("getUsefA")
  geno02 <- matrix(0, nrow = 6, ncol = 5,
                   dimnames = list(paste0("ind", 1:6), paste0("rs", 1:5)))
  geno02[1, 1] <- 1  # heterozygous call, disallowed when het_to_na = FALSE
  Gmat <- diag(6); dimnames(Gmat) <- list(rownames(geno02), rownames(geno02))
  gm <- data.frame(SNP = paste0("rs", 1:5), CHR = "1", cM = seq(0, 4, by = 1))
  fx <- setNames(rnorm(5), paste0("rs", 1:5))
  expect_error(
    usefulness_criterion(
      parent_ids = rownames(geno02)[1:4], gebv = setNames(rnorm(6), rownames(geno02)),
      variance_model = "simplemating", geno_matrix = geno02, G = Gmat,
      genetic_map = gm, snp_effects = fx, het_to_na = FALSE, verbose = FALSE
    ),
    "heterozygous"
  )
})

test_that("usefulness_criterion (simplemating): still errors on values outside 0/1/2/NA regardless of het_to_na", {
  skip_if_not_installed("SimpleMating")
  skip_if_simplemating_too_old("getUsefA")
  geno_bad <- matrix(0, nrow = 6, ncol = 5,
                     dimnames = list(paste0("ind", 1:6), paste0("rs", 1:5)))
  geno_bad[1, 1] <- 9  # not a recognised dosage value under any encoding
  Gmat <- diag(6); dimnames(Gmat) <- list(rownames(geno_bad), rownames(geno_bad))
  gm <- data.frame(SNP = paste0("rs", 1:5), CHR = "1", cM = seq(0, 4, by = 1))
  fx <- setNames(rnorm(5), paste0("rs", 1:5))
  expect_error(
    usefulness_criterion(
      parent_ids = rownames(geno_bad)[1:4], gebv = setNames(rnorm(6), rownames(geno_bad)),
      variance_model = "simplemating", geno_matrix = geno_bad, G = Gmat,
      genetic_map = gm, snp_effects = fx, verbose = FALSE
    ),
    "0/1/2"
  )
})

test_that("usefulness_criterion (simplemating): errors when neither genetic_map nor ld_matrix supplied", {
  skip_if_not_installed("SimpleMating")
  skip_if_simplemating_too_old("getUsefA")
  geno02 <- matrix(sample(c(0, 2), 30, replace = TRUE), nrow = 6, ncol = 5,
                   dimnames = list(paste0("ind", 1:6), paste0("rs", 1:5)))
  Gmat <- diag(6); dimnames(Gmat) <- list(rownames(geno02), rownames(geno02))
  fx <- setNames(rnorm(5), paste0("rs", 1:5))
  expect_error(
    usefulness_criterion(
      parent_ids = rownames(geno02)[1:4], gebv = setNames(rnorm(6), rownames(geno02)),
      variance_model = "simplemating", geno_matrix = geno02, G = Gmat,
      snp_effects = fx, verbose = FALSE
    ),
    "genetic_map"
  )
})
