library(testthat)
library(HapBlockR)

# ==============================================================================
# select_parents_by_family()
# ==============================================================================

# 4 families x 3 members each = 12 candidates. Scores hand-constructed so the
# family ranking and within-family ranking are both unambiguous:
#   FamA: 12, 11, 10   (top-3 mean = 11)   <- best family
#   FamB:  9,  8,  7   (top-3 mean = 8)
#   FamC:  6,  5,  4   (top-3 mean = 5)
#   FamD:  3,  2,  1   (top-3 mean = 2)    <- worst family
.fs_score <- setNames(
  c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
  c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3")
)
.fs_family <- setNames(
  c(rep("FamA", 3), rep("FamB", 3), rep("FamC", 3), rep("FamD", 3)),
  names(.fs_score)
)

# Converts an individuals x blocks character matrix into the shape
# extract_haplotypes() actually returns: a named list, one element per
# block, each a named character vector of one allele string per
# individual -- so tests can build a `haplotypes` fixture the same way a
# real pipeline would receive it, without hand-writing extract_haplotypes()
# internals.
.fs_make_haplotypes <- function(bh_matrix) {
  stats::setNames(
    lapply(colnames(bh_matrix), function(blk) {
      stats::setNames(bh_matrix[, blk], rownames(bh_matrix))
    }),
    colnames(bh_matrix)
  )
}

test_that("select_parents_by_family(): unnamed score errors", {
  expect_error(
    select_parents_by_family(unname(.fs_score), .fs_family, 2, 2),
    "score"
  )
})

test_that("select_parents_by_family(): unnamed family errors", {
  expect_error(
    select_parents_by_family(.fs_score, unname(.fs_family), 2, 2),
    "family"
  )
})

test_that("select_parents_by_family(): n_families/n_per_family < 1 error", {
  expect_error(select_parents_by_family(.fs_score, .fs_family, 0, 2), "n_families")
  expect_error(select_parents_by_family(.fs_score, .fs_family, 2, 0), "n_per_family")
})

test_that("select_parents_by_family(): ranks families by top-k mean and picks the best n_families", {
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2)
  expect_equal(res$family_ranking$family[res$family_ranking$rank == 1L], "FamA")
  expect_equal(res$family_ranking$family[res$family_ranking$rank == 2L], "FamB")
  expect_setequal(unique(res$by_family$family), c("FamA", "FamB"))
})

test_that("select_parents_by_family(): selects the top n_per_family lines within each chosen family", {
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2)
  expect_setequal(res$selected, c("A1", "A2", "B1", "B2"))
  expect_equal(nrow(res$by_family), 4L)
})

test_that("select_parents_by_family(): selected has length n_families * n_per_family in the unconstrained case", {
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 3, n_per_family = 2)
  expect_equal(length(res$selected), 6L)
})

test_that("select_parents_by_family(): family_ranking$topk_mean matches hand-computed values", {
  # n_per_family = 2 here (not 3): with n_per_family = 3, every family in
  # this fixture would have exactly 3 == 3 eligible members and be excluded
  # entirely by the new family-size eligibility rule (size <= quota). rank_k
  # = 3 independently controls the ranking sample size, so topk_mean is
  # still each family's full top-3 mean.
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 4, n_per_family = 2,
                                  rank_k = 3)
  fr <- res$family_ranking
  expect_equal(fr$topk_mean[fr$family == "FamA"], 11)
  expect_equal(fr$topk_mean[fr$family == "FamB"], 8)
  expect_equal(fr$topk_mean[fr$family == "FamC"], 5)
  expect_equal(fr$topk_mean[fr$family == "FamD"], 2)
})

test_that("select_parents_by_family(): n_families exceeding available families warns and uses all", {
  expect_warning(
    res <- select_parents_by_family(.fs_score, .fs_family, n_families = 10, n_per_family = 2),
    "n_families"
  )
  expect_equal(res$n_families, 4L)
  expect_setequal(unique(res$by_family$family), c("FamA", "FamB", "FamC", "FamD"))
})

test_that("select_parents_by_family(): n_per_family exceeding every family's size now errors via the eligibility rule (superseded 'warns and takes all members' behaviour)", {
  # Previously (before the family-size eligibility rule), requesting a
  # quota larger than a family's size truncated to the family's full
  # membership with a warning. That code path is now unreachable for a
  # scalar n_per_family: a family whose size does not exceed the quota is
  # excluded outright before ranking, so every one of the 4 size-3 families
  # here is excluded when n_per_family = 10, leaving nothing to select.
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1, n_per_family = 10,
                             verbose = FALSE),
    "quota's worth"
  )
})

test_that("select_parents_by_family(): min_sel_value floor excludes low-merit candidates before ranking", {
  # Keep only score >= 6 (drops C2, C3, D1-D3 -- 5 of 12 candidates). FamD
  # has no eligible members left at all (dropped before the eligibility
  # rule even runs); FamC has exactly 1 eligible member left (C1) which,
  # with n_per_family = 2, is now ALSO excluded by the family-size
  # eligibility rule (1 <= 2 -- no genuine within-group selection possible).
  # Only FamA and FamB (3 eligible members each) survive -- fewer than the
  # 4 requested, so this also exercises the n_families-shortfall warning.
  # NOTE: expect_warning()/expect_message() return the captured CONDITION
  # object, not the wrapped expression's value -- the assignment must
  # happen INSIDE the wrapped expression, not around the whole
  # expect_warning() call, or `res` ends up being the warning itself.
  expect_warning(
    res <- select_parents_by_family(.fs_score, .fs_family, n_families = 4, n_per_family = 2,
                                    min_sel_value = 6, min_sel_mode = "value", verbose = FALSE),
    "n_families"
  )
  expect_true(all(res$by_family$score >= 6))
  expect_equal(res$n_families, 2L)
  expect_setequal(unique(res$by_family$family), c("FamA", "FamB"))
  expect_equal(res$excluded_groups, "FamC")
})

test_that("select_parents_by_family(): min_sel_value with no candidates clearing the floor errors", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                             min_sel_value = 100, min_sel_mode = "value"),
    "min_sel_value"
  )
})

# -- ensure_haplotype_diversity ------------------------------------------------

test_that("select_parents_by_family(): ensure_haplotype_diversity = TRUE without value_matrix errors", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 1,
                             ensure_haplotype_diversity = TRUE),
    "value_matrix"
  )
})

test_that("select_parents_by_family(): ensure_haplotype_diversity avoids dominant-block collisions when possible", {
  # A1 (FamA's best) and B1 (FamB's best) share the SAME dominant block
  # (blk1). A2 (FamA's #2) has a DIFFERENT dominant block (blk2). With
  # diversity mode on, FamB's pick should be pushed to avoid colliding with
  # FamA's pick on blk1 where an alternative exists in FamB itself.
  vmat <- matrix(0, nrow = 12, ncol = 3,
                dimnames = list(names(.fs_score), c("blk1", "blk2", "blk3")))
  vmat["A1", "blk1"] <- 10; vmat["A1", "blk2"] <- 1; vmat["A1", "blk3"] <- 1
  vmat["B1", "blk1"] <- 10; vmat["B1", "blk2"] <- 1; vmat["B1", "blk3"] <- 1
  vmat["B2", "blk1"] <- 1;  vmat["B2", "blk2"] <- 10; vmat["B2", "blk3"] <- 1

  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 1,
                                  ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                                  diversity_method = "dominant_block")
  expect_equal(res$by_family$individual[res$by_family$family == "FamA"], "A1")
  # FamB's single slot should have skipped B1 (dominant block blk1, already
  # claimed by A1) in favour of B2 (dominant block blk2, free) -- score order
  # inside FamB is B1 > B2 > B3, so this only happens because of the
  # diversity adjustment, not because B2 outscored B1.
  expect_equal(res$by_family$individual[res$by_family$family == "FamB"], "B2")
  expect_false(any(res$by_family$collision))
})

test_that("select_parents_by_family(): forced collision is flagged when no alternative exists", {
  # Every candidate in both families has the SAME dominant block -- no
  # avoiding a collision once both families need to place a pick there.
  vmat <- matrix(0, nrow = 12, ncol = 2,
                dimnames = list(names(.fs_score), c("blk1", "blk2")))
  vmat[, "blk1"] <- 10
  vmat[, "blk2"] <- 1

  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 1,
                                  ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                                  diversity_method = "dominant_block")
  expect_equal(res$by_family$individual, c("A1", "B1"))  # merit rank still wins
  expect_equal(res$by_family$collision, c(FALSE, TRUE))  # FamA claims blk1 first
})

test_that("select_parents_by_family(): haplotypes (auto-derived allele matrix) lets a dominant-block collision resolve by allele", {
  # A1 and B1 both dominant at blk1, but carry DIFFERENT alleles there --
  # with `haplotypes` supplied (the extract_haplotypes()-shaped list), B1
  # should still be accepted (no forced collision), since the underlying
  # favourable variant actually differs. No manual matrix construction --
  # the function derives the individuals x blocks allele matrix itself.
  vmat <- matrix(0, nrow = 12, ncol = 2,
                dimnames = list(names(.fs_score), c("blk1", "blk2")))
  vmat[, "blk1"] <- 10
  vmat[, "blk2"] <- 1

  bh <- matrix("ref", nrow = 12, ncol = 2,
              dimnames = list(names(.fs_score), c("blk1", "blk2")))
  bh["A1", "blk1"] <- "hapX"
  bh["B1", "blk1"] <- "hapY"   # different allele from A1 at the shared block
  haps <- .fs_make_haplotypes(bh)

  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 1,
                                  ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                                  haplotypes = haps, diversity_method = "dominant_block")
  expect_equal(res$by_family$individual, c("A1", "B1"))
  expect_equal(res$by_family$collision, c(FALSE, FALSE))  # allele differs -> not a real collision
})

test_that("select_parents_by_family(): value_matrix missing eligible candidates errors", {
  vmat <- matrix(0, nrow = 6, ncol = 2,
                dimnames = list(names(.fs_score)[1:6], c("blk1", "blk2")))
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 4, n_per_family = 1,
                             ensure_haplotype_diversity = TRUE, value_matrix = vmat),
    "value_matrix"
  )
})

test_that("select_parents_by_family(): haplotypes missing a needed block errors", {
  vmat <- matrix(0, nrow = 12, ncol = 2,
                dimnames = list(names(.fs_score), c("blk1", "blk2")))
  bh <- matrix("ref", nrow = 12, ncol = 1,
              dimnames = list(names(.fs_score), c("blk1")))  # missing blk2
  haps <- .fs_make_haplotypes(bh)                            # -> only "blk1"
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 1,
                             ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                             haplotypes = haps, diversity_method = "dominant_block"),
    "haplotypes"
  )
})

test_that("select_parents_by_family(): haplotypes missing an eligible individual in a block errors", {
  vmat <- matrix(0, nrow = 12, ncol = 2,
                dimnames = list(names(.fs_score), c("blk1", "blk2")))
  bh <- matrix("ref", nrow = 12, ncol = 2,
              dimnames = list(names(.fs_score), c("blk1", "blk2")))
  haps <- .fs_make_haplotypes(bh)
  haps[["blk1"]] <- haps[["blk1"]][names(haps[["blk1"]]) != "A1"]  # drop A1
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 1,
                             ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                             haplotypes = haps, diversity_method = "dominant_block"),
    "haplotypes"
  )
})

# -- Backward-compatibility regression check --------------------------------

test_that("select_parents_by_family(): default shrinkage-corrected ranking reproduces the same selected set as legacy topk_mean ranking on the balanced fixture (equal family sizes)", {
  res_shrunk <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2)
  res_legacy <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                         family_rank_method = "topk_mean")
  expect_equal(res_shrunk$selected, res_legacy$selected)
  expect_true(all(!is.na(res_shrunk$family_ranking$shrinkage_weight)))
  expect_true(all(is.na(res_legacy$family_ranking$shrinkage_weight)))
})

test_that("select_parents_by_family(): group_by = 'family' (default) without a family argument errors", {
  expect_error(
    select_parents_by_family(.fs_score, n_families = 2, n_per_family = 2),
    "family"
  )
})

# -- Shrinkage-corrected family ranking (variance components) ---------------

test_that(".estimate_family_variance_components(): matches hand-computed values on the balanced fixture", {
  # grand_mean = 6.5; fam_means = 11,8,5,2; fam_ns = 3,3,3,3 (a=4, N=12)
  # SSB = 3*[(11-6.5)^2+(8-6.5)^2+(5-6.5)^2+(2-6.5)^2] = 135; MSB = 45
  # SSW = 4*2 = 8 (each family: 12,11,10 around mean 11 -> devs 1,0,-1); MSW = 1
  # n0 = (12 - 4*9/12)/3 = 3; tau2 = (45-1)/3 = 44/3
  vc <- HapBlockR:::.estimate_family_variance_components(.fs_score, .fs_family,
                                                          method = "anova", verbose = FALSE)
  expect_true(vc$ok)
  expect_equal(vc$sigma2, 1)
  expect_equal(vc$tau2, 44 / 3, tolerance = 1e-8)
  expect_equal(vc$n0, 3)
})

test_that(".estimate_family_variance_components(): not enough structure to estimate returns ok = FALSE", {
  score2 <- setNames(c(10, 5), c("X1", "X2"))
  fam2   <- setNames(c("FamA", "FamB"), c("X1", "X2"))
  vc <- HapBlockR:::.estimate_family_variance_components(score2, fam2,
                                                          method = "anova", verbose = FALSE)
  expect_false(vc$ok)
})

test_that(".rank_families(): shrunk_topk_mean preserves rank order and matches hand-computed values on the balanced fixture", {
  # target = mean(11,8,5,2) = 6.5; w = 3*(44/3) / (3*(44/3)+1) = 44/45 for
  # every family (equal size -> equal weight, a uniform affine transform of
  # topk_mean, so rank order is provably unaffected by shrinkage here).
  fr <- HapBlockR:::.rank_families(.fs_score, .fs_family, k = 3,
                                   method = "shrunk_topk_mean",
                                   variance_method = "anova",
                                   bias_correction = "none",
                                   G_fam = NULL, verbose = FALSE)
  fr <- fr[order(fr$family), ]
  expect_equal(fr$rank_score[fr$family == "FamA"], 10.9, tolerance = 1e-6)
  expect_equal(fr$rank_score[fr$family == "FamB"], 119.5 / 15, tolerance = 1e-6)
  expect_equal(fr$rank_score[fr$family == "FamC"], 75.5 / 15, tolerance = 1e-6)
  expect_equal(fr$rank_score[fr$family == "FamD"], 2.1, tolerance = 1e-6)
  expect_equal(fr$shrinkage_weight, rep(44 / 45, 4), tolerance = 1e-8)
  ord <- fr$family[order(-fr$rank_score)]
  expect_equal(ord, c("FamA", "FamB", "FamC", "FamD"))
})

test_that(".rank_families(): falls back to unshrunk ranking (with message) when variance components can't be estimated", {
  score2 <- setNames(c(10, 5), c("X1", "X2"))
  fam2   <- setNames(c("FamA", "FamB"), c("X1", "X2"))
  expect_message(
    fr <- HapBlockR:::.rank_families(score2, fam2, k = 1, method = "shrunk_topk_mean",
                                     variance_method = "anova", bias_correction = "none",
                                     G_fam = NULL, verbose = TRUE),
    "falling back"
  )
  expect_true(all(is.na(fr$shrinkage_weight)))
  expect_equal(fr$rank_score, fr$topk_mean)
})

test_that(".rank_families(): method = 'topk_mean' never shrinks", {
  fr <- HapBlockR:::.rank_families(.fs_score, .fs_family, k = 3, method = "topk_mean",
                                   verbose = FALSE)
  expect_true(all(is.na(fr$shrinkage_weight)))
  expect_equal(fr$rank_score, fr$topk_mean)
})

# -- rank_k, decoupled from n_per_family -------------------------------------

# FamE has one huge outlier and two very weak members; FamF is flat and
# moderate throughout -- engineered so which family ranks first flips
# between rank_k = 1 (FamE wins on its single best member) and rank_k = 2
# (FamF wins once FamE's weak second member drags its mean down).
.rk_score  <- setNames(c(100, -100, -100, 10, 10, 10),
                       c("E1", "E2", "E3", "F1", "F2", "F3"))
.rk_family <- setNames(c(rep("FamE", 3), rep("FamF", 3)), names(.rk_score))

test_that("select_parents_by_family(): rank_k decouples the ranking sample size from n_per_family", {
  res1 <- select_parents_by_family(.rk_score, .rk_family, n_families = 1, n_per_family = 1,
                                   rank_k = 1, family_rank_method = "topk_mean")
  expect_equal(res1$family_ranking$family[res1$family_ranking$rank == 1L], "FamE")

  res2 <- select_parents_by_family(.rk_score, .rk_family, n_families = 1, n_per_family = 1,
                                   rank_k = 2, family_rank_method = "topk_mean")
  expect_equal(res2$family_ranking$family[res2$family_ranking$rank == 1L], "FamF")
})

test_that("select_parents_by_family(): rank_k defaults to n_per_family when not supplied", {
  res <- select_parents_by_family(.rk_score, .rk_family, n_families = 1, n_per_family = 1,
                                  family_rank_method = "topk_mean")
  expect_equal(res$rank_k, 1L)
  expect_equal(res$family_ranking$family[res$family_ranking$rank == 1L], "FamE")
})

test_that("select_parents_by_family(): rank_k < 1 errors", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2, rank_k = 0),
    "rank_k"
  )
})

# -- diversity_method = "coverage_gain" --------------------------------------

test_that("select_parents_by_family(): diversity_method = 'coverage_gain' looks past a redundant dominant-block match", {
  # B1 shares A1's exact value profile (fully redundant -> zero marginal
  # coverage gain); B2 has a different, genuinely additive profile. A
  # dominant-block check alone would just flag a forced collision and take
  # B1 anyway (same dominant block, no allele info); coverage_gain instead
  # looks past B1 to B2, which actually adds coverage.
  # A0 is a padding member added ONLY so FamA's eligible size (2) exceeds
  # its n_per_family quota (1) and survives the family-size eligibility
  # rule (a single-member family requesting exactly 1 pick is, by
  # definition, no genuine selection at all, and would otherwise be
  # excluded outright). A0's score/value profile is deliberately
  # uncompetitive (lowest score, all-zero value row) so it never affects
  # which individual FamA actually picks.
  vmat <- matrix(0, nrow = 4, ncol = 2,
                dimnames = list(c("A1", "A0", "B1", "B2"), c("blk1", "blk2")))
  vmat["A1", "blk1"] <- 10; vmat["A1", "blk2"] <- 0
  vmat["B1", "blk1"] <- 10; vmat["B1", "blk2"] <- 0
  vmat["B2", "blk1"] <- 0;  vmat["B2", "blk2"] <- 5

  score_cg  <- setNames(c(10, 1, 9, 8), c("A1", "A0", "B1", "B2"))
  family_cg <- setNames(c("FamA", "FamA", "FamB", "FamB"), c("A1", "A0", "B1", "B2"))

  res <- select_parents_by_family(score_cg, family_cg, n_families = 2, n_per_family = 1,
                                  ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                                  diversity_method = "coverage_gain", verbose = FALSE)
  expect_equal(res$by_family$individual[res$by_family$family == "FamA"], "A1")
  expect_equal(res$by_family$individual[res$by_family$family == "FamB"], "B2")
  expect_false(any(res$by_family$collision))
})

# -- Uneven per-family/per-group quotas --------------------------------------

test_that(".resolve_n_per_family(): flat integer expands to every chosen group", {
  q <- HapBlockR:::.resolve_n_per_family(2L, c("FamA", "FamB"))
  expect_equal(unname(q), c(2L, 2L))
  expect_equal(names(q), c("FamA", "FamB"))
})

test_that(".resolve_n_per_family(): named vector is honoured per group", {
  q <- HapBlockR:::.resolve_n_per_family(c(FamA = 3L, FamB = 1L), c("FamA", "FamB"))
  expect_equal(unname(q[c("FamA", "FamB")]), c(3L, 1L))
})

test_that(".resolve_n_per_family(): unnamed multi-value vector errors", {
  expect_error(HapBlockR:::.resolve_n_per_family(c(3L, 1L), c("FamA", "FamB")), "NAMED")
})

test_that(".resolve_n_per_family(): missing an entry for a chosen group errors", {
  expect_error(
    HapBlockR:::.resolve_n_per_family(c(FamA = 3L), c("FamA", "FamB")),
    "missing an entry"
  )
})

test_that("select_parents_by_family(): uneven n_per_family (named vector) is honoured end-to-end", {
  # rank_k = 1 explicitly: with a named-vector n_per_family, the family-size
  # eligibility rule's threshold falls back to rank_k (n_per_family is only
  # defined for whichever families end up chosen, not known yet at
  # exclusion time). rank_k's default (max(n_per_family) = 3) would exclude
  # every size-3 family in this fixture outright; rank_k = 1 keeps all of
  # them eligible while leaving the family ordering (by top-1 score, still
  # FamA > FamB > FamC > FamD) unchanged.
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, rank_k = 1,
                                  n_per_family = c(FamA = 3L, FamB = 1L))
  expect_equal(sum(res$by_family$family == "FamA"), 3L)
  expect_equal(sum(res$by_family$family == "FamB"), 1L)
  expect_setequal(res$selected, c("A1", "A2", "A3", "B1"))
})

# -- Relatedness control (within_group_target_degree) ------------------------

# W is the top scorer; W-X are close (0.9), W-Y are the most distant pair in
# the whole family (-0.3 -- the unique minimum), everything else is 0.1 --
# engineered so the gain-end reference (top-2 by score) is unambiguously
# {W, X} and the diversity-end reference (maximin) is unambiguously {W, Y}.
.wg_ids <- c("W", "X", "Y", "Z")
G_wg <- diag(4)
dimnames(G_wg) <- list(.wg_ids, .wg_ids)
G_wg["W", "X"] <- G_wg["X", "W"] <- 0.9
G_wg["W", "Y"] <- G_wg["Y", "W"] <- -0.3
G_wg["W", "Z"] <- G_wg["Z", "W"] <- 0.1
G_wg["X", "Y"] <- G_wg["Y", "X"] <- 0.1
G_wg["X", "Z"] <- G_wg["Z", "X"] <- 0.1
G_wg["Y", "Z"] <- G_wg["Z", "Y"] <- 0.1

test_that(".rank_within_group_by_relatedness(): interpolates between the gain-end and diversity-end references", {
  r0 <- HapBlockR:::.rank_within_group_by_relatedness(.wg_ids, G_wg, n_take = 2, target_degree = 0)
  expect_equal(r0$gain_end, 0.9)
  expect_equal(r0$diversity_end, -0.3)
  expect_equal(r0$ceiling, 0.9)
  expect_setequal(r0$picked, c("W", "X"))

  r90 <- HapBlockR:::.rank_within_group_by_relatedness(.wg_ids, G_wg, n_take = 2, target_degree = 90)
  expect_equal(r90$ceiling, -0.3, tolerance = 1e-8)
  expect_setequal(r90$picked, c("W", "Y"))
})

test_that("select_parents_by_family(): within_group_target_degree shifts the within-group pick toward diversity", {
  score_wg  <- setNames(c(10, 9, 8, 7), .wg_ids)
  family_wg <- setNames(rep("FamW", 4), .wg_ids)

  res0 <- select_parents_by_family(score_wg, family_wg, n_families = 1, n_per_family = 2,
                                   G = G_wg, within_group_target_degree = 0, verbose = FALSE)
  expect_setequal(res0$selected, c("W", "X"))

  res90 <- select_parents_by_family(score_wg, family_wg, n_families = 1, n_per_family = 2,
                                    G = G_wg, within_group_target_degree = 90, verbose = FALSE)
  expect_setequal(res90$selected, c("W", "Y"))
  expect_equal(res90$mean_relationship, -0.3, tolerance = 1e-8)
})

test_that("select_parents_by_family(): within_group_target_degree out of [0, 90] errors", {
  G_id <- diag(12); dimnames(G_id) <- list(names(.fs_score), names(.fs_score))
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                             G = G_id, within_group_target_degree = 100),
    "within_group_target_degree"
  )
})

test_that("select_parents_by_family(): within_group_target_degree without G errors", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                             within_group_target_degree = 30),
    "requires G"
  )
})

# -- Unconditional relatedness diagnostics ------------------------------------

test_that("select_parents_by_family(): supplying G alone (no within_group_target_degree) still reports mean_relationship diagnostics", {
  G_id <- diag(12); dimnames(G_id) <- list(names(.fs_score), names(.fs_score))
  G_id["A1", "A2"] <- G_id["A2", "A1"] <- 0.4
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2, G = G_id)
  expect_false(is.na(res$mean_relationship))
  expect_false(any(is.na(res$family_ranking$mean_relationship[res$family_ranking$selected])))
})

test_that("select_parents_by_family(): mean_relationship is NA when G is not supplied", {
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2)
  expect_true(is.na(res$mean_relationship))
})

# -- Grouping: pedigree family vs. genetic cluster ---------------------------

# P-Q and R-S are each tightly related (0.8); every cross-pair is 0 --
# an unambiguous two-cluster structure. Pedigree family labels deliberately
# CROSS-CUT the genetic clusters (FamX = {P, R}, FamY = {Q, S}) so the
# crosswalk test actually exercises independence between the two groupings.
.fc_ids <- c("P", "Q", "R", "S")
G_fc <- diag(4)
dimnames(G_fc) <- list(.fc_ids, .fc_ids)
G_fc["P", "Q"] <- G_fc["Q", "P"] <- 0.8
G_fc["R", "S"] <- G_fc["S", "R"] <- 0.8
score_fc  <- setNames(c(10, 9, 8, 7), .fc_ids)
family_fc <- setNames(c("FamX", "FamY", "FamX", "FamY"), .fc_ids)

test_that(".cluster_by_relationship(): recovers two well-separated clusters, labelled by descending mean score", {
  cl <- HapBlockR:::.cluster_by_relationship(G_fc, .fc_ids, n_clusters = 2,
                                             cluster_method = "ward.D2", score = score_fc)
  expect_equal(unname(cl["P"]), unname(cl["Q"]))
  expect_equal(unname(cl["R"]), unname(cl["S"]))
  expect_false(unname(cl["P"]) == unname(cl["R"]))
  expect_equal(unname(cl["P"]), "cluster_1")  # P,Q mean score 9.5 > R,S mean 7.5
  expect_equal(unname(cl["R"]), "cluster_2")
})

test_that(".cluster_by_relationship(): n_clusters = 1 puts everyone in cluster_1", {
  cl <- HapBlockR:::.cluster_by_relationship(G_fc, .fc_ids, n_clusters = 1,
                                             cluster_method = "ward.D2", score = score_fc)
  expect_true(all(cl == "cluster_1"))
})

test_that(".cluster_by_relationship(): n_clusters exceeding candidate count errors", {
  expect_error(
    HapBlockR:::.cluster_by_relationship(G_fc, .fc_ids, n_clusters = 5,
                                         cluster_method = "ward.D2", score = score_fc),
    "n_clusters"
  )
})

test_that("select_parents_by_family(): group_by = 'genetic_cluster' selects by cluster and still echoes pedigree family", {
  res <- select_parents_by_family(score_fc, family = family_fc, n_families = 2, n_per_family = 1,
                                  group_by = "genetic_cluster", G = G_fc, n_clusters = 2,
                                  verbose = FALSE)
  expect_equal(res$group_by, "genetic_cluster")
  expect_setequal(res$selected, c("P", "R"))  # top scorer of each cluster
  expect_equal(res$by_family$genetic_group[res$by_family$individual == "P"], "cluster_1")
  expect_equal(res$by_family$genetic_group[res$by_family$individual == "R"], "cluster_2")
  # Pedigree family is echoed for cross-referencing only -- unaffected by
  # the cluster-based grouping decision itself.
  expect_equal(res$by_family$family[res$by_family$individual == "P"], "FamX")
  expect_equal(res$by_family$family[res$by_family$individual == "R"], "FamX")
})

test_that("select_parents_by_family(): group_by = 'genetic_cluster' without G errors", {
  expect_error(
    select_parents_by_family(score_fc, n_families = 2, n_per_family = 1,
                             group_by = "genetic_cluster", n_clusters = 2),
    "relationship matrix"
  )
})

test_that("select_parents_by_family(): group_by = 'genetic_cluster' without n_clusters errors", {
  expect_error(
    select_parents_by_family(score_fc, n_families = 2, n_per_family = 1,
                             group_by = "genetic_cluster", G = G_fc),
    "n_clusters"
  )
})

# ==============================================================================
# Family-size eligibility rule
# ==============================================================================

test_that("select_parents_by_family(): a family whose size does not exceed n_per_family is excluded (scalar quota)", {
  # 4 families of exactly 3 members each (.fs_family). n_per_family = 3 ->
  # size (3) <= quota (3) for EVERY family -> nothing survives.
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1, n_per_family = 3,
                             verbose = FALSE),
    "quota's worth"
  )
})

test_that("select_parents_by_family(): excluded_groups reports exactly the undersized families, and messages (unless verbose = FALSE)", {
  # Give FamD a single extra, very-low-score member so it's the ONLY
  # family whose size (4) exceeds the n_per_family = 3 quota -- FamA/FamB/
  # FamC (size 3 each) are excluded, FamD (size 4) survives alone.
  score2  <- c(.fs_score, setNames(-100, "D4"))
  family2 <- c(.fs_family, setNames("FamD", "D4"))

  expect_message(
    res <- select_parents_by_family(score2, family2, n_families = 1, n_per_family = 3,
                                    verbose = TRUE),
    "Excluding"
  )
  expect_setequal(res$excluded_groups, c("FamA", "FamB", "FamC"))
  expect_equal(unique(res$by_family$family), "FamD")
  expect_setequal(res$selected, c("D1", "D2", "D3"))  # D4 (score -100) not competitive

  # verbose = FALSE suppresses the message but the exclusion itself is
  # unaffected -- excluded_groups is populated identically.
  res_quiet <- select_parents_by_family(score2, family2, n_families = 1, n_per_family = 3,
                                        verbose = FALSE)
  expect_setequal(res_quiet$excluded_groups, c("FamA", "FamB", "FamC"))
})

test_that("select_parents_by_family(): a family whose size exceeds n_per_family survives untouched", {
  # FamW has 4 members; n_per_family = 3 -> 4 > 3, survives (not excluded).
  score_w  <- setNames(c(10, 9, 8, 7), c("W1", "W2", "W3", "W4"))
  family_w <- setNames(rep("FamOnly", 4), names(score_w))
  res <- select_parents_by_family(score_w, family_w, n_families = 1, n_per_family = 3,
                                  verbose = FALSE)
  expect_equal(res$excluded_groups, character(0))
  expect_setequal(res$selected, c("W1", "W2", "W3"))
})

test_that("select_parents_by_family(): named-vector n_per_family exclusion threshold falls back to rank_k", {
  # FamA (size 3) and FamB (size 3): with n_per_family a NAMED vector,
  # n_per_family itself is only defined for whichever families end up
  # CHOSEN (not known yet at exclusion time), so rank_k stands in as the
  # threshold. rank_k = 2 here -> both families (size 3 > 2) survive, even
  # though FamA's own specific quota (3) would, by itself, equal its size.
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, rank_k = 2,
                                  n_per_family = c(FamA = 3L, FamB = 1L), verbose = FALSE)
  expect_equal(res$excluded_groups, character(0))
  expect_equal(sum(res$by_family$family == "FamA"), 3L)
  expect_equal(sum(res$by_family$family == "FamB"), 1L)
})

# ==============================================================================
# REML vs. ANOVA variance-component estimation
# ==============================================================================

test_that(".estimate_family_variance_components_reml(): returns NULL when lme4 is not installed", {
  skip_if(requireNamespace("lme4", quietly = TRUE),
         "lme4 is installed -- this covers the not-installed branch only")
  expect_null(HapBlockR:::.estimate_family_variance_components_reml(.fs_score, .fs_family))
})

test_that(".estimate_family_variance_components(): dispatcher falls back from reml to anova (with message) when lme4 is unavailable", {
  skip_if(requireNamespace("lme4", quietly = TRUE),
         "lme4 is installed -- this covers the not-installed fallback only")
  expect_message(
    vc <- HapBlockR:::.estimate_family_variance_components(.fs_score, .fs_family,
                                                            method = "reml", verbose = TRUE),
    "falling back"
  )
  expect_equal(vc$method_used, "anova")
  expect_true(vc$ok)
})

test_that(".estimate_family_variance_components(): REML matches the ANOVA method-of-moments estimate exactly on this BALANCED fixture (known equivalence for balanced one-way designs)", {
  skip_if_not_installed("lme4")
  vc <- HapBlockR:::.estimate_family_variance_components(.fs_score, .fs_family,
                                                          method = "reml", verbose = FALSE)
  expect_true(vc$ok)
  expect_equal(vc$method_used, "reml")
  # Hand-computed ANOVA values from the balanced fixture (see the
  # ANOVA-estimator test above): sigma2 = 1, tau2 = 44/3. REML and ANOVA
  # method-of-moments variance-component estimates coincide exactly for a
  # BALANCED one-way random-effects design (Searle, Casella & McCulloch
  # 1992) -- this fixture (4 families x 3 members each) is balanced.
  expect_equal(vc$sigma2, 1, tolerance = 1e-5)
  expect_equal(vc$tau2, 44 / 3, tolerance = 1e-5)
})

test_that("select_parents_by_family(): variance_method is echoed, and variance_method_used reports which estimator actually ran", {
  res_anova <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                        variance_method = "anova", verbose = FALSE)
  expect_equal(res_anova$variance_method, "anova")
  expect_equal(res_anova$variance_method_used, "anova")
  # Regression check for the return-list attribute-capture fix: these
  # fields must be populated (not NULL/lost across the fam_tab reorder).
  expect_false(is.null(res_anova$tau2))
  expect_false(is.null(res_anova$sigma2))
  expect_equal(res_anova$tau2, 44 / 3, tolerance = 1e-8)
  expect_equal(res_anova$sigma2, 1, tolerance = 1e-8)
})

# ==============================================================================
# Finite-population order-statistics bias correction
# ==============================================================================

test_that(".expected_order_stat_std_normal(): matches known exact closed-form values for small n", {
  # E[max of 2 iid N(0,1)] = 1/sqrt(pi) (David & Nagaraja, Order
  # Statistics). By symmetry of the standard normal, E[min of 2] = -that.
  expect_equal(HapBlockR:::.expected_order_stat_std_normal(2, 2), 1 / sqrt(pi), tolerance = 1e-6)
  expect_equal(HapBlockR:::.expected_order_stat_std_normal(2, 1), -1 / sqrt(pi), tolerance = 1e-6)
  # E[max of 3 iid N(0,1)] = 3/(2*sqrt(pi)) (same reference).
  expect_equal(HapBlockR:::.expected_order_stat_std_normal(3, 3), 3 / (2 * sqrt(pi)), tolerance = 1e-6)
  # Sum of all n order-statistic expectations must be exactly 0 (they are
  # just a relabelling of n draws each with mean 0).
  s3 <- sum(vapply(1:3, function(i) HapBlockR:::.expected_order_stat_std_normal(3, i), numeric(1)))
  expect_equal(s3, 0, tolerance = 1e-6)
})

test_that(".expected_topk_mean_std_normal(): top-1 of n equals E[max]; k >= n returns exactly 0", {
  expect_equal(HapBlockR:::.expected_topk_mean_std_normal(2, 1), 1 / sqrt(pi), tolerance = 1e-6)
  expect_equal(HapBlockR:::.expected_topk_mean_std_normal(3, 1), 3 / (2 * sqrt(pi)), tolerance = 1e-6)
  # Taking every member is not a selection at all -- no bias possible.
  expect_equal(HapBlockR:::.expected_topk_mean_std_normal(5, 5), 0)
  expect_equal(HapBlockR:::.expected_topk_mean_std_normal(3, 4), 0)  # k > n clamped to k = n
})

test_that(".bias_correct_topk_mean(): matches hand-computed values from the exact order-statistic constants", {
  fam_tab <- data.frame(family = c("FamSmall", "FamBig"), topk_mean = c(9, 9),
                        n_members = c(2, 3), stringsAsFactors = FALSE)
  res <- HapBlockR:::.bias_correct_topk_mean(fam_tab, k = 1, sigma2 = 1)
  expect_equal(res$topk_mean_bias[res$family == "FamSmall"], 1 / sqrt(pi), tolerance = 1e-6)
  expect_equal(res$topk_mean_bias[res$family == "FamBig"], 3 / (2 * sqrt(pi)), tolerance = 1e-6)
  expect_equal(res$topk_mean_corrected, res$topk_mean - res$topk_mean_bias)
  # Same k, larger n_members -> strictly MORE selection-differential bias
  # (more candidates to cherry-pick the top-1 from), even though the raw
  # topk_mean tied at 9 for both families.
  expect_true(res$topk_mean_bias[res$family == "FamBig"] >
             res$topk_mean_bias[res$family == "FamSmall"])
  expect_true(res$topk_mean_corrected[res$family == "FamBig"] <
             res$topk_mean_corrected[res$family == "FamSmall"])
})

test_that(".bias_correct_topk_mean(): degenerate sigma2 (NA/<=0) leaves topk_mean unchanged", {
  fam_tab <- data.frame(family = "X", topk_mean = 5, n_members = 3, stringsAsFactors = FALSE)
  for (bad_sigma2 in list(NA_real_, 0, -1)) {
    res <- HapBlockR:::.bias_correct_topk_mean(fam_tab, k = 1, sigma2 = bad_sigma2)
    expect_equal(res$topk_mean_bias, 0)
    expect_equal(res$topk_mean_corrected, 5)
  }
})

test_that(".bias_correct_topk_mean(): k >= n_members gives exactly zero bias (taking the whole family)", {
  fam_tab <- data.frame(family = "Y", topk_mean = 7, n_members = 2, stringsAsFactors = FALSE)
  res <- HapBlockR:::.bias_correct_topk_mean(fam_tab, k = 2, sigma2 = 1)
  expect_equal(res$topk_mean_bias, 0, tolerance = 1e-8)
  expect_equal(res$topk_mean_corrected, 7, tolerance = 1e-8)
})

test_that("select_parents_by_family(): bias_correction is echoed, and 'none' leaves topk_mean_bias at 0 / topk_mean_corrected == topk_mean", {
  res_none <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                       bias_correction = "none", variance_method = "anova",
                                       verbose = FALSE)
  expect_equal(res_none$bias_correction, "none")
  expect_true(all(res_none$family_ranking$topk_mean_bias == 0))
  expect_equal(res_none$family_ranking$topk_mean_corrected, res_none$family_ranking$topk_mean)

  res_oc <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                     bias_correction = "order_stats", variance_method = "anova",
                                     verbose = FALSE)
  expect_equal(res_oc$bias_correction, "order_stats")
  # All 4 families here have equal size (3) and equal k (rank_k defaults
  # to n_per_family = 2) -> the bias is the SAME positive constant for
  # every family (sigma2 = 1 > 0, k < n), so it changes no ranking even
  # though it shifts every corrected value down.
  expect_true(all(res_oc$family_ranking$topk_mean_bias > 0))
  expect_equal(res_oc$family_ranking$topk_mean_corrected,
              res_oc$family_ranking$topk_mean - res_oc$family_ranking$topk_mean_bias)
})

# ==============================================================================
# Genomic-relationship-informed family shrinkage
# ==============================================================================

test_that(".aggregate_relationship_to_family(): matches hand-computed group-mean-block values", {
  G <- matrix(c(1, 0.5, 0.2,
               0.5, 1, 0.3,
               0.2, 0.3, 1), nrow = 3, byrow = TRUE,
             dimnames = list(c("a1", "a2", "b1"), c("a1", "a2", "b1")))
  family <- setNames(c("FamA", "FamA", "FamB"), c("a1", "a2", "b1"))
  G_fam <- HapBlockR:::.aggregate_relationship_to_family(G, family)
  # FamA-FamA block = mean(1, 0.5, 0.5, 1) = 0.75
  expect_equal(G_fam["FamA", "FamA"], 0.75, tolerance = 1e-8)
  # FamA-FamB block = mean(0.2, 0.3) = 0.25
  expect_equal(G_fam["FamA", "FamB"], 0.25, tolerance = 1e-8)
  expect_equal(G_fam["FamB", "FamA"], 0.25, tolerance = 1e-8)
  # FamB-FamB block = mean(1) = 1
  expect_equal(G_fam["FamB", "FamB"], 1, tolerance = 1e-8)
})

test_that(".aggregate_relationship_to_family(): returns NULL when G does not cover every individual", {
  G <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("a1", "a2"), c("a1", "a2")))
  family <- setNames(c("FamA", "FamA", "FamB"), c("a1", "a2", "b1"))  # b1 missing from G
  expect_null(HapBlockR:::.aggregate_relationship_to_family(G, family))
})

test_that(".rank_families_gblup(): matches a hand-solved Henderson mixed-model-equations system (G_fam = I)", {
  # a = 2 groups: FamA (n=3, y=10), FamB (n=2, y=6), tau2 = 4, sigma2 = 2,
  # G_fam = I (mutually unrelated groups). Hand-solved (3x3) MME system:
  #   mu_hat = 234/29, a_A = 48/29, a_B = -48/29
  #   rank_score_A = mu_hat + a_A = 282/29 ~= 9.7241379
  #   rank_score_B = mu_hat + a_B = 186/29 ~= 6.4137931
  fam_tab <- data.frame(family = c("FamA", "FamB"), n_members = c(3, 2),
                        topk_mean_corrected = c(10, 6), stringsAsFactors = FALSE)
  G_fam <- diag(2); dimnames(G_fam) <- list(c("FamA", "FamB"), c("FamA", "FamB"))

  res <- HapBlockR:::.rank_families_gblup(fam_tab, tau2 = 4, sigma2 = 2, G_fam = G_fam)
  expect_equal(res$rank_score[res$family == "FamA"], 282 / 29, tolerance = 1e-6)
  expect_equal(res$rank_score[res$family == "FamB"], 186 / 29, tolerance = 1e-6)
})

test_that(".rank_families_gblup(): returns NULL when G_fam does not cover every group", {
  fam_tab <- data.frame(family = c("FamA", "FamB"), n_members = c(3, 2),
                        topk_mean_corrected = c(10, 6), stringsAsFactors = FALSE)
  G_fam <- matrix(1, 1, 1, dimnames = list("FamA", "FamA"))  # FamB missing
  expect_null(HapBlockR:::.rank_families_gblup(fam_tab, tau2 = 4, sigma2 = 2, G_fam = G_fam))
})

test_that("select_parents_by_family(): use_family_relationship = TRUE with full G coverage triggers relationship_informed = TRUE", {
  G_id <- diag(12); dimnames(G_id) <- list(names(.fs_score), names(.fs_score))
  G_id["A1", "A2"] <- G_id["A2", "A1"] <- 0.3  # a little real structure

  res_rel <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                      variance_method = "anova", G = G_id,
                                      use_family_relationship = TRUE, verbose = FALSE)
  expect_true(res_rel$relationship_informed)
  # Under GBLUP, no single scalar shrinkage weight applies.
  expect_true(all(is.na(res_rel$family_ranking$shrinkage_weight)))

  res_norel <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                        variance_method = "anova", G = G_id,
                                        use_family_relationship = FALSE, verbose = FALSE)
  expect_false(res_norel$relationship_informed)
  expect_false(res_norel$use_family_relationship)
  expect_true(all(!is.na(res_norel$family_ranking$shrinkage_weight)))
})

test_that("select_parents_by_family(): use_family_relationship = TRUE without G falls back to i.i.d. shrinkage (relationship_informed = FALSE)", {
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                                  variance_method = "anova", use_family_relationship = TRUE,
                                  verbose = FALSE)
  expect_false(res$relationship_informed)
  expect_true(res$use_family_relationship)  # echoes the argument as requested...
  expect_true(all(!is.na(res$family_ranking$shrinkage_weight)))  # ...but no G -> plain w_f used
})

test_that("select_parents_by_family(): use_family_relationship must be a single logical", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                             use_family_relationship = c(TRUE, FALSE)),
    "use_family_relationship"
  )
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 2, n_per_family = 2,
                             use_family_relationship = "yes"),
    "use_family_relationship"
  )
})

# ==============================================================================
# family_select_mode: argument requirements
# ==============================================================================

test_that("select_parents_by_family(): family_select_mode = 'count' requires n_per_family", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1),
    "n_per_family is required"
  )
})

test_that("select_parents_by_family(): family_select_mode = 'percentage' requires a valid pct_per_family", {
  base_args <- list(score = .fs_score, family = .fs_family, n_families = 1,
                    family_select_mode = "percentage")
  expect_error(do.call(select_parents_by_family, base_args), "pct_per_family")
  expect_error(do.call(select_parents_by_family, c(base_args, list(pct_per_family = 0))),
              "pct_per_family")
  expect_error(do.call(select_parents_by_family, c(base_args, list(pct_per_family = 101))),
              "pct_per_family")
  expect_error(do.call(select_parents_by_family, c(base_args, list(pct_per_family = c(10, 20)))),
              "pct_per_family")
})

test_that("select_parents_by_family(): family_select_mode = 'sd_threshold' requires a valid sd_threshold", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "sd_threshold"),
    "sd_threshold"
  )
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "sd_threshold", sd_threshold = c(1, 2)),
    "sd_threshold"
  )
})

test_that("select_parents_by_family(): family_select_mode = 'check_relative' requires check_margin_pct and exactly one of check_id/check_value", {
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "check_relative",
                             check_id = "A1"),
    "check_margin_pct"
  )
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "check_relative",
                             check_margin_pct = 10),
    "Exactly one of check_id or check_value"
  )
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "check_relative",
                             check_margin_pct = 10, check_id = "A1", check_value = 10),
    "Exactly one of check_id or check_value"
  )
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "check_relative",
                             check_margin_pct = 10, check_id = "NotAScoredLine"),
    "check_id"
  )
})

test_that("select_parents_by_family(): a leftover argument from an inactive family_select_mode is reported, not silently ignored", {
  expect_message(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "sd_threshold", sd_threshold = -10,
                             n_per_family = 2, verbose = TRUE),
    "ignoring supplied argument.*n_per_family"
  )
  expect_message(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "count", n_per_family = 2,
                             pct_per_family = 50, verbose = TRUE),
    "ignoring supplied argument.*pct_per_family"
  )
  # The mode's OWN argument(s) never trigger this message (checked via
  # capture_messages() + grepl rather than expect_no_message(), since a
  # REML-unavailable fallback message may legitimately also fire here and
  # would otherwise make this assertion fragile across environments).
  msgs <- testthat::capture_messages(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1, n_per_family = 2,
                             variance_method = "anova", verbose = TRUE)
  )
  expect_false(any(grepl("ignoring supplied argument", msgs)))
})

test_that("select_parents_by_family(): a leftover diversity-adjustment argument is reported when ensure_haplotype_diversity = FALSE", {
  vmat <- matrix(0, nrow = 12, ncol = 1,
                dimnames = list(names(.fs_score), "blk1"))
  expect_message(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1, n_per_family = 2,
                             ensure_haplotype_diversity = FALSE, value_matrix = vmat,
                             verbose = TRUE),
    "ignoring supplied argument.*value_matrix"
  )
})

# ==============================================================================
# family_select_mode = "percentage"
# ==============================================================================

test_that("select_parents_by_family(): percentage mode selects ceiling(pct/100 * family_size) per chosen family", {
  # Every family in .fs_family has exactly 3 members; pct_per_family = 50 ->
  # ceiling(0.5*3) = 2 per family -> 3 > 2, none excluded.
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2,
                                  family_select_mode = "percentage", pct_per_family = 50,
                                  variance_method = "anova", verbose = FALSE)
  expect_equal(res$family_select_mode, "percentage")
  expect_equal(res$pct_per_family, 50)
  expect_equal(res$excluded_groups, character(0))
  expect_setequal(res$selected, c("A1", "A2", "B1", "B2"))
})

test_that("select_parents_by_family(): percentage mode's ceiling rounding can inflate the effective share for small families", {
  # pct_per_family = 34% of FamA's 3 members -> ceiling(1.02) = 2, an
  # effective ~67%, not the nominal 34% -- documented, unavoidable rounding
  # behaviour for small family sizes.
  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                                  family_select_mode = "percentage", pct_per_family = 34,
                                  variance_method = "anova", verbose = FALSE)
  expect_setequal(res$selected, c("A1", "A2"))
})

test_that("select_parents_by_family(): percentage mode excludes a family when its own resolved take would equal its full size", {
  # pct_per_family = 100 -> ceiling(1.0*3) = 3 == family size for EVERY
  # family in this fixture -> nothing survives (mirrors the count-mode
  # 100%-quota degenerate case).
  expect_error(
    select_parents_by_family(.fs_score, .fs_family, n_families = 1,
                             family_select_mode = "percentage", pct_per_family = 100,
                             verbose = FALSE),
    "quota's worth"
  )
})

test_that("select_parents_by_family(): percentage mode reports excluded_groups and messages (unless verbose = FALSE), while a larger surviving family is still selected from", {
  # FamSmall (3 members) and FamBig (5 members). pct_per_family = 75 ->
  # ceiling(0.75*3) = 3 == FamSmall's size (excluded); ceiling(0.75*5) = 4
  # < FamBig's size (survives) -- so only SOME families are excluded here,
  # unlike the 100%/90% cases above where every family in the shared
  # fixture is small enough to be excluded and the call errors instead.
  score_mix  <- setNames(c(10, 9, 8, 7, 6, 5, 4, 3),
                        c("S1", "S2", "S3", "B1", "B2", "B3", "B4", "B5"))
  family_mix <- setNames(c(rep("FamSmall", 3), rep("FamBig", 5)), names(score_mix))

  expect_message(
    res <- select_parents_by_family(score_mix, family_mix, n_families = 1,
                                    family_select_mode = "percentage", pct_per_family = 75,
                                    verbose = TRUE),
    "Excluding"
  )
  expect_equal(res$excluded_groups, "FamSmall")
  expect_setequal(unique(res$by_family$family), "FamBig")
  expect_setequal(res$selected, c("B1", "B2", "B3", "B4"))  # ceiling(0.75*5) = 4
})

test_that("select_parents_by_family(): at a realistic percentage (10%), only single-member families are excluded -- 2/3/10-member families are untouched", {
  # Exclusion under 'percentage' reduces to pct_per_family > 100*(n-1)/n
  # for a group of size n. At pct_per_family = 10, that threshold is only
  # crossed at n = 1 (threshold 0, so ANY positive percentage excludes a
  # lone member) -- n = 2 (threshold 50), n = 3 (threshold 66.7), and
  # n = 10 (threshold 90) are all far from being crossed at 10%, unlike
  # count mode where a fixed absolute quota can easily exceed a small
  # family's size regardless of what fraction that represents.
  score_sizes <- setNames(
    c(100, 90, 89, 80, 79, 78, seq(70, 61, by = -1)),
    c("O1", "T1", "T2", "H1", "H2", "H3", paste0("Z", 1:10))
  )
  family_sizes <- setNames(
    c("FamOne", rep("FamTwo", 2), rep("FamThree", 3), rep("FamTen", 10)),
    names(score_sizes)
  )

  res <- select_parents_by_family(score_sizes, family_sizes, n_families = 3,
                                  family_select_mode = "percentage", pct_per_family = 10,
                                  verbose = FALSE)
  expect_equal(res$excluded_groups, "FamOne")
  expect_setequal(unique(res$by_family$family), c("FamTwo", "FamThree", "FamTen"))
  # ceiling(0.1*2)=1, ceiling(0.1*3)=1, ceiling(0.1*10)=1 -- one line each.
  expect_equal(nrow(res$by_family), 3L)
})

# ==============================================================================
# family_select_mode = "sd_threshold"
# ==============================================================================

test_that("select_parents_by_family(): sd_threshold mode selects every member clearing population_mean + k*sd, auto-backfilling zero-contributing families", {
  # FamP = {10, 8}, FamQ = {6, 4}. Population = {10,8,6,4}; mean = 7 exactly
  # -> with sd_threshold = 0, the cutoff is exactly 7 regardless of sd,
  # avoiding any floating-point sd computation in the hand-check. FamP:
  # both 10 and 8 clear 7. FamQ: neither 6 nor 4 clears 7 -> contributes 0.
  score_pq  <- setNames(c(10, 8, 6, 4), c("P1", "P2", "Q1", "Q2"))
  family_pq <- setNames(c("FamP", "FamP", "FamQ", "FamQ"), names(score_pq))

  expect_warning(
    res <- select_parents_by_family(score_pq, family_pq, n_families = 2,
                                    family_select_mode = "sd_threshold", sd_threshold = 0,
                                    variance_method = "anova", verbose = TRUE),
    "Only 1 family"
  )
  expect_equal(res$family_select_mode, "sd_threshold")
  expect_equal(res$sd_threshold, 0)
  expect_equal(res$n_families, 1L)
  expect_equal(res$zero_selected_groups, "FamQ")
  expect_equal(res$excluded_groups, character(0))  # no size-based rule for this mode
  expect_setequal(res$selected, c("P1", "P2"))
})

test_that("select_parents_by_family(): sd_threshold mode auto-backfills past a mid-ranked, zero-contributing family to a lower-ranked one that still qualifies", {
  # FamA = {A1=30, A2=26} (topk_mean, k=2 -> 28); FamB = {B1=16, B2=14}
  # (topk_mean 15); FamC = {C1=20, C2=2} (topk_mean 11). Ranked by RAW
  # topk_mean (family_rank_method = "topk_mean", no shrinkage) so ranking
  # is decoupled from per-member threshold-clearing: FamA > FamB > FamC by
  # rank, even though FamB (ranked #2) contributes ZERO qualifiers while
  # FamC (ranked #3) contributes one.
  #
  # Population = {30,26,16,14,20,2}; sum = 108; mean = 18 exactly (sd_threshold
  # = 0 again sidesteps needing to hand-compute sd). FamA: 30,26 both >= 18.
  # FamB: 16,14 both < 18 -> zero. FamC: 20 >= 18, 2 < 18 -> one qualifier.
  score_bf  <- setNames(c(30, 26, 16, 14, 20, 2),
                       c("A1", "A2", "B1", "B2", "C1", "C2"))
  family_bf <- setNames(c("FamA", "FamA", "FamB", "FamB", "FamC", "FamC"),
                       names(score_bf))

  res <- select_parents_by_family(score_bf, family_bf, n_families = 2, rank_k = 2,
                                  family_rank_method = "topk_mean",
                                  family_select_mode = "sd_threshold", sd_threshold = 0,
                                  verbose = FALSE)
  expect_equal(res$n_families, 2L)  # exactly enough found -- no shortfall
  expect_equal(res$zero_selected_groups, "FamB")
  expect_setequal(res$selected, c("A1", "A2", "C1"))
  expect_equal(res$by_family$family_rank[res$by_family$individual == "A1"], 1L)
  expect_equal(res$by_family$family_rank[res$by_family$individual == "C1"], 2L)
  expect_false("FamB" %in% res$by_family$family)
})

test_that("select_parents_by_family(): sd_threshold mode errors when nobody anywhere clears the bar", {
  score_pq  <- setNames(c(10, 8, 6, 4), c("P1", "P2", "Q1", "Q2"))
  family_pq <- setNames(c("FamP", "FamP", "FamQ", "FamQ"), names(score_pq))
  expect_error(
    select_parents_by_family(score_pq, family_pq, n_families = 1,
                             family_select_mode = "sd_threshold", sd_threshold = 100,
                             verbose = FALSE),
    "No family/group has any member clearing"
  )
})

# ==============================================================================
# family_select_mode = "check_relative"
# ==============================================================================

test_that("select_parents_by_family(): check_relative mode (check_id) keeps everyone at least check_margin_pct above the check, auto-backfilling a family with nobody clearing it", {
  # check = "CHK" (score 100, no family entry -> ineligible as a candidate,
  # used purely as a reference value). margin = 10% -> cutoff = 110.
  # FamL: L1=121 (>=110), L2=115 (>=110), L3=90 (<110) -> 2 qualifiers.
  # FamM: M1=95, M2=80 -> both < 110 -> 0 qualifiers.
  score_ck  <- setNames(c(121, 115, 90, 95, 80, 100),
                       c("L1", "L2", "L3", "M1", "M2", "CHK"))
  family_ck <- setNames(c("FamL", "FamL", "FamL", "FamM", "FamM"),
                       c("L1", "L2", "L3", "M1", "M2"))

  expect_warning(
    res <- select_parents_by_family(score_ck, family_ck, n_families = 2,
                                    family_select_mode = "check_relative",
                                    check_id = "CHK", check_margin_pct = 10,
                                    family_rank_method = "topk_mean",
                                    verbose = TRUE),
    "Only 1 family"
  )
  expect_equal(res$family_select_mode, "check_relative")
  expect_equal(res$check_id, "CHK")
  expect_equal(res$check_margin_pct, 10)
  expect_equal(res$zero_selected_groups, "FamM")
  expect_setequal(res$selected, c("L1", "L2"))
})

test_that("select_parents_by_family(): check_relative mode (check_value) gives the same result as an equal check_id lookup", {
  score_ck  <- setNames(c(121, 115, 90, 95, 80),
                       c("L1", "L2", "L3", "M1", "M2"))
  family_ck <- setNames(c("FamL", "FamL", "FamL", "FamM", "FamM"), names(score_ck))

  res <- select_parents_by_family(score_ck, family_ck, n_families = 1,
                                  family_select_mode = "check_relative",
                                  check_value = 100, check_margin_pct = 10,
                                  family_rank_method = "topk_mean", verbose = FALSE)
  expect_equal(res$check_value, 100)
  expect_null(res$check_id)
  expect_setequal(res$selected, c("L1", "L2"))
})

test_that("select_parents_by_family(): check_relative margin is an inclusive floor -- lines beating it by more are still included", {
  # check_value = 100. FamOnly = {L1=125 (+25%), L2=118 (+18%), L3=90 (-10%)}.
  score_one  <- setNames(c(125, 118, 90), c("L1", "L2", "L3"))
  family_one <- setNames(rep("FamOnly", 3), names(score_one))

  res10 <- select_parents_by_family(score_one, family_one, n_families = 1,
                                    family_select_mode = "check_relative",
                                    check_value = 100, check_margin_pct = 10,
                                    verbose = FALSE)
  expect_setequal(res10$selected, c("L1", "L2"))  # both +25% and +18% clear +10%

  res20 <- select_parents_by_family(score_one, family_one, n_families = 1,
                                    family_select_mode = "check_relative",
                                    check_value = 100, check_margin_pct = 20,
                                    verbose = FALSE)
  expect_setequal(res20$selected, "L1")  # only +25% clears +20%; +18% no longer does
})

# ==============================================================================
# family_select_mode x diversity/relatedness interactions (previously
# deferred -- these are the "untested new-mode/diversity interaction combos"
# flagged by the critical review as its own gap, closed here).
# ==============================================================================

test_that("select_parents_by_family(): family_select_mode = 'percentage' combines correctly with ensure_haplotype_diversity (cross-group block avoidance still respects each family's own take_f quota)", {
  # FamA/FamB (3 members each, pct_per_family = 34 -> take_f = ceiling(0.34*3)
  # = 2 per family -- matches the plain percentage-mode test's own math).
  # 4 dominant blocks -- FamA's full membership (all 3, since Pass 1 walks
  # every member, not just the 2 actually taken) claims blk1/blk2/blk3,
  # leaving blk4 free for FamB's best member (B1) to avoid a collision, but
  # FamB's #2 pick (B2) has nowhere left to go (blk1, already claimed) --
  # exactly the "some collision unavoidable" pattern already covered for
  # count/quota = 1, now confirmed to hold for a multi-member percentage
  # quota too.
  vmat <- matrix(0, nrow = 12, ncol = 4,
                dimnames = list(names(.fs_score), c("blk1", "blk2", "blk3", "blk4")))
  vmat["A1", "blk1"] <- 10; vmat["A2", "blk2"] <- 10; vmat["A3", "blk3"] <- 10
  vmat["B1", "blk4"] <- 10; vmat["B2", "blk1"] <- 10; vmat["B3", "blk2"] <- 10

  res <- select_parents_by_family(.fs_score, .fs_family, n_families = 2,
                                  family_select_mode = "percentage", pct_per_family = 34,
                                  ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                                  diversity_method = "dominant_block", verbose = FALSE)
  expect_equal(res$family_select_mode, "percentage")
  expect_equal(res$pct_per_family, 34)
  expect_equal(res$by_family$individual[res$by_family$family == "FamA"], c("A1", "A2"))
  expect_equal(res$by_family$individual[res$by_family$family == "FamB"], c("B1", "B2"))
  expect_equal(res$by_family$collision[res$by_family$individual == "A1"], FALSE)
  expect_equal(res$by_family$collision[res$by_family$individual == "A2"], FALSE)
  expect_equal(res$by_family$collision[res$by_family$individual == "B1"], FALSE)
  expect_equal(res$by_family$collision[res$by_family$individual == "B2"], TRUE)
})

test_that("select_parents_by_family(): family_select_mode = 'sd_threshold' + ensure_haplotype_diversity reorders (but never drops) qualifying members", {
  # Threshold modes take EVERY qualifying member -- diversity mode has no
  # fixed per-family quota left to trim, so its only possible effect here
  # is on ORDER (rank_within_family / which member is listed first), never
  # on the SET actually selected. sd_threshold = -10 is deliberately far
  # below the population mean so all 3 members clear it trivially.
  score_p  <- setNames(c(12, 11, 9), c("P1", "P2", "P3"))
  family_p <- setNames(rep("FamP", 3), names(score_p))
  vmat <- matrix(0, nrow = 3, ncol = 2,
                dimnames = list(names(score_p), c("blk1", "blk2")))
  vmat["P1", "blk1"] <- 10; vmat["P1", "blk2"] <- 1
  vmat["P2", "blk1"] <- 10; vmat["P2", "blk2"] <- 1   # same dominant block as P1
  vmat["P3", "blk1"] <- 1;  vmat["P3", "blk2"] <- 10  # different dominant block

  res_div <- select_parents_by_family(score_p, family_p, n_families = 1,
                                      family_select_mode = "sd_threshold", sd_threshold = -10,
                                      ensure_haplotype_diversity = TRUE, value_matrix = vmat,
                                      diversity_method = "dominant_block", verbose = FALSE)
  expect_setequal(res_div$selected, c("P1", "P2", "P3"))  # nobody dropped
  expect_equal(res_div$by_family$individual, c("P1", "P3", "P2"))  # P3 pulled ahead of P2
  expect_equal(res_div$by_family$collision, c(FALSE, FALSE, TRUE))

  res_nodiv <- select_parents_by_family(score_p, family_p, n_families = 1,
                                        family_select_mode = "sd_threshold", sd_threshold = -10,
                                        verbose = FALSE)
  expect_equal(res_nodiv$by_family$individual, c("P1", "P2", "P3"))  # plain score order
  expect_null(res_nodiv$by_family$collision)  # column only exists when diversity mode is on
})

test_that("select_parents_by_family(): family_select_mode = 'check_relative' + use_family_relationship = TRUE combines G-informed shrinkage with threshold auto-backfill", {
  # FamX (20, 18) clearly clears check_value = 15; FamY (12, 10) and FamZ
  # (8, 6) both fall short -- exercising GBLUP-informed family ranking AND
  # the auto-backfill walk at the same time (both FamY and FamZ get
  # skipped as zero-contributing before the ranking is exhausted, and only
  # FamX is left).
  score_g  <- setNames(c(20, 18, 12, 10, 8, 6),
                      c("X1", "X2", "Y1", "Y2", "Z1", "Z2"))
  family_g <- setNames(c("FamX", "FamX", "FamY", "FamY", "FamZ", "FamZ"),
                      names(score_g))
  G_id <- diag(6); dimnames(G_id) <- list(names(score_g), names(score_g))
  G_id["X1", "X2"] <- G_id["X2", "X1"] <- 0.2  # a little real relatedness structure

  expect_warning(
    res <- select_parents_by_family(score_g, family_g, n_families = 2,
                                    family_select_mode = "check_relative",
                                    check_value = 15, check_margin_pct = 0,
                                    variance_method = "anova", G = G_id,
                                    use_family_relationship = TRUE, verbose = FALSE),
    "family/group.*contributed"
  )
  expect_true(res$relationship_informed)
  expect_true(all(is.na(res$family_ranking$shrinkage_weight)))  # GBLUP: no scalar w_f
  expect_equal(res$n_families, 1L)
  expect_setequal(res$zero_selected_groups, c("FamY", "FamZ"))
  expect_setequal(res$selected, c("X1", "X2"))
})
