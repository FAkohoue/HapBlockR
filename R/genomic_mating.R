# ==============================================================================
# genomic_mating.R
#
# Usefulness Criterion (UC) / genomic mating: score CANDIDATE CROSSES (pairs
# of parents), not just individual parents, by predicted mid-parent value
# PLUS predicted genetic variance of the cross (Schnell & Utz 1975; Bernardo
# 2003 "criterion II"; Zhong & Jannink 2007; genomic-mating formalisation in
# Akdemir et al. 2019 and Allier et al. 2019).
#
# This is the second of six planned parent-selection / mate-allocation
# strategy extensions to HapBlockR (the first was the optional coancestry
# penalty added to select_parents_ga()). It complements truncation_selection()
# and select_parents_ga() -- both of which choose a SET of parents based on
# their own merit/coverage -- by instead ranking CROSSES (pairs) by their
# predicted potential to produce superior transgressive-segregant progeny,
# which is a different and complementary breeding decision (which parents to
# combine, not just which parents to keep).
# ==============================================================================


# -- Internal: truncation-selection intensity for selecting the top
# proportion p of a cross's progeny, with an optional finite-progeny-size
# correction.
#
# Infinite-population / asymptotic case (n_progeny = NULL, the default): the
# classical Falconer & Mackay (1996, Ch. 11) formula i = phi(z)/p, z =
# qnorm(1-p). This assumes an effectively infinite, normally-distributed
# progeny population -- a good approximation for large progeny numbers, but
# it OVERSTATES the true expected selection intensity for small crosses
# (typical biparental cross sizes of 50-300 progeny), since it is the n -> Inf
# limit of the expected top-k order-statistic mean.
#
# Finite-population case (n_progeny supplied): rather than relying on an
# approximate closed-form small-sample correction (several exist in the
# literature and differ in their assumptions), this computes the TRUE
# finite-sample expected order-statistic mean directly by Monte Carlo: draw
# n_progeny iid N(0,1) values, take the mean of the top round(p*n_progeny) of
# them, repeat n_sim times, and average. This is exact up to Monte Carlo
# error (which shrinks with n_sim) rather than depending on any particular
# closed-form asymptotic-correction formula. The caller's global RNG state is
# saved and restored, so this does not disturb reproducibility elsewhere in
# a user's script even when `seed` is supplied.
.selection_intensity <- function(p, n_progeny = NULL, n_sim = 20000L,
                                 seed = NULL) {
  if (is.null(p) || length(p) != 1L || is.na(p) || p <= 0 || p > 1)
    stop("selected_proportion must be a single value in (0, 1].", call. = FALSE)
  if (p >= 1) return(0)

  if (is.null(n_progeny)) {
    z <- stats::qnorm(1 - p)
    return(stats::dnorm(z) / p)
  }

  if (!is.numeric(n_progeny) || length(n_progeny) != 1L || is.na(n_progeny) ||
      n_progeny < 2L)
    stop("n_progeny must be a single integer >= 2.", call. = FALSE)
  if (!is.numeric(n_sim) || length(n_sim) != 1L || n_sim < 100L)
    stop("n_sim must be a single integer >= 100.", call. = FALSE)
  k <- max(1L, round(p * n_progeny))
  if (k >= n_progeny) return(0)  # keeping (almost) everyone: no intensity

  if (!is.null(seed)) {
    have_seed <- exists(".Random.seed", envir = .GlobalEnv)
    old_seed  <- if (have_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (have_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  n_progeny <- as.integer(round(n_progeny))
  n_sim     <- as.integer(round(n_sim))
  top_means <- vapply(seq_len(n_sim), function(i) {
    x <- sort(stats::rnorm(n_progeny), decreasing = TRUE)
    mean(x[seq_len(k)])
  }, numeric(1L))
  mean(top_means)
}

.augment_uc_uncertainty <- function(
    out,
    gebv_se,
    gebv_reliability,
    phasing_reliability,
    require_phasing_reliability,
    n_progeny,
    downside_quantile,
    min_reliability
) {
  parent_ids <- unique(c(out$parent1, out$parent2))
  .validate_named_metric <- function(x, label, lower = 0, upper = Inf) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || is.null(names(x)))
      stop(label, " must be a named numeric vector.", call. = FALSE)
    missing_ids <- setdiff(parent_ids, names(x))
    if (length(missing_ids))
      stop(label, " is missing ", length(missing_ids), " parent(s).",
           call. = FALSE)
    if (any(!is.finite(x[parent_ids])) ||
        any(x[parent_ids] < lower | x[parent_ids] > upper))
      stop(label, " values must be finite and in [", lower, ", ", upper,
           "].", call. = FALSE)
    x
  }
  gebv_se <- .validate_named_metric(gebv_se, "gebv_se", 0, Inf)
  gebv_reliability <- .validate_named_metric(
    gebv_reliability, "gebv_reliability", 0, 1
  )
  phasing_reliability <- .validate_named_metric(
    phasing_reliability, "phasing_reliability", 0, 1
  )
  if (length(downside_quantile) != 1L || !is.finite(downside_quantile) ||
      downside_quantile <= 0 || downside_quantile >= 0.5)
    stop("downside_quantile must be a single number in (0, 0.5).",
         call. = FALSE)
  if (length(min_reliability) != 1L || !is.finite(min_reliability) ||
      min_reliability < 0 || min_reliability > 1)
    stop("min_reliability must be a single number in [0, 1].",
         call. = FALSE)

  out$predicted_sd <- sqrt(pmax(out$predicted_variance, 0))
  out$downside_quantile <- downside_quantile
  out$downside_value <- out$mid_parent_gebv +
    stats::qnorm(downside_quantile) * out$predicted_sd

  if (is.null(gebv_se)) {
    out$mid_parent_SE <- NA_real_
  } else {
    out$mid_parent_SE <- sqrt(
      gebv_se[out$parent1]^2 + gebv_se[out$parent2]^2
    ) / 2
  }
  progeny_mean_se <- if (is.null(n_progeny)) {
    rep(0, nrow(out))
  } else {
    out$predicted_sd / sqrt(as.integer(n_progeny))
  }
  out$UC_SE <- sqrt(out$mid_parent_SE^2 + progeny_mean_se^2)
  out$UC_lower_95 <- out$UC - stats::qnorm(0.975) * out$UC_SE
  out$UC_upper_95 <- out$UC + stats::qnorm(0.975) * out$UC_SE

  if (is.null(gebv_reliability)) {
    out$prediction_reliability <- NA_real_
  } else {
    out$prediction_reliability <- pmin(
      gebv_reliability[out$parent1],
      gebv_reliability[out$parent2]
    )
  }
  if (is.null(phasing_reliability)) {
    out$phasing_reliability <- NA_real_
  } else {
    out$phasing_reliability <- pmin(
      phasing_reliability[out$parent1],
      phasing_reliability[out$parent2]
    )
  }
  out$cross_reliability <- if (isTRUE(require_phasing_reliability)) {
    pmin(out$prediction_reliability, out$phasing_reliability)
  } else if (!is.null(phasing_reliability)) {
    pmin(out$prediction_reliability, out$phasing_reliability)
  } else {
    out$prediction_reliability
  }
  out$min_reliability <- min_reliability
  out$recommendation_eligible <- is.finite(out$UC) &
    is.finite(out$cross_reliability) &
    out$cross_reliability >= min_reliability
  out$eligibility_reason <- ifelse(
    !is.finite(out$UC),
    "cross_not_scored",
    ifelse(
      !is.finite(out$cross_reliability),
      ifelse(
        isTRUE(require_phasing_reliability) &
          !is.finite(out$phasing_reliability),
        "phasing_reliability_not_supplied",
        "prediction_reliability_not_supplied"
      ),
      ifelse(
        out$cross_reliability < min_reliability,
        "below_minimum_reliability",
        "eligible"
      )
    )
  )
  out
}

# -- Internal: SNP IDs of one block, in the SAME order extract_haplotypes()'s
# C++ backend uses when it builds hap strings (position-ascending within the
# block's CHR / [start_bp, end_bp] window). Re-sorted explicitly here rather
# than trusting incoming row order in snp_info -- decompose_block_effects()
# relies on an implicit intersect()-based order instead, which is fragile if
# snp_info isn't already POS-sorted; this helper is deliberately independent
# of that logic.
.block_snp_order <- function(snp_info, chr, start_bp, end_bp) {
  idx <- which(snp_info$CHR == chr & snp_info$POS >= start_bp &
                 snp_info$POS <= end_bp)
  idx <- idx[order(snp_info$POS[idx])]
  snp_info$SNP[idx]
}

# -- Internal: additive effect of ONE phased gamete allele string (one
# haplotype copy of a block, digits 0/1 per SNP in position order) given
# per-SNP additive effects. Written fresh rather than reusing
# decompose_block_effects(), which silently drops every block when given
# phased "h1|h2"-style strings (its fixed-width dosage-string length check
# never matches a single-gamete 0/1 string, and it NULLs the block out with
# no warning). This helper instead errors loudly on any mismatch: a stop()
# here is preferable to a silently wrong UC estimate.
.phased_allele_effect <- function(allele_str, snp_ids, snp_effects) {
  digits <- strsplit(allele_str, "", fixed = TRUE)[[1]]
  if (length(digits) != length(snp_ids))
    stop("Phased allele string length (", length(digits), ") does not ",
         "match the number of SNPs located for this block in snp_info (",
         length(snp_ids), "). This usually means `haplotypes`, `snp_info`, ",
         "and `snp_effects` were not all derived from the same genotype ",
         "data / block definitions used to fit the model.", call. = FALSE)
  dose <- suppressWarnings(as.integer(digits))
  if (anyNA(dose))
    stop("Non-numeric character found in phased allele string '",
         allele_str, "'.", call. = FALSE)
  eff <- snp_effects[snp_ids]
  eff[is.na(eff)] <- 0
  sum(dose * eff)
}

# -- Internal: block_independent-mode predicted mean/variance contribution of
# one block to one candidate cross, from the two parents' local GEBV at that
# block (vi, vj). Treats the block as a single bi-allelic locus contrasting
# the two parents' block values and applies the standard single-locus
# biparental segregation-variance formula (variance among F2/RIL-type
# segregants of a cross between two homozygous parents differing by
# (vi - vj) at a locus is proportional to ((vi - vj) / 2)^2). This assumes
# each parent is close to homozygous at the block (true by construction for
# inbred-line/RIL breeding programs; a coarser proxy for partially
# heterozygous outbred parents).
.block_contrib_independent <- function(vi, vj, segregation_factor) {
  c(mean = (vi + vj) / 2,
    var  = segregation_factor * ((vi - vj) / 2)^2)
}

# -- Internal: phased-mode predicted mean/variance contribution of one block
# to one candidate cross, by exhaustively enumerating the 4 equally-likely
# gamete-pair combinations from parent i's two haplotype-allele effects
# (eff_i1, eff_i2) and parent j's (eff_j1, eff_j2). This is exact for a
# single block under Mendelian segregation with no further recombination
# within the block (consistent with this package's premise that an LD block
# is a low-recombination unit). stats::var() divides by (n-1) = 3 for the 4
# combinations; since these are the full set of equally-likely outcomes (not
# a sample), the true population variance divides by n = 4, hence the 3/4
# correction below.
.block_contrib_phased <- function(eff_i1, eff_i2, eff_j1, eff_j2) {
  combos <- c(eff_i1 + eff_j1, eff_i1 + eff_j2, eff_i2 + eff_j1, eff_i2 + eff_j2)
  c(mean = mean(combos), var = stats::var(combos) * 3 / 4)
}

# -- Internal: "linked" variance mode support. -------------------------------
#
# .block_contrib_phased() above -- and variance_model = "phased" as a whole
# -- treats every target block as an independent locus: parent i's
# transmitted allele at each block is drawn 50/50 from {eff_i1, eff_i2}
# INDEPENDENTLY of what was drawn at every other block, and the total
# variance is just the sum of each block's own variance. That is only
# correct if blocks are far enough apart (or on different chromosomes) that
# recombination fully un-links them between generations; for blocks close
# together on the same chromosome, the alleles transmitted at nearby blocks
# are correlated (linked), which changes the true progeny variance -- this
# is exactly the "linkage-aware" gap flagged against the block_independent/
# phased models (variance_model = "simplemating" already closes this gap by
# wrapping SimpleMating::getUsefA(), a verified external multi-locus
# covariance solver, but that requires the SimpleMating package, which is
# GitHub-only, not on CRAN).
#
# variance_model = "linked" closes the same gap natively, without any new
# external dependency, via first-principles Monte Carlo simulation of
# meiosis rather than re-deriving a closed-form multi-locus covariance
# formula by hand (deliberately not attempted -- see the comment above
# .run_simplemating_uc() explaining why hand-deriving Lehermeier et al.
# (2017)'s covariance construction without a working R interpreter to check
# it against would be exactly the kind of risk this package's development
# has been careful to avoid). Monte Carlo simulation of the actual
# biological process is a fundamentally different, much lower-risk approach
# than re-deriving a published closed-form result: every building block
# below is simple, standard, and independently checkable --
# Haldane's mapping function (recombination fraction from genetic distance;
# Haldane 1919) and a single-step Markov-chain crossover walk along a
# chromosome (the same "at each step, either keep going or cross over with
# probability r" model used throughout multi-locus linkage simulation) --
# and the whole thing is hand-verified in tests against two EXACT closed-form
# boundary cases: blocks far enough apart that r -> 0.5 (must recover the
# existing "phased" mode's independent-sum result) and blocks at zero
# genetic distance so r = 0 (must recover an exact 4-combo enumeration
# treating the whole multi-block region as one fused super-locus). See
# test-genomic-mating.R for both.
#
# Scope: like "phased", this models the immediate F1-style segregation
# variance of a single cross (one meiosis per parent, gametes combined) --
# NOT a full multi-generation RIL/DH population variance (that is a
# different, harder quantity requiring a repeated-selfing/doubling
# simulation on top of this; "simplemating" is the option for that via
# its Type/Generation arguments). "linked" therefore does NOT use `type`/
# `generation` at all (unlike "simplemating") -- only `genetic_map`.

# Haldane (1919) mapping function: converts a genetic distance in
# CENTIMORGANS to a recombination fraction in [0, 0.5]. Standard, exact,
# no-interference model.
.haldane_r <- function(cM_dist) {
  0.5 * (1 - exp(-2 * abs(cM_dist) / 100))
}

# Maps each target block to a representative genetic-map position (mean cM
# of the SNPs from `snp_info` that fall inside the block's [start_bp, end_bp]
# window on its CHR and are also present in `genetic_map`), then sorts
# blocks by (CHR, cM). Blocks with zero matched SNPs are dropped (with the
# count reported by the caller).
#
# @param block_coords Data frame with (at least) block_id, CHR, start_bp,
#   end_bp -- one row per unique target block (e.g. the first row per
#   block_id from infer_block_haplotypes()'s `dip` output, which already
#   carries these columns).
# @return Data frame block_id, CHR, cM, sorted by (CHR, cM); one row per
#   block that had >= 1 matched SNP.
.block_genetic_positions <- function(block_coords, snp_info, genetic_map) {
  gm <- genetic_map[!duplicated(genetic_map$SNP), , drop = FALSE]
  gm_by_snp <- stats::setNames(gm$cM, gm$SNP)

  rows <- lapply(seq_len(nrow(block_coords)), function(i) {
    bid <- block_coords$block_id[i]
    ch  <- block_coords$CHR[i]
    sb  <- block_coords$start_bp[i]
    eb  <- block_coords$end_bp[i]
    snp_in_block <- snp_info$SNP[snp_info$CHR == ch &
                                   snp_info$POS >= sb & snp_info$POS <= eb]
    cm_vals <- gm_by_snp[intersect(snp_in_block, names(gm_by_snp))]
    if (!length(cm_vals)) return(NULL)
    data.frame(block_id = bid, CHR = as.character(ch),
              cM = mean(cm_vals, na.rm = TRUE), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  if (is.null(out) || !nrow(out))
    return(data.frame(block_id = character(), CHR = character(),
                      cM = numeric(), stringsAsFactors = FALSE))
  out[order(out$CHR, out$cM), , drop = FALSE]
}

# Recombination fraction between every pair of ADJACENT blocks in a
# (CHR, cM)-sorted block order: 0.5 (independent assortment) across a
# chromosome boundary, Haldane's mapping function otherwise. Length =
# length(block_order) - 1.
.adjacent_r <- function(block_pos) {
  n <- nrow(block_pos)
  if (n < 2L) return(numeric(0))
  r <- numeric(n - 1L)
  for (i in seq_len(n - 1L)) {
    if (block_pos$CHR[i] != block_pos$CHR[i + 1L]) {
      r[i] <- 0.5
    } else {
      r[i] <- .haldane_r(block_pos$cM[i + 1L] - block_pos$cM[i])
    }
  }
  r
}

# Orders target blocks by PHYSICAL position (CHR, start_bp) instead of a
# genetic map -- used for variance_model = "linked" when ld_matrix (not
# genetic_map) is supplied, since there is then no genetic distance to sort
# by. Unlike .block_genetic_positions(), no SNP/genetic_map matching is
# needed here (block_coords already carries CHR/start_bp directly), so no
# blocks are ever dropped by this step.
#
# @param block_coords Data frame with (at least) block_id, CHR, start_bp.
# @return Data frame block_id, CHR, start_bp, sorted by (CHR, start_bp).
.block_physical_positions <- function(block_coords) {
  out <- block_coords[, c("block_id", "CHR", "start_bp"), drop = FALSE]
  out$CHR <- as.character(out$CHR)
  out[order(out$CHR, out$start_bp), , drop = FALSE]
}

# Recombination-fraction PROXY between every pair of ADJACENT blocks in a
# (CHR, start_bp)-sorted (physical, not genetic) block order, derived from
# pairwise SNP r^2 in `ld_matrix` -- the "linked" mode's fallback when no
# genetic_map is supplied. Mirrors the same "ld_matrix as a recombination-
# fraction proxy (1 - LD)" convention already documented for
# variance_model = "simplemating"'s own ld_matrix argument, landing on the
# form c = 0.5 * (1 - LD) so it matches Haldane's mapping function's exact
# [0, 0.5] range and boundary behaviour above: LD = 1 (SNPs always co-
# inherited in this population) -> c = 0 (treated as fully linked, same as
# zero genetic distance); LD = 0 -> c = 0.5 (independent assortment, same as
# an unlinked/far-apart pair). This is a monotonic, boundary-consistent
# PROXY, not a validated genetic-distance estimator (LD reflects population-
# level historical recombination + drift + selection, not the true physical
# recombination fraction for THIS specific cross) -- prefer a real
# genetic_map when you have one; only fall back to this when you don't.
# Chromosome boundaries and pairs with no matched SNPs on either side both
# conservatively default to c = 0.5 (independent assortment), exactly as
# .adjacent_r() does for chromosome boundaries.
#
# @param block_pos Data frame block_id, CHR (+ any other columns, ignored) --
#   already in the desired (physical) block order, e.g. from
#   .block_physical_positions() or a keep-filtered subset of it.
# @param snp_ids_by_block Named list, block_id -> character vector of SNP IDs
#   physically inside that block (from snp_info, independent of ld_matrix
#   coverage -- intersected against ld_matrix's own dimnames below).
# @param ld_matrix Square, SNP-dimnamed LD (r^2, in [0, 1]) matrix, e.g. from
#   compute_r2() with dimnames set to the SNP IDs used.
# @return Numeric vector, length nrow(block_pos) - 1.
.adjacent_r_from_ld <- function(block_pos, snp_ids_by_block, ld_matrix) {
  n <- nrow(block_pos)
  if (n < 2L) return(numeric(0))
  ld_ids <- rownames(ld_matrix)
  r <- numeric(n - 1L)
  for (i in seq_len(n - 1L)) {
    if (block_pos$CHR[i] != block_pos$CHR[i + 1L]) {
      r[i] <- 0.5
      next
    }
    snp_a <- intersect(snp_ids_by_block[[block_pos$block_id[i]]], ld_ids)
    snp_b <- intersect(snp_ids_by_block[[block_pos$block_id[i + 1L]]], ld_ids)
    if (!length(snp_a) || !length(snp_b)) {
      r[i] <- 0.5  # no ld_matrix coverage spanning this gap -- conservative
      next
    }
    ld_mean <- mean(ld_matrix[snp_a, snp_b, drop = FALSE], na.rm = TRUE)
    if (!is.finite(ld_mean)) ld_mean <- 0
    r[i] <- 0.5 * (1 - ld_mean)
  }
  r
}

# Simulates n_sim independent gametes across the ordered blocks for ONE
# parent, given that parent's per-block effect vectors for hap1 and hap2 (in
# the SAME block order as r_adjacent), and returns each replicate's total
# transmitted effect (vectorised over n_sim via a matrix of block-to-block
# "switch" draws, not a per-replicate R-level loop, for speed since this
# runs once per candidate cross).
.simulate_gamete_totals <- function(eff1_vec, eff2_vec, r_adjacent, n_sim) {
  n_blk <- length(eff1_vec)
  if (n_blk == 0L) return(rep(0, n_sim))
  start_hap1 <- stats::runif(n_sim) < 0.5
  if (n_blk == 1L) {
    return(ifelse(start_hap1, eff1_vec[1L], eff2_vec[1L]))
  }
  # switches[s, k] = TRUE means a crossover occurs between block k and k+1
  # in replicate s (k = 1 .. n_blk-1).
  switches <- matrix(stats::runif(n_sim * (n_blk - 1L)) < rep(r_adjacent, each = n_sim),
                     nrow = n_sim, ncol = n_blk - 1L)
  # is_hap1[s, k]: which haplotype (TRUE = hap1) replicate s carries at
  # block k, built by cumulatively XOR-ing the switch indicators along the
  # block order starting from start_hap1.
  is_hap1 <- matrix(NA, nrow = n_sim, ncol = n_blk)
  is_hap1[, 1L] <- start_hap1
  for (k in 2:n_blk) {
    is_hap1[, k] <- xor(is_hap1[, k - 1L], switches[, k - 1L])
  }
  eff1_mat <- matrix(eff1_vec, nrow = n_sim, ncol = n_blk, byrow = TRUE)
  eff2_mat <- matrix(eff2_vec, nrow = n_sim, ncol = n_blk, byrow = TRUE)
  rowSums(ifelse(is_hap1, eff1_mat, eff2_mat))
}

# One candidate cross's predicted mean/variance under variance_model =
# "linked": simulate n_sim progeny (one gamete from each parent, summed),
# using the SAME RNG-save/restore discipline as .selection_intensity()'s
# finite-population correction so this does not disturb the caller's global
# RNG state.
.block_contrib_linked_mc <- function(eff_i1_vec, eff_i2_vec, eff_j1_vec, eff_j2_vec,
                                     r_adjacent, n_sim, seed = NULL) {
  if (!is.null(seed)) {
    have_seed <- exists(".Random.seed", envir = .GlobalEnv)
    old_seed  <- if (have_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (have_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  gi <- .simulate_gamete_totals(eff_i1_vec, eff_i2_vec, r_adjacent, n_sim)
  gj <- .simulate_gamete_totals(eff_j1_vec, eff_j2_vec, r_adjacent, n_sim)
  progeny <- gi + gj
  c(mean = mean(progeny), var = stats::var(progeny))
}

# -- Internal: "simplemating" variance mode. Wraps SimpleMating::getUsefA()
# (Peixoto et al. 2024, Resende-Lab/SimpleMating), which builds a proper
# multi-locus Mendelian-sampling covariance matrix per chromosome from
# genetic-map recombination fractions (Haldane) or an LD-matrix proxy,
# following Lehermeier et al. (2017) -- genuinely linkage-aware across ALL
# supplied SNPs, unlike this file's own block_independent/phased modes,
# which only ever sum contributions within/across independently-treated LD
# blocks. Deliberately NOT reimplemented natively: getUsefA()'s C++-backed
# multi-locus covariance construction was verified against the package's
# current (2025) source before wrapping, and re-deriving it by hand without
# being able to execute R here would be exactly the kind of risk this
# package's own development has been careful to avoid.
#
# Real scope restriction: getUsefA() requires Markers coded strictly 0/2
# (fully homozygous calls) -- it is built for DH/RIL founder material, not
# heterozygous outbred parents. Heterozygous (dose = 1) calls are handled
# via `het_to_na` rather than either (a) blindly erroring out on any real
# breeding program's genotypes, which realistically almost always carry SOME
# residual heterozygosity even in advanced RIL material, or (b) silently
# ROUNDING dose = 1 to 0 or 2, which would fabricate a specific homozygous
# allele call the data does not actually support and corrupt the result.
# `het_to_na = TRUE` (default) instead treats heterozygous cells as MISSING
# (NA) -- a real, defensible statistical choice (equivalent to "we don't
# have a confident homozygous call here"), not a fabricated one, and NA was
# already a supported value for this argument before `het_to_na` existed
# (SimpleMating::getUsefA() already receives and handles NA in Markers).
# `het_to_na = FALSE` restores the original strict behaviour (error) for
# callers who want to be stopped rather than have any conversion happen.
.run_simplemating_uc <- function(cross_pairs, geno_matrix, het_to_na,
                                 snp_effects, G, genetic_map, ld_matrix,
                                 selected_proportion, type, generation,
                                 n_threads, verbose) {
  if (!requireNamespace("SimpleMating", quietly = TRUE))
    stop("SimpleMating is required for variance_model = 'simplemating'. ",
         "Install with: remotes::install_github('Resende-Lab/SimpleMating')",
         call. = FALSE)
  # Same version-gap check as .run_simplemating_ocs() in R/ocs.R: getUsefA()
  # was added in SimpleMating 0.2.x (alongside GOCS()); older 0.1.x releases
  # only expose buildCrosses()/selectCrosses()/evoluteCandidates()/etc. and
  # do not export getUsefA() at all. Checked explicitly for a clear,
  # actionable error rather than letting the .Call()/argument-matching
  # failure surface as a confusing "not an exported object" error -- this
  # exact situation has already happened once in practice.
  if (!"getUsefA" %in% getNamespaceExports("SimpleMating"))
    stop("SimpleMating is installed (version ",
         as.character(utils::packageVersion("SimpleMating")), ") but does ",
         "not export getUsefA() -- this is an older release (pre-0.2.x) ",
         "that predates SimpleMating's getUsefA()/GOCS() API. Reinstall the ",
         "current version with: remotes::install_github(",
         "'Resende-Lab/SimpleMating', force = TRUE). If getUsefA() is still ",
         "missing after that, SimpleMating's API may have moved again -- ",
         "check https://github.com/Resende-Lab/SimpleMating for its current ",
         "exported functions.", call. = FALSE)
  if (is.null(geno_matrix))
    stop("variance_model = 'simplemating' requires geno_matrix (individuals ",
         "x SNPs, dosage-coded strictly 0/2/NA -- fully homozygous calls ",
         "only; SimpleMating::getUsefA() is built for DH/RIL founder ",
         "material, not heterozygous parents). Use variance_model = ",
         "'block_independent' or 'phased' for heterozygous parents.",
         call. = FALSE)
  if (is.null(G))
    stop("variance_model = 'simplemating' requires G (a relationship ",
         "matrix, dimnamed by individual ID): required by ",
         "SimpleMating::getUsefA()'s own API (its K argument), even though ",
         "it does not affect the variance calculation itself -- it is only ",
         "carried through for a downstream optimisation step.", call. = FALSE)
  if (is.null(genetic_map) && is.null(ld_matrix))
    stop("variance_model = 'simplemating' requires either genetic_map ",
         "(data frame: SNP, CHR, cM) or ld_matrix (a SNP x SNP LD matrix, ",
         "used as a recombination proxy when no genetic map is ",
         "available).", call. = FALSE)

  if (!is.matrix(geno_matrix)) geno_matrix <- as.matrix(geno_matrix)
  obs <- geno_matrix[!is.na(geno_matrix)]
  # Anything outside 0/1/2 is not a recognised dosage value at all (a real
  # data problem, e.g. a different encoding or stray sentinel/missing code)
  # -- always an error, regardless of het_to_na, since there is no
  # defensible automatic interpretation for it the way there is for 1.
  invalid <- obs[!(obs %in% c(0, 1, 2))]
  if (length(invalid))
    stop("geno_matrix must be dosage-coded 0/1/2/NA for variance_model = ",
         "'simplemating'. Found other value(s) (e.g. ",
         paste(utils::head(sort(unique(invalid)), 5L), collapse = ", "),
         ") that are not a recognised genotype dosage -- check geno_matrix's ",
         "encoding.", call. = FALSE)

  n_het <- sum(geno_matrix == 1, na.rm = TRUE)
  if (n_het > 0L) {
    if (isTRUE(het_to_na)) {
      n_obs <- length(obs)
      if (isTRUE(verbose))
        message("[usefulness_criterion] simplemating: ", n_het, " of ",
                n_obs, " non-missing genotype call(s) (",
                round(100 * n_het / n_obs, 2), "%) are heterozygous ",
                "(dosage = 1) -- SimpleMating::getUsefA() requires strictly ",
                "homozygous (0/2) calls, so these are treated as MISSING ",
                "(set to NA), not rounded to 0 or 2 (which would fabricate ",
                "a specific allele call the data does not support). NA is ",
                "handled by SimpleMating::getUsefA() the same as any other ",
                "missing call. Set het_to_na = FALSE to error instead of ",
                "converting.")
      geno_matrix[geno_matrix == 1] <- NA_real_
    } else {
      stop("geno_matrix contains ", n_het, " heterozygous (dosage = 1) ",
           "call(s) -- variance_model = 'simplemating' requires strictly ",
           "0/2/NA (fully homozygous DH/RIL-style calls); ",
           "SimpleMating::getUsefA() is built for DH/RIL founder material, ",
           "not heterozygous parents. Set het_to_na = TRUE (the default) ",
           "to automatically treat heterozygous calls as missing instead, ",
           "or use variance_model = 'block_independent' or 'phased' for ",
           "heterozygous parents.", call. = FALSE)
    }
  }

  if (is.null(names(snp_effects)))
    stop("snp_effects must be a named numeric vector (names = SNP IDs) ",
         "for variance_model = 'simplemating'.", call. = FALSE)

  snp_ids <- intersect(colnames(geno_matrix), names(snp_effects))
  if (!is.null(genetic_map)) {
    if (!all(c("SNP", "CHR", "cM") %in% names(genetic_map)))
      stop("genetic_map must have SNP, CHR, cM columns.", call. = FALSE)
    snp_ids <- intersect(snp_ids, genetic_map$SNP)
  }
  if (!is.null(ld_matrix))
    snp_ids <- intersect(snp_ids, rownames(ld_matrix))
  if (length(snp_ids) < 2L)
    stop("Fewer than 2 SNPs in common across geno_matrix, snp_effects, and ",
         "genetic_map/ld_matrix.", call. = FALSE)

  Markers <- geno_matrix[, snp_ids, drop = FALSE]
  addEff  <- snp_effects[snp_ids]

  MatePlan <- data.frame(Parent1 = cross_pairs$parent1,
                         Parent2 = cross_pairs$parent2,
                         stringsAsFactors = FALSE)

  # SimpleMating's own Map.In column ORDER differs depending on branch (a
  # quirk of its current source, verified directly rather than assumed):
  # genetic-map branch expects (CHR, cM, SNP); linkDes branch expects
  # (CHR, SNP). Handled here so HapBlockR users only ever supply one
  # consistent genetic_map format (SNP, CHR, cM columns, any order).
  if (!is.null(ld_matrix)) {
    if (!is.null(genetic_map)) {
      gm <- genetic_map[match(snp_ids, genetic_map$SNP), , drop = FALSE]
    } else {
      gm <- data.frame(SNP = snp_ids, CHR = "1", stringsAsFactors = FALSE)
    }
    Map.In  <- data.frame(CHR = gm$CHR, SNP = gm$SNP, stringsAsFactors = FALSE)
    linkDes <- ld_matrix[snp_ids, snp_ids, drop = FALSE]
  } else {
    gm <- genetic_map[match(snp_ids, genetic_map$SNP), , drop = FALSE]
    if (anyNA(gm$SNP))
      stop("Not every SNP in common between geno_matrix/snp_effects is ",
           "present in genetic_map.", call. = FALSE)
    Map.In  <- data.frame(CHR = gm$CHR, cM = gm$cM, SNP = gm$SNP,
                          stringsAsFactors = FALSE)
    linkDes <- NULL
  }

  cross_ids <- unique(c(cross_pairs$parent1, cross_pairs$parent2))
  missing_g <- setdiff(cross_ids, rownames(G))
  if (length(missing_g))
    stop(length(missing_g), " parent ID(s) in cross_pairs are missing from ",
         "G: ", paste(utils::head(missing_g, 10L), collapse = ", "),
         call. = FALSE)
  Kmat <- G[cross_ids, cross_ids, drop = FALSE]

  res <- tryCatch(
    SimpleMating::getUsefA(
      MatePlan = MatePlan, Markers = Markers, addEff = addEff, K = Kmat,
      Map.In = Map.In, linkDes = linkDes, propSel = selected_proportion,
      Type = type, Generation = generation, n_threads = n_threads,
      display_progress = isTRUE(verbose)
    ),
    error = function(e)
      stop("SimpleMating::getUsefA() failed: ", conditionMessage(e),
           "\nCheck that geno_matrix/snp_effects/genetic_map (or ld_matrix) ",
           "share consistent SNP IDs, and that geno_matrix is strictly ",
           "0/2/NA coded. If this looks like an argument-name mismatch, ",
           "SimpleMating's API may have changed since this wrapper was ",
           "written against its 2025 source -- check ?SimpleMating::getUsefA.",
           call. = FALSE)
  )

  out <- res[[1]]
  data.frame(
    parent1              = out$Parent1,
    parent2              = out$Parent2,
    mid_parent_gebv       = out$Mean,
    predicted_variance    = out$Variance,
    selection_intensity   = .selection_intensity(selected_proportion),
    UC                    = out$Usefulness,
    stringsAsFactors      = FALSE
  )
}


#' Usefulness Criterion (UC) / Genomic Mating: Rank Candidate Crosses
#'
#' Scores and ranks candidate two-parent crosses by the Usefulness Criterion
#' (Schnell & Utz 1975; Bernardo 2003; Zhong & Jannink 2007), i.e. by their
#' predicted potential to produce progeny in the top \code{selected_proportion}
#' of the cross -- not just by the parents' own individual merit. This
#' complements \code{\link{truncation_selection}} and
#' \code{\link{select_parents_ga}}: those two choose a SET of parents to
#' keep, while \code{usefulness_criterion()} instead ranks PAIRS of already-
#' chosen (or candidate) parents by how good a cross between them is expected
#' to be, which is the natural next breeding decision once a parent set has
#' been chosen (which crosses to actually make).
#'
#' @section The Usefulness Criterion:
#' For a cross between parents i and j, UC is defined as
#' \deqn{UC = \mu_{ij} + i_{sel} \sqrt{\sigma^2_{ij}}}
#' where \eqn{\mu_{ij}} is the predicted mid-parent value (mean of the
#' progeny distribution), \eqn{\sigma^2_{ij}} is the predicted genetic
#' variance of the cross's progeny, and \eqn{i_{sel}} is the standard
#' truncation-selection intensity for selecting the top
#' \code{selected_proportion} of that progeny (Falconer & Mackay 1996). UC
#' therefore favours crosses expected to *segregate* into superior progeny,
#' not just crosses between two already-good parents -- this is what lets it
#' identify transgressive-segregation potential that whole-genome GEBV alone
#' cannot see.
#'
#' @section Four variance-prediction modes:
#' \code{variance_model} lets you choose which assumptions to make about the
#' input data; this package does not force one over the other -- pick
#' whichever matches what you have:
#' \describe{
#'   \item{\code{"block_independent"}}{Works on UNPHASED data. Uses each
#'     parent's per-block local GEBV (\code{local_gebv}, e.g. from
#'     \code{\link{run_haplotype_prediction}}) and treats each block as a
#'     single bi-allelic locus contrasting the two parents' block values,
#'     via the standard biparental single-locus segregation-variance formula
#'     \eqn{segregation\_factor \times ((v_i - v_j)/2)^2}, summed across
#'     target blocks under an independent-assortment assumption. Assumes
#'     parents are close to homozygous at each block (exact for inbred-
#'     line/RIL programs; a coarser proxy otherwise). No extra dependency.}
#'   \item{\code{"phased"}}{Requires phased haplotypes. Reads each parent's
#'     two actual haplotype alleles per block via
#'     \code{\link{infer_block_haplotypes}}, computes each allele's own
#'     effect directly from per-SNP effects, and exactly enumerates the 4
#'     equally-likely gamete-pair combinations per block to get an exact
#'     within-block segregation mean/variance. Target blocks are still
#'     summed independently (no between-block linkage/recombination
#'     covariance term) -- this is a defensible block-level extension
#'     consistent with the package's "LD block = low-recombination unit"
#'     premise, but it is NOT a full multi-locus LD-aware variance model.
#'     No extra dependency.}
#'   \item{\code{"linked"}}{Same phased-haplotype input as \code{"phased"},
#'     plus EITHER a \code{genetic_map} OR an \code{ld_matrix} -- the
#'     linkage-aware upgrade of \code{"phased"} that closes its "target
#'     blocks summed independently" gap, without the external-package
#'     requirement of \code{"simplemating"}. \code{haplotypes} (phased) is
#'     required either way -- \code{ld_matrix} only relaxes the
#'     genetic-distance requirement, it does not substitute for phase
#'     information (population-level LD tells you nothing about which
#'     alleles a SPECIFIC individual's two chromosomes carry). Simulates
#'     \code{n_sim_linked} progeny by Monte Carlo: for each replicate, one
#'     gamete per parent is generated by a block-to-block crossover walk,
#'     with blocks ordered and recombination fractions between adjacent
#'     blocks obtained either (a) from \code{genetic_map} via Haldane's
#'     (1919) mapping function (preferred: a real, validated genetic
#'     distance), or (b) when no \code{genetic_map} is supplied, from
#'     \code{ld_matrix} instead: blocks are ordered by physical position
#'     (CHR, start_bp) and each adjacent pair's recombination fraction is
#'     PROXIED as \eqn{0.5 \times (1 - \overline{LD})} (mean pairwise SNP
#'     r\eqn{^2} between the two blocks' member SNPs) -- the same
#'     [0, 0.5] range and boundary behaviour as Haldane's mapping function,
#'     but a monotonic proxy rather than a validated genetic-distance
#'     estimator (LD reflects population history, not necessarily the true
#'     recombination fraction for a specific cross); prefer \code{genetic_map}
#'     when you have one. Independent assortment (\eqn{r = 0.5}) is used
#'     across chromosome boundaries either way. Same scope as
#'     \code{"phased"} -- single-cross F1-style segregation variance, NOT a
#'     multi-generation RIL/DH population variance (see
#'     \code{"simplemating"}'s \code{type}/\code{generation} for that
#'     distinct, harder quantity). No extra dependency (pure
#'     Monte Carlo, no external solver).}
#'   \item{\code{"simplemating"}}{The most rigorous option, and the one to
#'     prefer when your data fits its requirements. Wraps
#'     \code{SimpleMating::getUsefA()} (Peixoto et al. 2024,
#'     Resende-Lab/SimpleMating -- \code{Suggests}, not bundled), which
#'     builds a genuine multi-locus Mendelian-sampling covariance matrix per
#'     chromosome from a genetic map (Haldane-mapped recombination
#'     fractions) or an LD-matrix proxy when no map is available, following
#'     Lehermeier et al. (2017) -- linkage across ALL supplied SNPs, not
#'     just within/across independently-treated LD blocks like \code{"phased"}/
#'     \code{"linked"} above (which only model linkage BETWEEN target
#'     blocks, not within them). Real constraint: \code{SimpleMating::getUsefA()}
#'     requires \code{geno_matrix} coded strictly 0/2 (fully homozygous
#'     calls) -- built for DH/RIL founder material, not heterozygous outbred
#'     parents. By default (\code{het_to_na = TRUE}) heterozygous (dosage =
#'     1) calls are automatically treated as missing (set to \code{NA}) so
#'     this mode still runs on realistic 0/1/2-coded data rather than
#'     erroring outright -- see \code{het_to_na} below for the exact
#'     behaviour and how to opt back into strict error-on-heterozygous
#'     validation instead. Does not use
#'     \code{block_importance}/\code{local_gebv}/\code{haplotypes}/
#'     \code{n_progeny} at all -- see the parameter docs below for the
#'     arguments this mode actually needs (\code{geno_matrix},
#'     \code{het_to_na}, \code{G}, \code{genetic_map} or \code{ld_matrix},
#'     \code{type}, \code{generation}).}
#' }
#' @section Finite-population selection intensity:
#' The classical UC formula assumes an infinite, normally-distributed
#' progeny population. Real biparental crosses produce tens to a few hundred
#' progeny, and the TRUE expected selection intensity for a finite sample is
#' somewhat smaller than the asymptotic formula gives (order-statistics
#' theory: the expected mean of the top k of n draws converges to, but does
#' not equal, the top-p-quantile density ratio as n grows). Supplying
#' \code{n_progeny} switches to a Monte Carlo estimate of the exact
#' finite-sample intensity (see \code{n_sim}, \code{seed}) rather than an
#' approximate closed-form small-sample correction -- deliberately, since
#' several such closed-form corrections exist in the literature and they are
#' not all mutually consistent, whereas the finite-sample order-statistic
#' mean itself is estimable directly and exactly (up to Monte Carlo error) by
#' simulation. Leave \code{n_progeny = NULL} to use the standard asymptotic
#' formula (fine for large planned progeny numbers, and the pre-existing
#' default behaviour).
#'
#' In both variance modes, \code{predicted_variance} reflects segregation at the
#' target blocks in \code{block_importance} only, not the whole genome --
#' the same tractability trade-off \code{\link{select_parents_ga}} makes by
#' pre-filtering to a manageable set of top-ranked blocks. If
#' \code{block_importance} covers blocks explaining most of local-GEBV
#' variance (e.g. via \code{\link{select_top_blocks}}), this captures most of
#' the segregating variance that matters in practice, but is not exactly the
#' whole-genome quantity classical UC formulations assume. \code{mid_parent_gebv}
#' (the mean term), by contrast, always uses the full whole-genome \code{gebv}
#' you supply, since that is the best available estimate of the cross's mean.
#'
#' @section What this does and does not do:
#' This ranks pairs of parents for making a cross; it does not choose the
#' parent set itself (use \code{\link{truncation_selection}} or
#' \code{\link{select_parents_ga}} for that, optionally first), and it does
#' not allocate differential numbers of progeny/contributions across crosses
#' or constrain population-wide inbreeding the way true Optimal Contribution
#' Selection (OCS) does -- that is a separate, not-yet-implemented strategy
#' planned for this package. Treat \code{usefulness_criterion()} as a ranking
#' tool for choosing which of your candidate crosses to prioritise, not as a
#' full mate-allocation/OCS solver.
#'
#' @param parent_ids Character vector of candidate parent IDs to generate all
#'   pairwise crosses from via \code{\link{combn}}. Ignored if
#'   \code{cross_pairs} is supplied. Either this or \code{cross_pairs} is
#'   required.
#' @param cross_pairs Optional data frame or matrix with exactly two columns
#'   (parent1, parent2) giving a specific, possibly non-exhaustive, list of
#'   candidate crosses to score. Overrides \code{parent_ids} if both are
#'   given.
#' @param gebv Named numeric vector of whole-genome GEBV (one value per
#'   individual, names = individual IDs), used for the mid-parent mean term.
#'   Required.
#' @param selected_proportion Numeric in (0, 1]. Proportion of each cross's
#'   progeny you intend to keep, used to compute the selection intensity
#'   \eqn{i_{sel}} (default 0.1, i.e. keep the top 10\%). Smaller values give
#'   more weight to predicted variance (favour riskier, higher-upside
#'   crosses); values near 1 make UC converge to the mid-parent value alone.
#' @param n_progeny Optional integer, default \code{NULL}. If supplied,
#'   \eqn{i_{sel}} is computed for a FINITE progeny population of this size
#'   (by Monte Carlo -- see "Finite-population selection intensity" below)
#'   instead of the classical infinite-population asymptotic formula. Use
#'   this when your realistic cross size is small (tens to a few hundred
#'   progeny), where the asymptotic formula overstates \eqn{i_{sel}}.
#' @param n_sim Integer, default \code{20000L}. Number of Monte Carlo
#'   replicates used for the finite-population \eqn{i_{sel}} when
#'   \code{n_progeny} is supplied. Ignored otherwise. Larger values reduce
#'   Monte Carlo error at the cost of runtime.
#' @param seed Optional integer, default \code{NULL}. Seed for the
#'   finite-population Monte Carlo (ignored if \code{n_progeny} is
#'   \code{NULL}). The caller's global RNG state is saved and restored, so
#'   supplying this does not affect random-number generation elsewhere in
#'   your script.
#' @param variance_model Character, one of \code{"block_independent"}
#'   (default), \code{"phased"}, \code{"linked"}, or \code{"simplemating"}.
#'   See "Four variance-prediction modes" above.
#'   You choose which one to use based on what genotype data you have
#'   available; this package does not pick one for you.
#' @param block_importance Data frame with (at least) \code{block_id},
#'   \code{CHR}, \code{start_bp}, \code{end_bp} columns identifying which
#'   blocks to compute predicted variance over -- e.g.
#'   \code{run_haplotype_prediction()$block_importance}, optionally filtered
#'   first via \code{\link{select_top_blocks}} to a manageable, high-
#'   importance subset. Required for \code{"block_independent"}/
#'   \code{"phased"}/\code{"linked"}; unused (and not required) for
#'   \code{"simplemating"}.
#' @param block_ids Optional character vector to further restrict
#'   \code{block_importance} to specific \code{block_id} values before
#'   scoring (default \code{NULL} = use every block in \code{block_importance}
#'   as-is). Unused for \code{"simplemating"}.
#' @param local_gebv Numeric matrix of per-individual, per-block local GEBV
#'   (rows = individual IDs, columns = \code{block_id}), e.g.
#'   \code{run_haplotype_prediction()$local_gebv}. Required when
#'   \code{variance_model = "block_independent"}.
#' @param segregation_factor Numeric scalar, default \code{0.5}. Scales the
#'   single-locus biparental segregation-variance formula in
#'   \code{"block_independent"} mode (derived from standard F2/RIL additive
#'   segregation-variance theory: \eqn{0.5 \times ((v_i-v_j)/2)^2}). Exposed
#'   so you can adjust it for your population's mating design (e.g. a
#'   different generation of selfing/recombination) rather than being locked
#'   to the default assumption. Unused outside \code{"block_independent"}.
#' @param haplotypes Phased haplotype list as returned by
#'   \code{\link{extract_haplotypes}} (carrying its \code{block_info}
#'   attribute), with \code{phased = TRUE} blocks. Required when
#'   \code{variance_model} is \code{"phased"} or \code{"linked"}.
#' @param snp_info Data frame with (at least) \code{SNP}, \code{CHR},
#'   \code{POS} columns, matching what was used to build \code{haplotypes}
#'   and fit \code{snp_effects}. Required when \code{variance_model} is
#'   \code{"phased"} or \code{"linked"}. Unused for \code{"simplemating"}
#'   (use \code{genetic_map} instead, which carries the same SNP identity
#'   plus genetic, not physical, position).
#' @param snp_effects Named numeric vector of per-SNP additive marker
#'   effects (names = SNP IDs), e.g. from \code{\link{estimate_marker_effects}}
#'   or \code{\link{backsolve_snp_effects}}. Required when
#'   \code{variance_model} is \code{"phased"}, \code{"linked"}, or
#'   \code{"simplemating"}.
#' @param geno_matrix Numeric matrix, individuals x SNPs, dosage-coded
#'   0/1/2/NA. Required when \code{variance_model = "simplemating"}; unused
#'   otherwise. \code{SimpleMating::getUsefA()} itself requires strictly
#'   0/2/NA (fully homozygous DH/RIL-style calls) -- see \code{het_to_na} for
#'   how heterozygous (dosage = 1) calls are handled.
#' @param het_to_na Logical, default \code{TRUE}. \code{"simplemating"} only.
#'   If \code{TRUE} (default), heterozygous (dosage = 1) calls in
#'   \code{geno_matrix} are automatically treated as MISSING (set to
#'   \code{NA}) before calling \code{SimpleMating::getUsefA()}, which
#'   requires strictly homozygous 0/2/NA calls -- this is a defensible
#'   statistical choice (equivalent to "no confident homozygous call here"),
#'   NOT a fabricated one: heterozygous cells are never rounded to 0 or 2,
#'   which would invent a specific allele call the data does not support. If
#'   \code{FALSE}, \code{usefulness_criterion()} errors instead when any
#'   heterozygous call is found, leaving the decision to you rather than
#'   converting automatically. Unused otherwise.
#' @param G Dimnamed relationship matrix (row/column names = individual
#'   IDs), e.g. from \code{\link{compute_haplotype_grm}}. Required when
#'   \code{variance_model = "simplemating"} (passed through to
#'   \code{SimpleMating::getUsefA()}'s \code{K} argument -- it does not
#'   affect the variance calculation, only downstream optimisation
#'   metadata); unused otherwise.
#' @param genetic_map Data frame with \code{SNP}, \code{CHR}, \code{cM}
#'   columns (genetic, not physical, position). Preferred (real, validated
#'   distances) for \code{"linked"}; required unless \code{ld_matrix} is
#'   supplied instead. For \code{"simplemating"}, likewise required unless
#'   \code{ld_matrix} is supplied instead. Unused for
#'   \code{"block_independent"}/\code{"phased"}.
#' @param ld_matrix Optional square, SNP-dimnamed linkage-disequilibrium
#'   (r\eqn{^2}, in [0, 1]) matrix, e.g. from \code{\link{compute_r2}} with
#'   dimnames set to the SNP IDs used -- a recombination-fraction PROXY for
#'   when no \code{genetic_map} is available. Used by \code{"simplemating"}
#'   as \eqn{1 - LD} (passed through to \code{SimpleMating::getUsefA()}), and
#'   by \code{"linked"} as \eqn{0.5 \times (1 - \overline{LD})} between
#'   adjacent target blocks (mean pairwise SNP LD, blocks ordered by physical
#'   position when no \code{genetic_map} is given -- see "linked" in "Four
#'   variance-prediction modes" above for the full formula and its caveats).
#'   \strong{Does not relax the phased-\code{haplotypes} requirement for
#'   \code{"linked"}/\code{"phased"}} -- LD is a population-level statistic
#'   and cannot substitute for knowing which alleles a specific individual's
#'   two chromosomes actually carry. Unused for
#'   \code{"block_independent"}/\code{"phased"}.
#' @param type Character, one of \code{"RIL"} (default) or \code{"DH"}.
#'   \code{"simplemating"} only: the population type derived from each
#'   cross (recombinant inbred line vs. doubled haploid), which changes the
#'   Mendelian-sampling covariance formula. Unused otherwise -- including by
#'   \code{"linked"}, which models single-cross F1-style segregation
#'   variance rather than a multi-generation RIL/DH population variance
#'   (see "Four variance-prediction modes" above).
#' @param generation Integer, default \code{1L}. \code{"simplemating"} only:
#'   the generation at which DH lines are generated or RILs extracted
#'   (Lehermeier et al. 2017); values \code{>= 10} are treated as an
#'   effectively infinite generation. Unused otherwise.
#' @param n_threads Integer, default \code{1L}. \code{"simplemating"} only:
#'   threads used internally by \code{SimpleMating::getUsefA()}'s C++
#'   backend. Unused otherwise.
#' @param n_sim_linked Integer, default \code{2000L}. Number of Monte Carlo
#'   progeny simulated PER CANDIDATE CROSS for \code{variance_model =
#'   "linked"}. Ignored otherwise. Kept separate from \code{n_sim} (which is
#'   for the finite-population selection-intensity correction, computed
#'   ONCE regardless of cross count) because \code{"linked"} pays this cost
#'   once per candidate cross -- a large value here scales with the number
#'   of crosses being scored, unlike \code{n_sim}. \code{seed}, if supplied,
#'   is reused for both Monte Carlo procedures.
#' @param gebv_se Optional named numeric vector of GEBV standard errors.
#'   When supplied, the function propagates parental uncertainty to
#'   \code{mid_parent_SE} and the 95 percent UC interval.
#' @param gebv_reliability Optional named numeric vector in [0, 1], for
#'   example the \code{reliability} column returned by
#'   \code{\link{run_haplotype_prediction}}. Cross reliability is the
#'   conservative minimum of its two parental reliabilities.
#' @param phasing_reliability Optional named numeric vector in [0, 1].
#'   For \code{variance_model = "phased"} or \code{"linked"} -- both build
#'   progeny variance from phased haplotype blocks and depend equally on
#'   phasing accuracy -- each cross uses the conservative minimum of
#'   parental prediction and phasing reliability. Missing phasing
#'   reliability therefore cannot pass the recommendation gate in either
#'   mode.
#' @param downside_quantile Numeric in (0, 0.5). Progeny-distribution
#'   quantile reported as \code{downside_value}. Default \code{0.10}.
#' @param min_reliability Numeric in [0, 1]. A cross is marked
#'   \code{recommendation_eligible} only when its two-parent reliability
#'   meets this threshold. Missing reliability never passes the gate.
#'   Default \code{0.30}.
#' @param verbose Logical, default \code{TRUE}. Print progress and a summary
#'   of how many candidate crosses could not be scored (e.g. due to missing
#'   genotype or phase data at target blocks).
#'
#' @return A data frame, one row per candidate cross, sorted by descending
#'   UC, with columns:
#'   \describe{
#'     \item{\code{parent1}, \code{parent2}}{The two parent IDs.}
#'     \item{\code{mid_parent_gebv}}{Mean of the two parents' whole-genome
#'       \code{gebv}.}
#'     \item{\code{predicted_variance}}{Predicted genetic variance of the
#'       cross's progeny at the target blocks (see Details above for scope).
#'       \code{NA} if the cross could not be scored.}
#'     \item{\code{selection_intensity}}{The \eqn{i_{sel}} value used (same
#'       for every row, since it depends only on \code{selected_proportion}).}
#'     \item{\code{UC}}{\code{mid_parent_gebv + selection_intensity *
#'       sqrt(predicted_variance)}. \code{NA} if the cross could not be
#'       scored.}
#'     \item{\code{rank}}{Rank by descending UC (\code{NA} rows sort last).}
#'     \item{\code{downside_value}}{The requested lower progeny-distribution
#'       quantile, showing downside risk alongside expected selected gain.}
#'     \item{\code{mid_parent_SE}, \code{UC_SE}, \code{UC_lower_95},
#'       \code{UC_upper_95}}{Propagated uncertainty when \code{gebv_se} is
#'       supplied; the finite-progeny contribution is included when
#'       \code{n_progeny} is supplied.}
#'     \item{\code{prediction_reliability},
#'       \code{phasing_reliability}, \code{cross_reliability},
#'       \code{recommendation_eligible}, \code{eligibility_reason}}{Explicit
#'       reliability gate for promoting a ranked cross to a recommendation.}
#'   }
#'
#' @references
#' Schnell, F.W. & Utz, H.F. (1975). F1-Leistung und Elternwahl in der
#' Zuchtung von Selbstbefruchtern. \emph{Ber. Arbeitstagung Arbeitsgemeinschaft
#' Saatzuchtleiter}.
#'
#' Bernardo, R. (2003). Parental selection, number of breeding populations,
#' and size of each population in inbred development. \emph{Theoretical and
#' Applied Genetics}, 107, 1252-1256.
#'
#' Zhong, S. & Jannink, J.-L. (2007). Using quantitative trait loci results
#' to discriminate among crosses on the basis of their progeny mean and
#' variance. \emph{Genetics}, 177, 567-576.
#'
#' Haldane, J.B.S. (1919). The combination of linkage values and the
#' calculation of distances between the loci of linked factors.
#' \emph{Journal of Genetics}, 8, 299-309. (Mapping function used by
#' \code{variance_model = "linked"} to convert genetic distance in cM to
#' recombination fraction.)
#'
#' Falconer, D.S. & Mackay, T.F.C. (1996). \emph{Introduction to Quantitative
#' Genetics}, 4th ed. Longman.
#'
#' Akdemir, D., Beavis, W., Fritsche-Neto, R., Singh, A.K. & Isidro-Sanchez,
#' J. (2019). Multi-objective optimized genomic breeding strategies for
#' sustainable food improvement. \emph{Heredity}, 122, 672-683.
#'
#' Allier, A., Lehermeier, C., Charcosset, A., Moreau, L. & Teyssedre, S.
#' (2019). Improving short- and long-term genetic gain by accounting for
#' within-family variance in optimal cross-selection. \emph{Frontiers in
#' Genetics}, 10, 1006.
#'
#' @seealso \code{\link{truncation_selection}}, \code{\link{select_parents_ga}}
#'   for choosing the parent set itself; \code{\link{select_top_blocks}} for
#'   building a manageable \code{block_importance} subset.
#'
#' @export
usefulness_criterion <- function(
    parent_ids          = NULL,
    cross_pairs          = NULL,
    gebv,
    selected_proportion  = 0.1,
    n_progeny             = NULL,
    n_sim                 = 20000L,
    seed                  = NULL,
    variance_model        = c("block_independent", "phased", "linked", "simplemating"),
    block_importance      = NULL,
    block_ids             = NULL,
    local_gebv            = NULL,
    segregation_factor    = 0.5,
    haplotypes            = NULL,
    snp_info              = NULL,
    snp_effects           = NULL,
    geno_matrix            = NULL,
    het_to_na               = TRUE,
    G                       = NULL,
    genetic_map             = NULL,
    ld_matrix                = NULL,
    type                      = c("RIL", "DH"),
    generation                = 1L,
    n_threads                 = 1L,
    n_sim_linked              = 2000L,
    gebv_se                    = NULL,
    gebv_reliability           = NULL,
    phasing_reliability        = NULL,
    downside_quantile          = 0.10,
    min_reliability            = 0.30,
    verbose               = TRUE
) {
  variance_model <- match.arg(variance_model)
  type <- match.arg(type)

  if (is.null(cross_pairs)) {
    if (is.null(parent_ids) || length(unique(parent_ids)) < 2L)
      stop("Supply at least 2 unique parent_ids to generate pairwise crosses ",
           "from, or supply cross_pairs directly.", call. = FALSE)
    cross_pairs <- t(utils::combn(unique(parent_ids), 2L))
  } else {
    cross_pairs <- as.matrix(cross_pairs)
    if (ncol(cross_pairs) != 2L)
      stop("cross_pairs must have exactly 2 columns (parent1, parent2).",
           call. = FALSE)
  }
  cross_pairs <- as.data.frame(cross_pairs, stringsAsFactors = FALSE)
  names(cross_pairs) <- c("parent1", "parent2")
  cross_pairs$parent1 <- as.character(cross_pairs$parent1)
  cross_pairs$parent2 <- as.character(cross_pairs$parent2)

  if (missing(gebv) || is.null(names(gebv)))
    stop("gebv must be a named numeric vector (names = individual IDs).",
         call. = FALSE)
  missing_gebv <- setdiff(unique(c(cross_pairs$parent1, cross_pairs$parent2)),
                           names(gebv))
  if (length(missing_gebv))
    stop(length(missing_gebv), " parent ID(s) in the cross list have no ",
         "gebv entry: ",
         paste(utils::head(missing_gebv, 10L), collapse = ", "),
         call. = FALSE)

  if (variance_model == "simplemating") {
    out <- .run_simplemating_uc(cross_pairs, geno_matrix = geno_matrix,
                                het_to_na = het_to_na,
                                snp_effects = snp_effects, G = G,
                                genetic_map = genetic_map,
                                ld_matrix = ld_matrix,
                                selected_proportion = selected_proportion,
                                type = type, generation = generation,
                                n_threads = n_threads, verbose = verbose)
    ord <- order(-out$UC, na.last = TRUE)
    out <- out[ord, , drop = FALSE]
    out <- .augment_uc_uncertainty(
      out, gebv_se, gebv_reliability, phasing_reliability, FALSE,
      n_progeny,
      downside_quantile, min_reliability
    )
    out$rank <- seq_len(nrow(out))
    rownames(out) <- NULL
    return(out)
  }

  if (is.null(block_importance) || !is.data.frame(block_importance) ||
      !all(c("block_id", "CHR") %in% names(block_importance)))
    stop("block_importance must be a data frame with at least block_id and ",
         "CHR columns (e.g. run_haplotype_prediction()$block_importance, ",
         "optionally filtered via select_top_blocks()).", call. = FALSE)
  if (!is.null(block_ids))
    block_importance <- block_importance[block_importance$block_id %in% block_ids, , drop = FALSE]
  if (!nrow(block_importance))
    stop("No blocks remain in block_importance after filtering by block_ids.",
         call. = FALSE)

  i_sel <- .selection_intensity(selected_proportion, n_progeny = n_progeny,
                                n_sim = n_sim, seed = seed)
  if (isTRUE(verbose) && !is.null(n_progeny))
    message("[usefulness_criterion] Finite-population selection intensity ",
            "(n_progeny = ", n_progeny, "): i_sel = ", round(i_sel, 4),
            " (Monte Carlo, n_sim = ", n_sim, ").")

  if (variance_model == "block_independent") {
    if (is.null(local_gebv))
      stop("variance_model = 'block_independent' requires local_gebv ",
           "(e.g. run_haplotype_prediction()$local_gebv).", call. = FALSE)
    bids <- intersect(block_importance$block_id, colnames(local_gebv))
    if (!length(bids))
      stop("None of block_importance's block_id values are columns of ",
           "local_gebv.", call. = FALSE)
    if (isTRUE(verbose))
      message("[usefulness_criterion] block_independent mode: ", length(bids),
              " target block(s), ", nrow(cross_pairs), " candidate cross(es).")

    res <- lapply(seq_len(nrow(cross_pairs)), function(r) {
      p1 <- cross_pairs$parent1[r]; p2 <- cross_pairs$parent2[r]
      if (!(p1 %in% rownames(local_gebv)) || !(p2 %in% rownames(local_gebv)))
        return(c(var = NA_real_, block_mid = NA_real_))
      var_total <- 0; mid_total <- 0
      for (b in bids) {
        contrib <- .block_contrib_independent(local_gebv[p1, b],
                                                local_gebv[p2, b],
                                                segregation_factor)
        var_total <- var_total + unname(contrib["var"])
        mid_total <- mid_total + unname(contrib["mean"])
      }
      c(var = var_total, block_mid = mid_total)
    })

  } else if (variance_model == "phased") {
    if (is.null(haplotypes) || is.null(snp_info) || is.null(snp_effects))
      stop("variance_model = 'phased' requires haplotypes (phased ",
           "extract_haplotypes() output), snp_info, and snp_effects.",
           call. = FALSE)
    if (!all(c("SNP", "CHR", "POS") %in% names(snp_info)))
      stop("snp_info must have SNP, CHR, POS columns.", call. = FALSE)

    dip <- infer_block_haplotypes(haplotypes, resolve_unphased = FALSE)
    if (!nrow(dip))
      stop("infer_block_haplotypes(haplotypes) returned no rows -- check ",
           "that `haplotypes` is phased output from extract_haplotypes().",
           call. = FALSE)
    dip <- dip[dip$block_id %in% block_importance$block_id, , drop = FALSE]
    if (!nrow(dip))
      stop("None of block_importance's block_id values are present in ",
           "haplotypes.", call. = FALSE)

    if (isTRUE(verbose))
      message("[usefulness_criterion] phased mode: ",
              length(unique(dip$block_id)), " target block(s), ",
              nrow(cross_pairs), " candidate cross(es).")

    block_ids_use <- unique(dip$block_id)

    res <- lapply(seq_len(nrow(cross_pairs)), function(r) {
      p1 <- cross_pairs$parent1[r]; p2 <- cross_pairs$parent2[r]
      var_total <- 0; mid_total <- 0; any_block <- FALSE
      for (bid in block_ids_use) {
        row1 <- dip[dip$id == p1 & dip$block_id == bid, ]
        row2 <- dip[dip$id == p2 & dip$block_id == bid, ]
        if (!nrow(row1) || !nrow(row2)) next
        if (anyNA(c(row1$hap1[1], row1$hap2[1], row2$hap1[1], row2$hap2[1])))
          next
        snp_ids <- .block_snp_order(snp_info, row1$CHR[1], row1$start_bp[1],
                                     row1$end_bp[1])
        if (!length(snp_ids)) next
        eff_i1 <- .phased_allele_effect(row1$hap1[1], snp_ids, snp_effects)
        eff_i2 <- .phased_allele_effect(row1$hap2[1], snp_ids, snp_effects)
        eff_j1 <- .phased_allele_effect(row2$hap1[1], snp_ids, snp_effects)
        eff_j2 <- .phased_allele_effect(row2$hap2[1], snp_ids, snp_effects)
        contrib <- .block_contrib_phased(eff_i1, eff_i2, eff_j1, eff_j2)
        var_total <- var_total + unname(contrib["var"])
        mid_total <- mid_total + unname(contrib["mean"])
        any_block <- TRUE
      }
      if (!any_block) return(c(var = NA_real_, block_mid = NA_real_))
      c(var = var_total, block_mid = mid_total)
    })

  } else {
    # linked -- same phased-haplotype data as "phased" mode, PLUS a genetic
    # map, so blocks are no longer treated as independent: recombination
    # fractions between blocks (Haldane, from genetic distance) correlate
    # which haplotype is transmitted at nearby blocks within the same
    # parent's gamete. See the long comment above .block_contrib_linked_mc()
    # for the full derivation and its "Scope" note (this models a single
    # cross's F1-style segregation variance, like "phased" -- NOT a
    # multi-generation RIL/DH population variance, which is what
    # "simplemating"'s Type/Generation model instead).
    if (is.null(haplotypes) || is.null(snp_info) || is.null(snp_effects))
      stop("variance_model = 'linked' requires haplotypes (phased ",
           "extract_haplotypes() output), snp_info, and snp_effects (same ",
           "as variance_model = 'phased').", call. = FALSE)
    if (is.null(genetic_map) && is.null(ld_matrix))
      stop("variance_model = 'linked' additionally requires genetic_map ",
           "(data frame: SNP, CHR, cM) or ld_matrix (a SNP x SNP LD matrix, ",
           "used as a recombination-fraction PROXY when no genetic map is ",
           "available -- see 'ld_matrix' in ?usefulness_criterion) -- this ",
           "is what distinguishes 'linked' from 'phased', which assumes ",
           "every target block segregates independently of every other. ",
           "Note: 'linked' always additionally requires phased `haplotypes` ",
           "regardless of which of these two you supply -- ld_matrix only ",
           "relaxes the genetic-map requirement, not the phasing ",
           "requirement checked just above.", call. = FALSE)
    if (!all(c("SNP", "CHR", "POS") %in% names(snp_info)))
      stop("snp_info must have SNP, CHR, POS columns.", call. = FALSE)
    if (!is.null(genetic_map) && !all(c("SNP", "CHR", "cM") %in% names(genetic_map)))
      stop("genetic_map must have SNP, CHR, cM columns.", call. = FALSE)
    if (!is.null(ld_matrix) &&
        (is.null(rownames(ld_matrix)) || is.null(colnames(ld_matrix)) ||
         !identical(rownames(ld_matrix), colnames(ld_matrix))))
      stop("ld_matrix must be a square matrix with identical, non-NULL row ",
           "and column names (SNP IDs).", call. = FALSE)
    if (!is.numeric(n_sim_linked) || length(n_sim_linked) != 1L || n_sim_linked < 100L)
      stop("n_sim_linked must be a single integer >= 100.", call. = FALSE)

    dip <- infer_block_haplotypes(haplotypes, resolve_unphased = FALSE)
    if (!nrow(dip))
      stop("infer_block_haplotypes(haplotypes) returned no rows -- check ",
           "that `haplotypes` is phased output from extract_haplotypes().",
           call. = FALSE)
    dip <- dip[dip$block_id %in% block_importance$block_id, , drop = FALSE]
    if (!nrow(dip))
      stop("None of block_importance's block_id values are present in ",
           "haplotypes.", call. = FALSE)

    block_coords <- dip[!duplicated(dip$block_id),
                        c("block_id", "CHR", "start_bp", "end_bp"), drop = FALSE]

    # Genetic-map path (preferred: real cM distances via Haldane's mapping
    # function) when genetic_map is supplied; otherwise fall back to
    # physical block order (CHR, start_bp) with recombination fractions
    # PROXIED from ld_matrix -- see .adjacent_r_from_ld()'s header comment
    # for the proxy formula and its caveats. Both paths produce a block_pos
    # with a block_id column in the resolved order, so everything past this
    # point (snp_ids_by_block construction, the per-cross simulation loop)
    # is shared code, branching only on `using_genetic_map` where the two
    # paths actually differ (position source, recombination-fraction source).
    using_genetic_map <- !is.null(genetic_map)
    if (using_genetic_map) {
      block_pos <- .block_genetic_positions(block_coords, snp_info, genetic_map)
      n_dropped <- nrow(block_coords) - nrow(block_pos)
      if (!nrow(block_pos))
        stop("None of the target blocks have any SNP in common with ",
             "genetic_map -- cannot place any block on the genetic map.",
             call. = FALSE)
      if (n_dropped && isTRUE(verbose))
        message("[usefulness_criterion] linked mode: ", n_dropped, " target ",
                "block(s) dropped (no SNP in common with genetic_map).")
    } else {
      block_pos <- .block_physical_positions(block_coords)
      n_dropped <- 0L  # physical position is always known; nothing to drop
    }

    block_ids_use <- block_pos$block_id   # already sorted (CHR, cM or start_bp)
    snp_ids_by_block <- stats::setNames(
      lapply(seq_len(nrow(block_pos)), function(i) {
        bc <- block_coords[block_coords$block_id == block_pos$block_id[i], , drop = FALSE][1L, ]
        .block_snp_order(snp_info, bc$CHR[1], bc$start_bp[1], bc$end_bp[1])
      }),
      block_pos$block_id
    )
    r_adjacent <- if (using_genetic_map)
      .adjacent_r(block_pos)
    else
      .adjacent_r_from_ld(block_pos, snp_ids_by_block, ld_matrix)

    if (isTRUE(verbose))
      message("[usefulness_criterion] linked mode: ", length(block_ids_use),
              " target block(s) placed on the ",
              if (using_genetic_map) "genetic map" else
                "physical map (ld_matrix-proxied recombination -- see ?usefulness_criterion)",
              ", ", nrow(cross_pairs), " candidate cross(es), n_sim_linked = ",
              n_sim_linked, ".")

    res <- lapply(seq_len(nrow(cross_pairs)), function(r) {
      p1 <- cross_pairs$parent1[r]; p2 <- cross_pairs$parent2[r]
      eff_i1_vec <- numeric(0); eff_i2_vec <- numeric(0)
      eff_j1_vec <- numeric(0); eff_j2_vec <- numeric(0)
      r_use <- numeric(0)
      keep  <- logical(length(block_ids_use))
      for (k in seq_along(block_ids_use)) {
        bid <- block_ids_use[k]
        row1 <- dip[dip$id == p1 & dip$block_id == bid, ]
        row2 <- dip[dip$id == p2 & dip$block_id == bid, ]
        ok <- nrow(row1) > 0L && nrow(row2) > 0L &&
          !anyNA(c(row1$hap1[1], row1$hap2[1], row2$hap1[1], row2$hap2[1])) &&
          length(snp_ids_by_block[[bid]]) > 0L
        keep[k] <- ok
        if (!ok) next
        snp_ids <- snp_ids_by_block[[bid]]
        eff_i1_vec <- c(eff_i1_vec, .phased_allele_effect(row1$hap1[1], snp_ids, snp_effects))
        eff_i2_vec <- c(eff_i2_vec, .phased_allele_effect(row1$hap2[1], snp_ids, snp_effects))
        eff_j1_vec <- c(eff_j1_vec, .phased_allele_effect(row2$hap1[1], snp_ids, snp_effects))
        eff_j2_vec <- c(eff_j2_vec, .phased_allele_effect(row2$hap2[1], snp_ids, snp_effects))
      }
      if (!length(eff_i1_vec)) return(c(var = NA_real_, block_mid = NA_real_))
      # r_adjacent must be re-subset to only the ADJACENT PAIRS both of
      # whose blocks survived filtering above; blocks skipped for missing
      # data break adjacency, so recombination fractions are recomputed
      # from block_pos$cM (or, in the ld_matrix path, re-proxied from
      # ld_matrix) directly on the surviving block subset rather than
      # naively sub-indexing the original r_adjacent vector.
      surviving_pos <- block_pos[keep, , drop = FALSE]
      r_use <- if (using_genetic_map)
        .adjacent_r(surviving_pos)
      else
        .adjacent_r_from_ld(surviving_pos, snp_ids_by_block, ld_matrix)
      contrib <- .block_contrib_linked_mc(eff_i1_vec, eff_i2_vec,
                                          eff_j1_vec, eff_j2_vec,
                                          r_use, n_sim_linked, seed = seed)
      c(var = unname(contrib["var"]), block_mid = unname(contrib["mean"]))
    })
  }

  res_mat <- do.call(rbind, res)
  n_skipped <- sum(is.na(res_mat[, "var"]))
  if (n_skipped && isTRUE(verbose))
    message("[usefulness_criterion] ", n_skipped, " of ", nrow(cross_pairs),
            " candidate cross(es) could not be scored (missing genotype/",
            "phase data at every target block) and are NA.")

  mid_parent_gebv <- unname((gebv[cross_pairs$parent1] + gebv[cross_pairs$parent2]) / 2)
  pred_var <- as.numeric(res_mat[, "var"])

  out <- data.frame(
    parent1              = cross_pairs$parent1,
    parent2              = cross_pairs$parent2,
    mid_parent_gebv       = mid_parent_gebv,
    predicted_variance    = pred_var,
    selection_intensity   = i_sel,
    UC                    = mid_parent_gebv + i_sel * sqrt(pmax(pred_var, 0)),
    stringsAsFactors      = FALSE
  )
  ord <- order(-out$UC, na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  out <- .augment_uc_uncertainty(
    out, gebv_se, gebv_reliability, phasing_reliability,
    # "linked" depends on phased-haplotype accuracy just as much as "phased"
    # (both build progeny variance from phased haplotype blocks), so both
    # get the strict phasing-reliability gate; "block_independent" and
    # "simplemating" do not use phased haplotypes for progeny variance.
    variance_model %in% c("phased", "linked"), n_progeny,
    downside_quantile, min_reliability
  )
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}
