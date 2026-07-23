# ==============================================================================
# family_selection.R
#
# Family-quota parent selection -- a third parent-shortlist strategy,
# alongside truncation_selection() (Section 5.1-style, single-score ranking)
# and select_parents_ga() (Section 5.2-style, block-coverage search). Where
# those two ignore family/pedigree structure entirely (truncation_selection()
# by construction; select_parents_ga() unless its optional coancestry_weight/
# target_degree penalty is switched on), this function puts group structure
# at the centre of the selection rule itself: choose the n_families
# best-performing groups first, then the n_per_family best lines within each
# -- the natural, familiar shape of a real breeding-program shortlist ("take
# our best few families/groups, and the best few lines out of each").
#
# "Group" here means one of two things, chosen via group_by:
#   - "family" (default): the caller's own pedigree/cross-ID labels
#     (`family` argument) -- unchanged from the original design.
#   - "genetic_cluster": groups are derived directly from a relationship
#     matrix `G` via hierarchical clustering (Ward's minimum-variance
#     linkage by default) into `n_clusters` genetically coherent groups,
#     rather than trusting pedigree labels to be a faithful proxy for actual
#     relatedness (they aren't always -- mislabelling, incomplete records,
#     or loosely-applied family IDs can all decouple "same family label"
#     from "actually closely related"). Everything downstream (ranking,
#     quota-filling, diagnostics) operates identically regardless of which
#     mode produced the grouping.
# Both a pedigree `family` label and a genetic-cluster label are reported
# per individual whenever available (`by_family$family` /
# `by_family$genetic_group`), independent of which one actually drove the
# selection decision (`$group_by` echoes that) -- so a caller using genetic
# clusters for the actual quota logic can still see each pick's pedigree
# family for cross-referencing, and vice versa.
#
# Within that structure, several further refinements are layered on:
#
#   - Relatedness control: whenever `G` is supplied, the realised mean
#     pairwise relationship of the final shortlist (overall and per group)
#     is always reported, regardless of whether anything was asked to
#     optimise for it (mirrors select_parents_ga()'s own
#     $mean_relationship). On top of that, `within_group_target_degree`
#     (0-90, the SAME convention select_parents_ga()/select_parents_ocs()
#     use) turns this into an actual dial: filling a group's quota stops
#     being pure score-order truncation and instead balances score against
#     a relatedness ceiling, interpolated between a genuinely constructed
#     high-score reference group and select_core_collection()'s proven
#     diversity-maximising one, restricted to that group's own members.
#
#   - Shrinkage-corrected family ranking (the new default -- see NEWS for
#     the BREAKING note): families are no longer ranked by a raw sample
#     mean of their own top-k scores. A family's raw top-k mean carries two
#     DIFFERENT, independently-corrected biases:
#       (1) it is an unreliable estimator of the family's true quality when
#           the family is small -- exactly the kind of small-sample noise
#           classical quantitative genetics corrects for via shrinkage /
#           BLUP-style family evaluation. Variance components (true
#           between-family variance tau2 and within-family/residual
#           variance sigma2) are estimated from every eligible member's raw
#           score, via either REML (variance_method = "reml", the new
#           default -- lme4::lmer(score ~ 1 + (1|family)), more efficient
#           than method-of-moments under the unbalanced family sizes this
#           function always has to deal with) or the classical one-way
#           random-effects ANOVA method-of-moments estimator
#           (variance_method = "anova"; Searle, Casella & McCulloch 1992;
#           also the automatic, non-silent fallback whenever REML is
#           requested but lme4 is not installed or the model fails to fit).
#           Each family's top-k mean is then shrunk toward the across-family
#           mean by the resulting empirical-Bayes/BLUP reliability weight
#           w_f = n_f*tau2/(n_f*tau2+sigma2) (Robinson 1991) -- a family of 1
#           is pulled hard toward the pack; a 20-member family's genuine
#           top-k average is trusted almost fully.
#       (2) separately, "mean of a family's own best k" is a systematically
#           UPWARD-biased estimator of that family's true mean, purely from
#           having picked winners out of a finite pool -- present even for a
#           family whose individual scores carry zero estimation error of
#           their own. bias_correction = "order_stats" (the new default)
#           removes this via a finite-population order-statistics
#           correction computed by direct numerical integration (see
#           .expected_topk_mean_std_normal()), BEFORE shrinkage is applied.
#           bias_correction = "none" restores the old, uncorrected topk_mean.
#     On top of both corrections, when use_family_relationship = TRUE (the
#     default) and a relationship matrix G is supplied, plain i.i.d.
#     shrinkage toward the grand mean is replaced by a genomic-relationship-
#     informed (GBLUP-style) family-effect BLUP (see .rank_families_gblup())
#     that lets genetically related families borrow strength from EACH
#     OTHER too, not only from the grand mean -- a strict generalisation of
#     the scalar w_f formula (which it reduces to exactly when G_fam = I).
#     Falls back to the historical, unshrunk ranking (with a message, not
#     silently) whenever variance components can't be estimated at all, and
#     to plain i.i.d. shrinkage whenever G doesn't cover every eligible
#     family.
#
#   - Family-size eligibility rule (new): a family/group is excluded
#     entirely -- before ranking, before it can count toward n_families --
#     whenever its eligible membership does not exceed the number of lines
#     that would actually be taken from it (n_per_family in the common
#     scalar-quota case; rank_k as the practical stand-in for the uneven,
#     named-vector-quota case, where n_per_family is only defined for
#     whichever families end up chosen). Below that size there is no
#     genuine "select the best of" decision happening at all -- every
#     eligible member would be taken regardless of ranking -- so keeping
#     such a family in the pool would let it count toward n_families for
#     free, and would feed a topk_mean with zero real selection content
#     into the shrinkage/bias-correction machinery above.
#
#   - rank_k, decoupled from n_per_family: the family-ranking criterion's
#     sample size (how many of a family's own top members its ranking mean
#     is computed from) no longer has to equal the actual take-quota.
#
#   - Two haplotype/coverage diversity methods for the optional
#     ensure_haplotype_diversity adjustment: the original "dominant_block"
#     heuristic (each individual's single highest-value target block), or
#     the more rigorous "coverage_gain" (reuses select_parents_ga()'s own
#     .block_best_values() to ask whether a candidate's addition actually
#     increases total coverage across every target block, given everyone
#     already claimed so far -- not just whether one single block collides).
#
#   - Uneven per-family/per-group quotas: n_per_family may be a single
#     integer (flat quota, unchanged default) or a named integer vector
#     (explicit override per chosen family/group).
#
#   - family_select_mode (new): four independent, opt-in ways to decide
#     how many/which of a chosen family's members are actually taken.
#     "count" (default) is the original fixed-number-per-family quota
#     above. "percentage" takes a fixed SHARE of each family's own eligible
#     size (pct_per_family, rounded up via ceiling() so an included family
#     never contributes zero) -- still a quota (fixed count once resolved),
#     so it keeps its own version of the family-size eligibility rule
#     above, computed per family from its own size: excluded whenever
#     ceiling(pct_per_family/100 * group_size) == group_size. This is NOT
#     generally as fragile as count mode's version of the rule -- solving
#     that equality shows it triggers whenever pct_per_family exceeds
#     100*(n-1)/n for a group of size n, a threshold that climbs toward
#     100 as n grows. At a realistic rate like 10%, that threshold is only
#     crossed at n = 1 (a lone member always gets excluded, since ANY
#     positive percentage of 1 rounds up to taking it regardless of
#     merit) -- a 2-member, 10-member, or 100-member group is untouched by
#     the rule at 10%. It only starts excluding larger groups too at high
#     percentages (e.g. 90% excludes every group of size 9 or smaller,
#     because "90% of 9" already rounds up to all of them) -- a
#     percentage-dependent, narrowing edge case, not a blanket parallel to
#     count mode's much easier-to-trigger version. "sd_threshold" and
#     "check_relative" are different in kind: they are merit BARS, not
#     quotas -- every member of a chosen family clearing the bar is taken,
#     however many that turns out to be (zero to all), so there is no
#     size-based eligibility rule for them at all (a 1-member family can
#     meaningfully pass or fail a bar, unlike being forced to take its only
#     member regardless of merit under a fixed quota). "sd_threshold" is a
#     signed number of SDs from the population mean of ALL eligible
#     candidates (positive = an elite bar above the mean); "check_relative"
#     is a minimum percentage margin over a reference check's own score
#     (check_id, looked up from score, or a fixed check_value) -- e.g.
#     check_margin_pct = 10 keeps everyone at least 10% above the check,
#     including lines that clear it by 20% or 30%. Because a chosen family
#     can genuinely contribute zero members under either threshold mode,
#     family selection auto-backfills: any ranked family contributing
#     nothing is skipped (reported in $zero_selected_groups) and the next-
#     ranked family is tried instead, until n_families worth of
#     CONTRIBUTING groups are assembled or the ranking is exhausted.
#
# An optional second layer (ensure_haplotype_diversity = TRUE) addresses a
# specific risk the quota structure alone does not solve: picking the single
# best line from group A and the single best line from group B by score
# alone says nothing about whether those two lines are genomically
# complementary. When switched on, within-group selection is adjusted
# (greedily, group by group in rank order) to prefer lines that extend
# coverage of favourable target-block value beyond what has already been
# claimed by higher-ranked groups' picks, always keeping score as the
# tie-breaker so the constraint never silently discards a group's genuinely
# best line for a marginal diversity gain -- collisions are flagged, not
# hidden.
# ==============================================================================


# -- Internal: relationship matrix -> exact Euclidean-equivalent distance
# matrix, via the identity D_ij = G_ii + G_jj - 2*G_ij (the same conversion
# select_core_collection() uses -- not an approximation). Shared here rather
# than exported.
.relationship_to_distance <- function(G) {
  dg <- diag(G)
  D <- outer(dg, dg, "+") - 2 * G
  D[D < 0] <- 0
  diag(D) <- 0
  D
}

# -- Internal: hierarchical clustering of `ids` into `n_clusters` genetic
# groups from relationship matrix G, via hclust() on the distance matrix
# above (Ward's minimum-variance linkage, "ward.D2", by default --
# cluster_method is passed straight to hclust() for callers who want a
# different linkage criterion). Base R only (stats::hclust/cutree/as.dist),
# no new package dependency. Cluster labels are relabelled 1..n_clusters by
# DECREASING mean `score`, so "cluster_1" is consistently the
# strongest-performing genetic group rather than an arbitrary hclust/cutree
# merge-order artefact -- purely a labelling convenience, `score` plays no
# role in which individuals end up in which cluster (that is G/distance
# alone). Returns a named character vector (names = ids).
.cluster_by_relationship <- function(G, ids, n_clusters, cluster_method, score) {
  n_clusters <- as.integer(n_clusters)
  if (is.na(n_clusters) || n_clusters < 1L)
    stop("n_clusters must be >= 1.", call. = FALSE)
  if (n_clusters > length(ids))
    stop("n_clusters (", n_clusters, ") exceeds the number of eligible ",
         "candidates (", length(ids), ").", call. = FALSE)
  if (n_clusters == 1L)
    return(stats::setNames(rep("cluster_1", length(ids)), ids))

  Gs <- G[ids, ids, drop = FALSE]
  D  <- .relationship_to_distance(Gs)
  if (max(D) <= .Machine$double.eps)
    stop("All pairwise distances derived from G are ~0 among eligible ",
         "candidates (everyone is genetically identical under G) -- cannot ",
         "form meaningful genetic clusters.", call. = FALSE)

  hc  <- stats::hclust(stats::as.dist(D), method = cluster_method)
  raw <- stats::cutree(hc, k = n_clusters)
  raw <- raw[ids]  # defensive realignment to `ids` order

  cl_means <- tapply(score[ids], raw, mean)
  ord      <- order(-cl_means)
  relabel  <- stats::setNames(paste0("cluster_", seq_along(ord)),
                              names(cl_means)[ord])
  stats::setNames(unname(relabel[as.character(raw)]), ids)
}

# -- Internal: for every family/group level present in `score`'s names, the
# mean of that group's own top-k scores (k = rank_k, or fewer if the group
# has fewer members) -- the raw ranking criterion .rank_families() shrinks
# (or uses as-is under family_rank_method = "topk_mean"). Returns a data
# frame with one row per group: family, topk_mean, n_members.
.family_topk_mean <- function(score, family, k) {
  fam_levels <- unique(family)
  out <- lapply(fam_levels, function(f) {
    members <- names(family)[family == f]
    members <- intersect(members, names(score))
    s <- sort(score[members], decreasing = TRUE)
    kk <- min(k, length(s))
    data.frame(
      family     = f,
      topk_mean  = if (kk > 0L) mean(s[seq_len(kk)]) else NA_real_,
      n_members  = length(s),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# -- Internal: standard one-way random-effects ANOVA method-of-moments
# variance-components estimator (Searle, Casella & McCulloch 1992,
# "Variance Components"), used by .rank_families()'s shrinkage-corrected
# ranking as the variance_method = "anova" option (and as the automatic
# fallback whenever variance_method = "reml" cannot be used -- see
# .estimate_family_variance_components_reml() and the dispatcher
# .estimate_family_variance_components() below). Estimated from EVERY
# eligible individual's raw `score` (not from any already-selected/
# truncated statistic like topk_mean, which would bias the estimate),
# decomposed by whichever grouping is active. Returns tau2 (estimated true
# between-group variance, floored at 0), sigma2 (within-group/residual
# variance), and ok = FALSE when there is not enough structure to estimate
# them at all (fewer than 2 groups, or as many groups as individuals).
# Known limitation motivating variance_method = "reml": method-of-moments
# loses efficiency precisely when group sizes are unbalanced -- which they
# always are here, since uneven family size is the entire reason shrinkage
# is needed in the first place.
.estimate_family_variance_components_anova <- function(score, family) {
  fam        <- as.character(family)
  fam_levels <- unique(fam)
  a <- length(fam_levels)
  N <- length(score)
  if (a < 2L || N <= a)
    return(list(tau2 = NA_real_, sigma2 = NA_real_, ok = FALSE))

  grand_mean <- mean(score)
  fam_means  <- tapply(score, fam, mean)
  fam_ns     <- tapply(score, fam, length)

  SSB <- sum(fam_ns * (fam_means - grand_mean)^2)
  SSW <- sum((score - fam_means[fam])^2)
  MSB <- SSB / (a - 1)
  MSW <- SSW / (N - a)
  n0  <- (N - sum(fam_ns^2) / N) / (a - 1)

  sigma2 <- MSW
  tau2   <- if (n0 > sqrt(.Machine$double.eps)) max(0, (MSB - MSW) / n0) else 0

  list(tau2 = tau2, sigma2 = sigma2, ok = TRUE, n0 = n0, a = a, N = N)
}

# -- Internal: REML one-way random-effects variance-components estimator,
# via lme4::lmer(score ~ 1 + (1 | family)). More statistically efficient
# than .estimate_family_variance_components_anova() under unbalanced group
# sizes (see that function's own comment for why this matters here
# specifically). Returns NULL -- not an error -- whenever lme4 is not
# installed, there isn't enough structure to fit the model, or the fit
# itself fails; the dispatcher below falls back to the ANOVA estimator in
# every such case, with a message, not silently. A singular fit (tau2
# estimated at the 0 boundary) is NOT treated as a failure here -- it is a
# valid, informative REML result (essentially no real between-family signal
# once residual noise is accounted for), returned as tau2 = 0 exactly like
# the ANOVA path's own degenerate branch.
.estimate_family_variance_components_reml <- function(score, family) {
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)

  fam        <- as.character(family)
  fam_levels <- unique(fam)
  a <- length(fam_levels)
  N <- length(score)
  if (a < 2L || N <= a) return(NULL)

  df <- data.frame(score = as.numeric(score), family = fam,
                   stringsAsFactors = FALSE)
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(score ~ 1 + (1 | family), data = df, REML = TRUE,
                control = lme4::lmerControl(check.conv.singular = "ignore"))
    )),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(fit)) return(NULL)

  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  if (is.null(vc)) return(NULL)
  tau2_row   <- vc[vc$grp == "family", ]
  sigma2_row <- vc[vc$grp == "Residual", ]
  if (nrow(tau2_row) != 1L || nrow(sigma2_row) != 1L) return(NULL)

  singular <- tryCatch(isTRUE(lme4::isSingular(fit)), error = function(e) NA)

  list(tau2 = max(0, tau2_row$vcov), sigma2 = sigma2_row$vcov, ok = TRUE,
       n0 = NA_real_, a = a, N = N, singular = singular)
}

# -- Internal: dispatches variance-component estimation to REML or ANOVA
# per `method`, with a transparent, non-silent fallback from REML to ANOVA
# whenever REML is unavailable (lme4 not installed) or fails to fit. Tags
# the result with $method_used so callers/diagnostics can see which
# estimator actually produced tau2/sigma2, independent of what was
# requested.
.estimate_family_variance_components <- function(score, family, method, verbose) {
  if (method == "reml") {
    vc <- .estimate_family_variance_components_reml(score, family)
    if (!is.null(vc)) {
      vc$method_used <- "reml"
      return(vc)
    }
    if (isTRUE(verbose))
      message("[select_parents_by_family] REML variance-component ",
              "estimation unavailable (lme4 not installed, not enough ",
              "family/group structure, or the model failed to fit) -- ",
              "falling back to the ANOVA method-of-moments estimator.")
  }
  vc <- .estimate_family_variance_components_anova(score, family)
  vc$method_used <- "anova"
  vc
}

# -- Internal: E[X_(i)], the expected value of the i-th order statistic
# (ascending, 1-indexed) of n i.i.d. standard normal draws, via direct
# numerical integration of its known density
#   f_(i)(x) = n!/((i-1)!(n-i)!) * Phi(x)^(i-1) * (1-Phi(x))^(n-i) * phi(x)
# rather than a tabulated closed-form approximation (e.g. Burrows 1972) --
# see .expected_topk_mean_std_normal() for why. Evaluated in log-space
# (log_coef + (i-1)*log(Phi(x)) + (n-i)*log(1-Phi(x)) + log(phi(x)), then
# exponentiated) because Phi(x)^(i-1) or (1-Phi(x))^(n-i) can underflow to
# exactly 0 in ordinary scale well before the true integrand is negligible,
# for large (i, n-i). Integration is bounded to [-12, 12] rather than
# literal +/-Inf -- ample range for standard-normal tail mass at any n
# realistic for a breeding-program family/group size, and more numerically
# robust for stats::integrate() than open infinite bounds.
.expected_order_stat_std_normal <- function(n, i) {
  log_coef <- lgamma(n + 1) - lgamma(i) - lgamma(n - i + 1)
  integrand <- function(x) {
    log_dens <- log_coef +
      (i - 1) * stats::pnorm(x, log.p = TRUE) +
      (n - i) * stats::pnorm(x, lower.tail = FALSE, log.p = TRUE) +
      stats::dnorm(x, log = TRUE)
    x * exp(log_dens)
  }
  val <- tryCatch(
    stats::integrate(integrand, lower = -12, upper = 12,
                     rel.tol = 1e-8, stop.on.error = FALSE)$value,
    error = function(e) NA_real_
  )
  val
}

# -- Internal: E[mean of the top k order statistics out of n i.i.d.
# standard normal draws] -- the expected STANDARDISED selection
# differential a finite group of size n produces when its own best k
# members are averaged, from order-statistics sampling alone. Exactly 0
# when k >= n (taking every member is not a selection at all -- consistent
# with the new family-size exclusion rule in select_parents_by_family(),
# which removes groups in exactly that regime before this is ever called
# for them). By location-scale invariance of the normal distribution, a
# group with true within-group SD sqrt(sigma2) has an expected topk_mean
# inflated by sqrt(sigma2) * this value ABOVE its own true mean -- see
# .bias_correct_topk_mean(). This computes the same underlying quantity a
# tabulated closed-form approximation like Burrows (1972) targets, but
# directly from first principles via numerical integration of the order
# statistic's own density, rather than a transcribed/memorised coefficient
# -- deliberately chosen so its correctness can be checked by reasoning
# about the integral itself. Memoised per (n, k) pair for the life of the R
# session (a pure mathematical constant, independent of any data -- many
# families/groups in the same run, or across runs, often share a size).
.expected_topk_mean_std_normal <- local({
  cache <- new.env(parent = emptyenv())
  function(n, k) {
    n <- as.integer(n); k <- as.integer(min(k, n))
    if (k <= 0L || n <= 0L || k >= n) return(0)
    key <- paste0(n, "_", k)
    cached <- cache[[key]]
    if (!is.null(cached)) return(cached)
    idx  <- (n - k + 1L):n
    vals <- vapply(idx, function(i) .expected_order_stat_std_normal(n, i), numeric(1))
    out  <- if (any(!is.finite(vals))) NA_real_ else mean(vals)
    cache[[key]] <- out
    out
  }
})

# -- Internal: bias-corrected topk_mean, subtracting the expected
# order-statistics selection-differential inflation
# (sqrt(sigma2) * .expected_topk_mean_std_normal(n_members, k)) from each
# group's raw topk_mean, BEFORE any shrinkage is applied. This targets a
# DIFFERENT bias than shrinkage: shrinkage corrects for how much a group's
# mean can be TRUSTED given its size (a reliability question); this
# corrects for the fact that "mean of a group's own best k" systematically
# OVERSTATES that group's true mean, purely from having picked winners out
# of a finite pool -- present even for a group whose individual scores
# carry zero estimation error of their own. Requires sigma2 (see
# .estimate_family_variance_components()); returns topk_mean unchanged
# (bias = 0) when sigma2 is not finite or not positive, since there is
# nothing to correct in that degenerate case.
.bias_correct_topk_mean <- function(fam_tab, k, sigma2) {
  if (!isTRUE(is.finite(sigma2)) || sigma2 <= 0) {
    fam_tab$topk_mean_bias      <- rep(0, nrow(fam_tab))
    fam_tab$topk_mean_corrected <- fam_tab$topk_mean
    return(fam_tab)
  }
  bias <- vapply(fam_tab$n_members, function(n_f) {
    e <- sqrt(sigma2) * .expected_topk_mean_std_normal(n_f, k)
    if (is.finite(e)) e else 0
  }, numeric(1))
  fam_tab$topk_mean_bias      <- bias
  fam_tab$topk_mean_corrected <- fam_tab$topk_mean - bias
  fam_tab
}

# -- Internal: aggregate an individual-level relationship matrix G into a
# family/group-level relationship matrix -- G_fam[f, f2] = mean of every
# pairwise G entry between a member of f and a member of f2 (including
# f == f2's own within-group block, individual self-relationships and all
# -- the simplest, uniformly-defined aggregation, avoiding a separate
# special case for singleton groups). This is the same idea behind building
# a numerator relationship matrix among sires/families for classical
# General Combining Ability estimation from a full pedigree, applied here
# directly from genomic G instead. Returns NULL (caller falls back to
# independent/i.i.d. shrinkage) whenever G does not cover every eligible
# individual in `family`.
.aggregate_relationship_to_family <- function(G, family) {
  ids <- names(family)
  fam <- as.character(family)
  if (!length(ids) || !all(ids %in% rownames(G)) || !all(ids %in% colnames(G)))
    return(NULL)
  fam_levels <- unique(fam)
  a <- length(fam_levels)
  members <- lapply(fam_levels, function(f) ids[fam == f])
  names(members) <- fam_levels
  G_fam <- matrix(NA_real_, a, a, dimnames = list(fam_levels, fam_levels))
  for (i in seq_len(a)) {
    for (j in i:a) {
      block <- G[members[[i]], members[[j]], drop = FALSE]
      val <- mean(block)
      G_fam[i, j] <- G_fam[j, i] <- val
    }
  }
  G_fam
}

# -- Internal: genomic-relationship-informed family-effect BLUP, replacing
# the plain i.i.d. shrinkage weight w_f = n_f*tau2/(n_f*tau2+sigma2) with a
# proper GLS/BLUP solution that lets genetically related groups borrow
# strength from EACH OTHER, not only from the grand mean. Model:
#   y_f = mu + a_f + eps_f,  a ~ MVN(0, tau2 * G_fam),  eps_f ~ N(0, sigma2/n_f)
# (y_f = each group's bias-corrected topk_mean; the sigma2/n_f sampling-
# variance approximation for eps_f matches the one the plain i.i.d. formula
# already relies on -- this generalises it from a scalar to a full
# G_fam-structured case, rather than introducing a new assumption). Solved
# directly via Henderson's mixed-model equations (one small, (a+1)x(a+1)
# linear system, a = number of candidate groups) -- exact given tau2,
# sigma2 and G_fam, not an iterative approximation. Sanity check: when
# G_fam = I (no relatedness structure between groups at all), this reduces
# algebraically to mu_hat + w_f*(y_f - mu_hat) with w_f the same i.i.d.
# formula used elsewhere in this file -- i.e. it is a strict generalisation
# of the existing shrinkage, not a different, unrelated method. G_fam is
# ridge-regularised before inversion if it is (near-)singular -- small
# group counts make this easy to hit -- mirroring the `bend` convention
# used elsewhere in this package (e.g. compute_haplotype_grm()).
# Returns NULL (caller falls back to i.i.d. shrinkage) if G_fam does not
# cover every candidate group, or the linear system cannot be solved.
.rank_families_gblup <- function(fam_tab, tau2, sigma2, G_fam,
                                 y_col = "topk_mean_corrected") {
  fams <- fam_tab$family
  if (is.null(G_fam) || !all(fams %in% rownames(G_fam)) ||
      !all(fams %in% colnames(G_fam)))
    return(NULL)

  Gf <- G_fam[fams, fams, drop = FALSE]
  eig_min <- tryCatch(min(eigen(Gf, symmetric = TRUE, only.values = TRUE)$values),
                      error = function(e) NA_real_)
  if (!is.finite(eig_min) || eig_min < 1e-6)
    diag(Gf) <- diag(Gf) + (1e-6 - min(0, eig_min))

  Gf_inv <- tryCatch(solve(Gf), error = function(e) NULL)
  if (is.null(Gf_inv)) return(NULL)

  y   <- fam_tab[[y_col]]
  a   <- nrow(fam_tab)
  n_f <- fam_tab$n_members
  R_inv <- diag(n_f / sigma2, a)
  ones  <- rep(1, a)

  C11 <- sum(n_f) / sigma2
  C12 <- t(ones) %*% R_inv
  C21 <- t(C12)
  C22 <- R_inv + Gf_inv / tau2
  rhs1 <- sum(R_inv %*% y)
  rhs2 <- R_inv %*% y

  M   <- rbind(cbind(C11, C12), cbind(C21, C22))
  rhs <- c(rhs1, rhs2)
  sol <- tryCatch(solve(M, rhs), error = function(e) NULL)
  if (is.null(sol) || length(sol) != a + 1L) return(NULL)

  mu_hat <- sol[1]
  a_hat  <- sol[-1]
  fam_tab$rank_score <- mu_hat + a_hat
  fam_tab
}

# -- Internal: family/group ranking criterion. method = "shrunk_topk_mean"
# (default): raw topk_mean (see .family_topk_mean()) is (optionally
# bias-corrected, then) shrunk toward the across-group mean of topk_mean by
# an empirical-Bayes/BLUP reliability weight w_f = n_f*tau2 / (n_f*tau2 +
# sigma2) (e.g. Robinson 1991, "That BLUP Is a Good Thing") -- a group's
# ranking reflects how much its OWN DATA can be trusted, not just its raw
# point estimate -- UNLESS a genomic-relationship-informed alternative is
# available (G_fam supplied and covers every candidate group), in which
# case .rank_families_gblup() replaces the scalar w_f formula with a full
# GLS/BLUP solution that also lets related groups borrow strength from each
# other (see its own comment; algebraically a strict generalisation of the
# scalar formula, not a different method). tau2/sigma2 are estimated from
# every eligible member's raw score (not from topk_mean itself, to avoid a
# circular/biased estimate from an already-selected statistic) via either
# REML or ANOVA method-of-moments, per `variance_method` -- see
# .estimate_family_variance_components(). method = "topk_mean": the
# historical, unshrunk ranking (bias_correction/variance_method/G_fam are
# all ignored in that case). Falls back to unshrunk ranking (message, not
# silently) whenever variance components can't be estimated, or come out
# fully degenerate on the residual side.
.rank_families <- function(score, family, k, method, variance_method,
                           bias_correction, G_fam, verbose) {
  fam_tab <- .family_topk_mean(score, family, k)

  if (method == "topk_mean" || nrow(fam_tab) < 2L) {
    fam_tab$rank_score       <- fam_tab$topk_mean
    fam_tab$shrinkage_weight <- NA_real_
    return(fam_tab)
  }

  vc     <- .estimate_family_variance_components(score, family, variance_method, verbose)
  target <- mean(fam_tab$topk_mean, na.rm = TRUE)

  if (!isTRUE(vc$ok)) {
    if (isTRUE(verbose))
      message("[select_parents_by_family] Not enough family/group structure ",
              "to estimate shrinkage (need >= 2 groups with combined ",
              "membership exceeding the group count) -- falling back to ",
              "unshrunk topk_mean ranking.")
    fam_tab$rank_score       <- fam_tab$topk_mean
    fam_tab$shrinkage_weight <- NA_real_
    return(fam_tab)
  }

  fam_tab <- if (bias_correction == "order_stats") {
    .bias_correct_topk_mean(fam_tab, k, vc$sigma2)
  } else {
    fam_tab$topk_mean_bias      <- 0
    fam_tab$topk_mean_corrected <- fam_tab$topk_mean
    fam_tab
  }

  gblup_res <- NULL
  if (!is.null(G_fam) && vc$sigma2 > sqrt(.Machine$double.eps) && vc$tau2 > 0) {
    gblup_res <- .rank_families_gblup(fam_tab, vc$tau2, vc$sigma2, G_fam,
                                      y_col = "topk_mean_corrected")
  }

  if (!is.null(gblup_res)) {
    fam_tab <- gblup_res
    fam_tab$shrinkage_weight <- NA_real_  # not a single scalar under GBLUP
  } else {
    w <- if (vc$sigma2 <= sqrt(.Machine$double.eps)) {
      # Degenerate: essentially no within-group variance at all -- every
      # group mean is fully trustworthy regardless of size.
      rep(1, nrow(fam_tab))
    } else {
      fam_tab$n_members * vc$tau2 / (fam_tab$n_members * vc$tau2 + vc$sigma2)
    }
    fam_tab$shrinkage_weight <- w
    fam_tab$rank_score       <- target + w * (fam_tab$topk_mean_corrected - target)
  }

  attr(fam_tab, "tau2")                  <- vc$tau2
  attr(fam_tab, "sigma2")                <- vc$sigma2
  attr(fam_tab, "variance_method_used")  <- vc$method_used
  attr(fam_tab, "relationship_informed") <- !is.null(gblup_res)
  fam_tab
}

# -- Internal: resolve n_per_family (a single integer, or a named integer
# vector keyed by chosen family/group labels) into a named integer vector
# covering exactly `chosen_fams`.
.resolve_n_per_family <- function(n_per_family, chosen_fams) {
  # Dispatch on whether n_per_family is NAMED, not on its length -- a
  # length-1 NAMED vector (e.g. c(FamA = 3L)) means "a per-family quota
  # that happens to specify only one family so far", not "apply 3 to every
  # chosen family"; only a genuinely unnamed value is a flat/scalar quota.
  if (is.null(names(n_per_family))) {
    if (length(n_per_family) != 1L)
      stop("n_per_family, when supplying more than one value, must be a ",
           "NAMED integer vector (names = family/group labels) -- see ",
           "?select_parents_by_family's n_per_family argument.", call. = FALSE)
    return(stats::setNames(rep(as.integer(n_per_family), length(chosen_fams)),
                           chosen_fams))
  }
  missing_q <- setdiff(chosen_fams, names(n_per_family))
  if (length(missing_q))
    stop("n_per_family is missing an entry for selected family/group(s): ",
         paste(utils::head(missing_q, 10), collapse = ", "),
         if (length(missing_q) > 10) ", ..." else "", call. = FALSE)
  stats::setNames(as.integer(n_per_family[chosen_fams]), chosen_fams)
}

# -- Internal: dominant target block per individual -- the column of
# value_matrix at which that individual's own value is largest. Ties broken
# by column order (first max). Returns a named character vector (names =
# rownames(value_matrix)).
.dominant_block <- function(value_matrix) {
  db <- apply(value_matrix, 1L, function(row) colnames(value_matrix)[which.max(row)])
  stats::setNames(as.character(db), rownames(value_matrix))
}

# -- Internal: build an individuals x blocks_needed character matrix of each
# individual's raw allele/dosage string at each block, directly from
# extract_haplotypes()'s own return value. See ?select_parents_by_family's
# `haplotypes` argument for the full rationale.
.haplotypes_to_matrix <- function(haplotypes, blocks_needed, individuals_needed) {
  missing_blk <- setdiff(blocks_needed, names(haplotypes))
  if (length(missing_blk))
    stop(length(missing_blk), " value_matrix block(s) are not present in ",
         "haplotypes (the extract_haplotypes() output) -- value_matrix and ",
         "haplotypes must come from the same block-detection run: ",
         paste(utils::head(missing_blk, 10), collapse = ", "),
         if (length(missing_blk) > 10) ", ..." else "", call. = FALSE)

  out <- matrix(NA_character_, nrow = length(individuals_needed),
               ncol = length(blocks_needed),
               dimnames = list(individuals_needed, blocks_needed))
  for (blk in blocks_needed) {
    hs <- haplotypes[[blk]]
    missing_ind <- setdiff(individuals_needed, names(hs))
    if (length(missing_ind))
      stop(length(missing_ind), " eligible candidate(s) are missing from ",
           "haplotypes[[\"", blk, "\"]]: ",
           paste(utils::head(missing_ind, 10), collapse = ", "),
           if (length(missing_ind) > 10) ", ..." else "", call. = FALSE)
    out[, blk] <- unname(hs[individuals_needed])
  }
  out
}

# -- Internal: marginal block-coverage gain of adding `cand` to an
# already-claimed set `claimed_ids`, reusing select_parents_ga()'s own
# .block_best_values() so "coverage" means the same thing here as it does
# there. A single-individual value (like the "selfing"/"OPV" strategy in
# select_parents_ga() -- no specific pairing is modelled at the shortlist-
# construction stage; that is what usefulness_criterion()/
# select_parents_ocs() do downstream, once a shortlist exists).
.marginal_coverage_gain <- function(value_matrix, claimed_ids, cand, block_weights) {
  before <- if (length(claimed_ids))
    sum(.block_best_values(value_matrix, claimed_ids, strategy = "OPV") * block_weights)
  else 0
  after <- sum(.block_best_values(value_matrix, c(claimed_ids, cand), strategy = "OPV") *
              block_weights)
  after - before
}

# -- Internal: within-group pick order, trading score against relatedness
# via the SAME target_degree convention select_parents_ga()/
# select_parents_ocs() use (0 = max gain/score, 90 = max diversity),
# restricted to one group's own members and quota size. `members_by_score`
# is the group's FULL preference order (already accounting for any
# cross-group coverage/diversity adjustment -- see the main loop); this
# function walks that same order and additionally filters for relatedness,
# so cross-group diversity intent (if any) takes priority and relatedness is
# layered on top of it, in that disclosed order -- not jointly co-optimised.
#
# Gain end: the group's own top-n_take members from `members_by_score` --
# its "natural", no-relatedness-constraint pick -- own mean pairwise
# relationship. Diversity end: select_core_collection(strategy = "maximin")
# restricted to this group's own members -- proven 2-approximation
# (Gonzalez 1985), built ignoring score/coverage entirely. Ceiling
# interpolated between them. Slots are filled by walking
# `members_by_score` in order, accepting a candidate if the running mean
# relationship of picks-so-far would stay at or below the ceiling,
# deferring (not discarding) otherwise; if the quota still cannot be filled
# from acceptances alone, remaining slots are filled from the deferred pool
# in their original order -- this dial never drops a group's best remaining
# member purely to satisfy a relatedness preference, mirroring
# ensure_haplotype_diversity's own never-discard guarantee.
.rank_within_group_by_relatedness <- function(members_by_score, G, n_take, target_degree) {
  n_avail <- length(members_by_score)
  n_take  <- min(n_take, n_avail)
  if (n_take < 1L)
    return(list(picked = character(0), ceiling = NA_real_,
               gain_end = NA_real_, diversity_end = NA_real_))
  if (n_take < 2L || is.null(G))
    return(list(picked = members_by_score[seq_len(n_take)], ceiling = NA_real_,
               gain_end = NA_real_, diversity_end = NA_real_))

  Gs <- G[members_by_score, members_by_score, drop = FALSE]
  gain_ids <- members_by_score[seq_len(n_take)]
  gain_end <- .mean_pairwise_relationship(Gs, gain_ids)

  div_res <- select_core_collection(Gs, n_core = n_take, type = "relationship",
                                    strategy = "maximin", verbose = FALSE)
  diversity_end <- .mean_pairwise_relationship(Gs, div_res$selected)

  span    <- max(0, gain_end - diversity_end)
  ceiling <- gain_end - (target_degree / 90) * span

  picked   <- character(0)
  deferred <- character(0)
  for (cand in members_by_score) {
    if (length(picked) >= n_take) break
    trial <- c(picked, cand)
    rel <- if (length(trial) >= 2L) .mean_pairwise_relationship(Gs, trial) else 0
    if (length(picked) < 1L || rel <= ceiling) {
      picked <- trial
    } else {
      deferred <- c(deferred, cand)
    }
  }
  if (length(picked) < n_take) {
    fill <- setdiff(deferred, picked)
    need <- n_take - length(picked)
    picked <- c(picked, fill[seq_len(min(need, length(fill)))])
  }
  list(picked = picked, ceiling = ceiling, gain_end = gain_end,
       diversity_end = diversity_end)
}


#' Family- or Genetic-Cluster-Quota Parent Selection (Best Groups, Then Best Lines Within Them)
#'
#' @description
#' A third parent-shortlist strategy, alongside \code{\link{truncation_selection}}
#' (rank every candidate by a single score) and \code{\link{select_parents_ga}}
#' (search for the set jointly covering the most target-block value). This
#' function instead mirrors the shape a breeding program's shortlist usually
#' already takes in practice: pick the \code{n_families} best-performing
#' groups first, then the \code{n_per_family} best individual lines within
#' each chosen group. A "group" is either a pedigree/cross-ID label you
#' supply (\code{group_by = "family"}, the default) or a genetically
#' data-derived cluster built directly from a relationship matrix
#' (\code{group_by = "genetic_cluster"}) -- see \emph{Grouping: pedigree
#' family vs. genetic cluster} below. Unlike \code{\link{truncation_selection}},
#' which ranks the whole candidate pool as one list and can, by chance, draw
#' heavily from a small number of related groups, this function makes group
#' representation an explicit part of the selection rule itself -- every
#' chosen group contributes exactly \code{n_per_family} lines (or fewer, if
#' it has fewer eligible members), regardless of how its individual lines
#' would have ranked in a single population-wide list.
#'
#' @details
#' Groups are ranked by a (by default, shrinkage-corrected) function of each
#' group's own top \code{rank_k} members (\code{rank_k} defaults to
#' \code{n_per_family} if not supplied separately) -- not by a single best
#' individual (which would favour one-hit-wonder groups) and not by the
#' group's whole-membership mean (which would penalise a large group for
#' having a long tail below its own best lines). See \emph{Shrinkage-corrected
#' family ranking} below for exactly what "shrinkage-corrected" means and
#' why. The \code{n_families} highest-ranked groups by this criterion are
#' kept, and within each, \code{n_per_family} members are selected by
#' \code{score} (subject to the optional relatedness and haplotype-diversity
#' adjustments below).
#'
#' @section Grouping: pedigree family vs. genetic cluster:
#' \code{group_by = "family"} (default) uses your own \code{family} labels
#' exactly as before. \code{group_by = "genetic_cluster"} instead builds
#' groups directly from \code{G} (a relationship matrix): hierarchical
#' clustering (Ward's minimum-variance linkage, \code{cluster_method =
#' "ward.D2"} by default -- the standard, well-justified choice for compact,
#' genetically coherent clusters) on the distance matrix implied by
#' \code{G} (the same exact identity \code{\link{select_core_collection}}
#' uses, not an approximation), cut into \code{n_clusters} groups. This
#' matters because pedigree labels are not always a faithful proxy for
#' actual relatedness -- incomplete records, mislabelling, or loosely
#' applied family IDs can all decouple "same family label" from "actually
#' closely related," and two DIFFERENT family labels can still correspond to
#' closely related individuals (e.g. a shared grandparent). Genetic
#' clustering sidesteps that by building groups directly from the genomic
#' data itself.
#'
#' Both labels are reported per selected individual whenever available,
#' independent of which one actually drove the selection decision (echoed
#' in \code{$group_by}): \code{by_family$family} holds your pedigree label
#' (from the \code{family} argument, \code{NA} if not supplied -- supplying
#' it alongside \code{group_by = "genetic_cluster"} is purely for
#' cross-referencing, it does not affect the selection rule), and
#' \code{by_family$genetic_group} holds the auto-generated genetic cluster
#' label (\code{"cluster_1"}, \code{"cluster_2"}, ..., relabelled by
#' decreasing mean \code{score} so \code{"cluster_1"} is consistently the
#' strongest-performing genetic group) -- populated whenever clustering was
#' computed, which is always true under \code{group_by = "genetic_cluster"}
#' and also true under \code{group_by = "family"} if you additionally supply
#' \code{G} and \code{n_clusters} purely as a diagnostic cross-check (no
#' effect on which lines are chosen in that case). \code{family_ranking}'s
#' own \code{family} column always holds whichever grouping was actually
#' ranked/selected on (pedigree labels or genetic-cluster labels, per
#' \code{group_by}) -- it does not attempt a group-level crosswalk, since a
#' pedigree family can span multiple genetic clusters and vice versa; the
#' individual-level crosswalk in \code{by_family} is where that
#' correspondence is meaningful.
#'
#' @section Shrinkage-corrected family ranking:
#' A group's raw top-\code{rank_k} mean score carries two different,
#' independently-corrected biases rather than being trusted at face value
#' regardless of group size.
#'
#' \strong{Reliability shrinkage.} A small group's raw top-\code{rank_k}
#' mean is an unreliable estimator of that group's true quality -- exactly
#' the kind of small-sample noise classical quantitative genetics corrects
#' for via shrinkage / BLUP-style family evaluation.
#' \code{family_rank_method = "shrunk_topk_mean"} (the default) implements
#' this: treating each group's mean score as a random effect
#' (\eqn{\text{score} = \mu + a_f + e}, \eqn{a_f \sim N(0, \tau^2)},
#' \eqn{e \sim N(0, \sigma^2)}), \eqn{\tau^2} (true between-group variance)
#' and \eqn{\sigma^2} (within-group/residual variance) are estimated from
#' EVERY eligible member's raw \code{score} (not from \code{topk_mean}
#' itself -- fitting variance components on an already-selected/truncated
#' statistic would bias the estimate) via \code{variance_method}: either
#' \code{"reml"} (the default -- \code{lme4::lmer(score ~ 1 + (1|family))};
#' more statistically efficient than method-of-moments under unbalanced
#' group sizes, which this function always has to deal with, since uneven
#' family size is the entire reason shrinkage is needed here in the first
#' place) or \code{"anova"} (the classical one-way random-effects ANOVA
#' method-of-moments estimator, Searle, Casella & McCulloch 1992 -- also
#' the automatic, non-silent fallback whenever \code{"reml"} is requested
#' but \pkg{lme4} is not installed, there isn't enough structure to fit the
#' model, or the fit fails). Each group's (optionally bias-corrected -- see
#' below) \code{topk_mean} is then shrunk toward the across-group mean by
#' the resulting empirical-Bayes/BLUP reliability weight
#' \eqn{w_f = n_f \tau^2 / (n_f \tau^2 + \sigma^2)} (e.g. Robinson 1991,
#' "That BLUP Is a Good Thing") -- a group of 1 member is pulled hard
#' toward the pack; a 20-member group's genuine top-k average is trusted
#' almost fully -- UNLESS genomic-relationship-informed shrinkage (below)
#' is available, in which case it replaces this scalar formula.
#' \code{family_rank_method = "topk_mean"} restores the historical,
#' unshrunk ranking (\code{variance_method}, \code{bias_correction} and
#' \code{use_family_relationship} are all ignored in that case). Honest
#' limit: this shrinks the \emph{selection statistic} (\code{topk_mean})
#' using reliability weights estimated from the \emph{full} within-group
#' distribution -- a well-justified, practical empirical-Bayes approach,
#' not a single unified formal model of \code{topk_mean} itself. Falls back
#' to unshrunk ranking with a message, not silently, whenever variance
#' components cannot be estimated at all (e.g. every eligible group has
#' exactly 1 member).
#'
#' \strong{Finite-population selection-bias correction.} Separately from
#' reliability, "mean of a group's own best \code{rank_k}" is a
#' systematically \emph{upward}-biased estimator of that group's true mean,
#' purely from having picked winners out of a finite pool -- present even
#' for a group whose individual scores carry zero estimation error of their
#' own, and shrinking alone does not remove it. \code{bias_correction =
#' "order_stats"} (the default) removes this before shrinkage is applied,
#' by subtracting each group's expected order-statistics selection
#' differential (\eqn{\sqrt{\sigma^2} \times} the expected mean of the top
#' \code{rank_k} order statistics of \code{n_members} i.i.d. standard normal
#' draws, computed by direct numerical integration of the order statistic's
#' own density rather than a tabulated closed-form approximation).
#' \code{bias_correction = "none"} restores the uncorrected \code{topk_mean}.
#' Requires a finite, positive \eqn{\sigma^2} from the variance-component
#' step above; has no effect (bias = 0) otherwise.
#'
#' \strong{Genomic-relationship-informed shrinkage.} By default
#' (\code{use_family_relationship = TRUE}), whenever \code{G} is also
#' supplied and covers every eligible group, the plain scalar \eqn{w_f}
#' shrinkage above is replaced by a full GBLUP-style family-effect BLUP,
#' solved via Henderson's mixed-model equations from a family-level
#' relationship matrix (each entry the mean pairwise \code{G} value between
#' two groups' members). This lets genetically related groups borrow
#' strength from EACH OTHER, not only from the grand mean -- a strict
#' generalisation of the scalar \eqn{w_f} formula, which it reduces to
#' exactly when groups are mutually unrelated. Falls back to the plain
#' i.i.d. \eqn{w_f} formula (message, not silently) whenever \code{G} does
#' not cover every eligible group, or \code{use_family_relationship =
#' FALSE}.
#'
#' @section family_select_mode: four ways to decide a family's take:
#' Once the \code{n_families} best-ranked groups are chosen, four
#' independent, mutually exclusive modes decide how many/which of a
#' group's own members are actually taken:
#' \describe{
#'   \item{\code{"count"} (default)}{\code{n_per_family} lines per group --
#'     a fixed number, identical to the original design.}
#'   \item{\code{"percentage"}}{\code{pct_per_family} percent of THAT
#'     group's own eligible size (rounded UP via \code{ceiling()}, so an
#'     included group never contributes zero -- for small groups this means
#'     the effective share taken can be noticeably higher than the nominal
#'     percentage, an unavoidable discretisation effect of rounding a
#'     fraction of a small integer up rather than down).}
#'   \item{\code{"sd_threshold"}}{Every group member whose \code{score} is
#'     at least \code{sd_threshold} standard deviations above (or, if
#'     negative, below) the MEAN of every eligible candidate in this run
#'     (population-referenced, not group-referenced) is taken -- however
#'     many that turns out to be.}
#'   \item{\code{"check_relative"}}{Every group member whose \code{score}
#'     is at least \code{check_margin_pct} percent above a reference
#'     check's score is taken -- e.g. \code{check_margin_pct = 10} keeps
#'     everyone at least 10\% above the check, including lines that clear
#'     it by 20\% or 30\%. The check's own score is either looked up via
#'     \code{check_id} (a candidate already present in \code{score}) or
#'     supplied directly as a fixed \code{check_value} (for a check that
#'     is not itself a scored candidate in this run) -- exactly one of the
#'     two is required.}
#' }
#' \code{"count"} and \code{"percentage"} are QUOTA rules: they always
#' produce a fixed number of lines (once resolved) and so keep a
#' family-size eligibility rule (below) that excludes a group outright
#' before ranking whenever taking that many/that share would mean taking
#' every member regardless of merit. \code{"sd_threshold"} and
#' \code{"check_relative"} are THRESHOLD rules instead: how many members
#' clear a merit bar is genuine information, not a degenerate case tied to
#' group size, so there is no size-based pre-exclusion for them -- but a
#' chosen group can then legitimately contribute zero members if nobody
#' clears the bar. When that happens, that group is dropped (reported in
#' \code{$zero_selected_groups}, with a message unless \code{verbose =
#' FALSE}) and the next-ranked eligible group is tried instead, so the
#' final shortlist still ends up with \code{n_families} worth of
#' CONTRIBUTING groups (auto-backfill) whenever enough ranked groups exist.
#'
#' @section Family-size eligibility rule (\code{"count"}/\code{"percentage"} only):
#' A family/group is excluded entirely -- before ranking, before it can
#' count toward \code{n_families} -- whenever its eligible membership does
#' not exceed the number of lines that would actually be taken from it:
#' under \code{"count"}, \code{n_per_family} itself when it is a single
#' scalar (the common, flat-quota case), or \code{rank_k} when
#' \code{n_per_family} is a named vector for uneven per-group quotas (since
#' \code{n_per_family} is then only defined for whichever groups end up
#' chosen, which is not yet known at exclusion time -- \code{rank_k}, "the
#' number of lines this group's ranking is based on," is used as the
#' practical stand-in); under \code{"percentage"}, that group's OWN
#' \code{ceiling(pct_per_family/100 * group_size)}. Below that size there
#' is no genuine "select the best of" decision for that group at all --
#' every eligible member would be taken regardless of ranking -- so
#' including it would let a group that was never really competing count
#' toward \code{n_families} for free, and would feed a \code{topk_mean}
#' with zero real selection content into the shrinkage/bias-correction
#' machinery above (whose entire premise is that some winnowing happened).
#' Under \code{"percentage"}, this is considerably narrower than under
#' \code{"count"}: the exclusion condition reduces to \code{pct_per_family
#' > 100*(n-1)/n} for a group of size \code{n}, a threshold that rises
#' toward 100 as \code{n} grows. At a realistic rate such as 10\%, only
#' single-member groups (\code{n = 1}) are ever excluded -- any group of 2
#' or more is untouched. It only starts excluding larger groups too at
#' high percentages (e.g. 90\% excludes every group of size 9 or smaller).
#' Excluded family/group labels are reported, not silently dropped, in
#' \code{$excluded_groups}; an informational message names them (unless
#' \code{verbose = FALSE}). An error is raised if nothing remains once
#' undersized groups are excluded. Does not apply under \code{"sd_threshold"}/
#' \code{"check_relative"} -- see \emph{family_select_mode} above.
#'
#' @section Relatedness control (optional):
#' Whenever \code{G} is supplied, the realised mean off-diagonal pairwise
#' relationship of the final \code{selected} set (\code{$mean_relationship})
#' and of each selected group's own picks (\code{family_ranking}'s
#' \code{mean_relationship} column) is always reported -- present regardless
#' of whether anything was asked to optimise for it, mirroring
#' \code{\link{select_parents_ga}}'s own \code{$mean_relationship}. On top of
#' that, \code{within_group_target_degree} (numeric in \code{[0, 90]}, same
#' convention as \code{\link{select_parents_ga}}/\code{\link{select_parents_ocs}}:
#' \code{0} = max gain, prioritising \code{score} and accepting more
#' relatedness; \code{90} = max diversity, minimising relatedness) turns
#' this into an actual dial: filling a group's \code{n_per_family} slots
#' stops being pure score-order truncation and instead balances score
#' against a relatedness ceiling, interpolated between the group's own
#' top-scoring feasible subset and \code{\link{select_core_collection}
#' (strategy = "maximin")} restricted to that group's own members (a proven
#' 2-approximation, Gonzalez 1985). Never drops a group's best remaining
#' member purely to satisfy the relatedness preference -- if the quota
#' cannot be filled without a ceiling violation, the remaining slots are
#' filled from score order anyway (see \code{\link{select_parents_ga}} for
#' the analogous guarantee in its own \code{target_degree} documentation).
#' When \code{ensure_haplotype_diversity} is ALSO active,
#' relatedness filtering is layered on top of the coverage/diversity-adjusted
#' preference order, in that disclosed order -- the two constraints are not
#' jointly co-optimised, consistent with this function's overall greedy,
#' honestly-scoped design (see \emph{Haplotype/coverage diversity adjustment}
#' below).
#'
#' @section Haplotype/coverage diversity adjustment (optional):
#' With \code{ensure_haplotype_diversity = FALSE} (the default), each chosen
#' group's \code{n_per_family} slots are filled by \code{score} alone.  With
#' it set to \code{TRUE}, selection proceeds group by group in rank order
#' (best group first), slot by slot within each (best-scoring remaining
#' member first), maintaining a registry of what has already been claimed by
#' earlier picks across \emph{all} groups processed so far. Two methods:
#' \describe{
#'   \item{\code{diversity_method = "coverage_gain"} (default)}{Reuses
#'     \code{\link{select_parents_ga}}'s own block-coverage machinery
#'     (\code{.block_best_values()}) to ask whether a candidate's addition
#'     actually increases total value-weighted coverage across every target
#'     block in \code{value_matrix}, given everyone already claimed so far
#'     -- not merely whether one single block collides. The first remaining
#'     member (in descending \code{score} order) whose addition clears a
#'     near-zero gain threshold is taken; if every remaining member is
#'     already fully redundant with what is claimed, the highest-scoring
#'     one is taken anyway and flagged \code{collision = TRUE}.}
#'   \item{\code{diversity_method = "dominant_block"}}{The original,
#'     coarser heuristic: each individual's \emph{dominant} target block --
#'     the single target block at which that line's own value in
#'     \code{value_matrix} is highest. The first remaining member (score
#'     order) whose dominant block is not already claimed is taken; if
#'     every remaining member's dominant block is claimed and
#'     \code{haplotypes} was supplied, an allele-level check against the
#'     specific claim(s) at that block can still accept a candidate whose
#'     actual allele differs; otherwise the highest-scoring remaining member
#'     is taken anyway, flagged \code{collision = TRUE}.}
#' }
#' Both are greedy, group-by-group heuristics, not a joint search over the
#' whole shortlist the way \code{\link{select_parents_ga}}'s fitness function
#' is -- intentionally simpler, at the cost of not necessarily finding the
#' globally best assignment of lines to slots.
#'
#' @param score Named numeric vector, e.g. whole-genome GEBV
#'   (\code{run_haplotype_prediction()$gebv}) or any other selection index.
#'   Names are individual IDs.
#' @param family Named character or factor vector giving each individual's
#'   pedigree/family group. Required when \code{group_by = "family"} (used
#'   as the actual grouping). Optional when \code{group_by =
#'   "genetic_cluster"} (purely echoed in \code{by_family$family} for
#'   cross-referencing; does not affect which lines are selected in that
#'   mode). Names must cover every name in \code{score} to be considered
#'   (individuals in \code{family} but not \code{score}, or vice versa, are
#'   simply not matched).
#' @param group_by One of \code{"family"} (default) or
#'   \code{"genetic_cluster"}. See \emph{Grouping: pedigree family vs.
#'   genetic cluster}.
#' @param n_families Integer. Number of top-ranked groups to select from.
#'   If fewer than \code{n_families} distinct groups have at least one
#'   eligible member, all of them are used and a warning is issued.
#' @param n_per_family Integer, or a NAMED integer vector (names = chosen
#'   family/group labels) for uneven quotas. Required, and only used, when
#'   \code{family_select_mode = "count"} (the default). Number of best
#'   lines to select from each chosen group. If a chosen group has fewer
#'   than its quota of eligible members, all of its members are taken and
#'   a warning is issued. Also the default ranking sample size
#'   (\code{rank_k}) when \code{rank_k} is not supplied separately and
#'   \code{family_select_mode = "count"} -- see \emph{Details}.
#' @param family_select_mode One of \code{"count"} (default),
#'   \code{"percentage"}, \code{"sd_threshold"}, or \code{"check_relative"}.
#'   See \emph{family_select_mode: four ways to decide a family's take}.
#' @param pct_per_family Numeric in \code{(0, 100]}. Required, and only
#'   used, when \code{family_select_mode = "percentage"}.
#' @param sd_threshold Numeric (signed). Required, and only used, when
#'   \code{family_select_mode = "sd_threshold"}.
#' @param check_id Character, a single individual ID present in
#'   \code{score}. Used, and mutually exclusive with \code{check_value},
#'   when \code{family_select_mode = "check_relative"}.
#' @param check_value Numeric, a single fixed benchmark score. Used, and
#'   mutually exclusive with \code{check_id}, when
#'   \code{family_select_mode = "check_relative"} (for a check that is not
#'   itself a scored candidate in this run).
#' @param check_margin_pct Numeric (signed), e.g. \code{10} for "at least
#'   10\% better than the check". Required, and only used, when
#'   \code{family_select_mode = "check_relative"}.
#' @param rank_k Integer, or \code{NULL} (default). The number of each
#'   group's own top members its ranking criterion (\code{topk_mean}/
#'   \code{shrunk_topk_mean}) is computed from, independent of the actual
#'   take (\code{n_per_family} or otherwise, per \code{family_select_mode}).
#'   Defaults to \code{n_per_family} (or its maximum, if \code{n_per_family}
#'   is a vector) when \code{family_select_mode = "count"} and not supplied
#'   separately -- the historical behaviour; defaults to \code{1} for the
#'   other three modes, which have no take-quota to borrow a default from.
#' @param family_rank_method One of \code{"shrunk_topk_mean"} (default) or
#'   \code{"topk_mean"}. See \emph{Shrinkage-corrected family ranking}.
#' @param variance_method One of \code{"reml"} (default) or \code{"anova"}.
#'   How \eqn{\tau^2}/\eqn{\sigma^2} are estimated for shrinkage. Ignored
#'   under \code{family_rank_method = "topk_mean"}. See \emph{Shrinkage-
#'   corrected family ranking}.
#' @param bias_correction One of \code{"order_stats"} (default) or
#'   \code{"none"}. Whether each group's \code{topk_mean} is corrected for
#'   finite-population selection-differential bias before shrinkage.
#'   Ignored under \code{family_rank_method = "topk_mean"}. See
#'   \emph{Shrinkage-corrected family ranking}.
#' @param use_family_relationship Logical, default \code{TRUE}. Whether
#'   plain i.i.d. shrinkage is replaced by genomic-relationship-informed
#'   (GBLUP-style) family-effect shrinkage when \code{G} is supplied and
#'   covers every eligible group. Ignored under \code{family_rank_method =
#'   "topk_mean"} or when \code{G} is \code{NULL}. See \emph{Shrinkage-
#'   corrected family ranking}.
#' @param min_sel_value,min_sel_mode Optional merit floor applied to
#'   \code{score} \emph{before} group ranking or within-group selection,
#'   exactly as in \code{\link{truncation_selection}}/\code{\link{select_parents_ga}}.
#'   Default \code{min_sel_value = NULL} applies no floor.
#' @param ensure_haplotype_diversity Logical, default \code{FALSE}. See
#'   \emph{Haplotype/coverage diversity adjustment (optional)}. Requires
#'   \code{value_matrix}.
#' @param diversity_method One of \code{"coverage_gain"} (default) or
#'   \code{"dominant_block"}. See \emph{Haplotype/coverage diversity
#'   adjustment (optional)}. Ignored unless \code{ensure_haplotype_diversity
#'   = TRUE}.
#' @param value_matrix Numeric matrix (individuals x target blocks), e.g.
#'   local GEBV (\code{run_haplotype_prediction()$local_gebv}, typically
#'   restricted to \code{\link{select_top_blocks}}'s output first) or
#'   haplotype allele dosage -- the same kind of input
#'   \code{\link{select_parents_ga}} takes. Row names = individual IDs, must
#'   cover every individual under consideration once \code{min_sel_value}
#'   and group membership have been applied. Required when
#'   \code{ensure_haplotype_diversity = TRUE}; ignored (may be left
#'   \code{NULL}) otherwise.
#' @param haplotypes Optional named list, the direct, unmodified return
#'   value of \code{\link{extract_haplotypes}} (one element per block, each
#'   a named character vector of one dosage/phased allele string per
#'   individual). Only used by \code{diversity_method = "dominant_block"}'s
#'   allele-level collision check -- see \emph{Haplotype/coverage diversity
#'   adjustment (optional)}. Ignored under \code{diversity_method =
#'   "coverage_gain"}.
#' @param block_weights Numeric vector, length \code{ncol(value_matrix)}, or
#'   \code{NULL} (default: equal weight 1 for every block). Used by
#'   \code{diversity_method = "coverage_gain"} exactly as
#'   \code{\link{select_parents_ga}}'s argument of the same name is.
#' @param G Relationship/kinship matrix (n x n, dimnames = individual IDs),
#'   e.g. \code{run_haplotype_prediction()$G} or
#'   \code{\link{compute_haplotype_grm}} output, or \code{NULL} (default).
#'   Required when \code{group_by = "genetic_cluster"} or
#'   \code{within_group_target_degree} is supplied; optional otherwise (in
#'   which case supplying it still unlocks the \code{$mean_relationship}
#'   diagnostics -- see \emph{Relatedness control (optional)} -- and,
#'   together with \code{n_clusters}, the optional \code{genetic_group}
#'   cross-reference under \code{group_by = "family"}).
#' @param within_group_target_degree Numeric in \code{[0, 90]}, or
#'   \code{NULL} (default -- no effect). See \emph{Relatedness control
#'   (optional)}. Requires \code{G}.
#' @param n_clusters Integer. Number of genetic clusters to create. Required
#'   when \code{group_by = "genetic_cluster"}; also enables the optional
#'   \code{genetic_group} diagnostic cross-reference under \code{group_by =
#'   "family"} when supplied alongside \code{G}.
#' @param cluster_method Character, default \code{"ward.D2"}. Linkage method
#'   passed to \code{\link[stats]{hclust}} when building genetic clusters.
#' @param verbose Logical, default \code{TRUE}. Print informational messages
#'   (shrinkage fallback, skipped diagnostics).
#'
#' @return Named list:
#' \describe{
#'   \item{\code{selected}}{Character vector of all selected individual IDs,
#'     ordered by group rank, then by selection order within each group.}
#'   \item{\code{by_family}}{Data frame, one row per selected individual:
#'     \code{family} (pedigree label, \code{NA} if not available -- see
#'     \emph{Grouping}), \code{genetic_group} (auto-generated cluster label,
#'     \code{NA} if not computed), \code{family_rank}, \code{individual},
#'     \code{score}, \code{rank_within_family}, and, when
#'     \code{ensure_haplotype_diversity = TRUE}, \code{dominant_block} (when
#'     \code{diversity_method = "dominant_block"}) and \code{collision}
#'     (logical -- see \emph{Haplotype/coverage diversity adjustment}).}
#'   \item{\code{family_ranking}}{Data frame, one row per group considered
#'     (after the family-size eligibility rule has removed undersized
#'     groups -- see \emph{Family-size eligibility rule}): \code{family}
#'     (whichever grouping was ranked -- pedigree or genetic cluster, per
#'     \code{group_by}), \code{topk_mean} (raw ranking criterion),
#'     \code{topk_mean_bias}/\code{topk_mean_corrected} (present when
#'     shrinkage was applied -- the finite-population selection-bias
#'     estimate subtracted, and the corrected value, respectively; bias is
#'     0 and corrected equals raw whenever \code{bias_correction = "none"}
#'     or no correction could be computed), \code{shrinkage_weight}
#'     (\code{NA} unless \code{family_rank_method = "shrunk_topk_mean"},
#'     shrinkage was actually applied, AND it used the plain scalar i.i.d.
#'     formula rather than genomic-relationship-informed shrinkage, which
#'     has no single scalar weight), \code{rank_score} (the value groups
#'     were actually ranked on), \code{n_members}, \code{rank},
#'     \code{selected} (logical), and \code{mean_relationship} (that
#'     group's own selected picks' realised mean relationship, \code{NA}
#'     unless \code{G} was supplied or the group was not selected).}
#'   \item{\code{cutoff}}{Numeric. The \code{min_sel_value} cutoff actually
#'     applied (\code{-Inf} when \code{min_sel_value = NULL}).}
#'   \item{\code{n_families}, \code{n_per_family}, \code{family_select_mode},
#'     \code{pct_per_family}, \code{sd_threshold}, \code{check_id},
#'     \code{check_value}, \code{check_margin_pct}, \code{rank_k},
#'     \code{family_rank_method}, \code{variance_method},
#'     \code{bias_correction}, \code{use_family_relationship},
#'     \code{group_by}, \code{n_clusters}, \code{cluster_method},
#'     \code{diversity_method}, \code{within_group_target_degree}}{Echo the
#'     corresponding arguments (resolved values where applicable; \code{NULL}
#'     for whichever of \code{pct_per_family}/\code{sd_threshold}/
#'     \code{check_id}/\code{check_value}/\code{check_margin_pct} were not
#'     relevant to the \code{family_select_mode} actually used).}
#'   \item{\code{variance_method_used}}{Character, \code{"reml"} or
#'     \code{"anova"} -- which estimator actually produced \code{tau2}/
#'     \code{sigma2}, which can differ from \code{variance_method} when
#'     \code{"reml"} was requested but fell back to \code{"anova"} (see
#'     \emph{Shrinkage-corrected family ranking}). \code{NA} when shrinkage
#'     was not applied at all (\code{family_rank_method = "topk_mean"}, or
#'     variance components could not be estimated).}
#'   \item{\code{relationship_informed}}{Logical. Whether the final ranking
#'     actually used genomic-relationship-informed (GBLUP-style) shrinkage
#'     rather than the plain i.i.d. formula (or no shrinkage at all).}
#'   \item{\code{tau2}, \code{sigma2}}{Numeric. The estimated between-group
#'     and within-group/residual variance components underlying shrinkage
#'     (\code{NULL} when shrinkage was not applied).}
#'   \item{\code{excluded_groups}}{Character vector of family/group labels
#'     removed by the family-size eligibility rule before ranking (empty if
#'     none were excluded, and always empty under \code{family_select_mode
#'     \%in\% c("sd_threshold", "check_relative")}, which has no such rule).
#'     See \emph{Family-size eligibility rule}.}
#'   \item{\code{zero_selected_groups}}{Character vector of ranked
#'     family/group labels that were tried but contributed zero qualifying
#'     members and were passed over via auto-backfill (always empty under
#'     \code{family_select_mode \%in\% c("count", "percentage")}, which
#'     cannot produce a zero-member take). See \emph{family_select_mode:
#'     four ways to decide a family's take}.}
#'   \item{\code{ensure_haplotype_diversity}}{Echoes the argument.}
#'   \item{\code{mean_relationship}}{Numeric. Realised mean off-diagonal
#'     pairwise relationship among all of \code{selected}. \code{NA} unless
#'     \code{G} was supplied.}
#' }
#'
#' @section When to reach for this instead of truncation_selection() or select_parents_ga():
#' Reach for this function specifically when your program's real shortlist
#' decision already has the shape "our best few groups, and our best few
#' lines out of each" -- a common, familiar structure in programs organised
#' around discrete biparental or half-sib families, or when you want a
#' data-derived genetic-cluster version of that same structure
#' (\code{group_by = "genetic_cluster"}) rather than trusting pedigree
#' labels alone -- and you want that structure enforced directly by the
#' selection rule, rather than checking group balance only as a diagnostic
#' \emph{after} running \code{\link{truncation_selection}} or
#' \code{\link{select_parents_ga}}.
#'
#' Skip it, and use \code{\link{truncation_selection}} or
#' \code{\link{select_parents_ga}} instead, if your program does not
#' organise around discrete groups in the first place, or if a strict
#' per-group quota is not actually a constraint you want enforced.
#'
#' @seealso \code{\link{truncation_selection}}, \code{\link{select_parents_ga}},
#'   \code{\link{select_core_collection}} (genetic-cluster grouping's
#'   distance conversion and \code{within_group_target_degree}'s diversity
#'   end both reuse it), \code{\link{usefulness_criterion}}
#'
#' @references
#' Searle, S.R., Casella, G. & McCulloch, C.E. (1992). \emph{Variance
#' Components}. Wiley.
#'
#' Robinson, G.K. (1991). That BLUP is a Good Thing: The Estimation of
#' Random Effects. \emph{Statistical Science}, 6(1), 15-32.
#'
#' Gonzalez, T.F. (1985). Clustering to minimize the maximum intercluster
#' distance. \emph{Theoretical Computer Science}, 38, 293-306.
#'
#' Bates, D., Machler, M., Bolker, B. & Walker, S. (2015). Fitting Linear
#' Mixed-Effects Models Using lme4. \emph{Journal of Statistical Software},
#' 67(1), 1-48. (\code{variance_method = "reml"}.)
#'
#' Falconer, D.S. & Mackay, T.F.C. (1996). \emph{Introduction to Quantitative
#' Genetics} (4th ed.). Longman. (Within-family selection --
#' \code{family_select_mode = "percentage"}/\code{"sd_threshold"}.)
#'
#' @examples
#' \dontrun{
#' haps <- extract_haplotypes(geno, snp_info, blocks)   # same call already
#'                                                       # used upstream
#' res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
#' top  <- select_top_blocks(res$block_importance, n = 15)
#' vmat <- res$local_gebv[, top$block_id, drop = FALSE]
#'
#' # Pedigree-family mode, with relatedness control and coverage-gain
#' # diversity adjustment.
#' fam_sel <- select_parents_by_family(
#'   score         = res$gebv,
#'   family        = family_id,      # named vector, same names as res$gebv
#'   n_families    = 6L,
#'   n_per_family  = 3L,
#'   ensure_haplotype_diversity = TRUE,
#'   value_matrix  = vmat,
#'   G             = res$G,
#'   within_group_target_degree = 30
#' )
#' fam_sel$selected
#' fam_sel$by_family
#'
#' # Genetic-cluster mode instead of pedigree family, with the pedigree
#' # label still echoed for cross-referencing.
#' clust_sel <- select_parents_by_family(
#'   score         = res$gebv,
#'   family        = family_id,          # echoed only, not used for grouping
#'   group_by      = "genetic_cluster",
#'   G             = res$G,
#'   n_clusters    = 8L,
#'   n_families    = 6L,
#'   n_per_family  = 3L
#' )
#' clust_sel$by_family[, c("individual", "family", "genetic_group")]
#' }
#'
#' @export
select_parents_by_family <- function(score,
                                     family = NULL,
                                     n_families,
                                     n_per_family = NULL,
                                     family_select_mode = c("count", "percentage",
                                                            "sd_threshold",
                                                            "check_relative"),
                                     pct_per_family = NULL,
                                     sd_threshold = NULL,
                                     check_id = NULL,
                                     check_value = NULL,
                                     check_margin_pct = NULL,
                                     min_sel_value = NULL,
                                     min_sel_mode  = c("value", "percentile",
                                                      "sd_below_mean"),
                                     ensure_haplotype_diversity = FALSE,
                                     value_matrix = NULL,
                                     haplotypes   = NULL,
                                     block_weights = NULL,
                                     G = NULL,
                                     within_group_target_degree = NULL,
                                     group_by = c("family", "genetic_cluster"),
                                     rank_k = NULL,
                                     family_rank_method = c("shrunk_topk_mean",
                                                            "topk_mean"),
                                     variance_method = c("reml", "anova"),
                                     bias_correction = c("order_stats", "none"),
                                     use_family_relationship = TRUE,
                                     diversity_method = c("coverage_gain",
                                                          "dominant_block"),
                                     n_clusters = NULL,
                                     cluster_method = "ward.D2",
                                     verbose = TRUE) {
  if (is.null(names(score)))
    stop("score must be a named numeric vector (names = individual IDs).",
         call. = FALSE)
  group_by            <- match.arg(group_by)
  family_rank_method  <- match.arg(family_rank_method)
  variance_method     <- match.arg(variance_method)
  bias_correction     <- match.arg(bias_correction)
  diversity_method    <- match.arg(diversity_method)
  family_select_mode  <- match.arg(family_select_mode)
  if (!is.logical(use_family_relationship) || length(use_family_relationship) != 1L)
    stop("use_family_relationship must be a single logical value.", call. = FALSE)

  if (group_by == "family") {
    if (is.null(family) || is.null(names(family)))
      stop("family must be a named vector (names = individual IDs, ",
           "matching score) when group_by = \"family\".", call. = FALSE)
  } else {
    if (is.null(G) || is.null(rownames(G)) || is.null(colnames(G)))
      stop("group_by = \"genetic_cluster\" requires G (a dimnamed ",
           "relationship matrix).", call. = FALSE)
    if (is.null(n_clusters))
      stop("group_by = \"genetic_cluster\" requires n_clusters.",
           call. = FALSE)
    if (!is.null(family) && is.null(names(family)))
      stop("family, when supplied alongside group_by = \"genetic_cluster\" ",
           "for cross-referencing, must be a named vector.", call. = FALSE)
  }

  n_families <- as.integer(n_families)
  if (n_families < 1L) stop("n_families must be >= 1.", call. = FALSE)
  if (!is.null(rank_k) && (!is.numeric(rank_k) || rank_k < 1L))
    stop("rank_k must be >= 1, or NULL.", call. = FALSE)

  # -- family_select_mode-specific argument requirements. Exactly one of
  # the four mode-specific input sets below is actually used; the other
  # three are ignored (not validated) regardless of whether they happen to
  # be supplied, mirroring diversity_method/min_sel_mode's own "declare a
  # mode, only that mode's inputs matter" convention elsewhere in this
  # package.
  check_score <- NULL
  if (family_select_mode == "count") {
    if (is.null(n_per_family))
      stop("n_per_family is required when family_select_mode = \"count\" ",
           "(the default).", call. = FALSE)
    if (any(as.integer(n_per_family) < 1L))
      stop("n_per_family must be >= 1.", call. = FALSE)
  } else if (family_select_mode == "percentage") {
    if (is.null(pct_per_family) || !is.numeric(pct_per_family) ||
        length(pct_per_family) != 1L || pct_per_family <= 0 || pct_per_family > 100)
      stop("pct_per_family must be a single numeric value in (0, 100] when ",
           "family_select_mode = \"percentage\".", call. = FALSE)
  } else if (family_select_mode == "sd_threshold") {
    if (is.null(sd_threshold) || !is.numeric(sd_threshold) || length(sd_threshold) != 1L)
      stop("sd_threshold must be a single numeric value (signed number of ",
           "SDs above -- or, if negative, below -- the population mean) ",
           "when family_select_mode = \"sd_threshold\".", call. = FALSE)
  } else {  # "check_relative"
    if (is.null(check_margin_pct) || !is.numeric(check_margin_pct) ||
        length(check_margin_pct) != 1L)
      stop("check_margin_pct must be a single numeric value (e.g. 10 for ",
           "\"at least 10% better than the check\") when family_select_mode ",
           "= \"check_relative\".", call. = FALSE)
    if (is.null(check_id) == is.null(check_value))
      stop("Exactly one of check_id or check_value must be supplied when ",
           "family_select_mode = \"check_relative\" (check_id looks the ",
           "check's own score up from `score`; check_value is a fixed ",
           "benchmark number for a check that is not itself a candidate in ",
           "`score`).", call. = FALSE)
    if (!is.null(check_id)) {
      if (!is.character(check_id) || length(check_id) != 1L || !(check_id %in% names(score)))
        stop("check_id must be a single individual ID present in `score`.",
             call. = FALSE)
      check_score <- unname(score[check_id])
    } else {
      if (!is.numeric(check_value) || length(check_value) != 1L)
        stop("check_value must be a single numeric value.", call. = FALSE)
      check_score <- check_value
    }
  }

  # -- Warn (not silently ignore) when arguments belonging to an INACTIVE
  # family_select_mode are supplied anyway -- e.g. a leftover n_per_family
  # from before switching to family_select_mode = "sd_threshold" would
  # otherwise have no visible effect and no indication it was disregarded.
  # Does not apply to the active mode's own argument(s), already
  # required/validated above.
  if (isTRUE(verbose)) {
    ignored_mode_args <- character(0)
    if (family_select_mode != "count" && !is.null(n_per_family))
      ignored_mode_args <- c(ignored_mode_args, "n_per_family")
    if (family_select_mode != "percentage" && !is.null(pct_per_family))
      ignored_mode_args <- c(ignored_mode_args, "pct_per_family")
    if (family_select_mode != "sd_threshold" && !is.null(sd_threshold))
      ignored_mode_args <- c(ignored_mode_args, "sd_threshold")
    if (family_select_mode != "check_relative") {
      if (!is.null(check_id)) ignored_mode_args <- c(ignored_mode_args, "check_id")
      if (!is.null(check_value)) ignored_mode_args <- c(ignored_mode_args, "check_value")
      if (!is.null(check_margin_pct)) ignored_mode_args <- c(ignored_mode_args, "check_margin_pct")
    }
    if (length(ignored_mode_args))
      message("[select_parents_by_family] family_select_mode = \"",
              family_select_mode, "\" -- ignoring supplied argument(s) not ",
              "used by this mode: ",
              paste(ignored_mode_args, collapse = ", "))
  }

  if (ensure_haplotype_diversity && is.null(value_matrix))
    stop("ensure_haplotype_diversity = TRUE requires value_matrix (an ",
         "individuals x target-blocks matrix -- see ?select_parents_by_family's ",
         "value_matrix argument).", call. = FALSE)

  # -- Same "report, don't hide" treatment for the diversity-adjustment
  # inputs: a leftover value_matrix/haplotypes/block_weights from a
  # previous call with ensure_haplotype_diversity = TRUE would otherwise
  # be silently unused with no indication.
  if (isTRUE(verbose) && !isTRUE(ensure_haplotype_diversity)) {
    ignored_div_args <- character(0)
    if (!is.null(value_matrix)) ignored_div_args <- c(ignored_div_args, "value_matrix")
    if (!is.null(haplotypes)) ignored_div_args <- c(ignored_div_args, "haplotypes")
    if (!is.null(block_weights)) ignored_div_args <- c(ignored_div_args, "block_weights")
    if (length(ignored_div_args))
      message("[select_parents_by_family] ensure_haplotype_diversity = FALSE ",
              "-- ignoring supplied argument(s) that only apply when it is ",
              "TRUE: ", paste(ignored_div_args, collapse = ", "))
  }

  if (!is.null(within_group_target_degree)) {
    if (!is.numeric(within_group_target_degree) ||
        within_group_target_degree < 0 || within_group_target_degree > 90)
      stop("within_group_target_degree must be a single numeric value in ",
           "[0, 90] (same convention as select_parents_ga()/",
           "select_parents_ocs()'s target_degree).", call. = FALSE)
    if (is.null(G))
      stop("within_group_target_degree requires G (a relationship/kinship ",
           "matrix with dimnames covering the eligible candidates).",
           call. = FALSE)
  }

  # -- Merit floor (optional), then restrict to individuals with a score
  # and (when group_by = "family") a family assignment; under
  # group_by = "genetic_cluster", eligibility instead requires coverage by G.
  floor_res <- .apply_merit_floor(score, min_sel_value, min_sel_mode,
                                  label = "candidate")
  cutoff <- floor_res$cutoff

  if (group_by == "family") {
    eligible <- intersect(floor_res$eligible, names(family))
    if (!length(eligible))
      stop("No individuals have both a score (clearing min_sel_value, if ",
           "set) and a family assignment.", call. = FALSE)
  } else {
    eligible <- intersect(floor_res$eligible, rownames(G))
    if (!length(eligible))
      stop("No individuals have both a score (clearing min_sel_value, if ",
           "set) and coverage in G, required for group_by = ",
           "\"genetic_cluster\".", call. = FALSE)
  }

  score <- floor_res$score[eligible]

  # -- Pedigree-family crosswalk (family_out): populated whenever `family`
  # was supplied, regardless of group_by.
  family_out <- NULL
  if (!is.null(family)) {
    fam_full <- if (is.factor(family)) as.character(family) else family
    family_out <- stats::setNames(rep(NA_character_, length(eligible)), eligible)
    have <- intersect(eligible, names(fam_full))
    family_out[have] <- fam_full[have]
  }

  # -- Genetic-cluster crosswalk (genetic_group_out): computed whenever
  # group_by = "genetic_cluster" (mandatory), or as an opt-in diagnostic
  # under group_by = "family" when G + n_clusters are BOTH also supplied.
  genetic_group_out <- NULL
  if (group_by == "genetic_cluster") {
    genetic_group_out <- .cluster_by_relationship(G, eligible, n_clusters,
                                                  cluster_method, score)
  } else if (!is.null(G) && !is.null(n_clusters)) {
    missing_g_diag <- setdiff(eligible, rownames(G))
    if (!length(missing_g_diag)) {
      genetic_group_out <- .cluster_by_relationship(G, eligible, n_clusters,
                                                    cluster_method, score)
    } else if (isTRUE(verbose)) {
      message("[select_parents_by_family] G does not cover every eligible ",
              "candidate -- skipping the optional genetic_group diagnostic.")
    }
  }

  # -- The ACTIVE grouping used for ranking/quota logic from here on.
  active_group <- if (group_by == "family") {
    fg <- if (is.factor(family)) as.character(family) else family
    fg[eligible]
  } else {
    genetic_group_out
  }

  # -- rank_k default: under family_select_mode = "count", n_per_family
  # (its max, if a vector) when not supplied separately -- resolved here
  # (earlier than the rest of this function needs it) because the
  # small-family exclusion rule immediately below needs it for the
  # uneven-quota (named-vector n_per_family) case. The other three modes
  # have no natural "take-quota" to borrow a default from (pct_per_family/
  # sd_threshold/check_margin_pct don't translate into a ranking sample
  # size the way a count does), so rank_k defaults to 1 (rank by each
  # family's own single best member) unless supplied explicitly.
  if (is.null(rank_k)) {
    rank_k <- if (family_select_mode == "count") {
      if (length(n_per_family) == 1L) as.integer(n_per_family) else
        max(as.integer(n_per_family))
    } else {
      1L
    }
  }

  # -- Exclude families/groups whose eligible membership does not exceed
  # the number of lines that would actually be taken from them -- ONLY
  # meaningful under the two QUOTA-style modes ("count", "percentage"),
  # where the take-count is fixed before any merit is considered. Below
  # that size there is no genuine "select the best of" decision for that
  # group at all -- every eligible member is taken regardless of ranking --
  # so including it would let a group that was never really competing
  # count toward n_families, and would feed a topk_mean with zero real
  # selection content into the shrinkage/bias-correction machinery below
  # (whose entire premise is that SOME winnowing happened).
  #
  # Under the two THRESHOLD-style modes ("sd_threshold", "check_relative"),
  # this size-based pre-exclusion does not apply at all: how many of a
  # family's members clear a merit bar is genuine information, not a
  # degenerate case tied to family size (a 1-member family can meaningfully
  # pass or fail a threshold, unlike being forced to take its only member
  # regardless of merit under a fixed quota). Whether a chosen family ends
  # up contributing zero members under these modes is instead handled by
  # backfill further below (see zero_selected_groups).
  grp_sizes <- table(active_group)
  too_small <- character(0)
  excl_threshold_desc <- NULL
  if (family_select_mode == "count") {
    excl_threshold <- if (length(n_per_family) == 1L) as.integer(n_per_family) else rank_k
    too_small <- names(grp_sizes)[grp_sizes <= excl_threshold]
    excl_threshold_desc <- paste0("<= ", excl_threshold, " eligible member(s)")
  } else if (family_select_mode == "percentage") {
    take_f    <- ceiling(pct_per_family / 100 * as.numeric(grp_sizes))
    too_small <- names(grp_sizes)[as.numeric(grp_sizes) <= take_f]
    excl_threshold_desc <- paste0("membership too small for a genuine top-",
                                  pct_per_family, "% take")
  }
  excluded_groups <- character(0)
  if (length(too_small)) {
    excluded_groups <- too_small
    if (isTRUE(verbose))
      message("[select_parents_by_family] Excluding ", length(too_small),
              " family/group(s) with ", excl_threshold_desc,
              " -- no genuine within-group selection is possible below ",
              "that size: ",
              paste(utils::head(too_small, 10), collapse = ", "),
              if (length(too_small) > 10) ", ..." else "")
    keep         <- names(active_group)[!(active_group %in% too_small)]
    eligible     <- intersect(eligible, keep)
    score        <- score[eligible]
    active_group <- active_group[eligible]
  }
  if (!length(eligible))
    stop("Every family/group has <= its selection quota's worth of ",
         "eligible members -- nothing left to rank once groups with no ",
         "genuine within-group selection possible are excluded. Lower ",
         "n_per_family/pct_per_family, or supply more/larger families/",
         "groups.", call. = FALSE)

  # -- Family-level relationship matrix (optional), for
  # genomic-relationship-informed shrinkage -- see
  # .rank_families_gblup()/.aggregate_relationship_to_family(). Falls back
  # to independent/i.i.d. shrinkage (message, not silently) whenever G does
  # not cover every remaining eligible candidate.
  G_fam <- NULL
  if (isTRUE(use_family_relationship) && !is.null(G)) {
    G_fam <- .aggregate_relationship_to_family(G, active_group)
    if (is.null(G_fam) && isTRUE(verbose))
      message("[select_parents_by_family] G does not cover every eligible ",
              "candidate -- skipping genomic-relationship-informed family ",
              "shrinkage (falling back to independent/i.i.d. shrinkage).")
  }

  if (isTRUE(ensure_haplotype_diversity)) {
    if (!is.matrix(value_matrix)) value_matrix <- as.matrix(value_matrix)
    if (is.null(rownames(value_matrix)))
      stop("value_matrix must have row names (candidate individual IDs).",
           call. = FALSE)
    if (is.null(colnames(value_matrix)))
      colnames(value_matrix) <- paste0("block_", seq_len(ncol(value_matrix)))
    missing_v <- setdiff(eligible, rownames(value_matrix))
    if (length(missing_v))
      stop(length(missing_v), " eligible candidate(s) are missing from ",
           "value_matrix, required because ensure_haplotype_diversity = ",
           "TRUE: ", paste(utils::head(missing_v, 10), collapse = ", "),
           if (length(missing_v) > 10) ", ..." else "", call. = FALSE)
    if (diversity_method == "coverage_gain") {
      if (is.null(block_weights)) {
        block_weights <- rep(1, ncol(value_matrix))
      } else if (length(block_weights) != ncol(value_matrix)) {
        stop("block_weights must have length ncol(value_matrix).", call. = FALSE)
      }
    }
  }

  if (!is.null(within_group_target_degree)) {
    missing_g <- setdiff(eligible, rownames(G))
    if (length(missing_g))
      stop(length(missing_g), " eligible candidate(s) are not present in ",
           "G's dimnames, required because within_group_target_degree is ",
           "set: ", paste(utils::head(missing_g, 10), collapse = ", "),
           if (length(missing_g) > 10) ", ..." else "", call. = FALSE)
  }

  # -- Automatically derive the allele-identity matrix from `haplotypes`
  # when supplied and diversity_method = "dominant_block".
  block_haplotypes <- NULL
  if (isTRUE(ensure_haplotype_diversity) && diversity_method == "dominant_block" &&
      !is.null(haplotypes)) {
    block_haplotypes <- .haplotypes_to_matrix(
      haplotypes, blocks_needed = colnames(value_matrix),
      individuals_needed = eligible
    )
  }

  # -- Rank groups (shrinkage-corrected by default -- see Details).
  fam_tab <- .rank_families(score, active_group, k = rank_k,
                            method = family_rank_method,
                            variance_method = variance_method,
                            bias_correction = bias_correction,
                            G_fam = G_fam, verbose = verbose)
  # Captured BEFORE the data.frame subsetting below: `[.data.frame` does not
  # reliably preserve custom attr()-attached attributes across a row
  # reorder, so these are pulled out into plain local variables here rather
  # than re-read off fam_tab later.
  tau2_out                 <- attr(fam_tab, "tau2")
  sigma2_out                <- attr(fam_tab, "sigma2")
  variance_method_used_out  <- attr(fam_tab, "variance_method_used")
  relationship_informed_out <- isTRUE(attr(fam_tab, "relationship_informed"))

  fam_tab <- fam_tab[order(-fam_tab$rank_score), , drop = FALSE]
  fam_tab$rank <- seq_len(nrow(fam_tab))

  if (nrow(fam_tab) < n_families) {
    warning("Only ", nrow(fam_tab), " eligible family/group(s) available; ",
            "using all of them (fewer than n_families = ", n_families,
            " requested).", call. = FALSE)
    n_families <- nrow(fam_tab)
  }

  # -- Choose the families/groups to draw from, and (for the two
  # THRESHOLD-style modes) how many/which of each one's members actually
  # clear the bar. "count"/"percentage" take a fixed top-n_families slate
  # off the ranking (unchanged from before); "sd_threshold"/"check_relative"
  # walk DOWN the ranking and auto-backfill: any ranked family that would
  # contribute zero members (nobody clears the threshold) is skipped, not
  # counted toward n_families, and reported in zero_selected_groups -- the
  # next-ranked family is tried instead, until n_families worth of
  # CONTRIBUTING groups are assembled or the ranking is exhausted.
  zero_selected_groups <- character(0)
  qualifying_members_by_fam <- list()
  quota <- NULL

  if (family_select_mode %in% c("sd_threshold", "check_relative")) {
    threshold_cutoff <- if (family_select_mode == "sd_threshold") {
      mean(score, na.rm = TRUE) + sd_threshold * stats::sd(score, na.rm = TRUE)
    } else {
      check_score * (1 + check_margin_pct / 100)
    }

    chosen_fams <- character(0)
    idx <- 1L
    while (length(chosen_fams) < n_families && idx <= nrow(fam_tab)) {
      cand_fam     <- fam_tab$family[idx]
      cand_members <- names(active_group)[active_group == cand_fam]
      qual         <- cand_members[score[cand_members] >= threshold_cutoff]
      if (length(qual)) {
        chosen_fams <- c(chosen_fams, cand_fam)
        qualifying_members_by_fam[[cand_fam]] <- qual
      } else {
        zero_selected_groups <- c(zero_selected_groups, cand_fam)
      }
      idx <- idx + 1L
    }
    if (!length(chosen_fams))
      stop("No family/group has any member clearing the ", family_select_mode,
           " threshold -- nothing to select. Lower sd_threshold/",
           "check_margin_pct, or supply more/different candidates.",
           call. = FALSE)
    if (isTRUE(verbose) && length(zero_selected_groups))
      message("[select_parents_by_family] ", length(zero_selected_groups),
              " ranked family/group(s) contributed no member clearing the ",
              family_select_mode, " threshold and were skipped (auto-",
              "backfilled from the next-ranked family): ",
              paste(utils::head(zero_selected_groups, 10), collapse = ", "),
              if (length(zero_selected_groups) > 10) ", ..." else "")
    if (length(chosen_fams) < n_families) {
      warning("Only ", length(chosen_fams), " family/group(s) contributed ",
              "at least one qualifying line under family_select_mode = '",
              family_select_mode, "' (", n_families, " requested); using ",
              "all of them.", call. = FALSE)
      n_families <- length(chosen_fams)
    }
    fam_tab$selected <- fam_tab$family %in% chosen_fams
  } else {
    chosen_fams <- fam_tab$family[seq_len(n_families)]
    fam_tab$selected <- fam_tab$family %in% chosen_fams
    quota <- if (family_select_mode == "count") {
      .resolve_n_per_family(n_per_family, chosen_fams)
    } else {  # "percentage"
      stats::setNames(ceiling(pct_per_family / 100 * as.numeric(grp_sizes[chosen_fams])),
                      chosen_fams)
    }
  }

  # -- Within-group selection, group by group in rank order.
  dom_block <- if (isTRUE(ensure_haplotype_diversity) &&
                  diversity_method == "dominant_block")
    .dominant_block(value_matrix) else NULL

  claimed <- list()                    # dominant_block registry (cross-group)
  claimed_gain_registry <- character(0) # coverage_gain registry (cross-group)

  by_family_rows <- vector("list", n_families)

  for (fi in seq_len(n_families)) {
    fam_id  <- chosen_fams[fi]
    if (family_select_mode %in% c("sd_threshold", "check_relative")) {
      members <- qualifying_members_by_fam[[fam_id]]
      members <- members[order(-score[members])]
      n_take  <- length(members)
    } else {
      n_take  <- unname(quota[fam_id])
      members <- names(active_group)[active_group == fam_id]
      members <- members[order(-score[members])]
    }

    if (length(members) < n_take)
      warning("Family/group '", fam_id, "' has only ", length(members),
              " eligible member(s); taking all of them (fewer than ",
              n_take, " requested).", call. = FALSE)

    # -- Pass 1: FULL ordered preference list of this group's members (not
    # just n_take of them), honouring the cross-group coverage/diversity
    # registry when ensure_haplotype_diversity = TRUE -- runs to completion
    # rather than stopping at n_take, so an optional relatedness filter
    # (Pass 2) has the complete, collision-aware priority order to work
    # from.
    ordered <- character(0)
    collision_flag <- logical(0)
    remaining <- members
    while (length(remaining) > 0L) {
      chosen_this_slot <- NA_character_
      is_collision <- FALSE

      if (!isTRUE(ensure_haplotype_diversity)) {
        chosen_this_slot <- remaining[1L]
      } else if (diversity_method == "dominant_block") {
        blk_of_remaining <- dom_block[remaining]
        free_idx <- which(!(blk_of_remaining %in% names(claimed)))
        if (length(free_idx)) {
          chosen_this_slot <- remaining[free_idx[1L]]
        } else if (!is.null(block_haplotypes)) {
          found <- FALSE
          for (cand in remaining) {
            blk <- unname(dom_block[cand])
            allele <- unname(block_haplotypes[cand, blk])
            if (!(allele %in% claimed[[blk]])) {
              chosen_this_slot <- cand; found <- TRUE; break
            }
          }
          if (!found) { chosen_this_slot <- remaining[1L]; is_collision <- TRUE }
        } else {
          chosen_this_slot <- remaining[1L]; is_collision <- TRUE
        }
      } else {
        # diversity_method == "coverage_gain": same score-priority-with-a-
        # veto structure as "dominant_block" above, but the veto criterion
        # is real marginal coverage gain rather than single-block identity.
        found <- FALSE
        for (cand in remaining) {
          gain <- .marginal_coverage_gain(value_matrix, claimed_gain_registry,
                                          cand, block_weights)
          if (gain > sqrt(.Machine$double.eps)) {
            chosen_this_slot <- cand; found <- TRUE; break
          }
        }
        if (!found) { chosen_this_slot <- remaining[1L]; is_collision <- TRUE }
      }

      ordered <- c(ordered, chosen_this_slot)
      collision_flag <- c(collision_flag, is_collision)
      remaining <- setdiff(remaining, chosen_this_slot)

      if (isTRUE(ensure_haplotype_diversity)) {
        if (diversity_method == "dominant_block") {
          blk <- unname(dom_block[chosen_this_slot])
          allele <- if (!is.null(block_haplotypes))
            unname(block_haplotypes[chosen_this_slot, blk]) else ""
          claimed[[blk]] <- unique(c(claimed[[blk]], allele))
        } else {
          claimed_gain_registry <- c(claimed_gain_registry, chosen_this_slot)
        }
      }
    }

    # -- Pass 2 (optional): relatedness-ceiling filter over the full,
    # collision-aware preference order from Pass 1.
    if (!is.null(within_group_target_degree)) {
      rel_res <- .rank_within_group_by_relatedness(ordered, G, n_take,
                                                    within_group_target_degree)
      picked <- rel_res$picked
      collision_flag_final <- collision_flag[match(picked, ordered)]
    } else {
      n_pick <- min(n_take, length(ordered))
      picked <- ordered[seq_len(n_pick)]
      collision_flag_final <- collision_flag[seq_len(n_pick)]
    }

    row <- data.frame(
      family_rank         = fi,
      individual           = picked,
      score                 = unname(score[picked]),
      rank_within_family     = seq_along(picked),
      stringsAsFactors        = FALSE
    )
    row$.active_group <- fam_id
    if (isTRUE(ensure_haplotype_diversity)) {
      if (diversity_method == "dominant_block")
        row$dominant_block <- unname(dom_block[picked])
      row$collision <- collision_flag_final
    }
    by_family_rows[[fi]] <- row
  }

  by_family <- do.call(rbind, by_family_rows)

  # -- Crosswalk columns: `family` always pedigree (NA if unavailable),
  # `genetic_group` always the auto-generated cluster label (NA if not
  # computed) -- independent of which one was active for selection.
  if (group_by == "family") {
    by_family$family <- by_family$.active_group
    by_family$genetic_group <- if (!is.null(genetic_group_out))
      unname(genetic_group_out[by_family$individual]) else NA_character_
  } else {
    by_family$genetic_group <- by_family$.active_group
    by_family$family <- if (!is.null(family_out))
      unname(family_out[by_family$individual]) else NA_character_
  }
  by_family$.active_group <- NULL
  # Reorder columns for readability: family/genetic_group first.
  front <- c("family", "genetic_group", "family_rank", "individual", "score",
            "rank_within_family")
  by_family <- by_family[, c(front, setdiff(names(by_family), front)), drop = FALSE]

  # -- Relatedness diagnostics (unconditional whenever G is available).
  mean_relationship <- NA_real_
  fam_tab$mean_relationship <- NA_real_
  if (!is.null(G)) {
    sel <- by_family$individual
    if (length(sel) >= 2L && all(sel %in% rownames(G)))
      mean_relationship <- .mean_pairwise_relationship(
        G[sel, sel, drop = FALSE], sel)
    # Per-group realised mean relationship, matched via each row's own
    # (pre-crosswalk) active-group id, recovered from by_family_rows
    # directly rather than the now-overwritten by_family$family/
    # genetic_group columns.
    for (i in seq_len(n_families)) {
      fam_id <- chosen_fams[i]
      picks  <- by_family_rows[[i]]$individual
      if (length(picks) >= 2L && all(picks %in% rownames(G))) {
        fam_tab$mean_relationship[fam_tab$family == fam_id] <-
          .mean_pairwise_relationship(G[picks, picks, drop = FALSE], picks)
      }
    }
  }

  list(
    selected                    = by_family$individual,
    by_family                   = by_family,
    family_ranking              = fam_tab,
    cutoff                      = cutoff,
    n_families                  = n_families,
    n_per_family                = n_per_family,
    family_select_mode          = family_select_mode,
    pct_per_family              = pct_per_family,
    sd_threshold                = sd_threshold,
    check_id                    = check_id,
    check_value                 = check_value,
    check_margin_pct            = check_margin_pct,
    rank_k                      = rank_k,
    family_rank_method          = family_rank_method,
    variance_method             = variance_method,
    variance_method_used        = variance_method_used_out,
    bias_correction              = bias_correction,
    use_family_relationship      = isTRUE(use_family_relationship),
    relationship_informed        = relationship_informed_out,
    tau2                        = tau2_out,
    sigma2                      = sigma2_out,
    excluded_groups              = excluded_groups,
    zero_selected_groups        = zero_selected_groups,
    group_by                    = group_by,
    n_clusters                  = n_clusters,
    cluster_method              = cluster_method,
    ensure_haplotype_diversity  = isTRUE(ensure_haplotype_diversity),
    diversity_method            = diversity_method,
    within_group_target_degree  = within_group_target_degree,
    mean_relationship           = mean_relationship
  )
}
