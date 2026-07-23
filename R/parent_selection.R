# ==============================================================================
# parent_selection.R
# Genetic-algorithm founder-parent selection and truncation-selection baseline.
#
# This closes the single biggest capability gap identified against HapSelect
# (see HapBlockR_vs_HapSelect_comparison.md, section 3.1): HapSelect turns
# block-importance/localGEBV output into an actual *combination* of founder
# parents via a GA search (local_gebv_parent_selection() / haplotype_
# parent_selection(), using the GA package), with real crossing-scheme
# constraints (no_selfing/selfing, OHS/OPV/Haploid_OHS). HapBlockR previously
# stopped at ranking individuals (score_favorable_haplotypes(),
# summarize_parent_haplotypes()) -- useful, but a breeder still had to
# manually pick a complementary set.
#
# IMPORTANT CAVEAT: HapSelect's own GA fitness function and exact constraint
# semantics are not available to HapBlockR's author (closed source beyond the
# published docs/diagram). This module implements HapBlockR's OWN fitness
# function, directly inspired by the formula shown in HapSelect's own
# documentation diagram (`sum_j max((localGEBV_j1 + localGEBV_j2) / 2)`,
# i.e. for each top-ranked block, reward the best complementary PAIR of
# founders), generalised to a `value_matrix` that can be local GEBV values,
# haplotype-allele dosages, or any other per-individual-per-block "favourable
# value" score. It is not a byte-for-byte port of HapSelect's algorithm.
#
# Five exported functions:
#
#   truncation_selection()
#     Trivial baseline: top n_founders individuals by a single genome-wide
#     score (whole-genome GEBV, stacking index, etc.). Used standalone and
#     as the "TS" comparison arm for the GA-vs-TS forward simulation in
#     forward_simulation.R.
#
#   select_parents_ga()
#     GA::ga(type = "binary") search over founder-parent combinations
#     maximising coverage of favourable per-block values, subject to a
#     crossing-scheme `strategy`.
#
#   plot_parent_selection_pca()
#     PCA of a genomic/haplotype relationship matrix (or, via
#     `feature_matrix`, the target-block value_matrix GA actually searched),
#     coloured by which selection method (GA, TS, both, neither) each
#     individual belongs to -- mirrors HapSelect's parent-set PCA plot. Only
#     ever uses PC1/PC2, unscaled eigenvectors -- a quick 2-axis visual
#     cross-check, not a rigorous clustering.
#
#   cluster_selection_groups()
#     A more rigorous version of the same question "do these selections
#     represent the population's real genetic structure?": retains as many
#     leading PCs as needed to reach a user-set cumulative-variance
#     threshold (default 95%) -- not just PC1/PC2 -- properly variance-
#     scaled for Euclidean-distance clustering, runs hierarchical (Ward) or
#     k-means clustering on that subspace into `n_clusters` groups, and
#     cross-tabulates an arbitrary number of named selection strategies
#     (GA, TS, OCS, core-collection, or any other ID vector) against
#     cluster membership so representation gaps show up as numbers, not
#     just visual impressions from a 2D scatter.
#
#   plot_selection_clusters()
#     Companion PC1/PC2 plot for cluster_selection_groups() output, coloured
#     by cluster membership.
# ==============================================================================


# -- Internal: eligibility floor on a whole-genome merit score ---------------
# min_sel_value = NULL -> no filtering, every finite-score individual is
# eligible. Otherwise interpreted per min_sel_mode:
#   "value"          absolute cutoff: keep score >= min_sel_value
#   "percentile"     min_sel_value in (0, 1]: keep the top min_sel_value
#                    fraction by score (e.g. 0.6 = top 60%); cutoff computed
#                    via quantile(score, 1 - min_sel_value)
#   "sd_below_mean"  keep score >= mean(score) - min_sel_value * sd(score)
#                    (min_sel_value = number of SDs below the mean)
# Returns list(eligible = character vector of names clearing the floor,
# cutoff = the numeric cutoff actually applied (-Inf if min_sel_value=NULL),
# score = the input score vector restricted to finite values).
.apply_merit_floor <- function(score, min_sel_value, min_sel_mode,
                               label = "candidate") {
  if (is.null(names(score)))
    stop("score must be a named numeric vector (names = individual IDs).",
         call. = FALSE)
  ok    <- is.finite(score)
  score <- score[ok]

  if (is.null(min_sel_value))
    return(list(eligible = names(score), cutoff = -Inf, score = score))

  min_sel_mode <- match.arg(min_sel_mode,
                            c("value", "percentile", "sd_below_mean"))
  cutoff <- switch(
    min_sel_mode,
    "value" = min_sel_value,
    "percentile" = {
      if (min_sel_value <= 0 || min_sel_value > 1)
        stop("min_sel_value must be in (0, 1] when min_sel_mode = ",
             "'percentile' (fraction of candidates to keep from the top).",
             call. = FALSE)
      stats::quantile(score, probs = 1 - min_sel_value, na.rm = TRUE,
                      names = FALSE)
    },
    "sd_below_mean" = mean(score, na.rm = TRUE) -
      min_sel_value * stats::sd(score, na.rm = TRUE)
  )

  eligible <- names(score)[score >= cutoff]
  if (!length(eligible))
    stop("No ", label, "s clear the min_sel_value floor (cutoff = ",
         round(cutoff, 4), ", mode = '", min_sel_mode, "'). Lower ",
         "min_sel_value.", call. = FALSE)

  list(eligible = eligible, cutoff = cutoff, score = score)
}


#' Truncation Selection (Top-n by a Single Score)
#'
#' @description
#' The simplest possible parent-selection rule and the standard baseline
#' breeders compare any smarter method against: rank individuals by a single
#' genome-wide score and keep the top \code{n_founders}. Used here as the
#' "TS" comparison arm against \code{\link{select_parents_ga}} (see
#' \code{\link{ga_vs_ts_simulation}}), and standalone whenever a plain
#' best-GEBV shortlist is all that's needed.
#'
#' @param score Named numeric vector, e.g. whole-genome GEBV
#'   (\code{run_haplotype_prediction()$gebv}) or a stacking index
#'   (\code{score_favorable_haplotypes()$stacking_index}). Names are
#'   individual IDs.
#' @param n_founders Integer. Number of individuals to select.
#' @param min_sel_value Numeric or \code{NULL} (default \code{NULL} = no
#'   floor, every candidate with a finite \code{score} is eligible).
#'   Excludes candidates below a merit floor \emph{before} ranking, so that
#'   a plain top-\code{n_founders} call never has to hit deeper into the
#'   population than your program's own quality bar. Interpreted according
#'   to \code{min_sel_mode}. Named \code{min_sel_value} rather than
#'   \code{min_selection_index} because \code{score} need not be a
#'   selection index -- it can be any single genome-wide value.
#' @param min_sel_mode One of \code{"value"} (default), \code{"percentile"},
#'   \code{"sd_below_mean"}. Only used when \code{min_sel_value} is not
#'   \code{NULL}:
#'   \describe{
#'     \item{\code{"value"}}{\code{min_sel_value} is an absolute cutoff on
#'       \code{score} itself (same units/scale as \code{score} -- requires
#'       knowing that scale in advance).}
#'     \item{\code{"percentile"}}{\code{min_sel_value} in \code{(0, 1]} is
#'       the fraction of candidates to keep from the top, e.g. \code{0.6}
#'       keeps the top 60\% by \code{score}. Self-scaling -- no need to know
#'       \code{score}'s units.}
#'     \item{\code{"sd_below_mean"}}{\code{min_sel_value} is the number of
#'       standard deviations below \code{mean(score)} the cutoff sits, e.g.
#'       \code{1} keeps everyone within/above one SD of the mean. Also
#'       self-scaling.}
#'   }
#'
#' @return Named list:
#' \describe{
#'   \item{\code{selected}}{Character vector of the top \code{n_founders}
#'     individual IDs (after the \code{min_sel_value} floor, if set), sorted
#'     by descending \code{score}.}
#'   \item{\code{score}}{The corresponding scores, same order.}
#'   \item{\code{cutoff}}{Numeric. The merit-floor cutoff actually applied
#'     (\code{-Inf} when \code{min_sel_value = NULL}).}
#' }
#'
#' @section What this does and does not do:
#' This function ranks by a single number and takes the top
#' \code{n_founders} -- nothing more. It has no notion of which haplotype
#' blocks each individual carries, no complementarity check between the
#' individuals it selects (two selected parents could carry identical
#' favourable haplotypes and none of the ones another candidate has), and no
#' relatedness/coancestry management. That simplicity is the point: it is the
#' baseline every smarter method -- above all \code{\link{select_parents_ga}}
#' -- needs to outperform to justify its added complexity. See
#' \code{\link{select_parents_ga}}'s \emph{Choosing between this function and
#' truncation_selection()} section for a full comparison and a decision
#' guide, and the \emph{From Local GEBV to a Crossing Decision} vignette for
#' a worked example running both side by side.
#'
#' @seealso \code{\link{select_parents_ga}}
#'
#' @examples
#' \dontrun{
#' res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
#' ts  <- truncation_selection(res$gebv, n_founders = 20)
#' ts$selected
#' }
#'
#' @export
truncation_selection <- function(score, n_founders,
                                 min_sel_value = NULL,
                                 min_sel_mode  = c("value", "percentile",
                                                  "sd_below_mean")) {
  if (is.null(names(score)))
    stop("score must be a named numeric vector (names = individual IDs).",
         call. = FALSE)
  n_founders <- as.integer(n_founders)
  if (n_founders < 1L)
    stop("n_founders must be >= 1.", call. = FALSE)

  floor_res <- .apply_merit_floor(score, min_sel_value, min_sel_mode,
                                  label = "candidate")
  score <- floor_res$score[floor_res$eligible]

  if (length(score) < n_founders)
    warning("Only ", length(score), " individuals clear the current floor ",
            "(min_sel_value", if (!is.null(min_sel_value))
              paste0(" = ", min_sel_value, ", cutoff = ",
                     round(floor_res$cutoff, 4)) else " = NULL",
            "); returning all of them (fewer than n_founders = ", n_founders,
            " requested).", call. = FALSE)

  ord  <- order(score, decreasing = TRUE)
  keep <- ord[seq_len(min(n_founders, length(score)))]

  list(selected = names(score)[keep], score = unname(score[keep]),
       cutoff   = floor_res$cutoff)
}


# -- Internal: per-block "best achievable value" for a chosen founder subset -
# strategy determines how a block's value is realised from the chosen set:
#   "no_selfing", "OHS"          -> mean of the TWO largest values (distinct
#                                    parents each contribute one copy)
#   "selfing", "Haploid_OHS"     -> the single largest value (a parent can be
#                                    "selfed"/contribute both copies itself,
#                                    which always dominates any distinct pair
#                                    average, since the largest value is always
#                                    >= the average of itself and anything
#                                    smaller)
#   "OPV"                        -> the single largest value (population-value
#                                    framing: any one founder carrying it is
#                                    enough, no pairing needed) -- computed
#                                    identically to the self-allowed case
#                                    above, since both reduce to "is the best
#                                    value in the set achievable by the set".
# NOTE ON THE SELFING/OPV EQUIVALENCE: this is a direct mathematical
# consequence of the max-pair-average formulation (documented in
# select_parents_ga()'s help), not a shortcut taken for convenience --
# allowing self-pairing can never do worse than any distinct pair, since the
# best single value is always >= the mean of itself and a smaller value.
.block_best_values <- function(value_matrix, chosen, strategy) {
  sub <- value_matrix[chosen, , drop = FALSE]
  self_allowed <- strategy %in% c("selfing", "Haploid_OHS", "OPV")
  if (self_allowed || nrow(sub) < 2L) {
    apply(sub, 2L, max)
  } else {
    apply(sub, 2L, function(x) {
      s <- sort(x, decreasing = TRUE)
      mean(s[1:2])
    })
  }
}

# -- Internal: for reporting, which founder(s) realise each block's best value
.block_best_contributors <- function(value_matrix, chosen, strategy) {
  sub <- value_matrix[chosen, , drop = FALSE]
  self_allowed <- strategy %in% c("selfing", "Haploid_OHS", "OPV")
  n_blk <- ncol(sub)
  c1 <- character(n_blk); c2 <- character(n_blk)
  for (j in seq_len(n_blk)) {
    v <- sub[, j]
    ord <- order(v, decreasing = TRUE)
    c1[j] <- rownames(sub)[ord[1L]]
    c2[j] <- if (self_allowed || length(ord) < 2L) c1[j] else rownames(sub)[ord[2L]]
  }
  data.frame(contributor_1 = c1, contributor_2 = c2, stringsAsFactors = FALSE)
}


# -- Internal: mean OFF-DIAGONAL pairwise relationship among a chosen set,
# from a relationship/kinship matrix G (e.g. compute_haplotype_grm() or
# run_haplotype_prediction()$G). Excludes each individual's relationship to
# itself (the diagonal, i.e. 1+F_i under a VanRaden-style GRM) so the penalty
# reflects relatedness BETWEEN chosen parents, not their own inbreeding.
# Returns 0 for a set of size < 2 (no pairs to average).
.mean_pairwise_relationship <- function(G, chosen_ids) {
  if (length(chosen_ids) < 2L) return(0)
  sub <- G[chosen_ids, chosen_ids, drop = FALSE]
  mean(sub[upper.tri(sub)], na.rm = TRUE)
}

# -- Internal: greedy, feasible high-coverage founder group -------------------
# Builds an actual, real founder group of size n_founders by repeatedly
# adding whichever remaining candidate most improves total block coverage
# (sum(.block_best_values(...) * block_weights)), reusing .block_best_values()
# so the chosen `strategy` (no_selfing/selfing/OHS/OPV/Haploid_OHS) is
# respected automatically -- no separate feasibility layer needed, since
# strategy only changes HOW a block's value is computed from a chosen set,
# not which sets are valid.
#
# This is not an arbitrary heuristic: sum(.block_best_values(S)*block_weights)
# is a monotone submodular function of S under a fixed-cardinality constraint
# (adding a candidate to a smaller set can only help at least as much as
# adding it to a superset of that set, since best_j(S) is a max/mean-of-top-2
# that can only go up as S grows), so a plain greedy build is provably within
# a known factor (1 - 1/e, roughly 63%) of the true best-achievable coverage
# for a set of this size -- a real, feasible reference point, not merely a
# fast guess. This is what .calibrate_merit_scale() uses as the coverage
# term's upper reference instead of an unreachable per-block best-value sum
# (see its own comments for why that naive ceiling is biased).
#
# Optional merit_weight/merit_score: when supplied (merit_weight > 0), each
# step's marginal gain also includes merit_weight * mean(merit_score[trial]),
# matching the fitness function's own merit term exactly. This is used by
# .calibrate_relatedness_ceiling() to build target_degree's gain-end
# reference so it reflects what the GA would ACTUALLY converge to with no
# relatedness constraint when merit is also active -- see its own comments.
# NOTE: adding the merit term breaks the strict monotone-submodularity
# argument above (mean() does not have the same diminishing-returns
# property sum() does as the set grows), so with merit_weight > 0 this is a
# practical greedy hill-climb, not a formally guaranteed-near-optimal one --
# still a real, feasible group, just without the (1-1/e) proof. With the
# default merit_weight = 0 this function is byte-for-byte the original,
# provably-near-optimal coverage-only greedy (used as such by
# .calibrate_merit_scale()/suggest_merit_weight()).
#
# Returns a character vector of n_founders individual IDs (or fewer if
# n_cand < n_founders).
.greedy_coverage_group <- function(value_matrix, block_weights, strategy,
                                   n_founders, candidate_ids = NULL,
                                   merit_weight = 0, merit_score = NULL) {
  ids <- if (is.null(candidate_ids)) rownames(value_matrix) else candidate_ids
  n_founders <- min(n_founders, length(ids))
  chosen <- character(0)
  remaining <- ids
  use_merit <- merit_weight > 0 && !is.null(merit_score)
  for (step in seq_len(n_founders)) {
    gains <- vapply(remaining, function(cand) {
      trial <- c(chosen, cand)
      obj <- sum(.block_best_values(value_matrix, trial, strategy) * block_weights)
      if (use_merit) obj <- obj + merit_weight * mean(merit_score[trial])
      obj
    }, numeric(1))
    best_cand <- remaining[which.max(gains)]
    chosen <- c(chosen, best_cand)
    remaining <- setdiff(remaining, best_cand)
  }
  chosen
}

# -- Internal: trimmed low-reference group for a single ranking score ---------
# Sorts `ranking` ascending, drops the `trim` most extreme (lowest) values as
# a defence against a single outlier candidate (e.g. a data-entry error)
# artificially deflating the floor and inflating the estimated span, then
# returns the NEXT `n` IDs above that -- a real, valid group (never an
# unreachable construct), just deliberately not the most extreme possible
# one. `trim` should already be size-adapted by the caller (see
# .calibrate_merit_scale()).
.trimmed_low_group <- function(ranking, n, trim) {
  ord     <- names(sort(ranking))              # ascending
  n_avail <- length(ord)
  n       <- min(n, n_avail)
  if (n < 1L) return(character(0))
  # Shrink trim, rather than n, if there isn't room for both -- always
  # return the requested group size when the pool can supply it at all.
  trim <- max(0L, min(trim, n_avail - n))
  ord[(trim + 1L):(trim + n)]
}

# -- Internal: core scale-estimation shared by suggest_merit_weight() and
# select_parents_ga()'s merit_priority argument. Estimates how big a swing
# the merit term and the block-coverage term can realistically produce for a
# founder group of size n_founders, from THIS data, and returns the ratio
# needed to make merit_priority (0-100) a fair "how much do I care about
# merit vs. coverage" dial -- see ?suggest_merit_weight for the full
# definition and its documented, honestly-stated limits (matching realistic
# best-vs-worst SPREAD is one reasonable definition of "comparable scale",
# not the only possible one; distribution shape beyond that spread is not
# accounted for).
#
# value_matrix/merit_score/block_weights/strategy/n_founders are assumed
# already filtered/aligned to the SAME candidate pool the real search will
# use (the caller's job -- see select_parents_ga()'s merit_priority wiring).
.calibrate_merit_scale <- function(value_matrix, merit_score, block_weights,
                                   strategy, n_founders) {
  if (any(block_weights < 0))
    stop("block_weights must be non-negative for merit/coverage scale ",
         "calibration (suggest_merit_weight()/merit_priority) -- a negative ",
         "weight makes 'coverage span' ill-defined.", call. = FALSE)

  ids <- rownames(value_matrix)
  n_cand <- length(ids)
  # Adaptive trim: drop a small, size-scaled number of extreme low points
  # before taking the "worst realistic group" reference; skip trimming
  # entirely on small pools, where percentile-like trimming isn't meaningful
  # (see ?suggest_merit_weight).
  trim <- if (n_cand < 20L) 0L else max(1L, round(0.05 * n_cand))

  # -- Merit term: ceiling = mean of the actual best n_founders by
  # merit_score (exact, already a mean over multiple individuals so not a
  # single-point outlier risk); floor = mean of a trimmed low-scoring group.
  ms <- merit_score[ids]
  merit_ceiling_ids <- names(sort(ms, decreasing = TRUE))[seq_len(min(n_founders, n_cand))]
  merit_ceiling <- mean(ms[merit_ceiling_ids])
  merit_floor_ids <- .trimmed_low_group(ms, n_founders, trim)
  merit_floor <- mean(ms[merit_floor_ids])
  merit_span <- merit_ceiling - merit_floor

  # -- Coverage term: ceiling/floor/span, shared with target_degree's
  # penalty-scale reference -- see .coverage_span_reference().
  cov <- .coverage_span_reference(value_matrix, block_weights, strategy,
                                  n_founders, trim)

  ok <- is.finite(merit_span) && is.finite(cov$span) &&
    merit_span > sqrt(.Machine$double.eps)

  list(
    merit_ceiling = merit_ceiling, merit_floor = merit_floor,
    merit_span = merit_span,
    coverage_ceiling = cov$ceiling, coverage_floor = cov$floor,
    coverage_span = cov$span,
    scale_factor = if (ok) cov$span / merit_span else NA_real_,
    trim_n = trim, n_candidates = n_cand, ok = ok
  )
}

# -- Internal: coverage-term ceiling/floor/span, shared by
# .calibrate_merit_scale() (merit_priority) and .calibrate_relatedness_ceiling()
# (target_degree)'s penalty-scale reference. Ceiling = greedy-feasible
# high-coverage group (real, achievable, near-optimal by submodularity --
# see .greedy_coverage_group()); floor = a real, trimmed low-coverage group,
# ranked by each candidate's own row-max value (cheapest reasonable
# single-individual summary, same convention already used by
# select_parents_ga()'s top_candidates prefilter). `trim` is already
# size-adapted by the caller (see .calibrate_merit_scale()).
.coverage_span_reference <- function(value_matrix, block_weights, strategy,
                                     n_founders, trim) {
  ids <- rownames(value_matrix)
  ceiling_ids <- .greedy_coverage_group(value_matrix, block_weights, strategy,
                                        n_founders, candidate_ids = ids)
  coverage_ceiling <- sum(.block_best_values(value_matrix, ceiling_ids, strategy) *
                           block_weights)
  row_best <- apply(value_matrix, 1L, max)
  floor_ids <- .trimmed_low_group(row_best, n_founders, trim)
  coverage_floor <- sum(.block_best_values(value_matrix, floor_ids, strategy) *
                        block_weights)
  # Greedy is near-optimal (submodularity), not guaranteed-optimal, so in a
  # pathological edge case the trimmed "low row_best" floor group could in
  # principle score marginally higher than the greedy ceiling group -- clamp
  # at 0 rather than propagate a negative span into a negative scale factor.
  coverage_span <- max(0, coverage_ceiling - coverage_floor)
  list(ceiling = coverage_ceiling, floor = coverage_floor, span = coverage_span,
       ceiling_ids = ceiling_ids, floor_ids = floor_ids)
}


# -- Internal: relatedness-ceiling calibration for target_degree,
# select_parents_ga()'s easy relatedness-vs-gain dial. Mirrors
# select_parents_ocs()'s engine = "optisel" target_degree convention exactly
# (0 = the unconstrained max-gain solution's own relatedness; 90 = the
# minimum relatedness achievable at all, ignoring merit). Both frontier
# endpoints reuse existing, already-proven building blocks rather than a new
# unguaranteed heuristic:
#   - gain end:      .greedy_coverage_group(), extended to also account for
#                     merit_weight when merit is active in this same run (so
#                     degree 0 reflects what the GA would ACTUALLY converge
#                     to with no relatedness constraint, not a coverage-only
#                     proxy) -- exactly matching the "optisel" engine's own
#                     "max.Merit, no kinship constraint" extreme. Its OWN
#                     mean pairwise relationship is then measured directly
#                     (not assumed).
#   - diversity end: select_core_collection(strategy = "maximin"), a proven
#                     2-approximation for maximizing the minimum pairwise
#                     distance (Gonzalez 1985), built ignoring merit
#                     entirely -- exactly matching the "optisel" engine's own
#                     "min.Kin, ignoring merit" extreme. Its own mean
#                     pairwise relationship is the floor.
# Neither endpoint is a guaranteed EXACT optimum for its own criterion (the
# gain-end greedy has no guarantee once merit_weight > 0 breaks
# submodularity; maximin's guarantee covers minimum pairwise distance, not
# directly the mean), so the resulting span is clamped at 0 rather than
# allowed to invert.
#
# The penalty-scale reference (coverage_span) is the SAME quantity
# suggest_merit_weight()/merit_priority already compute from this candidate
# pool, so the caller can scale the relatedness-ceiling penalty to dominate
# any possible coverage/merit gain from violating it, without a new
# user-tunable constant (see .run_ga_once()'s fitness_fn).
.calibrate_relatedness_ceiling <- function(value_matrix, G, block_weights,
                                           strategy, n_founders,
                                           target_degree,
                                           merit_weight = 0,
                                           merit_score = NULL) {
  ids <- rownames(value_matrix)
  Gs  <- G[ids, ids, drop = FALSE]
  n_cand <- length(ids)
  trim <- if (n_cand < 20L) 0L else max(1L, round(0.05 * n_cand))

  gain_ids <- .greedy_coverage_group(value_matrix, block_weights, strategy,
                                     n_founders, candidate_ids = ids,
                                     merit_weight = merit_weight,
                                     merit_score = merit_score)
  gain_end <- .mean_pairwise_relationship(Gs, gain_ids)

  div_res <- select_core_collection(Gs, n_core = min(n_founders, n_cand),
                                    type = "relationship", strategy = "maximin",
                                    verbose = FALSE)
  diversity_end <- .mean_pairwise_relationship(Gs, div_res$selected)

  span <- max(0, gain_end - diversity_end)
  ceiling <- gain_end - (target_degree / 90) * span

  cov <- .coverage_span_reference(value_matrix, block_weights, strategy,
                                  n_founders, trim)

  list(
    gain_end = gain_end, diversity_end = diversity_end, span = span,
    ceiling = ceiling, coverage_span = cov$span,
    gain_ids = gain_ids, diversity_ids = div_res$selected
  )
}


#' Suggest a Starting merit_weight for select_parents_ga()
#'
#' @description
#' \code{\link{select_parents_ga}}'s \code{merit_weight} and
#' \code{coancestry_weight} arguments have no universal correct value: the
#' block-coverage term and the merit term live on different, problem-specific
#' scales, so a raw multiplier that works for one dataset can be meaningless
#' for another. This function estimates a sensible starting point directly
#' from your own data, in two ways: call it with \code{merit_priority} left
#' \code{NULL} to see the raw diagnostic numbers (the realistic spread of
#' each term, and the ratio between them), or supply \code{merit_priority}
#' (0-100, "how much do you care about merit vs. coverage") to also get a
#' literal \code{merit_weight} value ready to pass straight into
#' \code{\link{select_parents_ga}}. This is the same calculation
#' \code{select_parents_ga()}'s own \code{merit_priority} argument uses
#' internally -- calling this function first just lets you see the numbers
#' before committing to them.
#'
#' @details
#' What "comparable scale" means here, precisely: for a founder group of
#' size \code{n_founders}, this function estimates the realistic
#' \strong{best-vs-worst achievable spread} of the block-coverage sum, and
#' the realistic best-vs-worst achievable spread of mean \code{merit_score},
#' \emph{for this specific dataset}. \code{merit_priority = 100} sets
#' \code{merit_weight} so that merit's spread becomes comparable in
#' magnitude to coverage's spread; \code{merit_priority = 0} is identical to
#' \code{merit_weight = 0} (no merit term at all); values in between scale
#' linearly. This is one reasonable, explicitly-stated definition of
#' "comparable" -- not the only possible one (matching standard deviation
#' instead of spread, for instance, would give a different number) -- see
#' \emph{What this does not solve} below.
#'
#' The "best achievable" coverage reference is not a naive per-block sum of
#' each block's own maximum value across all candidates (which is usually
#' \emph{unreachable} by any single real group, since different blocks'
#' maxima often belong to different individuals, and would bias the estimate
#' toward an inflated ceiling). Instead, it is the coverage of an actual,
#' feasible founder group built by a greedy search: repeatedly add whichever
#' remaining candidate most improves total coverage. This is not an
#' arbitrary heuristic -- the coverage function is monotone submodular under
#' a fixed group-size constraint, so greedy is mathematically guaranteed to
#' land within a known factor (\eqn{1 - 1/e}, about 63\%) of the true best
#' achievable coverage, using a real, reachable group. The "worst
#' achievable" references (for both merit and coverage) are built from a
#' \emph{trimmed} low-ranked group -- the bottom \code{n_founders}
#' candidates after first setting aside a small, data-size-scaled number of
#' the most extreme low values -- so that a single outlier candidate (e.g. a
#' data-entry error) cannot single-handedly deflate the floor and distort
#' the estimated spread. Trimming is skipped entirely on small candidate
#' pools (under ~20), where it would not be meaningful.
#'
#' @section What this does not solve:
#' Two limits are inherent to \emph{any} scale-matching approach, not
#' specific to the method used here, and cannot be resolved by more
#' engineering -- they are documented rather than hidden:
#' \describe{
#'   \item{Matching spread is a choice, not a universal truth}{A different,
#'     equally defensible definition of "comparable" (e.g. matching standard
#'     deviation across many realistic groups, rather than the best-vs-worst
#'     achievable span) would produce a different scale factor. Treat
#'     \code{merit_priority}'s suggestion as a well-reasoned starting point
#'     to inspect and adjust, not a uniquely correct answer -- the literal
#'     \code{merit_weight} argument remains available in
#'     \code{\link{select_parents_ga}} for full manual control.}
#'   \item{A single span number does not capture distribution shape}{If
#'     \code{merit_score} or the block-coverage values are unusually shaped
#'     (e.g. strongly bimodal), the dial's practical effect may not feel
#'     perfectly linear across its 0-100 range even though the underlying
#'     calculation is exact for what it measures.}
#' }
#' This calibration also only weighs merit against coverage; if
#' \code{coancestry_weight} is also active in your \code{select_parents_ga()}
#' call, its effect is held fixed rather than jointly recalibrated -- use
#' \code{\link{select_parents_pareto}}'s sweep to explore that trade-off
#' separately, as already recommended for tuning \code{coancestry_weight} on
#' its own.
#'
#' @param value_matrix Numeric matrix (individuals x blocks), identical in
#'   shape/meaning to \code{\link{select_parents_ga}}'s own argument --
#'   ideally the exact same, already-filtered matrix you are about to pass
#'   to that call (after any \code{min_sel_value}/\code{top_candidates}
#'   filtering), so the calibration reflects the real candidate pool the GA
#'   will search.
#' @param merit_score Named numeric vector, whole-genome merit -- same as
#'   \code{\link{select_parents_ga}}'s argument of the same name. Must cover
#'   every individual in \code{value_matrix}.
#' @param n_founders Integer. Same as \code{\link{select_parents_ga}}'s
#'   argument of the same name -- the founder group size to calibrate for.
#' @param strategy One of \code{"no_selfing"} (default), \code{"selfing"},
#'   \code{"OHS"}, \code{"OPV"}, \code{"Haploid_OHS"} -- must match the
#'   \code{strategy} you intend to run \code{select_parents_ga()} with, since
#'   it changes how a block's achievable value is computed.
#' @param block_weights Numeric vector, length \code{ncol(value_matrix)}, or
#'   \code{NULL} (default: equal weight 1) -- same as
#'   \code{\link{select_parents_ga}}'s argument of the same name. Must be
#'   non-negative.
#' @param merit_priority Numeric in \code{[0, 100]}, or \code{NULL} (default).
#'   \code{NULL} returns only the diagnostic spread/scale numbers, with
#'   \code{suggested_merit_weight = NULL}. A number computes
#'   \code{suggested_merit_weight} too -- \code{0} is always equivalent to
#'   \code{merit_weight = 0}; \code{100} sets merit's spread comparable to
#'   coverage's spread; values between scale linearly.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{merit_span}, \code{coverage_span}}{The estimated realistic
#'     best-vs-worst achievable spread of each term for a group of size
#'     \code{n_founders}, on their own native scales.}
#'   \item{\code{merit_ceiling}, \code{merit_floor}, \code{coverage_ceiling},
#'     \code{coverage_floor}}{The four reference values \code{*_span} is
#'     computed from -- inspect these directly if you want to sanity-check
#'     the calibration by hand.}
#'   \item{\code{scale_factor}}{\code{coverage_span / merit_span}. The
#'     multiplier that would make merit's spread exactly equal to coverage's
#'     spread if used as \code{merit_weight} directly (equivalent to
#'     \code{merit_priority = 100}). \code{NA} if \code{merit_span} was too
#'     small to divide by safely (see \code{ok}).}
#'   \item{\code{merit_priority}}{Echoes the argument.}
#'   \item{\code{suggested_merit_weight}}{\code{merit_priority / 100 *
#'     scale_factor}, ready to pass to \code{\link{select_parents_ga}}'s
#'     \code{merit_weight} argument. \code{NULL} if \code{merit_priority}
#'     was \code{NULL}.}
#'   \item{\code{ok}}{Logical. \code{FALSE} if \code{merit_span} was too
#'     close to zero (e.g. \code{merit_score} has almost no spread among
#'     these candidates) to safely divide by -- \code{suggested_merit_weight}
#'     is \code{NULL} in that case regardless of \code{merit_priority}, with
#'     a warning, rather than returning a wild or infinite number.}
#'   \item{\code{trim_n}}{Integer. How many extreme low points were set
#'     aside per term before building the "worst achievable" reference
#'     group (\code{0} on small candidate pools, where trimming is skipped).}
#' }
#'
#' @seealso \code{\link{select_parents_ga}}'s \emph{Merit-weighted fitness
#'   (GA+TS hybrid, optional)} section for the fitness function this feeds
#'   into, and \code{\link{select_parents_pareto}} for exploring the
#'   \code{coancestry_weight} trade-off the same way.
#'
#' @examples
#' \dontrun{
#' res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
#' top  <- select_top_blocks(res$block_importance, n = 15)
#' vmat <- res$local_gebv[, top$block_id, drop = FALSE]
#' cal  <- suggest_merit_weight(vmat, res$gebv, n_founders = 20,
#'                              merit_priority = 50)
#' cal$merit_span; cal$coverage_span; cal$suggested_merit_weight
#' ga_out <- select_parents_ga(vmat, n_founders = 20, merit_score = res$gebv,
#'                             merit_weight = cal$suggested_merit_weight)
#' }
#'
#' @export
suggest_merit_weight <- function(value_matrix, merit_score, n_founders,
                                 strategy = c("no_selfing", "selfing", "OHS",
                                             "OPV", "Haploid_OHS"),
                                 block_weights = NULL, merit_priority = NULL) {
  strategy <- match.arg(strategy)
  if (!is.matrix(value_matrix)) value_matrix <- as.matrix(value_matrix)
  if (is.null(rownames(value_matrix)))
    stop("value_matrix must have row names (candidate individual IDs).",
         call. = FALSE)
  if (is.null(names(merit_score)))
    stop("merit_score must be a named numeric vector (names = individual ",
         "IDs).", call. = FALSE)
  missing_m <- setdiff(rownames(value_matrix), names(merit_score))
  if (length(missing_m))
    stop(length(missing_m), " value_matrix candidate(s) have no merit_score ",
         "entry: ", paste(utils::head(missing_m, 10), collapse = ", "),
         if (length(missing_m) > 10) ", ..." else "", call. = FALSE)

  n_founders <- as.integer(n_founders)
  if (n_founders < 1L) stop("n_founders must be >= 1.", call. = FALSE)
  if (n_founders > nrow(value_matrix))
    stop("n_founders (", n_founders, ") exceeds the number of candidates in ",
         "value_matrix (", nrow(value_matrix), ").", call. = FALSE)

  if (is.null(block_weights)) {
    block_weights <- rep(1, ncol(value_matrix))
  } else if (length(block_weights) != ncol(value_matrix)) {
    stop("block_weights must have length ncol(value_matrix).", call. = FALSE)
  }

  if (!is.null(merit_priority) &&
      (!is.numeric(merit_priority) || merit_priority < 0 || merit_priority > 100))
    stop("merit_priority must be a single numeric value in [0, 100], or NULL.",
         call. = FALSE)

  cal <- .calibrate_merit_scale(value_matrix, merit_score, block_weights,
                                strategy, n_founders)

  suggested <- NULL
  if (!is.null(merit_priority)) {
    if (!cal$ok) {
      warning("[suggest_merit_weight] merit_score has too little spread ",
              "among these candidates to calibrate a scale factor safely; ",
              "suggested_merit_weight is NULL regardless of merit_priority. ",
              "Supply merit_weight directly if you still want a merit term.",
              call. = FALSE)
    } else {
      suggested <- (merit_priority / 100) * cal$scale_factor
    }
  }

  list(
    merit_span       = cal$merit_span,
    coverage_span     = cal$coverage_span,
    merit_ceiling      = cal$merit_ceiling,
    merit_floor         = cal$merit_floor,
    coverage_ceiling      = cal$coverage_ceiling,
    coverage_floor          = cal$coverage_floor,
    scale_factor              = cal$scale_factor,
    merit_priority             = merit_priority,
    suggested_merit_weight      = suggested,
    ok                            = cal$ok,
    trim_n                         = cal$trim_n
  )
}


# -- Internal: run GA::ga() once for a given value_matrix/strategy/weights.
# Factored out of select_parents_ga() so it can be called n_reps times with
# different seeds for the replication/stability check documented there.
# Returns the selected founders, fitness, per-block contributor table, the
# raw ga_fit object, a convergence flag (ga_fit@iter < maxiter means
# GA::ga()'s own run-generations-without-improvement criterion triggered --
# a real plateau; running the full maxiter without stopping early means the
# search may not have converged), and the realised mean pairwise
# relationship of the chosen set (NA if G was not supplied).
#
# relatedness_ceiling/relatedness_penalty_coef implement target_degree (the
# easy alternative to coancestry_weight): a squared-violation penalty,
# applied only when the chosen set's mean relationship exceeds
# relatedness_ceiling, using the SAME idiom already used for the cardinality
# penalty above (k - n_founders)^2 * penalty_weight -- see
# .calibrate_relatedness_ceiling() for how both are derived. Mutually
# exclusive with coancestry_weight at the select_parents_ga() level, but
# both branches are harmless if somehow both were active (each computes
# `rel` independently).
.run_ga_once <- function(value_matrix, n_founders, strategy, block_weights,
                         popSize, maxiter, run, pmutation, pcrossover,
                         penalty_weight, seed, verbose,
                         G = NULL, coancestry_weight = 0,
                         merit_score = NULL, merit_weight = 0,
                         relatedness_ceiling = NULL,
                         relatedness_penalty_coef = 0) {
  n_cand <- nrow(value_matrix)
  ids    <- rownames(value_matrix)

  fitness_fn <- function(bits) {
    chosen_idx <- which(bits > 0.5)
    k <- length(chosen_idx)
    penalty <- (k - n_founders)^2 * penalty_weight
    if (k < 1L) return(-1e9)
    best_vals <- .block_best_values(value_matrix, chosen_idx, strategy)
    fit <- sum(best_vals * block_weights) - penalty
    if (!is.null(G) && coancestry_weight > 0 && k >= 2L) {
      rel <- .mean_pairwise_relationship(G, ids[chosen_idx])
      fit <- fit - coancestry_weight * rel
    }
    if (!is.null(G) && !is.null(relatedness_ceiling) && k >= 2L) {
      rel  <- .mean_pairwise_relationship(G, ids[chosen_idx])
      viol <- max(0, rel - relatedness_ceiling)
      fit  <- fit - relatedness_penalty_coef * viol^2
    }
    if (!is.null(merit_score) && merit_weight > 0) {
      mm <- mean(merit_score[ids[chosen_idx]])
      fit <- fit + merit_weight * mm
    }
    fit
  }

  if (!is.null(seed)) set.seed(seed)

  # Seed the initial population with valid-cardinality individuals (exactly
  # n_founders bits set) so the GA starts near the feasible region instead of
  # wasting early generations on the penalty term alone.
  suggested <- t(vapply(seq_len(min(popSize, 50L)), function(i) {
    bits <- integer(n_cand)
    bits[sample.int(n_cand, n_founders)] <- 1L
    bits
  }, integer(n_cand)))

  ga_fit <- GA::ga(
    type       = "binary",
    fitness    = fitness_fn,
    nBits      = n_cand,
    popSize    = popSize,
    maxiter    = maxiter,
    run        = run,
    pmutation  = pmutation,
    pcrossover = pcrossover,
    suggestions = suggested,
    monitor    = isTRUE(verbose),
    seed       = seed
  )

  best_bits  <- as.numeric(ga_fit@solution[1L, ])
  chosen_idx <- which(best_bits > 0.5)

  best_vals    <- .block_best_values(value_matrix, chosen_idx, strategy)
  contributors <- .block_best_contributors(value_matrix, chosen_idx, strategy)

  per_block <- data.frame(
    block_id   = colnames(value_matrix),
    best_value = as.numeric(best_vals),
    stringsAsFactors = FALSE
  )
  per_block <- cbind(per_block, contributors)

  mean_rel <- if (!is.null(G)) .mean_pairwise_relationship(G, ids[chosen_idx]) else NA_real_
  mean_mer <- if (!is.null(merit_score)) mean(merit_score[ids[chosen_idx]]) else NA_real_

  list(
    selected         = ids[chosen_idx],
    fitness          = sum(best_vals * block_weights),
    per_block        = per_block,
    ga_fit           = ga_fit,
    converged        = ga_fit@iter < maxiter,
    mean_relationship = mean_rel,
    mean_merit         = mean_mer
  )
}


#' Genetic-Algorithm Founder-Parent Selection
#'
#' @description
#' Searches for the set of \code{n_founders} individuals that jointly
#' maximises coverage of favourable per-block values (local GEBV, haplotype
#' allele dosage, or any comparable score) across a set of target blocks,
#' via \code{GA::ga(type = "binary")}. This is HapBlockR's answer to
#' HapSelect's \code{local_gebv_parent_selection()} /
#' \code{haplotype_parent_selection()} -- the single biggest capability gap
#' identified against HapSelect (see \code{vignette} / gap-analysis notes):
#' turning "here is what's good" (block importance, stacking scores) into
#' "here is the actual founder set", instead of leaving a breeder to eyeball
#' a ranked list.
#'
#' @details
#' The fitness function, for a candidate founder subset \eqn{S} of size
#' \code{n_founders}, is:
#' \deqn{\text{fitness}(S) = \sum_{j \in \text{blocks}} w_j \cdot \text{best}_j(S)}
#' where \eqn{\text{best}_j(S)} is the best value block \eqn{j} can achieve
#' from the founders in \eqn{S}, and \eqn{w_j} is that block's weight
#' (\code{block_weights}, default 1 for every block). How
#' \eqn{\text{best}_j(S)} is computed depends on \code{strategy} -- see the
#' \emph{Crossing-scheme strategies} section. This directly generalises the
#' pair-average formula shown in HapSelect's own documentation diagram
#' (\eqn{\sum_j \max((\text{localGEBV}_{j,1} + \text{localGEBV}_{j,2})/2)}),
#' but is HapBlockR's own implementation -- HapSelect's exact GA fitness
#' function and constraint semantics are not available beyond its published
#' documentation, so this is not a byte-for-byte port.
#'
#' A binary GA chromosome (one bit per candidate individual) is used rather
#' than a fixed-cardinality representation, with a quadratic penalty for
#' deviating from \code{n_founders} selected bits -- standard practice for
#' subset-selection problems with the \pkg{GA} package, since it has no
#' native hard-cardinality constraint.
#'
#' @section Crossing-scheme strategies:
#' \describe{
#'   \item{\code{"no_selfing"} (default)}{Each block's value is the mean of
#'     the \strong{two largest} values among the chosen founders -- two
#'     distinct parents must jointly contribute it. Mirrors HapSelect's
#'     localGEBV-mode \code{"no_selfing"}.}
#'   \item{\code{"selfing"}}{Each block's value is the \strong{single
#'     largest} value among the chosen founders -- a parent may be selfed to
#'     realise its own value alone, without needing a complementary partner.
#'     Mirrors HapSelect's localGEBV-mode \code{"selfing"}.}
#'   \item{\code{"OHS"}}{Same computation as \code{"no_selfing"}
#'     (two distinct parents required), named for haplotype-mode use where
#'     each parent contributes one haplotype copy.}
#'   \item{\code{"OPV"}}{Same computation as \code{"selfing"} (single best
#'     value) -- "population value" framing: any one founder carrying the
#'     favourable haplotype is enough, poolable across the whole founder set
#'     rather than requiring a specific complementary pairing.}
#'   \item{\code{"Haploid_OHS"}}{Same computation as \code{"selfing"} -- a
#'     single (heterozygous) parent may donate two non-homologous gametes.}
#' }
#' \strong{Why \code{"selfing"}, \code{"OPV"}, and \code{"Haploid_OHS"}
#' compute identically:} under the max-pair-average formulation, allowing a
#' founder to pair with itself can never do worse than any distinct pair,
#' because the largest value in a set is always \eqn{\geq} the mean of
#' itself and anything smaller. This is a direct mathematical consequence of
#' the fitness formula, not a simplification of convenience -- documented
#' here so it isn't mistaken for the three strategies being unimplemented.
#'
#' @param value_matrix Numeric matrix (individuals x blocks): local GEBV
#'   (\code{run_haplotype_prediction()$local_gebv}), haplotype allele dosage
#'   scaled to [0,1] (dosage / ploidy), or any comparable per-block
#'   favourable-value score. Row names = candidate individual IDs.
#'   \strong{Pre-filter to a manageable set of top-ranked blocks first} (e.g.
#'   via \code{\link{select_top_blocks}}) -- exactly as HapSelect's own
#'   documented workflow does before calling its parent-selection GA -- since
#'   the fitness function is evaluated \code{popSize * maxiter} times.
#' @param n_founders Integer. Target number of founders to select.
#' @param strategy One of \code{"no_selfing"} (default), \code{"selfing"},
#'   \code{"OHS"}, \code{"OPV"}, \code{"Haploid_OHS"}. See
#'   \emph{Crossing-scheme strategies}.
#' @param block_weights Numeric vector, length \code{ncol(value_matrix)}, or
#'   \code{NULL} (default: equal weight 1 for every block). E.g. pass the
#'   \code{var_scaled} column from \code{\link{select_top_blocks}} to weight
#'   higher-variance blocks more heavily.
#' @param top_candidates Integer or \code{NULL} (default). If supplied,
#'   restricts the GA's candidate pool to the top \code{top_candidates}
#'   individuals by their own row-max value in \code{value_matrix} before
#'   searching -- a heuristic prefilter to keep the binary chromosome length
#'   (and therefore search time) manageable for large populations. Set
#'   \code{NULL} to search the full population in \code{value_matrix}.
#' @param popSize,maxiter,run,pmutation,pcrossover Passed to
#'   \code{GA::ga()}. Defaults \code{100}, \code{200}, \code{50}, \code{0.1},
#'   \code{0.8} respectively.
#' @param penalty_weight Numeric. Quadratic penalty coefficient for deviating
#'   from exactly \code{n_founders} selected individuals. Default \code{10}
#'   times the largest single block weight, scaled so the penalty dominates
#'   the fitness for any deviation while not distorting the optimum for
#'   correctly-sized subsets. Increase if the GA returns solutions with the
#'   wrong count.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility. When
#'   \code{n_reps > 1}, this seeds replicate 1 and replicates
#'   \code{2..n_reps} use \code{seed + 1}, \code{seed + 2}, ... (deterministic
#'   and reproducible as a set, but each replicate explores a different
#'   starting point).
#' @param verbose Logical. Print \code{GA::ga()}'s iteration monitor for the
#'   first replicate only (silent for replicates \code{2..n_reps} to avoid
#'   flooding the console). Default \code{FALSE}.
#' @param merit_score Named numeric vector (e.g.
#'   \code{run_haplotype_prediction()$gebv}), or \code{NULL} (default).
#'   Whole-genome merit, used for up to two independent purposes depending on
#'   which of \code{min_sel_value}/\code{merit_weight} are set: (1) the
#'   \code{min_sel_value} eligibility floor below, a hard pre-search
#'   exclusion; and (2) when \code{merit_weight > 0}, a soft, additive term
#'   \emph{inside} the GA fitness function itself -- see \emph{Merit-weighted
#'   fitness (GA+TS hybrid, optional)}. The two are independent and
#'   commonly used together (a floor to exclude clearly ineligible
#'   candidates, plus a weight so merit keeps pulling on the search among
#'   those who remain), but either can be used alone. Required whenever
#'   \code{min_sel_value} is not \code{NULL}, or \code{merit_weight > 0}.
#' @param min_sel_value,min_sel_mode Merit floor applied to
#'   \code{value_matrix}'s candidate pool \emph{before} the GA searches, using
#'   \code{merit_score} to evaluate it. Same semantics as
#'   \code{\link{truncation_selection}}'s arguments of the same name; default
#'   \code{min_sel_value = NULL} applies no floor at all (the historical
#'   behaviour of this function). Existing from an earlier finding: without a
#'   floor, a genuinely poor overall performer who happens to uniquely carry
#'   one target block's favourable value will still be selected, purely to
#'   cover that block -- setting a floor here is the fix.
#' @param n_reps Integer >= 1. Default \code{5L}. Number of independent GA
#'   replicates to run (see \emph{GA rigour: replication and convergence}
#'   below); the replicate with the best fitness is returned as
#'   \code{selected}/\code{fitness}/\code{per_block}/\code{ga_fit}. Set to
#'   \code{1L} for the fastest, single-run legacy behaviour.
#' @param G Relationship/kinship matrix (n x n, dimnames = individual IDs
#'   covering every candidate in \code{value_matrix}), e.g.
#'   \code{run_haplotype_prediction()$G} or
#'   \code{\link{compute_haplotype_grm}} output, or \code{NULL} (default).
#'   Required whenever \code{coancestry_weight > 0} or \code{target_degree}
#'   is supplied; see \emph{Coancestry penalty (optional)} below. Unused (may
#'   be left \code{NULL}) otherwise.
#' @param coancestry_weight Numeric >= 0. Default \code{0} (no coancestry
#'   term at all -- the historical behaviour of this function, block coverage
#'   only). When positive, subtracts \code{coancestry_weight} times the mean
#'   off-diagonal pairwise relationship of the chosen set (from \code{G})
#'   from the fitness function, so the GA trades off block coverage against
#'   relatedness directly rather than leaving relatedness as an after-the-
#'   fact diagnostic. There is no universal default scale for this argument
#'   -- it must be tuned against your own \code{block_weights}/fitness scale
#'   (see \emph{Coancestry penalty (optional)}). Mutually exclusive with
#'   \code{target_degree} -- do not supply both.
#' @param target_degree Numeric in \code{[0, 90]}, or \code{NULL} (default --
#'   no effect, historical behaviour unchanged). The easy alternative to
#'   \code{coancestry_weight}: instead of a raw, unscaled penalty multiplier,
#'   state your relatedness-vs-gain preference on the SAME \code{[0, 90]}
#'   scale \code{\link{select_parents_ocs}} already uses (\code{0} = max
#'   gain, prioritising coverage/merit and accepting more relatedness;
#'   \code{90} = max diversity, minimising relatedness), and this function
#'   converts it into an internal relatedness-ceiling penalty from your own
#'   data before running the search -- see \emph{Relatedness ceiling
#'   (target_degree, optional)} below for exactly how, and its documented
#'   limits. Requires \code{G}. Mutually exclusive with an explicitly-supplied
#'   \code{coancestry_weight} (ambiguous otherwise) -- supplying both is an
#'   error.
#' @param merit_weight Numeric >= 0. Default \code{0} (no merit term at all
#'   -- the historical behaviour of this function, block coverage only, with
#'   \code{merit_score} usable purely as a pre-search eligibility floor via
#'   \code{min_sel_value}). When positive, \strong{adds}
#'   \code{merit_weight} times the chosen set's mean \code{merit_score} to
#'   the fitness function, so the GA rewards whole-genome merit directly
#'   inside the search itself, rather than only using it to exclude clearly
#'   ineligible candidates beforehand. See \emph{Merit-weighted fitness
#'   (GA+TS hybrid, optional)} below. Requires \code{merit_score}. Mutually
#'   exclusive with \code{merit_priority} -- do not supply both.
#' @param merit_priority Numeric in \code{[0, 100]}, or \code{NULL} (default
#'   -- no effect, historical behaviour unchanged). The easy alternative to
#'   \code{merit_weight}: instead of a raw, unscaled multiplier, state how
#'   much you care about merit vs. coverage as a plain percentage, and this
#'   function calibrates the right \code{merit_weight} from your own data
#'   before running the search -- see \code{\link{suggest_merit_weight}} for
#'   exactly what "comparable scale" means here and its documented limits.
#'   \code{0} is identical to \code{merit_weight = 0}; \code{100} scales
#'   merit's realistic influence on the fitness function to match coverage's;
#'   values between scale linearly. Requires \code{merit_score}. Mutually
#'   exclusive with an explicitly-supplied \code{merit_weight} (ambiguous
#'   otherwise) -- supplying both is an error.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{selected}}{Character vector of the selected individual IDs
#'     from the best-fitness replicate (length \code{n_founders}, or as close
#'     as the GA achieved -- check \code{length(selected) == n_founders}).}
#'   \item{\code{fitness}}{Numeric. Best fitness value found (raw block-value
#'     sum, after subtracting the cardinality penalty), best-fitness
#'     replicate.}
#'   \item{\code{per_block}}{Data frame for the best-fitness replicate:
#'     \code{block_id}, \code{best_value}, \code{contributor_1},
#'     \code{contributor_2} (equal to \code{contributor_1} under a
#'     self-allowed strategy).}
#'   \item{\code{strategy}}{Character, echoes the \code{strategy} argument.}
#'   \item{\code{ga_fit}}{The raw \code{GA::ga()} S4 result object for the
#'     best-fitness replicate, for convergence diagnostics (e.g.
#'     \code{plot(ga_fit)}).}
#'   \item{\code{converged}}{Logical. \code{TRUE} if the best-fitness
#'     replicate's GA stopped because fitness plateaued for \code{run}
#'     consecutive generations (\code{ga_fit@iter < maxiter}), rather than
#'     being cut off at \code{maxiter} without plateauing. \code{FALSE} means
#'     consider raising \code{maxiter}.}
#'   \item{\code{cutoff}}{Numeric. The \code{min_sel_value} cutoff actually
#'     applied to \code{merit_score} (\code{-Inf} when
#'     \code{min_sel_value = NULL}).}
#'   \item{\code{stability}}{List describing agreement across all
#'     \code{n_reps} replicates -- see \emph{GA rigour} below. Present even
#'     when \code{n_reps = 1} (trivially: every selected individual has
#'     \code{selection_freq = 1}). Also includes \code{$mean_relationship}
#'     (realised mean pairwise relationship of each replicate's chosen set,
#'     \code{NA} unless \code{G} was supplied) and \code{$mean_merit}
#'     (realised mean \code{merit_score} of each replicate's chosen set,
#'     \code{NA} unless \code{merit_score} was supplied).}
#'   \item{\code{mean_relationship}}{Numeric. The best-fitness replicate's
#'     realised mean off-diagonal pairwise relationship (from \code{G}) among
#'     the final \code{selected} set. \code{NA} unless \code{G} was supplied
#'     -- present regardless of whether \code{coancestry_weight > 0} or
#'     \code{target_degree} was set, so you can inspect relatedness even when
#'     the GA wasn't asked to optimise for it (compare against an
#'     unconstrained run to see the effect).}
#'   \item{\code{mean_merit}}{Numeric. The best-fitness replicate's realised
#'     mean \code{merit_score} among the final \code{selected} set. \code{NA}
#'     unless \code{merit_score} was supplied -- present regardless of
#'     whether \code{merit_weight > 0}, so you can inspect the chosen set's
#'     whole-genome merit even when the GA wasn't asked to optimise for it
#'     directly (compare against a \code{merit_weight = 0} run to see the
#'     term's effect).}
#'   \item{\code{merit_weight}}{Numeric. The merit weight actually used --
#'     either your explicit argument, or the value calibrated from
#'     \code{merit_priority} if that was supplied instead.}
#'   \item{\code{merit_priority}}{Echoes the \code{merit_priority} argument
#'     (\code{NULL} unless supplied).}
#'   \item{\code{coancestry_weight}}{Numeric. Echoes the \code{coancestry_weight}
#'     argument as actually used (\code{0} when \code{target_degree} was
#'     supplied instead -- the relatedness constraint is enforced via
#'     \code{relatedness_ceiling} in that case, not this argument).}
#'   \item{\code{target_degree}}{Echoes the \code{target_degree} argument
#'     (\code{NULL} unless supplied).}
#'   \item{\code{relatedness_ceiling}}{Numeric. The interpolated mean-
#'     relationship ceiling actually enforced when \code{target_degree} was
#'     supplied (\code{NULL} otherwise) -- see \emph{Relatedness ceiling
#'     (target_degree, optional)}. Compare against \code{mean_relationship}
#'     to check whether the search respected it.}
#' }
#'
#' @section Coancestry penalty (optional):
#' By default (\code{coancestry_weight = 0}) this function optimises block
#' coverage only, exactly as described in \emph{Details} -- relatedness among
#' the chosen founders plays no role in the search, which is why
#' \code{\link{plot_parent_selection_pca}} and a manual family-balance check
#' are recommended \emph{after} calling this function (see \emph{Choosing
#' between this function and truncation_selection()} below). Setting
#' \code{coancestry_weight > 0} (with \code{G} supplied) moves that check
#' \emph{into} the search itself: the fitness function becomes
#' \deqn{\text{fitness}(S) = \sum_j w_j \cdot \text{best}_j(S) -
#' \text{coancestry\_weight} \cdot \overline{G}(S)}
#' where \eqn{\overline{G}(S)} is the mean off-diagonal pairwise relationship
#' among the chosen set \eqn{S} (each individual's own diagonal
#' self-relationship/inbreeding term is excluded -- the penalty reflects
#' relatedness \emph{between} chosen parents, not their own inbreeding). This
#' is a direct, minimal extension of the existing block-coverage GA -- not a
#' full optimal-contribution-selection (OCS) formulation (see \emph{Choosing
#' between this function and truncation_selection()} for what a full OCS/
#' mate-allocation tool adds beyond this).
#'
#' \strong{Tuning \code{coancestry_weight}:} there is no universal correct
#' value -- \eqn{\overline{G}(S)} and the block-coverage sum are on different,
#' problem-specific scales (a VanRaden-style GRM is typically much smaller in
#' magnitude than a sum of local-GEBV values across many blocks). Start by
#' running with \code{coancestry_weight = 0} and inspecting \code{$fitness}'s
#' typical magnitude, then choose a weight that makes the coancestry term
#' large enough to visibly compete with it; increase further if the
#' resulting \code{selected} set is still too related for your program, and
#' compare \code{$mean_relationship} across a few candidate weights to see
#' the trade-off curve directly. An easier alternative to tuning this by
#' hand is \code{target_degree}, described next.
#'
#' @section Relatedness ceiling (target_degree, optional):
#' \code{coancestry_weight} has the same scale problem \code{merit_weight}
#' does (see \emph{Merit-weighted fitness} below): there is no universal
#' correct multiplier, because \eqn{\overline{G}(S)} and the block-coverage
#' sum live on different, dataset-specific scales. Rather than estimating a
#' matching scale for a weighted trade-off (the approach
#' \code{merit_priority} takes for the merit term), \code{target_degree}
#' sidesteps the scale problem altogether: it converts your \code{[0, 90]}
#' preference into a relatedness CEILING -- a bound in \code{G}'s own real
#' units, not an abstract multiplier -- interpolated between two genuinely
#' constructed reference groups for your actual candidate pool, exactly
#' mirroring how \code{\link{select_parents_ocs}}'s \code{engine = "optisel"}
#' interpolates its own kinship bound between two solved frontier extremes:
#' \itemize{
#'   \item \strong{Gain end (\code{target_degree = 0})}: the feasible,
#'     greedy-built high-coverage group (see \emph{Details}) -- extended to
#'     also account for \code{merit_weight} when merit is active in this same
#'     call, so this reference reflects what the search would actually
#'     converge to with no relatedness constraint, not a coverage-only proxy.
#'     Its own mean pairwise relationship becomes the ceiling at
#'     \code{target_degree = 0} (i.e. no additional constraint beyond what
#'     that group naturally has).
#'   \item \strong{Diversity end (\code{target_degree = 90})}:
#'     \code{\link{select_core_collection}(strategy = "maximin")} run on
#'     \code{G} for your candidate pool -- a proven 2-approximation for
#'     maximizing the minimum pairwise distance (Gonzalez 1985), built
#'     ignoring merit entirely. Its own mean pairwise relationship becomes
#'     the ceiling at \code{target_degree = 90}.
#' }
#' The ceiling at any \code{target_degree} in between is a linear
#' interpolation between these two real, feasible reference points. Inside
#' the search, exceeding the ceiling is penalised by a squared-violation term
#' scaled to dominate any possible coverage/merit gain from crossing it (the
#' same soft-constraint idiom already used for the cardinality penalty) --
#' so no separate penalty-strength argument is needed.
#'
#' \strong{Honest limits:} the gain-end reference is a plain greedy hill-
#' climb, not a formally guaranteed-near-optimal one, once
#' \code{merit_weight > 0} (adding a mean-based merit term breaks the strict
#' submodularity argument that gives the coverage-only greedy its proof) --
#' still a real, feasible group, just without that guarantee. Separately,
#' \code{select_core_collection}'s 2-approximation guarantee covers the
#' MINIMUM pairwise distance in its chosen group, not directly the MEAN
#' pairwise relationship used here as the diversity-end reference -- a
#' principled, reused building block, not a proof about the specific number
#' reported. Neither limitation is unique to this feature; they mirror the
#' same kind of honestly-stated approximation caveats already documented for
#' \code{merit_priority} (see \code{\link{suggest_merit_weight}}) and for
#' \code{select_parents_ocs}'s own \code{target_degree}.
#'
#' @section Merit-weighted fitness (GA+TS hybrid, optional):
#' By default (\code{merit_weight = 0}) this function's search rewards block
#' coverage only, exactly as described in \emph{Details} -- \code{merit_score},
#' if supplied at all, only ever acts as a pre-search eligibility floor via
#' \code{min_sel_value} (a hard yes/no cutoff), never as something the search
#' itself is rewarded for pursuing further. That matters because a floor
#' alone does not distinguish between a candidate that barely clears it and
#' one that clears it by a wide margin -- both are equally "eligible," and
#' among the eligible pool the search optimises purely for which target
#' blocks each candidate happens to carry good value at. On a real panel,
#' where per-block/local-GEBV estimates carry real estimation noise (they
#' are themselves statistical estimates, not ground truth), this can let the
#' search lean on a candidate whose apparent block coverage is partly or
#' wholly a noisy artefact, provided that candidate clears the floor at all.
#'
#' Setting \code{merit_weight > 0} closes this gap directly: the fitness
#' function becomes
#' \deqn{\text{fitness}(S) = \sum_j w_j \cdot \text{best}_j(S) +
#' \text{merit\_weight} \cdot \overline{\text{merit\_score}}(S)}
#' (combined with the coancestry penalty, when both are active, as a third
#' additive term). The search now has a direct, continuous incentive to
#' prefer higher-merit candidates \emph{throughout} the eligible pool, not
#' merely to avoid excluding low-merit ones -- pulling the result back
#' toward \code{\link{truncation_selection}}'s whole-genome-merit ranking
#' while still searching for joint block coverage, rather than treating
#' merit and coverage as two entirely separate stages (a floor, then an
#' unconstrained coverage search). This is the recommended way to combine
#' this function with a truncation-selection-style merit signal when you are
#' specifically concerned that block-coverage-only selection, run on noisy
#' per-block estimates, could otherwise select on noise: a nonzero
#' \code{merit_weight} keeps whole-genome merit -- typically a more stable,
#' better-estimated signal than any single block's local GEBV -- pulling on
#' every candidate's desirability throughout the search, not only at its
#' entry gate.
#'
#' \strong{Tuning \code{merit_weight}:} the same scale caveat as
#' \code{coancestry_weight} applies, and for the same reason --
#' \eqn{\overline{\text{merit\_score}}(S)} and the block-coverage sum are on
#' different, problem-specific scales, so there is no universal correct
#' value. Start by running with \code{merit_weight = 0} and inspecting
#' \code{$fitness}'s typical magnitude and \code{$mean_merit}'s typical
#' magnitude, then choose a weight that makes the merit term large enough to
#' visibly compete with the block-coverage sum; compare \code{$selected} and
#' \code{$mean_merit} across a few candidate weights (and against a plain
#' \code{\link{truncation_selection}} run on the same \code{merit_score}) to
#' see the effect directly, exactly as recommended for \code{coancestry_weight}
#' above. \code{min_sel_value} and \code{merit_weight} address different
#' failure modes and are not substitutes for each other: keep using
#' \code{min_sel_value} for a hard floor against candidates you never want
#' to consider at all, and add \code{merit_weight} on top of it when you also
#' want merit to keep influencing the choice among everyone who clears that
#' floor.
#'
#' \strong{An easier alternative to tuning \code{merit_weight} by hand:} set
#' \code{merit_priority} instead (a plain 0-100 "how much do I care about
#' merit vs. coverage" dial) and let this function calibrate the matching
#' \code{merit_weight} from your own data automatically -- see
#' \code{\link{suggest_merit_weight}} for exactly what it computes, why (a
#' greedy, provably near-optimal reachable high-coverage group as the
#' reference point, not an inflated unreachable one), and what it honestly
#' does not solve (there is no single universally "correct" notion of
#' comparable scale -- \code{merit_priority}'s suggestion is a well-reasoned
#' starting point, not the only defensible number). \code{merit_weight} and
#' \code{merit_priority} are mutually exclusive -- pick whichever style suits
#' you; you can also call \code{\link{suggest_merit_weight}} first to inspect
#' the calibration numbers, then pass a literal \code{merit_weight} yourself
#' if you want to fine-tune from there.
#'
#' @section GA rigour: replication and convergence:
#' A single GA run, taken at face value, tells you nothing about whether its
#' answer is a robust optimum or one of several near-equally-good solutions a
#' stochastic search happened to land on. This function addresses that
#' directly rather than leaving it to the caller:
#' \describe{
#'   \item{Convergence}{Every replicate's \code{ga_fit@iter} (generations
#'     actually run) is checked against \code{maxiter}. Stopping early means
#'     \pkg{GA}'s own \code{run}-generations-without-improvement criterion
#'     triggered -- a real plateau. Running the full \code{maxiter} without
#'     stopping early means the search may not have converged; \code{$converged}
#'     surfaces this instead of silently returning a possibly-unconverged
#'     result.}
#'   \item{Replication}{With \code{n_reps > 1} (default \code{5}), the search
#'     is repeated from \code{n_reps} different starting populations/seeds.
#'     \code{$stability$selection_freq} reports, for every candidate selected
#'     in at least one replicate, the fraction of replicates that selected
#'     them -- individuals at \code{1.0} are robustly supported regardless of
#'     the GA's random starting point; individuals selected in only one
#'     replicate out of several are borderline and worth a second look before
#'     committing to them. \code{$stability$fitness_range} shows how much the
#'     best achievable fitness varied across replicates -- a narrow range
#'     alongside high selection frequencies is the signature of a stable,
#'     trustworthy search.}
#' }
#' This is deliberately more rigorous by default than a single fixed-seed GA
#' run: HapSelect's own documented parent-selection GA (see \emph{Description})
#' does not report replication stability or a convergence flag at all. Set
#' \code{n_reps = 1} only once you have separately confirmed stability, or for
#' quick iteration during exploratory analysis.
#'
#' @section Choosing between this function and truncation_selection():
#' \code{\link{truncation_selection}} ranks by a single whole-genome score and
#' takes the top \code{n_founders}; this function instead searches for the
#' \emph{set} of \code{n_founders} that jointly covers the most favourable
#' value across a chosen set of target blocks (typically the output of
#' \code{\link{select_top_blocks}}). They answer genuinely different
#' questions and neither is a strict improvement on the other:
#' \describe{
#'   \item{\code{truncation_selection()} ignores}{which haplotype blocks each
#'     individual carries -- its top-N list can concentrate on the same
#'     favourable blocks while leaving others uncovered.}
#'   \item{\code{select_parents_ga()} ignores}{whole-genome merit by default,
#'     unless \code{merit_score} is supplied via \code{min_sel_value} (a hard
#'     pre-search floor), \code{merit_weight} (a soft, continuous term
#'     rewarded throughout the search -- see \emph{Merit-weighted fitness
#'     (GA+TS hybrid, optional)}), or both together. Without either, a
#'     genuinely poor-performing individual who happens to uniquely carry one
#'     target block's favourable value will still be selected, purely to
#'     cover that block -- and, even with only a floor and no weight, a
#'     candidate who barely clears the floor is treated identically to one
#'     who clears it by a wide margin. Setting \code{merit_weight > 0} is the
#'     tool's own built-in "GA+TS hybrid" mode: it keeps searching for joint
#'     block coverage while also rewarding whole-genome merit directly inside
#'     the search, which matters most when block-level/local-GEBV estimates
#'     carry real estimation noise and a coverage-only search risks selecting
#'     on that noise rather than on genuine signal.}
#'   \item{\strong{Neither function, by default}}{manages coancestry/
#'     inbreeding risk in the chosen set, assigns differential contributions,
#'     or decides who mates whom. This function's optional
#'     \code{coancestry_weight}/\code{target_degree}/\code{G} arguments (see
#'     \emph{Coancestry penalty (optional)} and \emph{Relatedness ceiling
#'     (target_degree, optional)}) add a relatedness term to the search
#'     itself, but still stop at a flat founder set -- no differential
#'     contributions, no mating list. Dedicated optimal-contribution-selection (OCS) and
#'     mate-allocation tools (e.g. AlphaMate, Gorjanc lab) solve all three
#'     problems jointly, typically maximising gain subject to an explicit
#'     coancestry constraint plus continuous contribution optimisation. A
#'     practical combined workflow: use this function (or
#'     \code{truncation_selection()}) to shortlist candidates worth
#'     considering, then hand that shortlist to a dedicated OCS tool for the
#'     coancestry-managed contribution and mating decision, rather than
#'     treating either function's output as a final mating plan.}
#' }
#' In practice, run both and compare with
#' \code{intersect(ga_out$selected, ts_out$selected)}: a high overlap means
#' truncation selection was already close to the GA's answer for this panel;
#' a low overlap means block coverage and whole-genome ranking disagree and
#' both founder sets are worth inspecting (including a family/diversity check
#' via \code{\link{plot_parent_selection_pca}}) before committing to either.
#' See the \emph{From Local GEBV to a Crossing Decision} vignette for a full
#' worked comparison.
#'
#' @references
#' Scrucca L (2013). GA: A Package for Genetic Algorithms in R.
#' \emph{Journal of Statistical Software} 53(4):1-37. \doi{10.18637/jss.v053.i04}
#'
#' @seealso \code{\link{truncation_selection}}, \code{\link{select_top_blocks}},
#'   \code{\link{plot_parent_selection_pca}}, \code{\link{ga_vs_ts_simulation}},
#'   \code{\link{suggest_merit_weight}}, \code{\link{select_parents_ocs}} (the
#'   \code{target_degree} convention this function reuses),
#'   \code{\link{select_core_collection}} (\code{target_degree}'s
#'   diversity-end reference)
#'
#' @examples
#' \dontrun{
#' res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
#' top  <- select_top_blocks(res$block_importance, n = 15)
#' vmat <- res$local_gebv[, top$block_id, drop = FALSE]
#' ga_out <- select_parents_ga(vmat, n_founders = 20, strategy = "no_selfing",
#'                             seed = 1)
#' ga_out$selected
#' }
#'
#' @export
select_parents_ga <- function(value_matrix,
                              n_founders,
                              strategy = c("no_selfing", "selfing", "OHS",
                                          "OPV", "Haploid_OHS"),
                              block_weights  = NULL,
                              top_candidates = NULL,
                              popSize    = 100L,
                              maxiter    = 200L,
                              run        = 50L,
                              pmutation  = 0.1,
                              pcrossover = 0.8,
                              penalty_weight = NULL,
                              seed       = NULL,
                              verbose    = FALSE,
                              merit_score   = NULL,
                              min_sel_value = NULL,
                              min_sel_mode  = c("value", "percentile",
                                               "sd_below_mean"),
                              n_reps     = 5L,
                              G = NULL,
                              coancestry_weight = 0,
                              merit_weight = 0,
                              merit_priority = NULL,
                              target_degree = NULL) {
  strategy <- match.arg(strategy)

  if (coancestry_weight > 0 && is.null(G))
    stop("coancestry_weight > 0 requires G (a relationship/kinship matrix ",
         "with dimnames covering value_matrix's candidates).", call. = FALSE)
  if (!is.null(G) && (is.null(rownames(G)) || is.null(colnames(G))))
    stop("G must have row and column names (individual IDs).", call. = FALSE)

  # -- target_degree: the easy 0-90 dial, mutually exclusive with an
  # explicitly-supplied coancestry_weight (ambiguous otherwise). Same
  # direction convention as select_parents_ocs()'s target_degree: 0 = max
  # gain (prioritises coverage/merit, accepting more relatedness), 90 = max
  # diversity (minimises relatedness). Converted into an internal
  # relatedness-ceiling penalty further down, AFTER merit_weight is fully
  # resolved and merit_score is validated to cover the final candidate pool
  # -- see .calibrate_relatedness_ceiling().
  if (!is.null(target_degree)) {
    if (!missing(coancestry_weight))
      stop("Supply exactly one of coancestry_weight or target_degree, not ",
           "both -- target_degree is converted into an internal ",
           "relatedness-ceiling penalty, so passing both is ambiguous.",
           call. = FALSE)
    if (!is.numeric(target_degree) || target_degree < 0 || target_degree > 90)
      stop("target_degree must be a single numeric value in [0, 90] (same ",
           "convention as select_parents_ocs()'s target_degree: 0 = max ",
           "gain, prioritising coverage/merit; 90 = max diversity, ",
           "minimising relatedness).", call. = FALSE)
    if (is.null(G))
      stop("target_degree requires G (a relationship/kinship matrix with ",
           "dimnames covering value_matrix's candidates).", call. = FALSE)
  }

  # -- merit_priority: the easy 0-100 dial, mutually exclusive with an
  # explicitly-supplied merit_weight (ambiguous otherwise -- see
  # ?suggest_merit_weight). Converted into an actual merit_weight value
  # further down, AFTER min_sel_value/top_candidates filtering, so the
  # calibration reflects the real candidate pool the GA will search.
  if (!is.null(merit_priority)) {
    if (!missing(merit_weight))
      stop("Supply exactly one of merit_weight or merit_priority, not both ",
           "-- merit_priority is converted into merit_weight internally, so ",
           "passing both is ambiguous. See ?suggest_merit_weight.",
           call. = FALSE)
    if (!is.numeric(merit_priority) || merit_priority < 0 || merit_priority > 100)
      stop("merit_priority must be a single numeric value in [0, 100].",
           call. = FALSE)
    if (is.null(merit_score))
      stop("merit_priority requires merit_score (a named whole-genome merit ",
           "vector) -- it is calibrated against merit_score's own spread. ",
           "See ?suggest_merit_weight.", call. = FALSE)
  }

  if (merit_weight > 0 && is.null(merit_score))
    stop("merit_weight > 0 requires merit_score (a named whole-genome merit ",
         "vector) -- it is the value the weighted merit term in the fitness ",
         "function is computed from. See ?select_parents_ga's 'Merit-",
         "weighted fitness (GA+TS hybrid, optional)' section.", call. = FALSE)

  if (!requireNamespace("GA", quietly = TRUE))
    stop("GA is required for select_parents_ga(). ",
         "Install with: install.packages('GA')", call. = FALSE)
  if (!is.matrix(value_matrix)) value_matrix <- as.matrix(value_matrix)
  if (is.null(rownames(value_matrix)))
    stop("value_matrix must have row names (candidate individual IDs).",
         call. = FALSE)
  if (is.null(colnames(value_matrix)))
    colnames(value_matrix) <- paste0("block_", seq_len(ncol(value_matrix)))

  # -- Merit floor (optional) -- filters the candidate POOL before anything
  # else, using merit_score (whole-genome value), not value_matrix itself.
  # See .apply_merit_floor() and ?select_parents_ga's min_sel_value docs.
  if (!is.null(min_sel_value) && is.null(merit_score))
    stop("min_sel_value requires merit_score (a named whole-genome value ",
         "vector) to evaluate the floor against.", call. = FALSE)

  cutoff <- -Inf
  if (!is.null(min_sel_value)) {
    missing_ids <- setdiff(rownames(value_matrix), names(merit_score))
    if (length(missing_ids))
      warning(length(missing_ids), " value_matrix candidate(s) have no ",
              "merit_score entry and are excluded from consideration: ",
              paste(utils::head(missing_ids, 10), collapse = ", "),
              if (length(missing_ids) > 10) ", ..." else "", call. = FALSE)
    floor_res <- .apply_merit_floor(merit_score, min_sel_value, min_sel_mode,
                                    label = "candidate")
    cutoff <- floor_res$cutoff
    keep   <- intersect(rownames(value_matrix), floor_res$eligible)
    if (!length(keep))
      stop("No value_matrix candidates clear the min_sel_value floor.",
           call. = FALSE)
    value_matrix <- value_matrix[keep, , drop = FALSE]
  }

  n_founders <- as.integer(n_founders)
  if (n_founders < 1L)
    stop("n_founders must be >= 1.", call. = FALSE)
  if (n_founders > nrow(value_matrix))
    stop("n_founders (", n_founders, ") exceeds the number of eligible ",
         "candidate individuals in value_matrix (", nrow(value_matrix),
         ").", call. = FALSE)

  # Mean-impute any NA (e.g. individuals missing from a singleton block)
  for (j in seq_len(ncol(value_matrix))) {
    na_j <- is.na(value_matrix[, j])
    if (any(na_j)) value_matrix[na_j, j] <- mean(value_matrix[, j], na.rm = TRUE)
  }

  if (is.null(block_weights)) {
    block_weights <- rep(1, ncol(value_matrix))
  } else if (length(block_weights) != ncol(value_matrix)) {
    stop("block_weights must have length ncol(value_matrix).", call. = FALSE)
  }

  # -- Optional candidate-pool prefilter (search-space size safeguard) -------
  if (!is.null(top_candidates)) {
    top_candidates <- as.integer(top_candidates)
    if (top_candidates < n_founders)
      stop("top_candidates must be >= n_founders.", call. = FALSE)
    row_best <- apply(value_matrix, 1L, max)
    keep_ind <- order(row_best, decreasing = TRUE)[
      seq_len(min(top_candidates, nrow(value_matrix)))]
    value_matrix <- value_matrix[keep_ind, , drop = FALSE]
  }

  # -- merit_priority -> merit_weight conversion, now that value_matrix
  # reflects the FINAL candidate pool (after min_sel_value/top_candidates
  # filtering) the GA will actually search -- see .calibrate_merit_scale()
  # and ?suggest_merit_weight for the full method and its documented limits.
  if (!is.null(merit_priority)) {
    # Same check .calibrate_merit_scale() itself has no way to detect on its
    # own (a missing name just silently becomes NA, corrupting the sort/mean
    # calculations below with no clear error) -- validated here, BEFORE
    # calibration runs, not only later for the merit_weight > 0 path.
    missing_p <- setdiff(rownames(value_matrix), names(merit_score))
    if (length(missing_p))
      stop(length(missing_p), " value_matrix candidate(s) have no ",
           "merit_score entry, required because merit_priority is set: ",
           paste(utils::head(missing_p, 10), collapse = ", "),
           if (length(missing_p) > 10) ", ..." else "", call. = FALSE)

    cal <- .calibrate_merit_scale(value_matrix, merit_score, block_weights,
                                  strategy, n_founders)
    if (!cal$ok) {
      warning("[select_parents_ga] merit_priority requested, but ",
              "merit_score has too little spread among the final candidate ",
              "pool to calibrate a scale factor safely; merit_weight stays ",
              "0 (no merit term). Supply merit_weight directly if you still ",
              "want one. See ?suggest_merit_weight.", call. = FALSE)
      merit_weight <- 0
    } else {
      merit_weight <- (merit_priority / 100) * cal$scale_factor
      if (isTRUE(verbose))
        message("[select_parents_ga] merit_priority = ", merit_priority,
                "% -> merit_weight = ", signif(merit_weight, 4),
                " (coverage_span = ", signif(cal$coverage_span, 4),
                ", merit_span = ", signif(cal$merit_span, 4), ").")
    }
  }

  if (is.null(penalty_weight))
    penalty_weight <- 10 * sum(block_weights) / max(1, ncol(value_matrix))

  n_reps <- as.integer(n_reps)
  if (n_reps < 1L)
    stop("n_reps must be >= 1.", call. = FALSE)

  if (!is.null(G)) {
    missing_g <- setdiff(rownames(value_matrix), rownames(G))
    if (length(missing_g))
      stop(length(missing_g), " value_matrix candidate(s) are not present ",
           "in G's dimnames, so their relationship to other candidates ",
           "cannot be evaluated: ",
           paste(utils::head(missing_g, 10), collapse = ", "),
           if (length(missing_g) > 10) ", ..." else "", call. = FALSE)
  }

  # merit_weight > 0 requires merit_score to cover every candidate still in
  # play after the floor/top_candidates filtering above -- unlike the
  # min_sel_value floor (which simply excludes anyone missing), a missing
  # entry here would silently propagate NA into the fitness function itself,
  # so this is a hard error rather than a warning-and-exclude.
  if (merit_weight > 0) {
    missing_m <- setdiff(rownames(value_matrix), names(merit_score))
    if (length(missing_m))
      stop(length(missing_m), " value_matrix candidate(s) have no ",
           "merit_score entry, required because merit_weight > 0: ",
           paste(utils::head(missing_m, 10), collapse = ", "),
           if (length(missing_m) > 10) ", ..." else "", call. = FALSE)
  }

  # -- target_degree -> relatedness-ceiling conversion, now that merit_weight
  # is fully resolved (from merit_priority if supplied) and merit_score is
  # known to cover every remaining candidate -- see
  # .calibrate_relatedness_ceiling() for the full method and
  # ?select_parents_ga's target_degree docs for the direction convention.
  relatedness_ceiling <- NULL
  relatedness_penalty_coef <- 0
  if (!is.null(target_degree)) {
    cal_rel <- .calibrate_relatedness_ceiling(
      value_matrix, G, block_weights, strategy, n_founders,
      target_degree, merit_weight = merit_weight, merit_score = merit_score
    )
    relatedness_ceiling <- cal_rel$ceiling
    # Scaled to dominate any possible coverage/merit gain from violating the
    # ceiling, without a new user-tunable constant -- NOT penalty_weight,
    # which is calibrated for the unrelated cardinality constraint on a
    # totally different scale (integer founder counts vs. a GRM's own
    # units) -- see ?select_parents_ga's target_degree docs.
    relatedness_penalty_coef <- 10 * max(cal_rel$coverage_span, sqrt(.Machine$double.eps))
    if (isTRUE(verbose))
      message("[select_parents_ga] target_degree = ", target_degree,
              " -> relatedness ceiling = ", signif(relatedness_ceiling, 4),
              " (gain end = ", signif(cal_rel$gain_end, 4),
              ", diversity end = ", signif(cal_rel$diversity_end, 4), ").")
  }

  # -- Replicated GA search (rigour: see ?select_parents_ga's "GA rigour"
  # section) -- n_reps independent runs, deterministic per-replicate seeds
  # derived from `seed` when supplied, so the SET of replicates is
  # reproducible even though each explores a different starting point.
  reps <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    this_seed <- if (!is.null(seed)) seed + i - 1L else NULL
    reps[[i]] <- .run_ga_once(
      value_matrix, n_founders, strategy, block_weights,
      popSize, maxiter, run, pmutation, pcrossover, penalty_weight,
      seed = this_seed, verbose = isTRUE(verbose) && i == 1L,
      G = G, coancestry_weight = coancestry_weight,
      merit_score = merit_score, merit_weight = merit_weight,
      relatedness_ceiling = relatedness_ceiling,
      relatedness_penalty_coef = relatedness_penalty_coef
    )
  }

  fitness_vals <- vapply(reps, function(r) r$fitness, numeric(1))
  best_i       <- which.max(fitness_vals)
  best         <- reps[[best_i]]

  if (length(best$selected) != n_founders) {
    warning("[select_parents_ga] Best replicate converged on ",
            length(best$selected), " founders, not the requested ",
            "n_founders = ", n_founders, ". Consider increasing ",
            "maxiter/penalty_weight, or accept this as the GA's best ",
            "trade-off.", call. = FALSE)
  }

  # -- Stability across replicates: per-individual selection frequency, and
  # how much the achievable fitness varied run to run. See "GA rigour".
  sel_lists <- lapply(reps, function(r) r$selected)
  sel_union <- unique(unlist(sel_lists))
  sel_freq  <- vapply(sel_union, function(id)
    mean(vapply(sel_lists, function(s) id %in% s, logical(1))), numeric(1))
  sel_freq  <- sort(sel_freq, decreasing = TRUE)

  stability <- list(
    n_reps             = n_reps,
    fitness_values     = fitness_vals,
    fitness_range      = range(fitness_vals),
    best_rep           = best_i,
    converged          = vapply(reps, function(r) r$converged, logical(1)),
    selection_freq     = sel_freq,
    mean_relationship  = vapply(reps, function(r) r$mean_relationship, numeric(1)),
    mean_merit         = vapply(reps, function(r) r$mean_merit, numeric(1))
  )

  list(
    selected          = best$selected,
    fitness           = best$fitness,
    per_block         = best$per_block,
    strategy          = strategy,
    ga_fit            = best$ga_fit,
    converged         = best$converged,
    cutoff            = cutoff,
    mean_relationship = best$mean_relationship,
    mean_merit        = best$mean_merit,
    merit_weight      = merit_weight,
    merit_priority    = merit_priority,
    coancestry_weight = coancestry_weight,
    target_degree     = target_degree,
    relatedness_ceiling = relatedness_ceiling,
    stability         = stability
  )
}


#' PCA of Parent-Selection Groups
#'
#' @description
#' Projects every genotyped individual onto the top two principal components,
#' coloured by which selection group they belong to (GA-selected,
#' truncation-selected, both, or neither) -- mirrors HapSelect's parent-set
#' PCA plot, which shows where each selection strategy's founders sit
#' relative to the whole population's diversity.
#'
#' Two different diversity spaces can be plotted, and they answer different
#' questions:
#' \itemize{
#'   \item \strong{\code{G} (default)}: genome-wide (or whole-haplotype-set)
#'     relationship-matrix PCA, via eigen-decomposition of \code{G} itself.
#'     Answers "do these selection groups also look diversity-covering
#'     across the whole genome?" -- a useful cross-check, but \code{G} is
#'     \emph{not} the space \code{\link{select_parents_ga}} actually
#'     searched.
#'   \item \strong{\code{feature_matrix}}: PCA (via \code{\link[stats]{prcomp}})
#'     of the individuals x target-blocks local-GEBV/haplotype-value matrix
#'     -- the literal \code{value_matrix} argument
#'     \code{\link{select_parents_ga}} optimised coverage over. Answers "does
#'     GA's selection actually spread out in the space it was asked to
#'     spread out in?" directly, at the cost of only reflecting the target
#'     blocks, not the rest of the genome.
#' }
#' Supply exactly one of the two; \code{G} remains the default so existing
#' calls are unaffected.
#'
#' @param G Genomic or haplotype relationship matrix (n x n), dimnames =
#'   individual IDs, e.g. \code{run_haplotype_prediction()$G} or
#'   \code{\link{compute_haplotype_grm}} output. Leave \code{NULL} (and
#'   supply \code{feature_matrix} instead) to plot in block-feature space
#'   instead of genome-wide diversity space.
#' @param ga_selected Character vector of individual IDs selected by
#'   \code{\link{select_parents_ga}} (its \code{$selected}).
#' @param ts_selected Character vector of individual IDs selected by
#'   \code{\link{truncation_selection}} (its \code{$selected}).
#' @param feature_matrix Individuals x target-blocks numeric matrix
#'   (dimnames rows = individual IDs), e.g. the same \code{value_matrix}
#'   passed to \code{\link{select_parents_ga}}
#'   (\code{res$local_gebv[, top_blocks$block_id, drop = FALSE]}). When
#'   supplied, PCA is computed on this matrix (centred, not scaled --
#'   local-GEBV columns are already on a common trait-value scale) instead
#'   of eigen-decomposing \code{G}. Leave \code{NULL} (the default) to use
#'   \code{G} instead. Supplying both \code{G} and \code{feature_matrix} is
#'   an error -- pick one diversity space per plot.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{select_parents_ga}}, \code{\link{truncation_selection}}
#'
#' @examples
#' \dontrun{
#' res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
#' top <- select_top_blocks(res$block_importance, n = 15)
#' vmat <- res$local_gebv[, top$block_id, drop = FALSE]
#' ga  <- select_parents_ga(vmat, n_founders = 20, seed = 1)
#' ts  <- truncation_selection(res$gebv, n_founders = 20)
#'
#' # Genome-wide diversity cross-check (default):
#' plot_parent_selection_pca(res$G, ga$selected, ts$selected)
#'
#' # The actual space GA searched (target-block local GEBV):
#' plot_parent_selection_pca(ga_selected = ga$selected, ts_selected = ts$selected,
#'                           feature_matrix = vmat)
#' }
#'
#' @export
plot_parent_selection_pca <- function(G = NULL, ga_selected, ts_selected,
                                      feature_matrix = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required: install.packages('ggplot2')", call. = FALSE)
  if (is.null(G) && is.null(feature_matrix))
    stop("Supply either G or feature_matrix (both are NULL).", call. = FALSE)
  if (!is.null(G) && !is.null(feature_matrix))
    stop("Supply exactly one of G or feature_matrix, not both.", call. = FALSE)

  if (!is.null(feature_matrix)) {
    if (!is.matrix(feature_matrix) || is.null(rownames(feature_matrix)))
      stop("feature_matrix must be a matrix with row names (individual IDs).",
           call. = FALSE)
    if (nrow(feature_matrix) < 3L || ncol(feature_matrix) < 2L)
      stop("feature_matrix must have >= 3 individuals (rows) and >= 2 ",
           "blocks (columns) for a 2-component PCA.", call. = FALSE)
    pr  <- stats::prcomp(feature_matrix, center = TRUE, scale. = FALSE)
    pc  <- pr$x[, 1:2, drop = FALSE]
    var_prop <- (pr$sdev^2 / sum(pr$sdev^2))[1:2]
    ids <- rownames(feature_matrix)
    space_title <- "Parent-selection groups in target-block feature space"
  } else {
    if (!is.matrix(G) || is.null(rownames(G)))
      stop("G must be a matrix with row/column names (individual IDs).",
           call. = FALSE)
    eig <- eigen(G, symmetric = TRUE)
    pc  <- eig$vectors[, 1:2, drop = FALSE]
    var_prop <- eig$values / sum(eig$values)
    ids <- rownames(G)
    space_title <- "Parent-selection groups in population diversity space"
  }

  grp <- ifelse(ids %in% ga_selected & ids %in% ts_selected, "Both",
         ifelse(ids %in% ga_selected, "GA-selected",
         ifelse(ids %in% ts_selected, "TS-selected", "Neither")))

  plot_df <- data.frame(
    id    = ids,
    PC1   = pc[, 1L],
    PC2   = pc[, 2L],
    group = factor(grp, levels = c("Neither", "TS-selected", "GA-selected", "Both")),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, colour = group)) +
    ggplot2::geom_point(data = plot_df[plot_df$group == "Neither", , drop = FALSE],
                        alpha = 0.35, size = 1.4) +
    ggplot2::geom_point(data = plot_df[plot_df$group != "Neither", , drop = FALSE],
                        size = 2.6, alpha = 0.9) +
    ggplot2::scale_colour_manual(
      values = c("Neither" = "grey70", "TS-selected" = "#0072B2",
                "GA-selected" = "#D55E00", "Both" = "#009E73"),
      name = "Selection group") +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_prop[1L]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_prop[2L]),
      title = space_title
    ) +
    ggplot2::theme_minimal()
}


#' Cluster Individuals on Retained Principal Components and Cross-Tabulate
#' Selection Coverage
#'
#' @description
#' \code{\link{plot_parent_selection_pca}} only ever looks at PC1/PC2, which
#' can be a small fraction of total genetic variance in a genetically
#' complex/structured population (e.g. 3 distinct clusters captured by only
#' ~26\% cumulative variance on PC1+PC2 combined). This function instead
#' retains as many leading principal components as needed to reach a
#' user-set cumulative-variance threshold (default 95\%), clusters
#' individuals in that higher-dimensional subspace (hierarchical or
#' k-means), and cross-tabulates an arbitrary number of named selection
#' strategies (GA, truncation, OCS, core-collection, or any other ID vector
#' you have) against cluster membership -- turning "does this method
#' represent every genetically distinct sub-group?" into a table of counts
#' and proportions per cluster, per strategy, instead of a visual impression
#' from a 2-axis scatter. Individuals may belong to more than one strategy
#' at once (e.g. selected by both GA and OCS) -- each strategy's
#' representation is counted independently, not as a single mutually
#' exclusive category.
#'
#' @param G Genomic or haplotype relationship matrix (n x n), dimnames =
#'   individual IDs. Same duality as \code{\link{plot_parent_selection_pca}}:
#'   leave \code{NULL} and supply \code{feature_matrix} instead to cluster
#'   in target-block feature space rather than genome-wide diversity space.
#'   \strong{Unlike} \code{plot_parent_selection_pca()}, the retained PC
#'   scores here are eigenvectors scaled by \code{sqrt(eigenvalue)} (the
#'   standard PCA-score convention), not raw unit-norm eigenvectors --
#'   necessary so that Euclidean distance across many retained PCs is
#'   properly variance-weighted for clustering; a 2-axis scatter plot does
#'   not need this because each axis is independently labelled with its own
#'   \% variance.
#' @param feature_matrix Individuals x target-blocks numeric matrix, e.g.
#'   the same \code{value_matrix} passed to \code{\link{select_parents_ga}}.
#'   Supply exactly one of \code{G}/\code{feature_matrix}.
#' @param groups Named list of character vectors, one per selection strategy
#'   to check, e.g. \code{list(GA = ga$selected, TS = ts$selected, OCS =
#'   ocs_parents, Core = core_res$selected)}. Names label that strategy's
#'   columns in the returned \code{table}; IDs not present in \code{G}/
#'   \code{feature_matrix} are ignored (with a message).
#' @param variance_threshold Numeric in (0, 1]. Retain the smallest number
#'   of leading PCs whose cumulative proportion of variance is \code{>=}
#'   this threshold (always at least 2, so \code{\link{plot_selection_clusters}}
#'   can still show a PC1/PC2 view of the retained subspace). Default
#'   \code{0.95}.
#' @param method \code{"hierarchical"} (default -- \code{\link[stats]{hclust}}
#'   with Ward's minimum-variance linkage (\code{"ward.D2"}), then
#'   \code{\link[stats]{cutree}} to \code{n_clusters} groups; deterministic,
#'   no seed needed) or \code{"kmeans"} (\code{\link[stats]{kmeans}},
#'   \code{nstart = 25}; needs \code{seed} for reproducibility).
#' @param n_clusters Integer \code{>= 2}. Number of clusters to cut the
#'   population into. \strong{Required} -- this package does not silently
#'   guess a number of clusters for you. If you don't already have a
#'   biological reason for a specific number (e.g. "this program has 3
#'   founder families"), a reasonable starting point is to try a small range
#'   (e.g. 2-6) and look at \code{cluster_fit} (the \code{hclust}/
#'   \code{kmeans} object returned) with standard diagnostics
#'   (\code{plot(cluster_fit)} for a dendrogram; \code{kmeans}'s
#'   \code{tot.withinss} across a range of \code{k} for an elbow plot) before
#'   settling on \code{n_clusters}.
#' @param seed Integer or \code{NULL}. Only used by \code{method = "kmeans"};
#'   ignored (with a message if supplied) for \code{"hierarchical"}, which is
#'   deterministic.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{cluster}}{Named factor (names = individual IDs), cluster
#'     membership, levels \code{"Cluster 1"} .. \code{"Cluster n_clusters"}.}
#'   \item{\code{table}}{Data frame, one row per cluster: \code{n_total}
#'     (population size in that cluster), \code{pct_population} (that
#'     cluster's \% share of the whole population), then for every name in
#'     \code{groups} a \code{n_<name>} (how many of that strategy's
#'     selections fall in this cluster) and \code{prop_<name>}
#'     (\code{n_<name> / n_total} -- what fraction of \emph{this cluster}
#'     that strategy selected) -- this is the actual answer to "which
#'     clusters does strategy X under-represent?".}
#'   \item{\code{n_pcs_retained}}{Integer, how many leading PCs met
#'     \code{variance_threshold}.}
#'   \item{\code{variance_explained}}{Numeric, actual cumulative proportion
#'     of variance captured by \code{n_pcs_retained} PCs (\code{>=}
#'     \code{variance_threshold}).}
#'   \item{\code{variance_threshold}, \code{method}}{Echoed back.}
#'   \item{\code{space}}{Character, \code{"population diversity space"} or
#'     \code{"target-block feature space"} depending on whether \code{G} or
#'     \code{feature_matrix} was supplied.}
#'   \item{\code{groups}}{Character vector, the \code{names(groups)} echoed
#'     back (the strategy labels used as \code{table} column suffixes).}
#'   \item{\code{pc_scores}}{The retained, variance-scaled PC score matrix
#'     (individuals x \code{n_pcs_retained}) actually clustered on -- input
#'     to \code{\link{plot_selection_clusters}}.}
#'   \item{\code{cluster_fit}}{The raw \code{hclust} or \code{kmeans} object,
#'     for diagnostics (dendrograms, within-cluster sum of squares, etc.).}
#' }
#'
#' @seealso \code{\link{plot_parent_selection_pca}},
#'   \code{\link{plot_selection_clusters}}, \code{\link{select_parents_ga}},
#'   \code{\link{truncation_selection}}
#'
#' @examples
#' \dontrun{
#' res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
#' ga  <- select_parents_ga(res$local_gebv, n_founders = 20, seed = 1)
#' ts  <- truncation_selection(res$gebv, n_founders = 20)
#' cl  <- cluster_selection_groups(res$G, groups = list(GA = ga$selected,
#'                                                      TS = ts$selected),
#'                                 variance_threshold = 0.95,
#'                                 method = "hierarchical", n_clusters = 3)
#' cl$table
#' }
#'
#' @export
cluster_selection_groups <- function(G = NULL, feature_matrix = NULL,
                                     groups,
                                     variance_threshold = 0.95,
                                     method = c("hierarchical", "kmeans"),
                                     n_clusters,
                                     seed = NULL,
                                     verbose = TRUE) {
  method <- match.arg(method)
  if (is.null(G) && is.null(feature_matrix))
    stop("Supply either G or feature_matrix (both are NULL).", call. = FALSE)
  if (!is.null(G) && !is.null(feature_matrix))
    stop("Supply exactly one of G or feature_matrix, not both.", call. = FALSE)
  if (missing(groups) || !is.list(groups) || length(groups) < 1L)
    stop("groups must be a non-empty named list of character vectors (one ",
         "per selection strategy to check), e.g. list(GA = ga$selected, ",
         "TS = ts$selected).", call. = FALSE)
  if (is.null(names(groups)) || any(!nzchar(names(groups))))
    stop("groups must be a NAMED list -- each name labels that selection ",
         "strategy's columns in the returned table.", call. = FALSE)
  if (anyDuplicated(names(groups)))
    stop("groups names must be unique.", call. = FALSE)
  if (!all(vapply(groups, is.character, logical(1L))))
    stop("Every element of groups must be a character vector of individual IDs.",
         call. = FALSE)
  if (!is.numeric(variance_threshold) || length(variance_threshold) != 1L ||
      is.na(variance_threshold) || variance_threshold <= 0 || variance_threshold > 1)
    stop("variance_threshold must be a single number in (0, 1].", call. = FALSE)
  if (missing(n_clusters) || !is.numeric(n_clusters) || length(n_clusters) != 1L ||
      is.na(n_clusters) || n_clusters < 2 || n_clusters != as.integer(n_clusters))
    stop("n_clusters must be a single integer >= 2.", call. = FALSE)
  n_clusters <- as.integer(n_clusters)

  if (!is.null(feature_matrix)) {
    if (!is.matrix(feature_matrix) || is.null(rownames(feature_matrix)))
      stop("feature_matrix must be a matrix with row names (individual IDs).",
           call. = FALSE)
    if (nrow(feature_matrix) < 3L || ncol(feature_matrix) < 2L)
      stop("feature_matrix must have >= 3 individuals (rows) and >= 2 ",
           "blocks (columns).", call. = FALSE)
    pr <- stats::prcomp(feature_matrix, center = TRUE, scale. = FALSE)
    # prcomp()'s $x is already the properly variance-scaled score matrix
    # (columns have variance == eigenvalue) -- no extra scaling needed here.
    pc_all <- pr$x
    var_prop_all <- (pr$sdev^2) / sum(pr$sdev^2)
    ids <- rownames(feature_matrix)
    space_label <- "target-block feature space"
  } else {
    if (!is.matrix(G) || is.null(rownames(G)))
      stop("G must be a matrix with row/column names (individual IDs).",
           call. = FALSE)
    eig <- eigen(G, symmetric = TRUE)
    # Unlike plot_parent_selection_pca() (raw unit-norm eigenvectors, fine
    # for a 2-axis scatter with independently labelled % variance), scale
    # eigenvectors by sqrt(eigenvalue) here so Euclidean distance across
    # potentially many retained PCs is properly variance-weighted --
    # otherwise a near-noise PC50 would contribute to cluster distances on
    # equal footing with PC1, defeating the point of a variance threshold.
    ev_pos <- pmax(eig$values, 0)
    pc_all <- sweep(eig$vectors, 2, sqrt(ev_pos), `*`)
    var_prop_all <- ev_pos / sum(ev_pos)
    ids <- rownames(G)
    space_label <- "population diversity space"
  }

  if (n_clusters > length(ids))
    stop("n_clusters (", n_clusters, ") cannot exceed the number of ",
         "individuals (", length(ids), ").", call. = FALSE)

  cum_var <- cumsum(var_prop_all)
  n_pcs <- which(cum_var >= variance_threshold)[1L]
  if (is.na(n_pcs)) n_pcs <- length(var_prop_all)
  n_pcs <- max(n_pcs, 2L)
  n_pcs <- min(n_pcs, ncol(pc_all))

  pc_scores <- pc_all[, seq_len(n_pcs), drop = FALSE]
  rownames(pc_scores) <- ids
  colnames(pc_scores) <- paste0("PC", seq_len(n_pcs))

  if (isTRUE(verbose))
    message("[cluster_selection_groups] Retaining ", n_pcs, " of ",
            ncol(pc_all), " PCs to reach >= ", round(100 * variance_threshold, 1),
            "% cumulative variance (", round(100 * cum_var[n_pcs], 1),
            "% actually captured) in ", space_label, ".")

  if (method == "hierarchical") {
    if (!is.null(seed) && isTRUE(verbose))
      message("[cluster_selection_groups] seed is ignored for method = ",
              "'hierarchical' (deterministic).")
    d  <- stats::dist(pc_scores)
    hc <- stats::hclust(d, method = "ward.D2")
    cl_int <- stats::cutree(hc, k = n_clusters)
    cluster_fit <- hc
  } else {
    if (!is.null(seed)) set.seed(seed)
    km <- stats::kmeans(pc_scores, centers = n_clusters, nstart = 25L)
    cl_int <- km$cluster
    cluster_fit <- km
  }
  cl_levels <- paste0("Cluster ", sort(unique(cl_int)))
  cl <- factor(paste0("Cluster ", cl_int), levels = cl_levels)
  names(cl) <- ids

  # Cross-tabulate: one row per cluster, n_total/pct_population, then
  # n_<name>/prop_<name> for every strategy in `groups`, computed
  # INDEPENDENTLY per strategy (not as a single mutually exclusive category)
  # -- an individual selected by two strategies at once (e.g. GA and OCS)
  # is counted in both, which is exactly what "does strategy X represent
  # every cluster" needs to answer honestly.
  cluster_counts <- table(cl)
  tab_df <- data.frame(
    Cluster = names(cluster_counts),
    n_total = as.integer(cluster_counts),
    stringsAsFactors = FALSE
  )
  tab_df$pct_population <- 100 * tab_df$n_total / sum(tab_df$n_total)

  for (gname in names(groups)) {
    sel_ids <- groups[[gname]]
    unknown <- setdiff(sel_ids, ids)
    if (length(unknown) && isTRUE(verbose))
      message("[cluster_selection_groups] ", length(unknown), " ID(s) in ",
              "groups[[\"", gname, "\"]] are not among the ", length(ids),
              " clustered individuals and will be ignored: ",
              paste(utils::head(unknown, 5L), collapse = ", "),
              if (length(unknown) > 5L) ", ..." else "")
    in_grp <- ids %in% sel_ids
    n_by_cluster <- vapply(levels(cl), function(cc) sum(in_grp & cl == cc),
                           integer(1L))
    n_col    <- paste0("n_", gname)
    prop_col <- paste0("prop_", gname)
    tab_df[[n_col]]    <- unname(n_by_cluster[tab_df$Cluster])
    tab_df[[prop_col]] <- ifelse(tab_df$n_total > 0,
                                 tab_df[[n_col]] / tab_df$n_total, NA_real_)
  }
  rownames(tab_df) <- NULL

  list(
    cluster            = cl,
    table              = tab_df,
    n_pcs_retained     = n_pcs,
    variance_explained = cum_var[n_pcs],
    variance_threshold = variance_threshold,
    method             = method,
    space              = space_label,
    groups             = names(groups),
    pc_scores          = pc_scores,
    cluster_fit        = cluster_fit
  )
}


#' Plot Individuals Coloured by Genetic Cluster
#'
#' @description
#' Companion plot to \code{\link{cluster_selection_groups}}: PC1 vs PC2 of
#' the same (variance-scaled) retained-PC subspace used for clustering,
#' coloured by cluster membership -- distinct from
#' \code{\link{plot_parent_selection_pca}}, which colours by GA/TS/Both/
#' Neither selection-group membership instead. Read them together: this
#' plot shows where the genetic clusters are; \code{cluster_selection_groups()
#' $table} shows how well each selection method represents each cluster.
#'
#' @param cluster_res List returned by \code{\link{cluster_selection_groups}}.
#' @param save_path Optional file path ending in \code{.pdf}. When supplied,
#'   the plot is also saved via \code{ggplot2::ggsave()} at 300 dpi.
#'   \code{NULL} (default) does not save to disk.
#' @param width,height Numeric, PDF dimensions in inches, used only when
#'   \code{save_path} is supplied. Defaults \code{7}, \code{5.5}.
#'
#' @return A \code{ggplot2} object (returned whether or not \code{save_path}
#'   is supplied).
#'
#' @seealso \code{\link{cluster_selection_groups}},
#'   \code{\link{plot_parent_selection_pca}}
#'
#' @examples
#' \dontrun{
#' cl <- cluster_selection_groups(res$G, ga_selected = ga$selected,
#'                                ts_selected = ts$selected, n_clusters = 3)
#' plot_selection_clusters(cl, save_path = "clusters.pdf")
#' }
#'
#' @export
plot_selection_clusters <- function(cluster_res, save_path = NULL,
                                    width = 7, height = 5.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required: install.packages('ggplot2')", call. = FALSE)
  if (!is.list(cluster_res) || is.null(cluster_res$pc_scores) ||
      is.null(cluster_res$cluster))
    stop("cluster_res must be the list returned by cluster_selection_groups().",
         call. = FALSE)
  if (ncol(cluster_res$pc_scores) < 2L)
    stop("cluster_res$pc_scores has fewer than 2 retained PCs -- cannot ",
         "plot PC1 vs PC2 (should not happen; cluster_selection_groups() ",
         "always retains >= 2 PCs).", call. = FALSE)

  plot_df <- data.frame(
    PC1     = cluster_res$pc_scores[, 1L],
    PC2     = cluster_res$pc_scores[, 2L],
    Cluster = cluster_res$cluster,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, colour = Cluster)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::labs(
      x = "PC1 (retained-PC subspace, variance-scaled)",
      y = "PC2 (retained-PC subspace, variance-scaled)",
      title = paste0("Genetic clusters in ", cluster_res$space),
      subtitle = sprintf(
        "%s clustering, k = %d, %d PCs retained (%.1f%% cumulative variance, threshold %.0f%%)",
        cluster_res$method, nlevels(cluster_res$cluster), cluster_res$n_pcs_retained,
        100 * cluster_res$variance_explained, 100 * cluster_res$variance_threshold
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey30"),
      axis.title    = ggplot2::element_text(size = 13, face = "bold"),
      axis.text     = ggplot2::element_text(size = 11),
      legend.title  = ggplot2::element_text(size = 11, face = "bold"),
      legend.text   = ggplot2::element_text(size = 10)
    )

  if (!is.null(save_path)) {
    if (!grepl("\\.pdf$", save_path, ignore.case = TRUE))
      stop("save_path must end in '.pdf'.", call. = FALSE)
    ggplot2::ggsave(save_path, p, width = width, height = height,
                    device = "pdf", dpi = 300, useDingbats = FALSE)
  }
  p
}
