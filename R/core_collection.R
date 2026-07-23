# ==============================================================================
# core_collection.R
#
# Core-collection / diversity-maximizing subset selection -- strategy 6 (the
# last) of six planned parent-selection strategy extensions. Every other
# selection tool in this package (truncation_selection(), select_parents_ga(),
# select_parents_ocs(), usefulness_criterion()) is ultimately built around a
# MERIT criterion, with diversity/relatedness folded in as a penalty or
# constraint. select_core_collection() inverts that: the primary objective
# is diversity itself -- select a subset that best REPRESENTS the genetic
# diversity of a panel, the way a genebank core collection, a training/
# reference panel, or a founder set for a new diversification programme
# needs. An optional merit floor is still available (reusing
# .apply_merit_floor() from parent_selection.R) so you can ask for "the most
# diverse subset among those that already clear a merit bar," rather than
# diversity in a vacuum.
#
# Implemented as the classical farthest-point / maximin greedy heuristic
# (Gonzalez 1985): start from the single most distant pair, then repeatedly
# add whichever remaining candidate is farthest from everything already
# selected. This has a proven 2-approximation guarantee for the maximin
# diversity objective, is the standard "M strategy" in the core-collection
# literature (Schoen & Brown 1993), and -- unlike this package's other
# statistical-model-heavy functions -- is a simple, exactly-verifiable
# greedy algorithm with no external dependency and no execution-testing
# risk. A mean-distance ("MD strategy") variant is also available.
# ==============================================================================

#' Core-Collection / Diversity-Maximizing Subset Selection
#'
#' Selects a subset of \code{n_core} individuals from a candidate panel that
#' best represents its genetic diversity -- the classical core-collection
#' problem (Schoen & Brown 1993; Gonzalez 1985), useful for genebank
#' curation, building a diverse training/reference panel, or choosing a
#' broad founder set for a new diversification programme. Unlike this
#' package's other parent-selection tools, diversity itself is the primary
#' objective here, not a penalty on top of a merit criterion -- though an
#' optional merit floor lets you restrict to "the most diverse subset among
#' candidates that already clear a merit bar" (see \code{min_sel_value}).
#'
#' @section Two greedy strategies:
#' \describe{
#'   \item{\code{"maximin"}}{(Default.) Farthest-point traversal (Gonzalez
#'     1985): start from the single most genetically distant pair, then
#'     repeatedly add whichever remaining candidate has the largest MINIMUM
#'     distance to everything already selected. Maximizes the smallest
#'     pairwise distance in the final set (guards against near-duplicate
#'     individuals slipping in), with a proven 2-approximation guarantee.
#'     The standard "M strategy" in the core-collection literature.}
#'   \item{\code{"mean_distance"}}{At each step, adds whichever remaining
#'     candidate most increases the selected set's MEAN pairwise distance
#'     (the "MD strategy"). Tends to spread the selection more evenly across
#'     the whole diversity space rather than prioritising the single most
#'     extreme outliers; more expensive to compute (recomputed per candidate
#'     per step) but still tractable for typical core-collection sizes.}
#' }
#' Neither is an exact solver for its objective (both are the standard
#' greedy heuristics used throughout the core-collection literature, not a
#' guaranteed-optimal search) -- see \code{\link{validate_crosses_exact}} if
#' you need a true optimum on a small-enough problem (that function targets
#' cross selection specifically, not this subset-diversity problem, but
#' shares the same "exact validation of a heuristic" philosophy).
#'
#' @param G Dimnamed relationship or distance matrix (row/column names =
#'   individual IDs), e.g. from \code{\link{compute_haplotype_grm}}.
#' @param n_core Integer. Size of the core collection to select.
#' @param type Character, one of \code{"relationship"} (default) or
#'   \code{"distance"}. If \code{"relationship"}, \code{G} is converted to a
#'   genetic distance matrix via \eqn{D_{ij} = G_{ii} + G_{jj} - 2G_{ij}}
#'   (the exact identity relating a Gram/relationship matrix to squared
#'   Euclidean distance in the space it represents -- not an approximation).
#'   If \code{"distance"}, \code{G} is used as a distance matrix directly.
#' @param strategy Character, one of \code{"maximin"} (default) or
#'   \code{"mean_distance"}. See Details above.
#' @param merit Optional named numeric vector (e.g. GEBV or a selection
#'   index), names = individual IDs. Only used together with
#'   \code{min_sel_value} to pre-filter the candidate pool; does not
#'   otherwise influence which individuals are chosen (this function
#'   optimizes diversity, not merit, among whichever candidates remain
#'   eligible).
#' @param min_sel_value,min_sel_mode Optional merit floor applied to
#'   \code{merit} before diversity selection, via the same
#'   \code{.apply_merit_floor()} logic used by
#'   \code{\link{truncation_selection}}/\code{\link{select_parents_ga}} --
#'   \code{min_sel_mode} one of \code{"value"}, \code{"percentile"},
#'   \code{"sd_below_mean"}. Both ignored if \code{merit} is \code{NULL}.
#' @param seed Optional integer. Currently only relevant for the
#'   degenerate \code{n_core = 1} case (no merit supplied), where the single
#'   selected individual is otherwise chosen at random; included for
#'   reproducibility and API consistency with this package's other
#'   selection functions.
#' @param verbose Logical, default \code{TRUE}.
#'
#' @return A list with \code{selected} (character vector of chosen
#'   individual IDs, in selection order), \code{n_core}, \code{strategy},
#'   \code{mean_distance} and \code{min_distance} (of the final selected
#'   set), and \code{trace} (data frame, one row per selection step:
#'   \code{step}, \code{added}, \code{criterion} -- the maximin distance or
#'   mean-distance value achieved at that step, useful for plotting how
#'   diversity accumulates as the core collection grows).
#'
#' @references
#' Schoen, D.J. & Brown, A.H.D. (1993). Conservation of allelic richness in
#' wild crop relatives is aided by assessment of genetic markers.
#' \emph{Proceedings of the National Academy of Sciences}, 90, 10623-10627.
#'
#' Gonzalez, T.F. (1985). Clustering to minimize the maximum intercluster
#' distance. \emph{Theoretical Computer Science}, 38, 293-306.
#'
#' @seealso \code{\link{compute_haplotype_grm}}, \code{\link{select_parents_ga}},
#'   \code{\link{validate_crosses_exact}}
#' @export
select_core_collection <- function(
    G,
    n_core,
    type          = c("relationship", "distance"),
    strategy      = c("maximin", "mean_distance"),
    merit         = NULL,
    min_sel_value = NULL,
    min_sel_mode  = c("value", "percentile", "sd_below_mean"),
    seed          = NULL,
    verbose       = TRUE
) {
  type         <- match.arg(type)
  strategy     <- match.arg(strategy)
  min_sel_mode <- match.arg(min_sel_mode)

  if (!is.matrix(G) || is.null(rownames(G)) || is.null(colnames(G)))
    stop("G must be a dimnamed matrix (row/column names = individual IDs).",
         call. = FALSE)
  if (!identical(rownames(G), colnames(G)))
    stop("G must have identical row and column names.", call. = FALSE)

  ids <- rownames(G)

  if (!is.null(merit)) {
    if (is.null(names(merit)))
      stop("merit must be a named numeric vector (names = individual IDs).",
           call. = FALSE)
    common <- intersect(ids, names(merit))
    if (!length(common))
      stop("No individuals in common between G and merit.", call. = FALSE)
    if (!is.null(min_sel_value)) {
      floor_res <- .apply_merit_floor(merit[common], min_sel_value,
                                      min_sel_mode, label = "candidate")
      # floor_res$eligible is already the filtered vector of eligible ID
      # names (names(score)[score >= cutoff]) -- it is NOT an index into
      # `common`. `common` itself is an unnamed character vector (from
      # intersect()), so the previous `common[floor_res$eligible]` performed
      # character-based indexing against a vector with no `names()`, which
      # silently returns NA for every element rather than an error -- caught
      # by test-core-collection.R's merit-floor test, not by any prior
      # execution (this branch was previously untested).
      common <- floor_res$eligible
      if (isTRUE(verbose))
        message("[select_core_collection] Merit floor: ", length(common),
                " of ", length(ids), " individual(s) eligible (cutoff = ",
                round(floor_res$cutoff, 4), ").")
    }
    ids <- common
  }

  if (length(ids) < 2L)
    stop("Fewer than 2 eligible individuals; cannot build a core collection.",
         call. = FALSE)
  if (is.null(n_core) || n_core < 1L)
    stop("n_core must be >= 1.", call. = FALSE)
  if (n_core > length(ids))
    stop("n_core (", n_core, ") exceeds the number of eligible individuals ",
         "(", length(ids), ").", call. = FALSE)

  Gs <- G[ids, ids, drop = FALSE]
  D <- if (type == "relationship") {
    dg <- diag(Gs)
    outer(dg, dg, "+") - 2 * Gs
  } else {
    Gs
  }
  D[D < 0] <- 0
  diag(D) <- 0

  if (max(D) <= .Machine$double.eps)
    stop("All pairwise distances are ~0 (every eligible individual is ",
         "genetically identical under G) -- cannot build a meaningful core ",
         "collection.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)
  n <- length(ids)

  if (n_core == 1L) {
    sel_idx <- if (!is.null(merit)) which.max(merit[ids]) else sample.int(n, 1L)
    return(list(
      selected      = ids[sel_idx],
      n_core        = 1L,
      strategy      = strategy,
      mean_distance = NA_real_,
      min_distance  = NA_real_,
      trace         = data.frame(step = 1L, added = ids[sel_idx],
                                 criterion = NA_real_, stringsAsFactors = FALSE)
    ))
  }

  first_pair <- which(D == max(D), arr.ind = TRUE)[1, ]
  sel <- unname(c(first_pair[1], first_pair[2]))
  trace <- data.frame(step = c(1L, 2L), added = ids[sel],
                      criterion = c(NA_real_, D[sel[1], sel[2]]),
                      stringsAsFactors = FALSE)

  if (strategy == "maximin") {
    min_dist_to_sel <- apply(D[, sel, drop = FALSE], 1, min)
    while (length(sel) < n_core) {
      min_dist_to_sel[sel] <- -Inf
      nxt <- which.max(min_dist_to_sel)
      trace <- rbind(trace, data.frame(step = length(sel) + 1L,
                                       added = ids[nxt],
                                       criterion = min_dist_to_sel[nxt],
                                       stringsAsFactors = FALSE))
      sel <- c(sel, nxt)
      min_dist_to_sel <- pmin(min_dist_to_sel, D[, nxt])
    }
  } else {
    while (length(sel) < n_core) {
      remaining <- setdiff(seq_len(n), sel)
      scores <- vapply(remaining, function(cand) {
        cand_set <- c(sel, cand)
        sub <- D[cand_set, cand_set, drop = FALSE]
        mean(sub[upper.tri(sub)])
      }, numeric(1L))
      best <- which.max(scores)
      nxt  <- remaining[best]
      trace <- rbind(trace, data.frame(step = length(sel) + 1L,
                                       added = ids[nxt],
                                       criterion = scores[best],
                                       stringsAsFactors = FALSE))
      sel <- c(sel, nxt)
    }
  }

  sub_final <- D[sel, sel, drop = FALSE]
  offdiag   <- sub_final[upper.tri(sub_final)]
  rownames(trace) <- NULL

  list(
    selected      = ids[sel],
    n_core        = n_core,
    strategy      = strategy,
    mean_distance = mean(offdiag),
    min_distance  = min(offdiag),
    trace         = trace
  )
}
