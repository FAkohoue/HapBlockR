# ==============================================================================
# pareto_selection.R
#
# Multi-objective (Pareto) selection -- strategy 4 of six planned parent-
# selection strategy extensions. Where select_parents_ocs()'s target_degree
# and select_parents_ga()'s coancestry_weight each collapse the merit-vs-
# diversity tradeoff into a SINGLE scalar dial you have to guess a value
# for, this file lets you see the actual tradeoff curve (the Pareto
# frontier: candidate sets/crosses where you cannot improve one objective
# without making another worse) and pick a point off it deliberately,
# instead of guessing a dial setting and hoping.
#
# Two pieces, deliberately layered:
#   pareto_front()          -- a general, dependency-free non-dominated-sort
#     utility (standard algorithm, NSGA-II-style crowding distance as a
#     secondary "how well spread out is the front" diagnostic). Works on
#     ANY data frame with objective columns -- candidate parents, candidate
#     crosses (e.g. usefulness_criterion() output), or anything else.
#   select_parents_pareto() -- sweeps select_parents_ga()'s coancestry_weight
#     across a grid, runs it at each value, and Pareto-filters the resulting
#     (merit, relatedness) points into an empirical frontier. Deliberately
#     built on top of the already-implemented, already-verified GA solver
#     rather than a new from-scratch multi-objective search algorithm (e.g.
#     NSGA-II) -- consistent with this package's general preference for
#     reusing validated machinery over hand-building new optimizers.
# ==============================================================================


#' Pareto (Non-Dominated) Front of a Set of Candidates
#'
#' Given a data frame of candidates (parents, crosses, or anything else) each
#' scored on two or more objectives, flags which candidates are
#' Pareto-optimal (non-dominated): no other candidate is at least as good on
#' every objective and strictly better on at least one. This is the general
#' tool underneath \code{\link{select_parents_pareto}}, but works on any
#' objective columns you give it -- e.g. \code{usefulness_criterion()}
#' output scored on \code{mid_parent_gebv} (maximize) and
#' \code{predicted_variance} (context-dependent), or a multi-trait selection
#' index's separate trait columns.
#'
#' @details
#' Uses the standard \eqn{O(n^2 \times k)} pairwise-dominance algorithm
#' (\eqn{n} = candidates, \eqn{k} = objectives) -- exact, not a heuristic,
#' but not intended for huge candidate sets (low hundreds is comfortable;
#' thousands will be slow). Among the non-dominated set, also computes the
#' NSGA-II crowding distance (Deb et al. 2002) per objective as a secondary
#' diagnostic: candidates near the extremes of the front get \code{Inf};
#' candidates in sparsely populated regions of the front get a larger value
#' than candidates crowded next to near-identical alternatives. This is
#' purely descriptive (which non-dominated points are most "distinct" from
#' their neighbours on the front) -- it does not change which points are
#' Pareto-optimal.
#'
#' @param data Data frame of candidates.
#' @param objectives Character vector of column names in \code{data} to
#'   treat as objectives.
#' @param directions Character vector, same length as \code{objectives} (or
#'   length 1, recycled), each \code{"max"} or \code{"min"}. Default
#'   \code{"max"} for every objective.
#'
#' @return \code{data} with two columns appended: \code{pareto_optimal}
#'   (logical) and \code{crowding_distance} (numeric, \code{NA} for
#'   dominated candidates), sorted with Pareto-optimal candidates first
#'   (by descending crowding distance among them).
#'
#' @references
#' Deb, K., Pratap, A., Agarwal, S. & Meyarivan, T. (2002). A fast and
#' elitist multiobjective genetic algorithm: NSGA-II. \emph{IEEE Transactions
#' on Evolutionary Computation}, 6, 182-197.
#'
#' @seealso \code{\link{select_parents_pareto}}, \code{\link{usefulness_criterion}}
#' @export
pareto_front <- function(data, objectives, directions = "max") {
  if (!is.data.frame(data))
    stop("data must be a data frame.", call. = FALSE)
  if (!all(objectives %in% names(data)))
    stop("objectives not found in data: ",
         paste(setdiff(objectives, names(data)), collapse = ", "),
         call. = FALSE)
  if (length(directions) == 1L)
    directions <- rep(directions, length(objectives))
  if (length(directions) != length(objectives))
    stop("directions must have length 1 or the same length as objectives.",
         call. = FALSE)
  if (!all(directions %in% c("max", "min")))
    stop("directions must be 'max' or 'min'.", call. = FALSE)

  n <- nrow(data)
  M <- as.matrix(data[, objectives, drop = FALSE])
  storage.mode(M) <- "double"
  for (j in seq_along(objectives))
    if (directions[j] == "min") M[, j] <- -M[, j]

  dominated <- rep(FALSE, n)
  if (n > 1L) {
    for (i in seq_len(n)) {
      for (k in seq_len(n)) {
        if (i == k) next
        if (all(M[k, ] >= M[i, ]) && any(M[k, ] > M[i, ])) {
          dominated[i] <- TRUE
          break
        }
      }
    }
  }

  out <- data
  out$pareto_optimal    <- !dominated
  out$crowding_distance <- NA_real_

  idx <- which(out$pareto_optimal)
  if (length(idx) >= 2L) {
    cd <- rep(0, length(idx))
    for (j in seq_along(objectives)) {
      vals <- M[idx, j]
      rng  <- diff(range(vals))
      if (rng < .Machine$double.eps) next
      ord <- order(vals)
      cd[ord[1]]           <- Inf
      cd[ord[length(ord)]] <- Inf
      if (length(ord) > 2L) {
        for (p in 2:(length(ord) - 1L))
          cd[ord[p]] <- cd[ord[p]] +
            (vals[ord[p + 1L]] - vals[ord[p - 1L]]) / rng
      }
    }
    out$crowding_distance[idx] <- cd
  }

  ord <- order(-out$pareto_optimal,
              -ifelse(is.na(out$crowding_distance), -Inf, out$crowding_distance))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Pareto Frontier of Parent Sets: Merit vs. Relatedness
#'
#' Sweeps \code{\link{select_parents_ga}}'s \code{coancestry_weight} across a
#' grid of values, runs the GA at each one, and Pareto-filters the resulting
#' (merit, relatedness) points into an empirical frontier. This displays the
#' realised gain-vs-diversity trade-off for the candidate population so that
#' the programme can apply a declared policy to select a frontier point.
#'
#' @details
#' The function uses \code{\link{select_parents_ga}} as the solver at each
#' grid point and \code{\link{pareto_front}} to retain the non-dominated
#' results. This weight-sweep formulation supplies complete GA diagnostics at
#' every point and removes any solution that is lower in merit and higher in
#' relatedness than another solution.
#'
#' Runtime is the sum of every grid point's GA run -- with the default 6
#' grid points and \code{n_reps = 3}, that is 18 GA searches. Lower
#' \code{n_reps} or narrow \code{coancestry_weights} for faster, coarser
#' sweeps; the defaults favour a broad first look over speed.
#'
#' @param value_matrix,n_founders,strategy,block_weights,top_candidates,popSize,maxiter,run,pmutation,pcrossover
#'   Passed through to \code{\link{select_parents_ga}} at every grid point;
#'   see its documentation. The strategy values include
#'   \code{"OHS"} (Optimal Haplotype Selection) and
#'   \code{"OPV"} (Optimal Population Value).
#' @param G Relationship/kinship matrix, required (used both for the
#'   coancestry penalty during each GA run and to report each resulting
#'   set's realised \code{mean_relationship} for the frontier, including at
#'   \code{coancestry_weight = 0}).
#' @param coancestry_weights Numeric vector of \code{coancestry_weight}
#'   values to sweep. Default \code{c(0, 0.25, 0.5, 1, 2, 4)} -- a broad
#'   first pass; problem-specific, since block-coverage and relationship
#'   scores have no common natural scale (same caveat as
#'   \code{coancestry_weight} itself in \code{?select_parents_ga}).
#' @param merit Optional named numeric vector (e.g. whole-genome GEBV or a
#'   selection index) used only to report each resulting set's
#'   \code{mean_merit} on the frontier and as the merit axis for Pareto
#'   filtering. If \code{NULL} (default), the GA's own block-coverage
#'   \code{fitness} is used as the merit axis instead -- a valid but less
#'   directly interpretable stand-in for whole-genome merit.
#' @param n_reps Integer, default \code{3L}. Passed to
#'   \code{\link{select_parents_ga}} at each grid point (lower than that
#'   function's own default of 5, to keep the sweep's total runtime
#'   reasonable).
#' @param seed Optional integer. If supplied, grid point \code{i} uses
#'   \code{seed + i - 1L} for reproducibility across the sweep.
#' @param verbose Logical, default \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{frontier}}{Data frame, one row per grid point, with
#'       \code{coancestry_weight}, \code{fitness}, \code{mean_relationship},
#'       \code{mean_merit} (\code{NA} if \code{merit} not supplied),
#'       \code{n_selected}, \code{converged}, \code{run_index} (row index
#'       into \code{runs}), and the \code{pareto_optimal}/
#'       \code{crowding_distance} columns from \code{\link{pareto_front}}.}
#'     \item{\code{runs}}{List of the full \code{\link{select_parents_ga}}
#'       return value at each grid point (in original sweep order, indexed
#'       by \code{frontier$run_index}) -- use this to get the actual
#'       \code{$selected} individual IDs for any frontier point you choose.}
#'   }
#'
#' @seealso \code{\link{select_parents_ga}}, \code{\link{pareto_front}},
#'   \code{\link{select_parents_ocs}}
#' @export
select_parents_pareto <- function(
    value_matrix, n_founders,
    strategy       = c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS"),
    block_weights  = NULL,
    top_candidates = NULL,
    G,
    coancestry_weights = c(0, 0.25, 0.5, 1, 2, 4),
    merit           = NULL,
    popSize         = 100L,
    maxiter         = 200L,
    run             = 50L,
    pmutation       = 0.1,
    pcrossover      = 0.8,
    n_reps          = 3L,
    seed            = NULL,
    verbose         = TRUE
) {
  strategy <- match.arg(strategy)
  if (missing(G) || is.null(G))
    stop("G is required (a relationship/kinship matrix) -- ",
         "select_parents_pareto() sweeps the merit-vs-relatedness tradeoff, ",
         "which needs G even at coancestry_weight = 0 (to report each ",
         "candidate set's realised mean_relationship for the frontier).",
         call. = FALSE)
  if (!length(coancestry_weights))
    stop("coancestry_weights must have at least one value.", call. = FALSE)

  if (isTRUE(verbose))
    message("[select_parents_pareto] Sweeping ", length(coancestry_weights),
            " coancestry_weight value(s): ",
            paste(coancestry_weights, collapse = ", "))

  runs <- vector("list", length(coancestry_weights))
  rows <- vector("list", length(coancestry_weights))

  for (i in seq_along(coancestry_weights)) {
    w <- coancestry_weights[i]
    if (isTRUE(verbose))
      message("[select_parents_pareto] coancestry_weight = ", w, " (", i,
              "/", length(coancestry_weights), ") ...")
    res <- select_parents_ga(
      value_matrix = value_matrix, n_founders = n_founders, strategy = strategy,
      block_weights = block_weights, top_candidates = top_candidates,
      popSize = popSize, maxiter = maxiter, run = run,
      pmutation = pmutation, pcrossover = pcrossover,
      seed = if (is.null(seed)) NULL else seed + i - 1L,
      verbose = FALSE, n_reps = n_reps, G = G, coancestry_weight = w
    )
    runs[[i]] <- res
    mean_merit <- if (!is.null(merit)) mean(merit[res$selected], na.rm = TRUE) else NA_real_
    rows[[i]] <- data.frame(
      coancestry_weight = w,
      fitness           = res$fitness,
      mean_relationship = res$mean_relationship,
      mean_merit        = mean_merit,
      n_selected        = length(res$selected),
      converged         = isTRUE(res$converged),
      stringsAsFactors  = FALSE
    )
  }

  frontier <- do.call(rbind, rows)
  frontier$run_index <- seq_len(nrow(frontier))

  obj_col  <- if (!is.null(merit)) "mean_merit" else "fitness"
  frontier <- pareto_front(frontier, objectives = c(obj_col, "mean_relationship"),
                           directions = c("max", "min"))

  if (isTRUE(verbose))
    message("[select_parents_pareto] ", sum(frontier$pareto_optimal), " of ",
            nrow(frontier), " sweep point(s) are Pareto-optimal (", obj_col,
            " vs. mean_relationship).")

  list(frontier = frontier, runs = runs)
}
