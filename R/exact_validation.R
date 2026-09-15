# ==============================================================================
# exact_validation.R
#
# Exact optimisation as a validation check -- strategy 5 of six planned
# parent-selection strategy extensions. Every mate-allocation tool this
# package uses or wraps (SimpleMating::selectCrosses()/GOCS(), AlphaMate,
# select_parents_ga()'s GA search) is a HEURISTIC: none of them is
# guaranteed to find the true optimum. This file provides a genuine EXACT
# solver -- binary integer linear programming via the lpSolve package -- for
# the same core cross-selection problem (choose n_cross crosses from a
# candidate list, maximizing a criterion, subject to a per-parent maximum
# contribution and an optional relatedness-based culling threshold), so you
# can sanity-check a heuristic plan against the true optimum on a
# small-enough candidate set.
#
# Deliberately scoped as a VALIDATION tool, not a production workflow
# replacement: integer programming does not scale the way a GA or a greedy
# heuristic does, so this is intended for tens to low hundreds of candidate
# crosses, not the thousands a full breeding program might generate. See
# max_vars below.
#
# Deliberate scope decision: a per-parent MINIMUM contribution (SimpleMating's
# min.cross) is NOT implemented as a hard constraint here. Enforcing "min_cross
# applies only to parents that end up used at all" exactly requires a
# conditional (big-M) formulation with an extra binary indicator variable per
# parent -- solvable, but a well-known source of subtle bugs if the big-M
# constant isn't chosen carefully, and this file's whole point is to be the
# trustworthy exact reference, not another thing to doubt. max_cross (a
# simple, unconditional per-parent upper bound) and the total n_cross count
# are both implemented as straightforward linear constraints instead.
# ==============================================================================

#' Exact (ILP) Cross Selection: A Validation Check for Heuristic Mating Plans
#'
#' Solves the cross-selection problem -- choose exactly \code{n_cross}
#' crosses from a candidate list, maximizing a criterion, subject to a
#' per-parent maximum contribution and an optional relatedness-based culling
#' threshold -- as a binary integer linear program via
#' \code{lpSolve::lp()}, guaranteeing the true optimum (not a heuristic
#' approximation). Use this to check how close a heuristic mating plan (from
#' \code{\link{select_parents_ocs}}, \code{SimpleMating::selectCrosses()}/
#' \code{GOCS()}, or AlphaMate) came to the best achievable plan under the
#' same constraints, on a small enough candidate set that solving exactly is
#' practical.
#'
#' @section What this is and is not:
#' This is a validation/sanity-check tool, not a replacement for
#' \code{\link{select_parents_ocs}} in normal use. Integer programming does
#' not scale to the size of a real candidate cross list (all pairwise
#' combinations of a large parent set can easily run into the tens of
#' thousands of candidate crosses) the way a GA or greedy heuristic does --
#' see \code{max_vars}. A per-parent MINIMUM contribution (SimpleMating's
#' \code{min.cross}) is deliberately not supported as a hard constraint; see
#' the source-level comment in \code{R/exact_validation.R} for why.
#'
#' @param data Data frame of candidate crosses (e.g.
#'   \code{\link{usefulness_criterion}} output, or a
#'   \code{SimpleMating::getUsefA()}/\code{selectCrosses()} table).
#' @param n_cross Integer. Exact number of crosses the plan must contain.
#' @param max_cross Optional integer. Maximum number of crosses any single
#'   parent can participate in. \code{NULL} (default) leaves this
#'   unconstrained -- the true optimum may then concentrate heavily on very
#'   few parents; supply the same value your heuristic plan used for an
#'   apples-to-apples comparison.
#' @param culling_pairwise_k Optional numeric. Candidate crosses with
#'   relatedness above this value are excluded before solving (matching
#'   \code{SimpleMating}/\code{select_parents_ocs()}'s culling convention).
#'   Requires \code{relatedness_col} (already in \code{data}) or \code{G}
#'   (to compute relatedness per cross on the fly).
#' @param parent1_col,parent2_col Character, default \code{"parent1"}/
#'   \code{"parent2"}. Column names in \code{data} identifying each cross's
#'   two parents.
#' @param criterion_col Character, default \code{"UC"}. Column in
#'   \code{data} to maximize (e.g. \code{"UC"} from
#'   \code{\link{usefulness_criterion}}, or \code{"mid_parent_gebv"}).
#' @param relatedness_col Optional character. Column in \code{data} giving
#'   each cross's relatedness value, used for \code{culling_pairwise_k}. If
#'   \code{NULL} and \code{G} is supplied, computed automatically as
#'   \code{G[parent1, parent2]} per row.
#' @param G Optional dimnamed relationship matrix, used to compute
#'   relatedness per cross when \code{relatedness_col} is not already in
#'   \code{data}.
#' @param heuristic_plan Optional data frame (same \code{parent1_col}/
#'   \code{parent2_col} convention) -- a heuristic mating plan to compare
#'   against the exact optimum. If supplied, the return value includes the
#'   heuristic plan's total criterion and its percentage gap below the exact
#'   optimum.
#' @param max_vars Integer, default \code{2000L}. Safety cap on the number
#'   of candidate crosses (after culling) the solver will attempt --
#'   integer programming is exact but can become slow well before this in
#'   the worst case; lower it if solving is too slow, or pre-filter
#'   \code{data} (e.g. a tighter \code{culling_pairwise_k}, or restrict to
#'   your top-ranked crosses by \code{criterion_col}) rather than raising it
#'   blindly.
#' @param verbose Logical, default \code{TRUE}.
#'
#' @return A list inheriting from \code{HapBlockR_exact_cross_validation}
#'   and \code{hapblockr_result}, with \code{exact_plan} (data frame: the optimal cross
#'   selection, a subset of \code{data}'s rows), \code{exact_objective}
#'   (the true optimal total criterion), \code{n_candidates} (candidate
#'   crosses considered after culling), \code{status} (lpSolve's solver
#'   status; \code{0} = optimal solution found), and, if
#'   \code{heuristic_plan} was supplied, \code{heuristic_objective} and
#'   \code{gap_pct} (the heuristic plan's percentage shortfall below the
#'   exact optimum).
#'
#' @seealso \code{\link{select_parents_ocs}}, \code{\link{usefulness_criterion}}
#' @export
validate_crosses_exact <- function(
    data,
    n_cross,
    max_cross           = NULL,
    culling_pairwise_k  = NULL,
    parent1_col         = "parent1",
    parent2_col         = "parent2",
    criterion_col       = "UC",
    relatedness_col     = NULL,
    G                   = NULL,
    heuristic_plan      = NULL,
    max_vars            = 2000L,
    verbose             = TRUE
) {
  result_call <- match.call()
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("lpSolve is required for validate_crosses_exact(). Install with: ",
         "install.packages('lpSolve')", call. = FALSE)
  if (!is.data.frame(data))
    stop("data must be a data frame.", call. = FALSE)
  data_input <- data

  req  <- c(parent1_col, parent2_col, criterion_col)
  miss <- setdiff(req, names(data))
  if (length(miss))
    stop("data is missing column(s): ", paste(miss, collapse = ", "),
         " (set parent1_col/parent2_col/criterion_col to match your data's ",
         "column names).", call. = FALSE)

  input_rows <- seq_len(nrow(data))
  missing_criterion <- is.na(data[[criterion_col]])
  excluded_records <- data.frame(
    input_row = input_rows[missing_criterion],
    reason = rep("missing_selection_criterion", sum(missing_criterion)),
    stringsAsFactors = FALSE
  )
  data <- data[!missing_criterion, , drop = FALSE]
  input_rows <- input_rows[!missing_criterion]
  if (!nrow(data))
    stop("No rows in data have a non-NA criterion_col value.", call. = FALSE)

  p1 <- as.character(data[[parent1_col]])
  p2 <- as.character(data[[parent2_col]])

  Kvec <- NULL
  if (!is.null(relatedness_col) && relatedness_col %in% names(data)) {
    Kvec <- data[[relatedness_col]]
  } else if (!is.null(G)) {
    miss_g <- setdiff(unique(c(p1, p2)), rownames(G))
    if (length(miss_g))
      stop(length(miss_g), " parent ID(s) missing from G: ",
           paste(utils::head(miss_g, 10L), collapse = ", "), call. = FALSE)
    Kvec <- vapply(seq_len(nrow(data)), function(i) G[p1[i], p2[i]], numeric(1L))
  } else if (!is.null(culling_pairwise_k)) {
    stop("culling_pairwise_k was supplied but no relatedness information ",
         "is available -- provide relatedness_col (a column already in ",
         "data) or G (a relationship matrix to compute it from).",
         call. = FALSE)
  }

  if (!is.null(culling_pairwise_k) && !is.null(Kvec)) {
    keep <- which(Kvec <= culling_pairwise_k)
    removed <- setdiff(seq_along(Kvec), keep)
    if (length(removed)) {
      excluded_records <- rbind(
        excluded_records,
        data.frame(
          input_row = input_rows[removed],
          reason = rep(
            "pairwise_relatedness_above_culling_threshold",
            length(removed)
          ),
          stringsAsFactors = FALSE
        )
      )
    }
    data <- data[keep, , drop = FALSE]
    input_rows <- input_rows[keep]
    p1 <- p1[keep]; p2 <- p2[keep]
  }

  n_vars <- nrow(data)
  if (!n_vars)
    stop("No candidate crosses remain after culling by culling_pairwise_k.",
         call. = FALSE)
  if (n_vars > max_vars)
    stop(n_vars, " candidate crosses exceed max_vars = ", max_vars, ". ",
         "This is an EXACT solver intended for small-scale validation, not ",
         "a primary large-N workflow tool -- pre-filter data (a tighter ",
         "culling_pairwise_k, or restrict to your top-ranked crosses by ",
         "criterion_col) rather than raising max_vars blindly.",
         call. = FALSE)
  if (n_cross > n_vars)
    stop("n_cross (", n_cross, ") exceeds the number of candidate crosses ",
         "available (", n_vars, ") after culling.", call. = FALSE)

  if (isTRUE(verbose)) {
    message("[validate_crosses_exact] ", n_vars, " candidate cross(es) ",
            "after culling; solving exactly for n_cross = ", n_cross,
            if (!is.null(max_cross)) paste0(", max_cross = ", max_cross)
            else " (max_cross unconstrained)", " ...")
  }

  Y <- as.numeric(data[[criterion_col]])
  parents <- unique(c(p1, p2))

  const.mat <- matrix(1, nrow = 1L, ncol = n_vars)
  const.dir <- "="
  const.rhs <- n_cross

  if (!is.null(max_cross)) {
    cap_mat <- t(vapply(parents, function(pp)
      as.numeric(p1 == pp | p2 == pp), numeric(n_vars)))
    const.mat <- rbind(const.mat, cap_mat)
    const.dir <- c(const.dir, rep("<=", length(parents)))
    const.rhs <- c(const.rhs, rep(max_cross, length(parents)))
  }

  sol <- tryCatch(
    lpSolve::lp(direction = "max", objective.in = Y, const.mat = const.mat,
               const.dir = const.dir, const.rhs = const.rhs, all.bin = TRUE),
    error = function(e)
      stop("lpSolve::lp() failed: ", conditionMessage(e), call. = FALSE)
  )

  if (sol$status != 0)
    stop("lpSolve::lp() did not find an optimal solution (status = ",
         sol$status, "; typically means the problem is infeasible under ",
         "the given n_cross/max_cross/culling_pairwise_k combination). ",
         "Try relaxing max_cross or culling_pairwise_k, or lowering ",
         "n_cross.", call. = FALSE)

  chosen <- which(round(sol$solution) == 1)
  exact_plan <- data[chosen, , drop = FALSE]
  rownames(exact_plan) <- NULL
  exact_objective <- sol$objval

  out <- list(exact_plan = exact_plan, exact_objective = exact_objective,
             n_candidates = n_vars, status = sol$status)

  if (!is.null(heuristic_plan)) {
    hp1 <- as.character(heuristic_plan[[parent1_col]])
    hp2 <- as.character(heuristic_plan[[parent2_col]])
    hkey <- ifelse(hp1 < hp2, paste(hp1, hp2), paste(hp2, hp1))
    dkey <- ifelse(p1 < p2, paste(p1, p2), paste(p2, p1))
    match_idx <- match(hkey, dkey)
    n_unmatched <- sum(is.na(match_idx))
    if (n_unmatched && isTRUE(verbose))
      message("[validate_crosses_exact] ", n_unmatched, " of ",
              length(hkey), " heuristic_plan cross(es) not found in the ",
              "(possibly culled) candidate list -- excluded from the ",
              "heuristic objective total.")
    heuristic_objective <- sum(Y[match_idx], na.rm = TRUE)
    gap_pct <- 100 * (exact_objective - heuristic_objective) /
      abs(exact_objective)
    out$heuristic_objective <- heuristic_objective
    out$gap_pct <- gap_pct
    if (isTRUE(verbose))
      message("[validate_crosses_exact] Exact optimum = ",
              round(exact_objective, 4), "; heuristic plan = ",
              round(heuristic_objective, 4), " (", round(gap_pct, 2),
              "% below exact).")
  } else if (isTRUE(verbose)) {
    message("[validate_crosses_exact] Exact optimum = ",
            round(exact_objective, 4), " across ", nrow(exact_plan),
            " cross(es).")
  }

  parent_load <- table(c(
    as.character(exact_plan[[parent1_col]]),
    as.character(exact_plan[[parent2_col]])
  ))
  contract_inputs <- list(candidate_crosses = data_input)
  if (!is.null(G)) contract_inputs$relationship_matrix <- G
  if (!is.null(heuristic_plan))
    contract_inputs$heuristic_plan <- heuristic_plan
  plan_score <- as.numeric(exact_plan[[criterion_col]])
  decision_table <- data.frame(
    parent1 = as.character(exact_plan[[parent1_col]]),
    parent2 = as.character(exact_plan[[parent2_col]]),
    score = plan_score,
    rank = rank(-plan_score, ties.method = "first"),
    selected = TRUE,
    stringsAsFactors = FALSE
  )
  uncertainty <- data.frame(
    exact_objective = exact_objective,
    heuristic_objective = if (is.null(out$heuristic_objective))
      NA_real_ else out$heuristic_objective,
    gap_pct = if (is.null(out$gap_pct)) NA_real_ else out$gap_pct,
    stringsAsFactors = FALSE
  )

  class(out) <- c("HapBlockR_exact_cross_validation", "list")
  .add_hapblockr_contract(
    result = out,
    method = "validate_crosses_exact",
    call = result_call,
    parameters = list(
      n_cross = n_cross,
      max_cross = max_cross,
      culling_pairwise_k = culling_pairwise_k,
      parent1_col = parent1_col,
      parent2_col = parent2_col,
      criterion_col = criterion_col,
      relatedness_col = relatedness_col,
      max_vars = max_vars
    ),
    sample_ids = parents,
    inputs = contract_inputs,
    transformations = c(
      "missing-criterion exclusion",
      if (!is.null(culling_pairwise_k))
        "pairwise-relatedness culling" else character(),
      "exact binary integer linear optimisation"
    ),
    quality_gates = c(
      solver_optimal = isTRUE(sol$status == 0L),
      exact_cross_count = nrow(exact_plan) == n_cross,
      parent_contribution_cap = is.null(max_cross) ||
        all(parent_load <= max_cross),
      selected_criteria_finite = all(is.finite(plan_score))
    ),
    excluded_records = excluded_records,
    decision_table = decision_table,
    uncertainty = uncertainty
  )
}
