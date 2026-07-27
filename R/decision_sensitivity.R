# ==============================================================================
# Threshold, population, and environment decision sensitivity
# ==============================================================================

.hb_decision_entities <- function(decision) {
  if ("id" %in% names(decision))
    return(as.character(decision$id))
  pair1 <- intersect(c("female", "parent1"), names(decision))
  pair2 <- intersect(c("male", "parent2"), names(decision))
  if (length(pair1) && length(pair2))
    return(paste(
      decision[[pair1[1L]]], decision[[pair2[1L]]], sep = " x "
    ))
  paste0("decision_", seq_len(nrow(decision)))
}

#' Assess Recommendation Stability Across Scenarios
#'
#' Compares validated result objects from threshold sweeps, leave-one-
#' environment-out analyses, leave-one-population-out analyses, or alternative
#' model assumptions. Scenarios are supplied explicitly so that the function
#' never refits or silently changes a model.
#'
#' @param results Named list of at least two \code{hapblockr_result} objects.
#' @param top_n Number of top-ranked rows treated as selected when no explicit
#'   selection column is present. Default \code{10}.
#' @param score_col Optional common numeric score column.
#' @param selected_col Optional common logical selection column.
#' @param minimum_selection_frequency Minimum scenario frequency required for
#'   a stable recommendation.
#' @param minimum_baseline_jaccard Minimum Jaccard overlap required between
#'   every scenario and the first named baseline scenario.
#' @param require_valid Logical. Refuse failed input result contracts.
#'
#' @return A \code{hapblockr_result} with item frequencies, pairwise Jaccard
#'   overlap, rank and score summaries, and scenario validation status.
#' @export
assess_decision_stability <- function(
    results,
    top_n = 10L,
    score_col = NULL,
    selected_col = NULL,
    minimum_selection_frequency = 0.50,
    minimum_baseline_jaccard = 0.50,
    require_valid = TRUE
) {
  result_call <- match.call()
  if (!is.list(results) || length(results) < 2L ||
      is.null(names(results)) || any(!nzchar(names(results))) ||
      anyDuplicated(names(results)))
    stop("results must be a uniquely named list of at least two result ",
         "objects.", call. = FALSE)
  valid_class <- vapply(
    results, inherits, logical(1L), what = "hapblockr_result"
  )
  if (!all(valid_class))
    stop("Every scenario must inherit from 'hapblockr_result'.",
         call. = FALSE)
  scenario_valid <- vapply(
    results,
    function(x) identical(
      x$result_contract$validation_status, "passed"
    ),
    logical(1L)
  )
  if (isTRUE(require_valid) && !all(scenario_valid))
    stop("Failed result contracts cannot enter a required-valid stability ",
         "analysis.", call. = FALSE)
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L)
    stop("top_n must be a positive integer.", call. = FALSE)
  for (value in c(minimum_selection_frequency,
                  minimum_baseline_jaccard)) {
    if (length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1)
      stop("Stability thresholds must be in [0, 1].", call. = FALSE)
  }

  scenario_rows <- list()
  selected_sets <- list()
  for (scenario in names(results)) {
    decision <- as.data.frame(
      results[[scenario]]$result_contract$decision_table,
      stringsAsFactors = FALSE
    )
    if (!nrow(decision))
      stop("Scenario '", scenario, "' has no decision rows.",
           call. = FALSE)
    entity <- .hb_decision_entities(decision)
    if (anyNA(entity) || any(!nzchar(entity)) || anyDuplicated(entity))
      stop("Scenario '", scenario, "' has missing or duplicate decision ",
           "entities.", call. = FALSE)
    selected_name <- if (!is.null(selected_col)) {
      selected_col
    } else {
      found <- intersect(
        c("recommendable", "recommendation_eligible", "selected", "feasible"),
        names(decision)
      )
      if (length(found)) found[1L] else NULL
    }
    score_name <- if (!is.null(score_col)) {
      score_col
    } else {
      found <- intersect(
        c("selection_index", "UC", "score", "prediction", "objective"),
        names(decision)
      )
      if (length(found)) found[1L] else {
        numeric_columns <- names(decision)[
          vapply(decision, is.numeric, logical(1L))
        ]
        setdiff(numeric_columns, "rank")[1L]
      }
    }
    if (!is.null(selected_name) && !selected_name %in% names(decision))
      stop("selected_col was not found in scenario '", scenario, "'.",
           call. = FALSE)
    if (is.null(score_name) || is.na(score_name) ||
        !score_name %in% names(decision) ||
        !is.numeric(decision[[score_name]]))
      stop("No common numeric decision score is available in scenario '",
           scenario, "'.", call. = FALSE)
    score <- decision[[score_name]]
    rank_value <- if ("rank" %in% names(decision) &&
                      is.numeric(decision$rank)) {
      decision$rank
    } else {
      rank(-score, ties.method = "min", na.last = "keep")
    }
    selected <- if (is.null(selected_name)) {
      rank_value <= min(as.integer(top_n), nrow(decision))
    } else {
      as.logical(decision[[selected_name]])
    }
    selected[is.na(selected)] <- FALSE
    selected_sets[[scenario]] <- entity[selected]
    scenario_rows[[scenario]] <- data.frame(
      scenario = scenario,
      entity = entity,
      score = score,
      rank = rank_value,
      selected = selected,
      scenario_valid = unname(scenario_valid[scenario]),
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, scenario_rows)
  rownames(long) <- NULL
  entities <- sort(unique(long$entity))
  stability <- do.call(rbind, lapply(entities, function(entity) {
    x <- long[long$entity == entity, , drop = FALSE]
    data.frame(
      entity = entity,
      scenarios_observed = nrow(x),
      selection_frequency = mean(x$selected),
      mean_score = mean(x$score, na.rm = TRUE),
      score_sd = if (sum(is.finite(x$score)) > 1L)
        stats::sd(x$score, na.rm = TRUE) else NA_real_,
      median_rank = stats::median(x$rank, na.rm = TRUE),
      worst_rank = max(x$rank, na.rm = TRUE),
      stable_recommendation =
        mean(x$selected) >= minimum_selection_frequency,
      stringsAsFactors = FALSE
    )
  }))
  rownames(stability) <- NULL

  scenario_pairs <- utils::combn(names(results), 2L)
  overlap <- do.call(rbind, lapply(seq_len(ncol(scenario_pairs)), function(i) {
    first <- scenario_pairs[1L, i]
    second <- scenario_pairs[2L, i]
    union_set <- union(selected_sets[[first]], selected_sets[[second]])
    intersection_set <- intersect(
      selected_sets[[first]], selected_sets[[second]]
    )
    data.frame(
      scenario1 = first,
      scenario2 = second,
      n_selected1 = length(selected_sets[[first]]),
      n_selected2 = length(selected_sets[[second]]),
      n_intersection = length(intersection_set),
      n_union = length(union_set),
      jaccard = if (length(union_set))
        length(intersection_set) / length(union_set) else 1,
      stringsAsFactors = FALSE
    )
  }))
  baseline <- names(results)[1L]
  baseline_overlap <- overlap$jaccard[
    overlap$scenario1 == baseline | overlap$scenario2 == baseline
  ]
  result <- list(
    stability = stability,
    pairwise_overlap = overlap,
    scenario_decisions = long,
    baseline = baseline
  )
  .add_hapblockr_contract(
    result = result,
    method = "assess_decision_stability",
    call = result_call,
    parameters = list(
      top_n = as.integer(top_n),
      score_col = score_col,
      selected_col = selected_col,
      minimum_selection_frequency = minimum_selection_frequency,
      minimum_baseline_jaccard = minimum_baseline_jaccard,
      require_valid = isTRUE(require_valid)
    ),
    sample_ids = entities,
    inputs = lapply(results, function(x) x$result_contract$input_hashes),
    transformations = c(
      "scenario decision extraction",
      "selection-frequency analysis",
      "pairwise Jaccard overlap"
    ),
    quality_gates = c(
      all_scenarios_valid = all(scenario_valid),
      baseline_overlap =
        all(baseline_overlap >= minimum_baseline_jaccard),
      stable_recommendation_available =
        any(stability$stable_recommendation)
    ),
    warnings = if (all(scenario_valid)) character() else
      paste(sum(!scenario_valid), "failed scenario(s) included"),
    decision_table = stability,
    uncertainty = overlap
  )
}
