# ==============================================================================
# Operational cross screening and mating-plan certification
# ==============================================================================

.hb_pair_key <- function(parent1, parent2, directional = FALSE) {
  parent1 <- as.character(parent1)
  parent2 <- as.character(parent2)
  if (isTRUE(directional)) {
    paste(parent1, parent2, sep = "\r")
  } else {
    paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = "\r")
  }
}

.hb_optional_status <- function(status, column, default) {
  if (is.null(column) || !column %in% names(status))
    rep(default, nrow(status))
  else
    status[[column]]
}

#' Screen Crosses for Operational Feasibility
#'
#' Screens directed crosses, where parent 1 is the female and parent 2 is the
#' male, against role, fertility, flowering, reciprocal, heterotic-group,
#' forbidden-pair, and quarantine rules. The function reports every failed
#' rule rather than stopping at the first failure.
#'
#' @param cross_pairs Two-column matrix or data frame of directed crosses.
#' @param candidate_status Candidate table with one row per parent.
#' @param id_col Candidate identifier column.
#' @param female_allowed_col,male_allowed_col Optional logical role columns.
#' @param fertility_col Optional logical fertility column.
#' @param flowering_start_col,flowering_end_col Optional numeric or Date
#'   flowering-window columns. Both must be supplied together.
#' @param heterotic_group_col Optional heterotic-group column.
#' @param require_different_heterotic_groups Logical.
#' @param quarantine_group_col Optional quarantine-group column.
#' @param quarantine_compatibility Optional named logical matrix whose rows
#'   are female quarantine groups and columns are male quarantine groups.
#' @param forbidden_pairs Optional two-column table. Pairs are treated as
#'   unordered unless it contains a logical \code{directional} column.
#' @param reciprocal_effects Optional table with \code{female}, \code{male},
#'   and logical \code{allowed} columns.
#' @param no_selfing Logical.
#'
#' @return A data frame with female, male, feasibility, and a semicolon-
#'   separated reason ledger.
#' @export
screen_candidate_crosses <- function(
    cross_pairs,
    candidate_status,
    id_col = "id",
    female_allowed_col = "female_allowed",
    male_allowed_col = "male_allowed",
    fertility_col = "fertile",
    flowering_start_col = "flowering_start",
    flowering_end_col = "flowering_end",
    heterotic_group_col = "heterotic_group",
    require_different_heterotic_groups = FALSE,
    quarantine_group_col = "quarantine_group",
    quarantine_compatibility = NULL,
    forbidden_pairs = NULL,
    reciprocal_effects = NULL,
    no_selfing = TRUE
) {
  pairs <- as.data.frame(cross_pairs, stringsAsFactors = FALSE)
  if (ncol(pairs) < 2L)
    stop("cross_pairs must contain at least two columns.", call. = FALSE)
  pairs <- pairs[, 1:2, drop = FALSE]
  names(pairs) <- c("female", "male")
  pairs$female <- as.character(pairs$female)
  pairs$male <- as.character(pairs$male)
  status <- as.data.frame(candidate_status, stringsAsFactors = FALSE)
  if (!id_col %in% names(status))
    stop("candidate_status is missing id_col.", call. = FALSE)
  status$id <- as.character(status[[id_col]])
  if (anyNA(status$id) || any(!nzchar(status$id)) ||
      anyDuplicated(status$id))
    stop("Candidate IDs must be non-missing, non-empty, and unique.",
         call. = FALSE)
  unknown <- setdiff(unique(c(pairs$female, pairs$male)), status$id)
  if (length(unknown))
    stop("cross_pairs contains unknown candidate ID(s): ",
         paste(head(unknown, 5L), collapse = ", "), call. = FALSE)

  female_status <- status[match(pairs$female, status$id), , drop = FALSE]
  male_status <- status[match(pairs$male, status$id), , drop = FALSE]
  reasons <- vector("list", nrow(pairs))
  add_reason <- function(which_rows, reason) {
    for (row in which(which_rows))
      reasons[[row]] <<- c(reasons[[row]], reason)
  }

  if (isTRUE(no_selfing))
    add_reason(pairs$female == pairs$male, "selfing_forbidden")

  female_allowed <- .hb_optional_status(
    female_status, female_allowed_col, TRUE
  )
  male_allowed <- .hb_optional_status(male_status, male_allowed_col, TRUE)
  if (!is.logical(female_allowed) || anyNA(female_allowed) ||
      !is.logical(male_allowed) || anyNA(male_allowed))
    stop("Role-allowed columns must be complete logical values.",
         call. = FALSE)
  add_reason(!female_allowed, "female_role_forbidden")
  add_reason(!male_allowed, "male_role_forbidden")

  female_fertile <- .hb_optional_status(female_status, fertility_col, TRUE)
  male_fertile <- .hb_optional_status(male_status, fertility_col, TRUE)
  if (!is.logical(female_fertile) || anyNA(female_fertile) ||
      !is.logical(male_fertile) || anyNA(male_fertile))
    stop("The fertility column must contain complete logical values.",
         call. = FALSE)
  add_reason(!female_fertile, "female_infertile")
  add_reason(!male_fertile, "male_infertile")

  has_flowering <- all(
    c(flowering_start_col, flowering_end_col) %in% names(status)
  )
  one_flowering <- any(
    c(flowering_start_col, flowering_end_col) %in% names(status)
  )
  if (one_flowering && !has_flowering)
    stop("Supply both flowering start and end columns, or neither.",
         call. = FALSE)
  flowering_overlap <- rep(NA, nrow(pairs))
  if (has_flowering) {
    fs <- female_status[[flowering_start_col]]
    fe <- female_status[[flowering_end_col]]
    ms <- male_status[[flowering_start_col]]
    me <- male_status[[flowering_end_col]]
    if (anyNA(c(fs, fe, ms, me)) || any(fs > fe) || any(ms > me))
      stop("Flowering windows must be complete and ordered start <= end.",
           call. = FALSE)
    flowering_overlap <- pmax(fs, ms) <= pmin(fe, me)
    add_reason(!flowering_overlap, "flowering_not_synchronised")
  }

  if (isTRUE(require_different_heterotic_groups)) {
    if (!heterotic_group_col %in% names(status))
      stop("A heterotic-group column is required.", call. = FALSE)
    female_group <- as.character(female_status[[heterotic_group_col]])
    male_group <- as.character(male_status[[heterotic_group_col]])
    if (anyNA(c(female_group, male_group)) ||
        any(!nzchar(c(female_group, male_group))))
      stop("Heterotic groups must be complete and non-empty.",
           call. = FALSE)
    add_reason(
      female_group == male_group,
      "same_heterotic_group_forbidden"
    )
  }

  if (!is.null(quarantine_compatibility)) {
    if (!quarantine_group_col %in% names(status))
      stop("A quarantine-group column is required.", call. = FALSE)
    Q <- as.matrix(quarantine_compatibility)
    if (!is.logical(Q) || is.null(rownames(Q)) || is.null(colnames(Q)) ||
        anyNA(Q))
      stop("quarantine_compatibility must be a complete named logical ",
           "matrix.", call. = FALSE)
    fg <- as.character(female_status[[quarantine_group_col]])
    mg <- as.character(male_status[[quarantine_group_col]])
    if (any(!fg %in% rownames(Q)) || any(!mg %in% colnames(Q)))
      stop("A candidate quarantine group is absent from the compatibility ",
           "matrix.", call. = FALSE)
    permitted <- Q[cbind(match(fg, rownames(Q)), match(mg, colnames(Q)))]
    add_reason(!permitted, "quarantine_incompatible")
  }

  if (!is.null(forbidden_pairs)) {
    forbidden <- as.data.frame(forbidden_pairs, stringsAsFactors = FALSE)
    if (ncol(forbidden) < 2L)
      stop("forbidden_pairs must contain at least two columns.",
           call. = FALSE)
    directional <- if ("directional" %in% names(forbidden))
      as.logical(forbidden$directional) else rep(FALSE, nrow(forbidden))
    if (anyNA(directional))
      stop("forbidden_pairs$directional must be complete logical values.",
           call. = FALSE)
    forbidden_keys <- c(
      .hb_pair_key(
        forbidden[[1L]][!directional],
        forbidden[[2L]][!directional],
        FALSE
      ),
      .hb_pair_key(
        forbidden[[1L]][directional],
        forbidden[[2L]][directional],
        TRUE
      )
    )
    pair_unordered <- .hb_pair_key(pairs$female, pairs$male, FALSE)
    pair_directed <- .hb_pair_key(pairs$female, pairs$male, TRUE)
    add_reason(
      pair_unordered %in% forbidden_keys |
        pair_directed %in% forbidden_keys,
      "forbidden_pair"
    )
  }

  if (!is.null(reciprocal_effects)) {
    reciprocal <- as.data.frame(reciprocal_effects,
                                stringsAsFactors = FALSE)
    if (!all(c("female", "male", "allowed") %in% names(reciprocal)) ||
        !is.logical(reciprocal$allowed) || anyNA(reciprocal$allowed))
      stop("reciprocal_effects must contain female, male, and complete ",
           "logical allowed columns.", call. = FALSE)
    key <- .hb_pair_key(reciprocal$female, reciprocal$male, TRUE)
    if (anyDuplicated(key))
      stop("reciprocal_effects contains duplicate directed pairs.",
           call. = FALSE)
    matched <- match(.hb_pair_key(pairs$female, pairs$male, TRUE), key)
    disallowed <- !is.na(matched) & !reciprocal$allowed[matched]
    disallowed[is.na(disallowed)] <- FALSE
    add_reason(disallowed, "reciprocal_cross_forbidden")
  }

  reason_text <- vapply(
    reasons,
    function(x) paste(unique(x), collapse = ";"),
    character(1L)
  )
  data.frame(
    female = pairs$female,
    male = pairs$male,
    feasible = !nzchar(reason_text),
    flowering_overlap = flowering_overlap,
    reasons = reason_text,
    stringsAsFactors = FALSE
  )
}

.hb_violation_row <- function(constraint, subject, observed, limit, detail) {
  data.frame(
    constraint = constraint,
    subject = as.character(subject),
    observed = as.numeric(observed),
    limit = as.numeric(limit),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

#' Certify an Operational Mating Plan
#'
#' Validates a complete plan against cross-level rules, family size, parent
#' capacities, required pairs, period-specific capacity, and subpopulation
#' contribution quotas. It returns a feasibility certificate and all binding
#' constraints.
#'
#' @param plan Data frame containing female, male, family size, and optionally
#'   period.
#' @param candidate_status Candidate status table passed to
#'   \code{\link{screen_candidate_crosses}}. Optional numeric columns
#'   \code{female_capacity}, \code{male_capacity}, \code{total_capacity},
#'   \code{seed_available}, and \code{pollen_available} are enforced.
#' @param female_col,male_col,family_size_col,period_col Plan column names.
#' @param min_family_size Minimum permitted family size.
#' @param required_pairs Optional two-column table of required unordered pairs.
#' @param capacity_by_period Optional table containing \code{id},
#'   \code{period}, and one or more of \code{female_capacity},
#'   \code{male_capacity}, and \code{total_capacity}.
#' @param subpopulation_col Candidate subpopulation column.
#' @param subpopulation_quotas Optional table with \code{subpopulation} and
#'   optional minimum or maximum contribution columns.
#' @param no_repeated_cross Logical.
#' @param ... Additional arguments passed to
#'   \code{\link{screen_candidate_crosses}}.
#'
#' @return A \code{hapblockr_result} with a row-level certificate,
#'   violations, capacity use, quota use, and binding constraints.
#' @export
certify_mating_plan <- function(
    plan,
    candidate_status,
    female_col = "female",
    male_col = "male",
    family_size_col = "family_size",
    period_col = "period",
    min_family_size = 1L,
    required_pairs = NULL,
    capacity_by_period = NULL,
    subpopulation_col = "subpopulation",
    subpopulation_quotas = NULL,
    no_repeated_cross = TRUE,
    ...
) {
  result_call <- match.call()
  mating <- as.data.frame(plan, stringsAsFactors = FALSE)
  required_columns <- c(female_col, male_col, family_size_col)
  if (!all(required_columns %in% names(mating)))
    stop("plan is missing required columns: ",
         paste(setdiff(required_columns, names(mating)), collapse = ", "),
         call. = FALSE)
  if (!nrow(mating))
    stop("plan must contain at least one cross.", call. = FALSE)
  mating$female <- as.character(mating[[female_col]])
  mating$male <- as.character(mating[[male_col]])
  mating$family_size <- as.numeric(mating[[family_size_col]])
  if (any(!is.finite(mating$family_size)) ||
      any(mating$family_size <= 0) ||
      any(mating$family_size != floor(mating$family_size)))
    stop("Family sizes must be positive whole numbers.", call. = FALSE)
  if (length(min_family_size) != 1L || !is.finite(min_family_size) ||
      min_family_size < 1)
    stop("min_family_size must be at least 1.", call. = FALSE)
  mating$period <- if (period_col %in% names(mating))
    as.character(mating[[period_col]]) else "all"
  if (anyNA(mating$period) || any(!nzchar(mating$period)))
    stop("Plan periods must be complete and non-empty.", call. = FALSE)

  screened <- screen_candidate_crosses(
    mating[c("female", "male")],
    candidate_status,
    ...
  )
  mating$cross_rules_passed <- screened$feasible
  mating$reasons <- screened$reasons
  violations <- list()
  add_violation <- function(x) violations[[length(violations) + 1L]] <<- x

  invalid_rows <- which(!screened$feasible)
  if (length(invalid_rows)) {
    for (row in invalid_rows)
      add_violation(.hb_violation_row(
        "cross_rule", row, 0, 1, screened$reasons[row]
      ))
  }
  small <- which(mating$family_size < min_family_size)
  if (length(small)) {
    for (row in small)
      add_violation(.hb_violation_row(
        "minimum_family_size", row, mating$family_size[row],
        min_family_size, "family_size_below_minimum"
      ))
    mating$reasons[small] <- paste0(
      ifelse(nzchar(mating$reasons[small]),
             paste0(mating$reasons[small], ";"), ""),
      "family_size_below_minimum"
    )
  }
  pair_key <- .hb_pair_key(mating$female, mating$male, FALSE)
  if (isTRUE(no_repeated_cross) && anyDuplicated(pair_key)) {
    repeated <- duplicated(pair_key) | duplicated(pair_key, fromLast = TRUE)
    for (row in which(repeated))
      add_violation(.hb_violation_row(
        "no_repeated_cross", row, 2, 1, "repeated_unordered_pair"
      ))
    mating$reasons[repeated] <- paste0(
      ifelse(nzchar(mating$reasons[repeated]),
             paste0(mating$reasons[repeated], ";"), ""),
      "repeated_unordered_pair"
    )
  }

  status <- as.data.frame(candidate_status, stringsAsFactors = FALSE)
  status$id <- as.character(status$id)
  capacity_columns <- intersect(
    c("female_capacity", "male_capacity", "total_capacity",
      "seed_available", "pollen_available"),
    names(status)
  )
  use <- data.frame(
    id = status$id,
    female_use = vapply(status$id, function(id)
      sum(mating$family_size[mating$female == id]), numeric(1L)),
    male_use = vapply(status$id, function(id)
      sum(mating$family_size[mating$male == id]), numeric(1L)),
    stringsAsFactors = FALSE
  )
  use$total_use <- use$female_use + use$male_use
  binding <- list()
  for (column in capacity_columns) {
    limit <- as.numeric(status[[column]])
    if (anyNA(limit) || any(limit < 0))
      stop(column, " must contain complete non-negative values.",
           call. = FALSE)
    observed <- switch(
      column,
      female_capacity = use$female_use,
      seed_available = use$female_use,
      male_capacity = use$male_use,
      pollen_available = use$male_use,
      total_capacity = use$total_use
    )
    exceeded <- observed > limit
    for (row in which(exceeded))
      add_violation(.hb_violation_row(
        column, status$id[row], observed[row], limit[row],
        "parent_capacity_exceeded"
      ))
    exact <- is.finite(limit) & abs(limit - observed) < 1e-8
    if (any(exact)) {
      binding[[length(binding) + 1L]] <- data.frame(
        constraint = column,
        subject = status$id[exact],
        observed = observed[exact],
        limit = limit[exact],
        slack = 0,
        stringsAsFactors = FALSE
      )
    }
    use[[column]] <- limit
    use[[paste0(column, "_slack")]] <- limit - observed
  }

  period_use <- data.frame()
  if (!is.null(capacity_by_period)) {
    period_capacity <- as.data.frame(capacity_by_period,
                                     stringsAsFactors = FALSE)
    if (!all(c("id", "period") %in% names(period_capacity)))
      stop("capacity_by_period must contain id and period.", call. = FALSE)
    period_capacity$id <- as.character(period_capacity$id)
    period_capacity$period <- as.character(period_capacity$period)
    if (anyDuplicated(paste(period_capacity$id, period_capacity$period)))
      stop("capacity_by_period contains duplicate id-period rows.",
           call. = FALSE)
    period_columns <- intersect(
      c("female_capacity", "male_capacity", "total_capacity"),
      names(period_capacity)
    )
    if (!length(period_columns))
      stop("capacity_by_period requires at least one capacity column.",
           call. = FALSE)
    period_use <- period_capacity
    period_use$female_use <- mapply(
      function(id, period) sum(mating$family_size[
        mating$female == id & mating$period == period
      ]),
      period_use$id, period_use$period
    )
    period_use$male_use <- mapply(
      function(id, period) sum(mating$family_size[
        mating$male == id & mating$period == period
      ]),
      period_use$id, period_use$period
    )
    period_use$total_use <- period_use$female_use + period_use$male_use
    for (column in period_columns) {
      limit <- as.numeric(period_use[[column]])
      if (anyNA(limit) || any(limit < 0))
        stop("Period capacities must be complete and non-negative.",
             call. = FALSE)
      observed <- switch(
        column,
        female_capacity = period_use$female_use,
        male_capacity = period_use$male_use,
        total_capacity = period_use$total_use
      )
      exceeded <- observed > limit
      for (row in which(exceeded))
        add_violation(.hb_violation_row(
          paste0("period_", column),
          paste(period_use$id[row], period_use$period[row], sep = "@"),
          observed[row], limit[row], "period_capacity_exceeded"
        ))
      period_use[[paste0(column, "_slack")]] <- limit - observed
      exact <- abs(limit - observed) < 1e-8
      if (any(exact)) {
        binding[[length(binding) + 1L]] <- data.frame(
          constraint = paste0("period_", column),
          subject = paste(
            period_use$id[exact], period_use$period[exact], sep = "@"
          ),
          observed = observed[exact],
          limit = limit[exact],
          slack = 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!is.null(required_pairs)) {
    required <- as.data.frame(required_pairs, stringsAsFactors = FALSE)
    if (ncol(required) < 2L)
      stop("required_pairs must contain at least two columns.",
           call. = FALSE)
    required_key <- unique(.hb_pair_key(required[[1L]], required[[2L]], FALSE))
    missing_required <- setdiff(required_key, pair_key)
    if (length(missing_required)) {
      for (key in missing_required)
        add_violation(.hb_violation_row(
          "required_pair", gsub("\r", " x ", key, fixed = TRUE),
          0, 1, "required_pair_absent"
        ))
    }
  }

  quota_use <- data.frame()
  if (!is.null(subpopulation_quotas)) {
    if (!subpopulation_col %in% names(status))
      stop("candidate_status is missing subpopulation_col.", call. = FALSE)
    quotas <- as.data.frame(subpopulation_quotas,
                            stringsAsFactors = FALSE)
    if (!"subpopulation" %in% names(quotas))
      stop("subpopulation_quotas must contain subpopulation.",
           call. = FALSE)
    contribution <- setNames(
      vapply(status$id, function(id)
        sum(mating$family_size[
          mating$female == id | mating$male == id
        ]), numeric(1L)),
      status$id
    )
    status_subpopulation <- as.character(status[[subpopulation_col]])
    quota_use <- data.frame(
      subpopulation = as.character(quotas$subpopulation),
      contribution = vapply(
        as.character(quotas$subpopulation),
        function(group) sum(contribution[
          status$id[status_subpopulation == group]
        ]),
        numeric(1L)
      ),
      stringsAsFactors = FALSE
    )
    for (bound in c("min_contribution", "max_contribution")) {
      if (!bound %in% names(quotas)) next
      limit <- as.numeric(quotas[[bound]])
      if (anyNA(limit) || any(limit < 0))
        stop(bound, " must be complete and non-negative.", call. = FALSE)
      failed <- if (bound == "min_contribution")
        quota_use$contribution < limit else quota_use$contribution > limit
      for (row in which(failed))
        add_violation(.hb_violation_row(
          bound, quota_use$subpopulation[row],
          quota_use$contribution[row], limit[row],
          "subpopulation_quota_failed"
        ))
      slack <- if (bound == "min_contribution")
        quota_use$contribution - limit else limit - quota_use$contribution
      quota_use[[bound]] <- limit
      quota_use[[paste0(bound, "_slack")]] <- slack
      exact <- abs(slack) < 1e-8
      if (any(exact)) {
        binding[[length(binding) + 1L]] <- data.frame(
          constraint = bound,
          subject = quota_use$subpopulation[exact],
          observed = quota_use$contribution[exact],
          limit = limit[exact],
          slack = 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  violation_table <- if (length(violations))
    do.call(rbind, violations) else data.frame(
      constraint = character(), subject = character(),
      observed = numeric(), limit = numeric(), detail = character(),
      stringsAsFactors = FALSE
    )
  binding_table <- if (length(binding))
    unique(do.call(rbind, binding)) else data.frame(
      constraint = character(), subject = character(),
      observed = numeric(), limit = numeric(), slack = numeric(),
      stringsAsFactors = FALSE
    )
  mating$feasible <- !nzchar(mating$reasons)
  feasible <- !nrow(violation_table)
  certificate <- data.frame(
    feasible = feasible,
    n_crosses = nrow(mating),
    total_families = sum(mating$family_size),
    n_violations = nrow(violation_table),
    n_binding_constraints = nrow(binding_table),
    stringsAsFactors = FALSE
  )
  result <- list(
    certificate = certificate,
    plan = mating,
    violations = violation_table,
    binding_constraints = binding_table,
    capacity_use = use,
    period_capacity_use = period_use,
    subpopulation_quota_use = quota_use
  )
  .add_hapblockr_contract(
    result = result,
    method = "certify_mating_plan",
    call = result_call,
    parameters = list(
      min_family_size = min_family_size,
      no_repeated_cross = isTRUE(no_repeated_cross)
    ),
    sample_ids = unique(c(mating$female, mating$male)),
    inputs = list(
      plan = mating[c("female", "male", "family_size", "period")],
      candidate_status = candidate_status,
      required_pairs = required_pairs,
      capacity_by_period = capacity_by_period,
      subpopulation_quotas = subpopulation_quotas
    ),
    transformations = "operational feasibility certification",
    quality_gates = c(
      cross_rules = all(screened$feasible),
      minimum_family_size = !length(small),
      no_repeated_cross =
        !isTRUE(no_repeated_cross) || !anyDuplicated(pair_key),
      no_capacity_violations =
        !any(grepl("capacity", violation_table$constraint)),
      required_pairs_present =
        !any(violation_table$constraint == "required_pair"),
      subpopulation_quotas =
        !any(grepl("contribution", violation_table$constraint))
    ),
    warnings = if (feasible) character() else
      paste(nrow(violation_table), "operational violation(s)"),
    excluded_records = violation_table,
    decision_table = certificate,
    uncertainty = binding_table
  )
}
