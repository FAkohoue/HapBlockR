# ==============================================================================
# Multivariate and environment-aware breeding models
# ==============================================================================

.hb_validate_covariance <- function(x, traits, name, positive_definite = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) != length(traits) ||
      ncol(x) != length(traits))
    stop(name, " must be a numeric ", length(traits), " x ",
         length(traits), " matrix.", call. = FALSE)
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      !setequal(rownames(x), traits) || !setequal(colnames(x), traits))
    stop(name, " must have row and column names matching the traits.",
         call. = FALSE)
  x <- x[traits, traits, drop = FALSE]
  if (any(!is.finite(x)) || max(abs(x - t(x))) > 1e-8)
    stop(name, " must be finite and symmetric.", call. = FALSE)
  eigenvalues <- eigen((x + t(x)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values
  tolerance <- 1e-8 * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -tolerance ||
      (positive_definite && min(eigenvalues) <= tolerance))
    stop(name, " must be ", if (positive_definite) "positive definite." else
           "positive semidefinite.", call. = FALSE)
  (x + t(x)) / 2
}

.hb_validate_kernel <- function(K, ids, name = "K") {
  K <- as.matrix(K)
  if (!is.numeric(K) || nrow(K) != ncol(K) ||
      is.null(rownames(K)) || is.null(colnames(K)))
    stop(name, " must be a square numeric matrix with dimnames.",
         call. = FALSE)
  if (!identical(rownames(K), colnames(K)))
    stop(name, " row and column names must be identical and ordered.",
         call. = FALSE)
  missing_ids <- setdiff(ids, rownames(K))
  if (length(missing_ids))
    stop(name, " is missing ", length(missing_ids), " required ID(s): ",
         paste(head(missing_ids, 5L), collapse = ", "), call. = FALSE)
  K <- K[ids, ids, drop = FALSE]
  if (any(!is.finite(K)) || max(abs(K - t(K))) > 1e-8)
    stop(name, " must be finite and symmetric.", call. = FALSE)
  eigenvalues <- eigen((K + t(K)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values
  tolerance <- 1e-8 * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -tolerance)
    stop(name, " must be positive semidefinite.", call. = FALSE)
  if (any(diag(K) <= 0))
    stop(name, " must have strictly positive diagonal entries.",
         call. = FALSE)
  (K + t(K)) / 2
}

.hb_validate_objective_vector <- function(x, traits, name) {
  if (!is.numeric(x) || is.complex(x) || is.null(names(x)) ||
      length(x) != length(traits) || anyDuplicated(names(x)) ||
      !setequal(names(x), traits) || any(!is.finite(x[traits]))) {
    stop(name, " must be a finite real numeric vector named by trait.",
         call. = FALSE)
  }
  x <- as.numeric(x[traits])
  names(x) <- traits
  if (any(x < 0)) {
    stop(
      name,
      " must contain non-negative favourable-direction magnitudes. ",
      "Use directions to declare traits that should decrease; signed ",
      "objectives are not accepted.",
      call. = FALSE
    )
  }
  x
}

.hb_validate_engine_control <- function(x, name) {
  if (!is.list(x)) {
    stop(name, " must be a named list.", call. = FALSE)
  }
  if (!length(x)) return(x)
  if (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x))) {
    stop(name, " must have unique, non-empty argument names.", call. = FALSE)
  }
  x
}

.hb_engine_vector <- function(x, traits, name) {
  if (!is.numeric(x) || is.complex(x) || is.null(names(x)) ||
      length(x) != length(traits) || anyDuplicated(names(x)) ||
      !setequal(names(x), traits) || any(!is.finite(x[traits]))) {
    stop(
      "DesiredGainR returned an invalid ", name,
      "; expected one finite named value per trait.",
      call. = FALSE
    )
  }
  out <- as.numeric(x[traits])
  names(out) <- traits
  out
}

.hb_check_engine_arguments <- function(control, engine, name, protected) {
  overridden <- intersect(names(control), protected)
  if (length(overridden)) {
    stop(name, " cannot override: ", paste(overridden, collapse = ", "),
         call. = FALSE)
  }
  # do.call() otherwise accepts partial argument names, including typos.
  unknown <- setdiff(names(control), names(formals(engine)))
  if (length(unknown)) {
    stop(name, " contains unknown argument(s): ",
         paste(unknown, collapse = ", "),
         ". Use exact DesiredGainR argument names.", call. = FALSE)
  }
  absent <- setdiff(protected, names(formals(engine)))
  if (length(absent)) {
    stop("The installed DesiredGainR API is incompatible; missing argument(s): ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.hb_engine_transformation <- function(engine, traits, direction_sign) {
  transform <- lapply(c("centre", "scale", "direction"), function(field) {
    .hb_engine_vector(engine$transformation[[field]], traits,
                      paste("transformation", field))
  })
  names(transform) <- c("centre", "scale", "direction")
  if (any(transform$scale <= 0)) {
    stop("DesiredGainR returned a non-positive transformation scale.",
         call. = FALSE)
  }
  if (!identical(unname(transform$direction), as.numeric(direction_sign))) {
    stop("DesiredGainR returned trait directions inconsistent with the ",
         "HapBlockR request.", call. = FALSE)
  }
  transform
}

.hb_engine_ranking <- function(engine, ids, score_col, n_select, exact) {
  decision <- as.data.frame(engine$ranked_geno)
  if (!all(c("id", score_col, "Selected") %in% names(decision))) {
    stop("DesiredGainR returned an incomplete ranking table.", call. = FALSE)
  }
  if (nrow(decision) != length(ids) || anyDuplicated(decision$id) ||
      !setequal(as.character(decision$id), ids)) {
    stop("DesiredGainR returned candidate IDs inconsistent with the input.",
         call. = FALSE)
  }
  score <- decision[[score_col]]
  selected <- decision$Selected
  if (!is.numeric(score) || is.complex(score) || any(!is.finite(score)) ||
      !is.logical(selected) || anyNA(selected) || sum(selected) < 1L ||
      sum(selected) > n_select || (exact && sum(selected) != n_select)) {
    stop("DesiredGainR returned invalid scores or selection decisions.",
         call. = FALSE)
  }
  names(decision)[names(decision) == score_col] <- "selection_index"
  decision$merit_score <- decision$selection_index
  decision
}

.hb_engine_positive_scalar <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(
      "DesiredGainR returned an invalid ", name, ".",
      call. = FALSE
    )
  }
  lower_ok <- if (allow_zero) x >= 0 else x > 0
  if (!lower_ok) {
    stop(
      "DesiredGainR returned an invalid ", name, ".",
      call. = FALSE
    )
  }
  as.numeric(x)
}

.hb_validate_symmetric_matrix <- function(x, traits, name) {
  x <- as.matrix(x)
  if (!is.numeric(x) || is.complex(x) ||
      nrow(x) != length(traits) || ncol(x) != length(traits) ||
      is.null(rownames(x)) || is.null(colnames(x)) ||
      !setequal(rownames(x), traits) || !setequal(colnames(x), traits)) {
    stop(name, " must be a named real numeric square matrix matching traits.",
         call. = FALSE)
  }
  x <- x[traits, traits, drop = FALSE]
  if (any(!is.finite(x)) || max(abs(x - t(x))) > 1e-8) {
    stop(name, " must be finite and symmetric.", call. = FALSE)
  }
  (x + t(x)) / 2
}

.hb_cov_to_parameters <- function(S) {
  L <- t(chol(S))
  out <- numeric(nrow(S) * (nrow(S) + 1L) / 2L)
  cursor <- 1L
  for (i in seq_len(nrow(S))) {
    for (j in seq_len(i)) {
      out[cursor] <- if (i == j) log(L[i, j]) else L[i, j]
      cursor <- cursor + 1L
    }
  }
  out
}

.hb_parameters_to_cov <- function(parameters, dimension) {
  L <- matrix(0, dimension, dimension)
  cursor <- 1L
  for (i in seq_len(dimension)) {
    for (j in seq_len(i)) {
      L[i, j] <- if (i == j) exp(parameters[cursor]) else parameters[cursor]
      cursor <- cursor + 1L
    }
  }
  tcrossprod(L)
}

.hb_gls_components <- function(V, X, y) {
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  V_inv <- chol2inv(R)
  XtVinvX <- crossprod(X, V_inv %*% X)
  RX <- tryCatch(chol(XtVinvX), error = function(e) NULL)
  if (is.null(RX)) return(NULL)
  XtVinvX_inv <- chol2inv(RX)
  beta <- XtVinvX_inv %*% crossprod(X, V_inv %*% y)
  residual <- y - X %*% beta
  list(
    R = R,
    V_inv = V_inv,
    XtVinvX_inv = XtVinvX_inv,
    beta = as.numeric(beta),
    residual = as.numeric(residual),
    logdet_V = 2 * sum(log(diag(R))),
    logdet_X = 2 * sum(log(diag(RX)))
  )
}

#' Construct an Economic or Desired-Gain Selection Index
#'
#' Constructs a Smith-Hazel economic index or a Pesek-Baker desired-gain
#' index, or delegates an optimised Desired-Gain Selection Index (DGSI) or
#' Quadratic Genomic Selection Index (QGSI) to \pkg{DesiredGainR}.
#' Trait directions are required; declared units are optional metadata.
#'
#' @param trait_values Numeric matrix or data frame with candidates in rows
#'   and traits in columns. Row names are required. DGSI uses internally
#'   modelled trait values; QGSI uses genomic estimated breeding values.
#'   The trait name \code{"id"} is reserved by the delegated interface.
#' @param genetic_cov Named genetic covariance matrix in original trait units.
#'   Used in the classical and DGSI coefficient calculations. For QGSI, this
#'   matrix is retained as upstream context; it is not substituted for Gamma.
#' @param phenotypic_cov Named positive-definite phenotypic or index-variable
#'   covariance matrix in original trait units. Used by classical methods and
#'   DGSI; retained as upstream context for QGSI.
#' @param economic_weights Named numeric vector of economic weights. Supply
#'   exactly one of \code{economic_weights} and \code{desired_gains}. Values
#'   must be non-negative importance magnitudes after trait orientation;
#'   \code{directions}, rather than a negative weight, declares decrease.
#' @param desired_gains Named numeric vector of non-negative desired-gain
#'   magnitudes in the favourable-direction trait space. For Pesek-Baker these
#'   are in original trait units. For DGSI these are in candidate standard
#'   deviations, regardless of \code{dgsi_control$scale_traits}; divide an
#'   original-unit target by its candidate standard deviation before passing it.
#' @param directions Named character vector containing \code{"increase"} or
#'   \code{"decrease"} for every trait.
#' @param units Optional named character vector of trait units.
#' @param method Index method. The classical methods are `"smith_hazel"` and
#'   `"pesek_baker"`. `"dgsi"` and `"qgsi"` are delegated to
#'   \pkg{DesiredGainR}. If `NULL`, the method is inferred from the supplied
#'   classical objective for backward compatibility.
#' @param n_select Number of candidates selected by DGSI or QGSI. DGSI may
#'   retain fewer when fewer candidates satisfy its eligibility thresholds.
#'   When omitted
#'   for either delegated method, the default is 10 percent of candidates,
#'   with a minimum of one. It does not apply to classical indices.
#' @param quadratic_weights Symmetric quadratic and cross-product weight matrix
#'   required by QGSI. It must describe the favourable-direction trait space
#'   obtained after applying \code{directions} and any requested scaling.
#'   QGSI linear economic weights must also refer to that analysis scale.
#' @param dgsi_control Named list of additional arguments passed to
#'   \code{DesiredGainR::run_dgsi()}. Names must match exactly; partial names
#'   and structural overrides are rejected. Reference and validation data
#'   contain the same original-unit trait columns as \code{trait_values}.
#' @param qgsi_control Named list of additional arguments passed to
#'   \code{DesiredGainR::run_qgsi()}. Structural arguments supplied by
#'   HapBlockR, including the candidate data, weights, directions and
#'   selection count, cannot be overridden. Names must match exactly. Supply
#'   \code{Gamma} in original trait units, or \code{relationship_matrix} with
#'   row and column names matching the reference IDs. Reference data use an
#'   \code{id} column and the original-unit trait columns. Without an explicit
#'   reference, candidate values form the reference population.
#' @param selection_intensity Optional positive selection intensity used to
#'   report the
#'   expected response vector for the classical Smith-Hazel and Pesek-Baker
#'   methods; the default for those methods is one. DGSI and QGSI derive normal
#'   selection intensity from the proportion actually selected under their
#'   requested rule, so this argument does not apply to delegated methods.
#'
#' @return A \code{hapblockr_result} containing candidate scores, objective
#'   information, covariance diagnostics, and method-appropriate response
#'   summaries. Classical and DGSI methods return linear coefficients. QGSI
#'   instead returns linear weights, the quadratic-weight matrix, and
#'   candidate-specific contributions because it has no single global
#'   coefficient vector. The unmodified delegated fit is in \code{engine_result}
#'   and can be passed to DesiredGainR comparison tools. DGSI's
#'   \code{coefficients$coefficient} remains on the engine's favourable-direction
#'   analysis scale. Use the named \code{coefficients_original_units} vector
#'   to combine original-unit marker, haplotype or block effects. Candidate
#'   scores equal \code{trait_values \%*\% coefficients_original_units +
#'   score_intercept}; the intercept accounts for reference centring.
#' @export
build_selection_index <- function(
    trait_values,
    genetic_cov,
    phenotypic_cov,
    economic_weights = NULL,
    desired_gains = NULL,
    directions,
    units = NULL,
    method = NULL,
    n_select = NULL,
    quadratic_weights = NULL,
    dgsi_control = list(),
    qgsi_control = list(),
    selection_intensity = NULL
) {
  result_call <- match.call()
  values <- as.matrix(trait_values)
  if (!is.numeric(values) || !nrow(values) || !ncol(values) ||
      is.null(rownames(values)) || is.null(colnames(values)))
    stop("trait_values must be a non-empty numeric matrix with row and ",
         "column names.", call. = FALSE)
  if (any(!is.finite(values)) || anyDuplicated(rownames(values)) ||
      anyDuplicated(colnames(values)) || anyNA(dimnames(values)[[1L]]) ||
      anyNA(dimnames(values)[[2L]]) || any(!nzchar(rownames(values))) ||
      any(!nzchar(colnames(values))))
    stop("trait_values must be finite with unique candidate and trait IDs.",
         call. = FALSE)
  traits <- colnames(values)
  if (is.null(method)) {
    method <- if (!is.null(economic_weights)) "smith_hazel" else
      "pesek_baker"
  }
  method <- match.arg(
    method,
    c("smith_hazel", "pesek_baker", "dgsi", "qgsi")
  )
  dgsi_control <- .hb_validate_engine_control(dgsi_control, "dgsi_control")
  qgsi_control <- .hb_validate_engine_control(qgsi_control, "qgsi_control")
  if (method != "dgsi" && length(dgsi_control)) {
    stop("dgsi_control applies only to method = 'dgsi'.", call. = FALSE)
  }
  if (method != "qgsi" && length(qgsi_control)) {
    stop("qgsi_control applies only to method = 'qgsi'.", call. = FALSE)
  }
  if (xor(is.null(economic_weights), is.null(desired_gains)) == FALSE)
    stop("Supply exactly one of economic_weights and desired_gains.",
         call. = FALSE)
  if (is.null(names(directions)) || !setequal(names(directions), traits) ||
      any(!directions[traits] %in% c("increase", "decrease")))
    stop("directions must be named by trait and contain only 'increase' or ",
         "'decrease'.", call. = FALSE)
  if (is.null(units)) {
    units <- stats::setNames(rep(NA_character_, length(traits)), traits)
  } else if (is.null(names(units)) || !setequal(names(units), traits) ||
             anyNA(units[traits]) || any(!nzchar(units[traits]))) {
    stop("When supplied, units must contain one named, non-empty value for ",
         "every trait.", call. = FALSE)
  }
  delegated <- method %in% c("dgsi", "qgsi")
  if (delegated && !is.null(selection_intensity)) {
    stop(
      "selection_intensity applies only to smith_hazel and pesek_baker. ",
      "DGSI and QGSI derive intensity from n_select.",
      call. = FALSE
    )
  }
  if (!delegated && !is.null(n_select)) {
    stop(
      "n_select applies only to dgsi and qgsi.",
      call. = FALSE
    )
  }
  if (delegated && is.null(n_select)) {
    n_select <- max(1L, floor(0.10 * nrow(values)))
  }
  if (delegated &&
      (!is.numeric(n_select) || length(n_select) != 1L ||
       !is.finite(n_select) || n_select < 1 ||
       n_select != as.integer(n_select) || n_select > nrow(values))) {
    stop(
      "n_select must be one positive integer no larger than the number ",
      "of candidates.",
      call. = FALSE
    )
  }
  if (!delegated && is.null(selection_intensity)) {
    selection_intensity <- 1
  }
  if (!delegated &&
      (!is.numeric(selection_intensity) ||
       length(selection_intensity) != 1L ||
       !is.finite(selection_intensity) || selection_intensity <= 0)) {
    stop(
      "selection_intensity must be one positive finite real number.",
      call. = FALSE
    )
  }
  if (delegated) n_select <- as.integer(n_select)
  if (delegated && "id" %in% traits) {
    stop("The trait name 'id' is reserved by the DesiredGainR interface; ",
         "rename this trait and its covariance/objective entries.",
         call. = FALSE)
  }

  G <- .hb_validate_covariance(genetic_cov, traits, "genetic_cov")
  P <- .hb_validate_covariance(
    phenotypic_cov, traits, "phenotypic_cov", positive_definite = TRUE
  )
  direction_sign <- ifelse(directions[traits] == "increase", 1, -1)
  D <- diag(direction_sign, nrow = length(traits))
  dimnames(D) <- list(traits, traits)
  values_oriented <- sweep(values, 2L, direction_sign, "*")
  G_oriented <- D %*% G %*% D
  P_oriented <- D %*% P %*% D

  if (method %in% c("dgsi", "qgsi")) {
    return(.hb_build_desiredgain_index(
      method = method,
      values = values,
      G = G,
      P = P,
      economic_weights = economic_weights,
      desired_gains = desired_gains,
      directions = directions,
      units = units,
      n_select = n_select,
      quadratic_weights = quadratic_weights,
      dgsi_control = dgsi_control,
      qgsi_control = qgsi_control,
      result_call = result_call
    ))
  }

  if (method == "smith_hazel") {
    if (is.null(economic_weights)) {
      stop("method = 'smith_hazel' requires economic_weights.",
           call. = FALSE)
    }
    objective <- .hb_validate_objective_vector(
      economic_weights, traits, "economic_weights"
    )
    coefficients <- solve(P_oriented, G_oriented %*% objective)
    index_type <- "smith_hazel"
    aggregate_variance <- as.numeric(
      crossprod(objective, G_oriented %*% objective)
    )
    index_objective_cov <- as.numeric(
      crossprod(coefficients, G_oriented %*% objective)
    )
  } else {
    if (is.null(desired_gains)) {
      stop("method = 'pesek_baker' requires desired_gains.",
           call. = FALSE)
    }
    objective <- .hb_validate_objective_vector(
      desired_gains, traits, "desired_gains"
    )
    middle <- G_oriented %*% solve(P_oriented, G_oriented)
    middle <- .hb_validate_covariance(
      middle, traits, "G P^-1 G", positive_definite = TRUE
    )
    coefficients <- solve(
      P_oriented,
      G_oriented %*% solve(middle, objective)
    )
    index_type <- "pesek_baker"
    aggregate_variance <- NA_real_
    index_objective_cov <- NA_real_
  }
  coefficients <- as.numeric(coefficients)
  names(coefficients) <- traits
  index_variance <- as.numeric(
    crossprod(coefficients, P_oriented %*% coefficients)
  )
  if (!is.finite(index_variance) || index_variance <= 0)
    stop("The resulting index has non-positive variance.", call. = FALSE)
  expected_response_oriented <- as.numeric(
    selection_intensity * G_oriented %*% coefficients /
      sqrt(index_variance)
  )
  expected_response <- expected_response_oriented * direction_sign
  accuracy <- if (index_type == "smith_hazel" &&
                  aggregate_variance > 0) {
    index_objective_cov / sqrt(index_variance * aggregate_variance)
  } else {
    NA_real_
  }
  scores <- as.numeric(values_oriented %*% coefficients)
  decision <- data.frame(
    id = rownames(values),
    selection_index = scores,
    rank = rank(-scores, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  decision <- decision[order(decision$rank, decision$id), , drop = FALSE]
  coefficient_table <- data.frame(
    trait = traits,
    unit = units[traits],
    direction = directions[traits],
    coefficient = coefficients,
    linear_weight = NA_real_,
    expected_response = expected_response,
    realised_response_sd_favourable = NA_real_,
    realised_response_sd_original = NA_real_,
    stringsAsFactors = FALSE
  )
  response_summary <- data.frame(
    trait = traits,
    direction = directions[traits],
    response = expected_response,
    response_kind = "expected_genetic_response",
    response_scale = "original trait units",
    stringsAsFactors = FALSE
  )
  result <- list(
    index_type = index_type,
    coefficients = coefficient_table,
    scores = decision,
    genetic_cov = G,
    phenotypic_cov = P,
    genetic_cor = stats::cov2cor(G),
    index_variance = index_variance,
    index_accuracy = accuracy,
    response_summary = response_summary
  )
  .add_hapblockr_contract(
    result = result,
    method = "build_selection_index",
    call = result_call,
    parameters = list(
      index_type = index_type,
      directions = directions[traits],
      units = units[traits],
      selection_intensity = selection_intensity
    ),
    sample_ids = rownames(values),
    inputs = list(
      trait_values = values,
      genetic_cov = G,
      phenotypic_cov = P,
      objective = objective
    ),
    transformations = c("trait direction orientation", index_type,
                        "selection index"),
    quality_gates = c(
      complete_trait_values = all(is.finite(values)),
      genetic_covariance_psd =
        min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) >= -1e-8,
      phenotypic_covariance_pd =
        min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) > 0,
      positive_index_variance = index_variance > 0
    ),
    decision_table = decision,
    uncertainty = coefficient_table
  )
}

.hb_build_desiredgain_index <- function(
    method,
    values,
    G,
    P,
    economic_weights,
    desired_gains,
    directions,
    units,
    n_select,
    quadratic_weights,
    dgsi_control,
    qgsi_control,
    result_call
) {
  if (!requireNamespace("DesiredGainR", quietly = TRUE)) {
    stop(
      "method = '", method,
      "' requires the DesiredGainR package (version 0.5.0 or later).",
      call. = FALSE
    )
  }
  desiredgainr_version <- utils::packageVersion("DesiredGainR")
  if (desiredgainr_version < numeric_version("0.5.0")) {
    stop(
      "DesiredGainR 0.5.0 or later is required; installed version is ",
      as.character(desiredgainr_version), ".",
      call. = FALSE
    )
  }
  traits <- colnames(values)
  init <- data.frame(
    id = rownames(values),
    stringsAsFactors = FALSE
  )
  candidate <- data.frame(
    id = rownames(values),
    values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  lower <- traits[directions[traits] == "decrease"]
  direction_sign <- ifelse(directions[traits] == "increase", 1, -1)

  if (method == "dgsi") {
    if (is.null(desired_gains) || !is.null(economic_weights)) {
      stop("method = 'dgsi' requires desired_gains and no economic_weights.",
           call. = FALSE)
    }
    desired_gains <- .hb_validate_objective_vector(
      desired_gains, traits, "desired_gains"
    )
    protected <- c(
      "init_data", "cand_data", "trait_cols", "dg", "P", "G",
      "id_col", "lower_is_better", "n_select"
    )
    .hb_check_engine_arguments(
      dgsi_control, DesiredGainR::run_dgsi, "dgsi_control", protected
    )
    engine <- do.call(
      DesiredGainR::run_dgsi,
      c(
        list(
          init_data = init,
          cand_data = candidate,
          trait_cols = traits,
          dg = desired_gains,
          P = P,
          G = G,
          id_col = "id",
          lower_is_better = lower,
          n_select = n_select
        ),
        dgsi_control
      )
    )
    decision <- .hb_engine_ranking(
      engine, rownames(values), "SelectionIndex", n_select, exact = FALSE
    )
    engine_transform <- .hb_engine_transformation(engine, traits, direction_sign)
    realised_favourable <- .hb_engine_vector(
      engine$realised_response, traits, "DGSI realised_response"
    )
    realised_original <- realised_favourable * direction_sign
    coefficients <- .hb_engine_vector(
      engine$coefficients, traits, "DGSI coefficients"
    )
    coefficients_original <- coefficients * engine_transform$direction /
      engine_transform$scale
    score_intercept <- -sum(engine_transform$centre * coefficients_original)
    theory <- engine$theoretical_response
    if (!is.list(theory)) {
      stop("DesiredGainR returned no DGSI theoretical_response.",
           call. = FALSE)
    }
    expected_original <- .hb_engine_vector(
      theory$original_units, traits,
      "DGSI theoretical_response$original_units"
    )
    engine_selection_intensity <- .hb_engine_positive_scalar(
      theory$selection_intensity, "DGSI selection intensity",
      allow_zero = TRUE
    )
    .hb_engine_positive_scalar(theory$index_sd, "DGSI index SD")
    coefficient_table <- data.frame(
      trait = traits,
      unit = units[traits],
      direction = directions[traits],
      coefficient = coefficients,
      linear_weight = NA_real_,
      expected_response = expected_original,
      realised_response_sd_favourable = realised_favourable,
      realised_response_sd_original = realised_original,
      stringsAsFactors = FALSE
    )
    response_summary <- rbind(
      data.frame(
        trait = traits,
        direction = directions[traits],
        response = expected_original,
        response_kind = "model_expected_genetic_response",
        response_scale = "original trait units",
        stringsAsFactors = FALSE
      ),
      data.frame(
        trait = traits,
        direction = directions[traits],
        response = realised_original,
        response_kind = "realised_selected_set_differential",
        response_scale = "trait standard deviations in original direction",
        stringsAsFactors = FALSE
      )
    )
    transformations <- c(
      "trait direction orientation",
      "DesiredGainR DGSI optimisation",
      paste("automatic replicate selection:", engine$optimism$selection_rule),
      "DesiredGainR theoretical transmitted response"
    )
  } else {
    if (is.null(economic_weights) || !is.null(desired_gains)) {
      stop(
        "method = 'qgsi' requires economic_weights and no desired_gains.",
        call. = FALSE
      )
    }
    if (is.null(quadratic_weights)) {
      stop("method = 'qgsi' requires quadratic_weights.", call. = FALSE)
    }
    economic_weights <- .hb_validate_objective_vector(
      economic_weights, traits, "economic_weights"
    )
    quadratic_weights <- .hb_validate_symmetric_matrix(
      quadratic_weights, traits, "quadratic_weights"
    )
    protected <- c(
      "init_data", "gebv_data", "trait_cols", "linear_weights", "W",
      "id_col", "lower_is_better", "n_select", "selection_proportion"
    )
    .hb_check_engine_arguments(
      qgsi_control, DesiredGainR::run_qgsi, "qgsi_control", protected
    )
    engine <- do.call(
      DesiredGainR::run_qgsi,
      c(
        list(
          init_data = init,
          gebv_data = candidate,
          trait_cols = traits,
          linear_weights = economic_weights,
          W = quadratic_weights,
          id_col = "id",
          lower_is_better = lower,
          n_select = n_select
        ),
        qgsi_control
      )
    )
    decision <- .hb_engine_ranking(
      engine, rownames(values), "QGSI", n_select, exact = TRUE
    )
    engine_transform <- .hb_engine_transformation(engine, traits, direction_sign)
    gain_table <- as.data.frame(engine$expected_gain_per_trait)
    if (!all(c("Trait", "Expected_Genetic_Gain") %in% names(gain_table)) ||
        nrow(gain_table) != length(traits) || anyDuplicated(gain_table$Trait) ||
        anyNA(match(traits, gain_table$Trait))) {
      stop("DesiredGainR returned an incomplete QGSI expected-gain table.",
           call. = FALSE)
    }
    expected_analysis <- as.numeric(
      gain_table$Expected_Genetic_Gain[match(traits, gain_table$Trait)]
    )
    names(expected_analysis) <- traits
    if (any(!is.finite(expected_analysis))) {
      stop("DesiredGainR returned non-finite QGSI expected gains.",
           call. = FALSE)
    }
    expected_original <- expected_analysis * engine_transform$scale *
      engine_transform$direction
    engine_selection_intensity <- .hb_engine_positive_scalar(
      engine$theoretical_parameters$selection_intensity,
      "QGSI selection intensity", allow_zero = TRUE
    )
    coefficient_table <- data.frame(
      trait = traits,
      unit = units[traits],
      direction = directions[traits],
      coefficient = NA_real_,
      linear_weight = as.numeric(economic_weights[traits]),
      expected_response = expected_original,
      realised_response_sd_favourable = NA_real_,
      realised_response_sd_original = NA_real_,
      stringsAsFactors = FALSE
    )
    response_summary <- data.frame(
      trait = traits,
      direction = directions[traits],
      response = expected_original,
      response_kind = "model_expected_genetic_gain",
      response_scale = "original trait units",
      stringsAsFactors = FALSE
    )
    transformations <- c(
      "trait direction orientation",
      "DesiredGainR QGSI scoring",
      "candidate-specific quadratic contributions",
      "DesiredGainR model-expected gain"
    )
  }
  decision <- decision[order(-decision$selection_index), , drop = FALSE]
  result <- list(
    index_type = method,
    coefficients = coefficient_table,
    scores = decision,
    genetic_cov = G,
    phenotypic_cov = P,
    genetic_cor = stats::cov2cor(G),
    response_summary = response_summary,
    coefficients_original_units = if (method == "dgsi") coefficients_original
      else NULL,
    score_intercept = if (method == "dgsi") score_intercept else NULL,
    scoring_transformation = engine_transform,
    selection_intensity = engine_selection_intensity,
    linear_weights = if (method == "qgsi") economic_weights else NULL,
    quadratic_weights = if (method == "qgsi") quadratic_weights else NULL,
    engine_result = engine,
    contribution_scope = if (method == "qgsi") {
      paste(
        "QGSI quadratic contributions are candidate-specific and are not",
        "interpreted as one global additive marker-effect vector."
      )
    } else {
      paste(
        "Use coefficients_original_units to combine original-unit marker,",
        "haplotype and block effects as a trait-weighted sum."
      )
    }
  )
  .add_hapblockr_contract(
    result = result,
    method = "build_selection_index",
    call = result_call,
    parameters = list(
      index_type = method,
      directions = directions[traits],
      units = units[traits],
      n_select = n_select,
      n_selected = sum(decision$Selected),
      selection_intensity = engine_selection_intensity
    ),
    sample_ids = rownames(values),
    inputs = list(
      trait_values = values,
      genetic_cov = G,
      phenotypic_cov = P,
      objective = if (method == "dgsi") desired_gains else economic_weights,
      quadratic_weights = quadratic_weights,
      engine_control = if (method == "dgsi") dgsi_control else qgsi_control
    ),
    external_tools = list(DesiredGainR = as.character(desiredgainr_version)),
    transformations = transformations,
    quality_gates = c(
      complete_trait_values = all(is.finite(values)),
      finite_scores = all(is.finite(decision$selection_index)),
      objective_magnitudes_non_negative =
        if (method == "dgsi") all(desired_gains >= 0) else
          all(economic_weights >= 0),
      automatic_replicate_choice = if (method == "dgsi")
        sum(engine$replicate_diagnostics$Chosen) == 1L else TRUE
    ),
    decision_table = decision,
    uncertainty = if (method == "dgsi")
      engine$replicate_diagnostics else engine$component_summary
  )
}

.hb_normal_selection_intensity <- function(n_select, n_candidates) {
  proportion <- n_select / n_candidates
  if (proportion >= 1) return(0)
  stats::dnorm(stats::qnorm(1 - proportion)) / proportion
}

#' Fit a Multivariate Genomic BLUP
#'
#' Fits a multivariate GBLUP by restricted maximum likelihood. Without a full
#' sampling covariance, the model uses
#' \eqn{V = Sigma_g \otimes K + Sigma_e \otimes R}, where \eqn{R} carries
#' relative precision. With a full sampling covariance \eqn{S}, the default is
#' \eqn{V = Sigma_g \otimes K + S}; an explicit option adds
#' \eqn{Sigma_e \otimes I}. Missing phenotype cells are allowed; all
#' candidates must be represented in the genomic relationship matrix.
#'
#' @param phenotypes Numeric candidate-by-trait matrix with dimnames, or a
#'   multi-trait result from \code{\link{prepare_breeding_targets}}.
#' @param K Named genomic relationship matrix.
#' @param genetic_cov Optional named genetic covariance matrix. Supply with
#'   \code{residual_cov} to fit with fixed covariance components.
#' @param residual_cov Optional named residual covariance matrix.
#' @param sampling_covariance_mode Treatment of a full sampling covariance
#'   supplied through \code{\link{prepare_breeding_targets}}. With
#'   \code{"sampling_only"} (default), the known matrix is the complete
#'   record-error covariance and no additional residual covariance is fitted.
#'   With \code{"sampling_plus_residual"}, an additional
#'   \eqn{Sigma_e \otimes I} nugget is fitted. The option has no effect when no
#'   full sampling covariance is present.
#' @param estimate_covariances Logical. Estimate covariance matrices by REML.
#' @param maxit Maximum optimiser iterations.
#' @param reltol Relative optimiser tolerance.
#' @param min_reliability Minimum reliability for a recommendation.
#'
#' @return A \code{hapblockr_result} with trait covariance diagnostics,
#'   multivariate predictions, fixed-effect-adjusted PEV, reliability, fitted
#'   means, the REML log-likelihood, and sampling-covariance provenance.
#'   Compare REML likelihoods only between models with the same fixed-effect
#'   design.
#' @export
fit_multitrait_gblup <- function(
    phenotypes,
    K,
    genetic_cov = NULL,
    residual_cov = NULL,
    sampling_covariance_mode = c(
      "sampling_only", "sampling_plus_residual"
    ),
    estimate_covariances = is.null(genetic_cov) && is.null(residual_cov),
    maxit = 200L,
    reltol = 1e-8,
    min_reliability = 0.30
) {
  result_call <- match.call()
  sampling_covariance_mode <- match.arg(sampling_covariance_mode)
  target_bundle <- .hb_unpack_model_targets(phenotypes)
  if (is.null(target_bundle)) {
    Y <- as.matrix(phenotypes)
    precision_matrix <- matrix(1, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    target_provenance <- NULL
    known_sampling_covariance <- NULL
  } else {
    target_rows <- phenotypes$targets
    ids_target <- unique(target_rows$id)
    traits_target <- unique(target_rows$trait)
    Y <- matrix(
      NA_real_, length(ids_target), length(traits_target),
      dimnames = list(ids_target, traits_target)
    )
    precision_matrix <- Y
    for (row in seq_len(nrow(target_rows))) {
      Y[target_rows$id[row], target_rows$trait[row]] <-
        target_rows$model_value[row]
      precision_matrix[target_rows$id[row], target_rows$trait[row]] <-
        target_rows$precision_weight[row]
    }
    target_provenance <- list(
      contract = target_bundle$contract,
      precision = target_bundle$precision_provenance,
      full_sampling_covariance = phenotypes$covariance
    )
    known_sampling_covariance <- NULL
    if (!is.null(phenotypes$covariance)) {
      pair_lookup <- stats::setNames(
        target_rows$record_key,
        paste(target_rows$id, target_rows$trait, sep = "::")
      )
      full_pairs <- paste(
        rep(ids_target, length(traits_target)),
        rep(traits_target, each = length(ids_target)),
        sep = "::"
      )
      full_keys <- unname(pair_lookup[full_pairs])
      observed_positions <- which(!is.na(full_keys))
      covariance_keys <- rownames(phenotypes$covariance)
      missing_covariance_keys <- setdiff(
        full_keys[observed_positions], covariance_keys
      )
      if (length(missing_covariance_keys)) {
        stop(
          "The sampling covariance is missing target record key(s): ",
          paste(head(missing_covariance_keys, 5L), collapse = ", "),
          call. = FALSE
        )
      }
      environment_token <- unique(ifelse(
        is.na(target_rows$environment),
        "ACROSS",
        target_rows$environment
      ))
      full_matrix_keys <- ifelse(
        is.na(full_keys),
        paste(full_pairs, environment_token, sep = "::"),
        full_keys
      )
      known_sampling_covariance <- matrix(
        0, length(full_keys), length(full_keys),
        dimnames = list(full_matrix_keys, full_matrix_keys)
      )
      known_sampling_covariance[
        observed_positions, observed_positions
      ] <- phenotypes$covariance[
        full_keys[observed_positions],
        full_keys[observed_positions],
        drop = FALSE
      ]

      covariance_precision <- 1 / diag(phenotypes$covariance)[
        match(target_rows$record_key, rownames(phenotypes$covariance))
      ]
      precision_group <- .hb_target_precision_group(target_rows)
      covariance_precision <- stats::ave(
        covariance_precision,
        precision_group,
        FUN = function(x) x / mean(x)
      )
      discrepancy <- max(
        abs(covariance_precision - target_rows$precision_weight)
      )
      discrepancy_scale <- max(
        1,
        max(abs(covariance_precision)),
        max(abs(target_rows$precision_weight))
      )
      if (discrepancy > 1e-6 * discrepancy_scale) {
        stop(
          "The diagonal of the full sampling covariance is inconsistent ",
          "with the target precision weights.",
          call. = FALSE
        )
      }
    }
  }
  if (!is.numeric(Y) || !nrow(Y) || ncol(Y) < 2L ||
      is.null(rownames(Y)) || is.null(colnames(Y)))
    stop("phenotypes must be a numeric matrix with candidate row names and ",
         "at least two named trait columns.", call. = FALSE)
  if (anyDuplicated(rownames(Y)) || anyDuplicated(colnames(Y)))
    stop("Phenotype candidate and trait identifiers must be unique.",
         call. = FALSE)
  ids <- rownames(Y)
  traits <- colnames(Y)
  K <- .hb_validate_kernel(K, ids)
  if (any(colSums(is.finite(Y)) < 3L))
    stop("Every trait requires at least three finite observations.",
         call. = FALSE)
  if (length(min_reliability) != 1L || !is.finite(min_reliability) ||
      min_reliability < 0 || min_reliability > 1)
    stop("min_reliability must be in [0, 1].", call. = FALSE)
  estimate_covariances <- isTRUE(estimate_covariances)
  use_residual_component <- is.null(known_sampling_covariance) ||
    sampling_covariance_mode == "sampling_plus_residual"
  if (use_residual_component &&
      xor(is.null(genetic_cov), is.null(residual_cov))) {
    stop("Supply both genetic_cov and residual_cov, or neither.",
         call. = FALSE)
  }
  if (!use_residual_component && !is.null(residual_cov)) {
    stop(
      "residual_cov must be NULL when sampling_covariance_mode = ",
      "'sampling_only'; the full sampling covariance is the record-error ",
      "covariance.",
      call. = FALSE
    )
  }

  n <- nrow(Y)
  q <- ncol(Y)
  y_full <- as.vector(Y)
  observed <- is.finite(y_full)
  y <- y_full[observed]
  precision_full <- as.vector(precision_matrix)
  precision_full[!is.finite(precision_full)] <- 1
  residual_scale <- diag(1 / sqrt(precision_full))
  if (!is.null(known_sampling_covariance)) {
    residual_scale <- diag(length(precision_full))
  }
  sampling_term <- known_sampling_covariance
  X_full <- kronecker(diag(q), matrix(1, nrow = n, ncol = 1L))
  X <- X_full[observed, , drop = FALSE]
  I_n <- diag(n)

  if (is.null(genetic_cov)) {
    empirical <- stats::cov(Y, use = "pairwise.complete.obs")
    empirical[!is.finite(empirical)] <- 0
    diag(empirical) <- pmax(diag(empirical), 1e-6)
    eig <- eigen((empirical + t(empirical)) / 2, symmetric = TRUE)
    empirical <- eig$vectors %*%
      diag(pmax(eig$values, 1e-6), q) %*% t(eig$vectors)
    dimnames(empirical) <- list(traits, traits)
    G_start <- if (use_residual_component) empirical * 0.5 else empirical
    E_start <- if (use_residual_component) empirical * 0.5 else NULL
  } else {
    G_start <- .hb_validate_covariance(
      genetic_cov, traits, "genetic_cov", positive_definite = TRUE
    )
    E_start <- if (use_residual_component) {
      .hb_validate_covariance(
        residual_cov, traits, "residual_cov", positive_definite = TRUE
      )
    } else {
      NULL
    }
  }

  covariance_parameter_count <- q * (q + 1L) / 2L
  fixed_effect_rank <- qr(X)$rank
  reml_residual_df <- length(y) - fixed_effect_rank
  if (reml_residual_df <= 0L) {
    stop("The fixed-effect design leaves no residual degrees of freedom.",
         call. = FALSE)
  }
  objective <- function(parameters) {
    G_t <- .hb_parameters_to_cov(
      parameters[seq_len(covariance_parameter_count)], q
    )
    V_full <- kronecker(G_t, K)
    if (use_residual_component) {
      E_t <- .hb_parameters_to_cov(
        parameters[-seq_len(covariance_parameter_count)], q
      )
      V_full <- V_full +
        residual_scale %*% kronecker(E_t, I_n) %*% residual_scale
    }
    if (!is.null(sampling_term)) V_full <- V_full + sampling_term
    V <- V_full[observed, observed, drop = FALSE]
    fit <- .hb_gls_components(V, X, y)
    if (is.null(fit)) return(.Machine$double.xmax / 100)
    reml_residual_df * log(2 * pi) +
      fit$logdet_V + fit$logdet_X +
      sum(fit$residual * (fit$V_inv %*% fit$residual))
  }

  start <- .hb_cov_to_parameters(G_start)
  if (use_residual_component) {
    start <- c(start, .hb_cov_to_parameters(E_start))
  }
  optimiser <- NULL
  if (estimate_covariances) {
    optimiser <- stats::optim(
      start,
      objective,
      method = "BFGS",
      control = list(maxit = as.integer(maxit), reltol = reltol)
    )
    if (optimiser$convergence != 0L || !is.finite(optimiser$value))
      stop("Multivariate REML optimisation failed: convergence code ",
           optimiser$convergence, ". ", optimiser$message, call. = FALSE)
    parameters <- optimiser$par
  } else {
    parameters <- start
  }
  split_at <- covariance_parameter_count
  Sigma_g <- .hb_parameters_to_cov(parameters[seq_len(split_at)], q)
  Sigma_e <- if (use_residual_component) {
    .hb_parameters_to_cov(parameters[-seq_len(split_at)], q)
  } else {
    NULL
  }
  dimnames(Sigma_g) <- list(traits, traits)
  if (!is.null(Sigma_e)) {
    dimnames(Sigma_e) <- list(traits, traits)
  }

  C <- kronecker(Sigma_g, K)
  V_full <- C
  if (use_residual_component) {
    V_full <- V_full +
      residual_scale %*% kronecker(Sigma_e, I_n) %*% residual_scale
  }
  if (!is.null(sampling_term)) V_full <- V_full + sampling_term
  V <- V_full[observed, observed, drop = FALSE]
  fit <- .hb_gls_components(V, X, y)
  if (is.null(fit))
    stop("The fitted multivariate covariance system is singular.",
         call. = FALSE)
  C_observed <- C[, observed, drop = FALSE]
  u <- as.numeric(C_observed %*% fit$V_inv %*% fit$residual)
  C_Vinv <- C_observed %*% fit$V_inv
  fixed_projection <- C_Vinv %*% X
  pev <- pmax(
    diag(C) -
      rowSums(C_Vinv * C_observed) +
      rowSums(
        (fixed_projection %*% fit$XtVinvX_inv) *
          fixed_projection
      ),
    0
  )
  fitted_mean <- as.numeric(X_full %*% fit$beta)
  prediction <- fitted_mean + u
  prior_variance <- rep(diag(Sigma_g), each = n) * rep(diag(K), q)
  reliability_raw <- 1 - pev / prior_variance
  reliability <- pmin(1, pmax(0, reliability_raw))
  decision <- data.frame(
    id = rep(ids, q),
    trait = rep(traits, each = n),
    observed = y_full,
    fitted_mean = fitted_mean,
    breeding_value = u,
    prediction = prediction,
    PEV = pev,
    reliability = reliability,
    min_reliability = min_reliability,
    recommendable = is.finite(reliability) &
      reliability >= min_reliability,
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    trait1 = rep(traits, each = q),
    trait2 = rep(traits, q),
    genetic_covariance = as.vector(Sigma_g),
    residual_covariance = if (is.null(Sigma_e)) {
      rep(NA_real_, q * q)
    } else {
      as.vector(Sigma_e)
    },
    genetic_correlation = as.vector(stats::cov2cor(Sigma_g)),
    residual_correlation = if (is.null(Sigma_e)) {
      rep(NA_real_, q * q)
    } else {
      as.vector(stats::cov2cor(Sigma_e))
    },
    stringsAsFactors = FALSE
  )
  if (!is.null(target_provenance)) {
    target_provenance$sampling_covariance_mode <-
      if (is.null(known_sampling_covariance)) {
        "none"
      } else {
        sampling_covariance_mode
      }
  }
  result <- list(
    predictions = decision,
    genetic_cov = Sigma_g,
    residual_cov = Sigma_e,
    genetic_cor = stats::cov2cor(Sigma_g),
    residual_cor = if (is.null(Sigma_e)) NULL else stats::cov2cor(Sigma_e),
    sampling_covariance = known_sampling_covariance,
    trait_means = stats::setNames(fit$beta, traits),
    optimiser = optimiser,
    log_likelihood = -0.5 * objective(parameters),
    fixed_effect_rank = fixed_effect_rank,
    reml_residual_df = reml_residual_df,
    sampling_covariance_mode = if (is.null(known_sampling_covariance)) {
      "none"
    } else {
      sampling_covariance_mode
    },
    n_observed = sum(observed),
    target_provenance = target_provenance
  )
  .add_hapblockr_contract(
    result = result,
    method = "fit_multitrait_gblup",
    call = result_call,
    parameters = list(
      estimate_covariances = estimate_covariances,
      sampling_covariance_mode = if (is.null(known_sampling_covariance)) {
        "none"
      } else {
        sampling_covariance_mode
      },
      maxit = as.integer(maxit),
      reltol = reltol,
      min_reliability = min_reliability
    ),
    sample_ids = ids,
    inputs = list(phenotypes = Y, genomic_relationship = K),
    transformations = c(
      if (is.null(known_sampling_covariance)) {
        "precision-weighted multivariate REML"
      } else if (use_residual_component) {
        "full sampling covariance plus residual-nugget multivariate REML"
      } else {
        "full sampling-covariance multivariate REML"
      },
      "multivariate genomic BLUP"
    ),
    quality_gates = c(
      covariance_optimisation_converged =
        !estimate_covariances || optimiser$convergence == 0L,
      genetic_covariance_pd =
        min(eigen(Sigma_g, symmetric = TRUE,
                  only.values = TRUE)$values) > 0,
      residual_covariance_pd = if (is.null(Sigma_e)) {
        TRUE
      } else {
        min(eigen(Sigma_e, symmetric = TRUE,
                  only.values = TRUE)$values) > 0
      },
      sampling_covariance_used = if (is.null(known_sampling_covariance)) {
        TRUE
      } else {
        any(abs(known_sampling_covariance) > 0)
      },
      reliability_bounded =
        all(reliability >= 0 & reliability <= 1, na.rm = TRUE)
    ),
    decision_table = decision,
    uncertainty = diagnostics
  )
}

#' Build an Environmental Similarity Kernel
#'
#' @param environmental_covariates Data frame containing one row per
#'   environment.
#' @param environment_col Column containing unique environment IDs.
#' @param covariate_cols Numeric environmental covariates. Defaults to every
#'   numeric column other than \code{environment_col}.
#' @param standardise Logical. Centre and scale covariates before calculating
#'   the linear kernel.
#'
#' @return A positive-semidefinite environment-by-environment kernel.
#' @export
build_environment_kernel <- function(
    environmental_covariates,
    environment_col = "environment",
    covariate_cols = NULL,
    standardise = TRUE
) {
  x <- as.data.frame(environmental_covariates, stringsAsFactors = FALSE)
  if (!environment_col %in% names(x))
    stop("environment_col was not found.", call. = FALSE)
  environments <- as.character(x[[environment_col]])
  if (anyNA(environments) || any(!nzchar(environments)) ||
      anyDuplicated(environments))
    stop("Environment IDs must be non-missing, non-empty, and unique.",
         call. = FALSE)
  if (is.null(covariate_cols)) {
    covariate_cols <- setdiff(
      names(x)[vapply(x, is.numeric, logical(1L))],
      environment_col
    )
  }
  if (!length(covariate_cols) || !all(covariate_cols %in% names(x)))
    stop("At least one valid numeric environmental covariate is required.",
         call. = FALSE)
  X <- as.matrix(x[, covariate_cols, drop = FALSE])
  if (!is.numeric(X) || any(!is.finite(X)))
    stop("Environmental covariates must be finite numeric values.",
         call. = FALSE)
  if (isTRUE(standardise)) {
    sd_x <- apply(X, 2L, stats::sd)
    if (any(!is.finite(sd_x) | sd_x == 0))
      stop("Every environmental covariate must vary when standardise = TRUE.",
           call. = FALSE)
    X <- scale(X)
  }
  K <- tcrossprod(X) / ncol(X)
  dimnames(K) <- list(environments, environments)
  attr(K, "covariates") <- covariate_cols
  attr(K, "standardised") <- isTRUE(standardise)
  K
}

#' Fit a Reaction-Norm Genomic GxE Model
#'
#' Fits a genomic main-effect plus reaction-norm model with record covariance
#' \eqn{sigma_g^2 K_g + sigma_{ge}^2(K_g \circ K_e) + sigma_e^2 I}.
#' Environment means are fitted as fixed effects by default.
#'
#' @param data Data frame with genotype, environment, and phenotype columns.
#' @param K Named genomic relationship matrix.
#' @param K_environment Optional named environment kernel. The identity matrix
#'   gives independent environment-specific deviations.
#' @param id_col,environment_col,phenotype_col Column names.
#' @param precision_col Optional positive precision-weight column. Prepared
#'   breeding targets supply their normalised precision automatically.
#' @param environment_fixed Fit an independent fixed mean per environment.
#' @param min_reliability Minimum reliability for recommendation.
#' @param maxit Maximum REML optimiser iterations.
#'
#' @return A \code{hapblockr_result} with variance components,
#'   environment-specific genomic predictions, fixed-effect-adjusted PEV,
#'   reliability, and the REML log-likelihood. Compare REML likelihoods only
#'   between models with the same fixed-effect design.
#' @export
fit_gxe_gblup <- function(
    data,
    K,
    K_environment = NULL,
    id_col = "id",
    environment_col = "environment",
    phenotype_col = "phenotype",
    precision_col = NULL,
    environment_fixed = TRUE,
    min_reliability = 0.30,
    maxit = 200L
) {
  result_call <- match.call()
  target_bundle <- .hb_unpack_model_targets(data, allow_environment = TRUE)
  if (!is.null(target_bundle)) {
    records <- as.data.frame(data$targets, stringsAsFactors = FALSE)
    if (length(unique(records$trait)) != 1L)
      stop("fit_gxe_gblup() accepts one prepared trait per call.",
           call. = FALSE)
    if (anyNA(records$environment))
      stop("GxE fitting requires an environment for every prepared target.",
           call. = FALSE)
    records$phenotype <- records$model_value
    records$precision <- records$precision_weight
    target_provenance <- list(
      contract = target_bundle$contract,
      precision = target_bundle$precision_provenance
    )
  } else {
    records <- as.data.frame(data, stringsAsFactors = FALSE)
    required <- c(id_col, environment_col, phenotype_col)
    if (!all(required %in% names(records)))
      stop("data is missing required columns: ",
           paste(setdiff(required, names(records)), collapse = ", "),
           call. = FALSE)
    records$id <- as.character(records[[id_col]])
    records$environment <- as.character(records[[environment_col]])
    records$phenotype <- as.numeric(records[[phenotype_col]])
    records$precision <- if (is.null(precision_col)) {
      rep(1, nrow(records))
    } else {
      as.numeric(records[[precision_col]])
    }
    target_provenance <- NULL
  }
  valid <- nzchar(records$id) & nzchar(records$environment) &
    is.finite(records$phenotype) & is.finite(records$precision) &
    records$precision > 0
  if (!all(valid))
    stop("Genotype IDs, environment IDs, and phenotypes must be complete and ",
         "finite.", call. = FALSE)
  ids <- unique(records$id)
  environments <- unique(records$environment)
  K <- .hb_validate_kernel(K, ids)
  if (is.null(K_environment)) {
    K_environment <- diag(length(environments))
    dimnames(K_environment) <- list(environments, environments)
  } else {
    K_environment <- .hb_validate_kernel(
      K_environment, environments, "K_environment"
    )
  }
  if (length(min_reliability) != 1L || !is.finite(min_reliability) ||
      min_reliability < 0 || min_reliability > 1)
    stop("min_reliability must be in [0, 1].", call. = FALSE)

  Kg_records <- K[records$id, records$id, drop = FALSE]
  Ke_records <- K_environment[
    records$environment, records$environment, drop = FALSE
  ]
  Kge_records <- Kg_records * Ke_records
  y <- records$phenotype
  records$precision <- records$precision / mean(records$precision)
  residual_shape <- diag(1 / records$precision)
  X <- if (isTRUE(environment_fixed)) {
    environment_factor <- factor(
      records$environment,
      levels = environments
    )
    design <- stats::model.matrix(~ environment_factor - 1L)
    colnames(design) <- environments
    design
  } else {
    matrix(1, nrow(records), 1L, dimnames = list(NULL, "intercept"))
  }
  fixed_effect_rank <- qr(X)$rank
  reml_residual_df <- length(y) - fixed_effect_rank
  if (reml_residual_df <= 0L) {
    stop("The fixed-effect design leaves no residual degrees of freedom.",
         call. = FALSE)
  }
  start_variance <- max(stats::var(y), 1e-6)
  objective <- function(log_variances) {
    variances <- exp(log_variances)
    V <- variances[1L] * Kg_records +
      variances[2L] * Kge_records +
      variances[3L] * residual_shape
    fit <- .hb_gls_components(V, X, y)
    if (is.null(fit)) return(.Machine$double.xmax / 100)
    reml_residual_df * log(2 * pi) +
      fit$logdet_V + fit$logdet_X +
      sum(fit$residual * (fit$V_inv %*% fit$residual))
  }
  optimiser <- stats::optim(
    log(c(start_variance * 0.4, start_variance * 0.3,
          start_variance * 0.3)),
    objective,
    method = "BFGS",
    control = list(maxit = as.integer(maxit), reltol = 1e-8)
  )
  if (optimiser$convergence != 0L || !is.finite(optimiser$value))
    stop("GxE REML optimisation failed: convergence code ",
         optimiser$convergence, ". ", optimiser$message, call. = FALSE)
  variances <- exp(optimiser$par)
  names(variances) <- c("genomic_main", "genomic_by_environment", "residual")
  V <- variances[1L] * Kg_records +
    variances[2L] * Kge_records +
    variances[3L] * residual_shape
  fit <- .hb_gls_components(V, X, y)
  if (is.null(fit))
    stop("The fitted GxE covariance system is singular.", call. = FALSE)

  grid <- expand.grid(
    id = ids,
    environment = environments,
    stringsAsFactors = FALSE
  )
  C_main <- variances[1L] * K[grid$id, records$id, drop = FALSE]
  C_gxe <- variances[2L] *
    K[grid$id, records$id, drop = FALSE] *
    K_environment[
      grid$environment, records$environment, drop = FALSE
    ]
  C_total <- C_main + C_gxe
  u_main <- as.numeric(C_main %*% fit$V_inv %*% fit$residual)
  u_gxe <- as.numeric(C_gxe %*% fit$V_inv %*% fit$residual)
  fixed <- if (isTRUE(environment_fixed)) {
    stats::setNames(fit$beta, colnames(X))[
      grid$environment
    ]
  } else {
    rep(fit$beta[1L], nrow(grid))
  }
  prior <- variances[1L] * diag(K)[match(grid$id, rownames(K))] +
    variances[2L] *
      diag(K)[match(grid$id, rownames(K))] *
      diag(K_environment)[
        match(grid$environment, rownames(K_environment))
      ]
  reduction <- rowSums((C_total %*% fit$V_inv) * C_total)
  fixed_projection <- C_total %*% fit$V_inv %*% X
  fixed_effect_correction <- rowSums(
    (fixed_projection %*% fit$XtVinvX_inv) * fixed_projection
  )
  pev <- pmax(prior - reduction + fixed_effect_correction, 0)
  reliability <- pmin(1, pmax(0, 1 - pev / prior))
  grid$environment_mean <- as.numeric(fixed)
  grid$genomic_main <- u_main
  grid$genomic_by_environment <- u_gxe
  grid$prediction <- grid$environment_mean + u_main + u_gxe
  grid$PEV <- pev
  grid$reliability <- reliability
  grid$min_reliability <- min_reliability
  grid$recommendable <- reliability >= min_reliability

  observed_predictions <- merge(
    records[c("id", "environment", "phenotype")],
    grid,
    by = c("id", "environment"),
    all.x = TRUE,
    sort = FALSE
  )
  environment_stability <- do.call(rbind, lapply(
    split(grid, grid$id),
    function(x) data.frame(
      id = x$id[1L],
      mean_prediction = mean(x$prediction),
      sd_across_environments = if (nrow(x) > 1L)
        stats::sd(x$prediction) else 0,
      worst_environment_prediction = min(x$prediction),
      min_reliability = min(x$reliability),
      recommendation_stable = all(x$recommendable),
      stringsAsFactors = FALSE
    )
  ))
  rownames(environment_stability) <- NULL
  result <- list(
    predictions = grid,
    observed_predictions = observed_predictions,
    environment_stability = environment_stability,
    variance_components = variances,
    optimiser = optimiser,
    log_likelihood = -0.5 * objective(optimiser$par),
    fixed_effect_rank = fixed_effect_rank,
    reml_residual_df = reml_residual_df,
    K_environment = K_environment,
    # genomic_main is mathematically constant within ID, but BLAS may return
    # values that differ at the final stored decimal between environment rows.
    # unique(grid[c("id", "genomic_main")]) can therefore retain duplicate
    # IDs on some R/BLAS combinations. Select the first row by ID explicitly.
    across_environment_predictions = {
      across <- grid[
        !duplicated(grid$id), c("id", "genomic_main"), drop = FALSE
      ]
      rownames(across) <- NULL
      across
    },
    target_provenance = target_provenance
  )
  .add_hapblockr_contract(
    result = result,
    method = "fit_gxe_gblup",
    call = result_call,
    parameters = list(
      environment_fixed = isTRUE(environment_fixed),
      min_reliability = min_reliability,
      maxit = as.integer(maxit)
    ),
    sample_ids = ids,
    inputs = list(
      phenotypes = records[
        c("id", "environment", "phenotype", "precision")
      ],
      genomic_relationship = K,
      environment_relationship = K_environment
    ),
    transformations = c(
      "precision-weighted reaction-norm REML",
      "environment-specific GBLUP",
      "across-environment genomic main effect"
    ),
    quality_gates = c(
      optimiser_converged = optimiser$convergence == 0L,
      positive_variance_components = all(variances > 0),
      complete_prediction_grid =
        nrow(grid) == length(ids) * length(environments),
      reliability_bounded = all(reliability >= 0 & reliability <= 1)
    ),
    decision_table = environment_stability,
    uncertainty = grid
  )
}
