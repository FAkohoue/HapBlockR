# ==============================================================================
# Common decision-result and provenance contract
# ==============================================================================

.normalise_quality_gates <- function(quality_gates) {
  if (is.null(quality_gates) || !length(quality_gates)) {
    return(data.frame(
      gate = character(), passed = logical(), detail = character(),
      stringsAsFactors = FALSE
    ))
  }
  if (is.logical(quality_gates)) {
    gate_names <- names(quality_gates)
    if (is.null(gate_names))
      gate_names <- paste0("gate_", seq_along(quality_gates))
    return(data.frame(
      gate = gate_names,
      passed = unname(quality_gates),
      detail = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(quality_gates) ||
      !all(c("gate", "passed") %in% names(quality_gates)))
    stop("quality_gates must be a named logical vector or a data frame ",
         "containing 'gate' and 'passed'.", call. = FALSE)
  quality_gates$gate <- as.character(quality_gates$gate)
  quality_gates$passed <- as.logical(quality_gates$passed)
  if (!"detail" %in% names(quality_gates))
    quality_gates$detail <- NA_character_
  quality_gates
}

.hash_result_inputs <- function(inputs) {
  if (is.null(inputs) || !length(inputs)) return(character())
  if (is.null(names(inputs)) || any(!nzchar(names(inputs))))
    stop("inputs must be a named list.", call. = FALSE)
  vapply(inputs, digest::digest, character(1L),
         algo = "sha256", serialize = TRUE)
}

.add_hapblockr_contract <- function(
    result,
    method,
    call,
    parameters = list(),
    seed = NULL,
    sample_ids = character(),
    variant_ids = character(),
    inputs = list(),
    transformations = character(),
    quality_gates = NULL,
    warnings = character(),
    fallbacks = character(),
    excluded_records = data.frame(),
    decision_table = data.frame(),
    uncertainty = data.frame(),
    external_tools = list()
) {
  if (!is.list(result))
    stop("A hapblockr_result must be based on a list.", call. = FALSE)
  if (length(method) != 1L || is.na(method) || !nzchar(method))
    stop("method must be one non-empty string.", call. = FALSE)

  sample_ids <- as.character(sample_ids)
  variant_ids <- as.character(variant_ids)
  gates <- .normalise_quality_gates(quality_gates)
  gate_pass <- !nrow(gates) || all(!is.na(gates$passed) & gates$passed)

  result$result_contract <- list(
    schema_version = "1.0.0",
    method = as.character(method),
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    call = paste(deparse(call, width.cutoff = 500L), collapse = " "),
    parameters = parameters,
    random_seed = if (is.null(seed)) NA_integer_ else as.integer(seed)[1L],
    identifiers = list(
      sample_ids = sample_ids,
      variant_ids = variant_ids
    ),
    input_hashes = .hash_result_inputs(inputs),
    transformations = as.character(transformations),
    software = list(
      R = R.version.string,
      HapBlockR = as.character(utils::packageVersion("HapBlockR")),
      external_tools = external_tools
    ),
    quality_gates = gates,
    warnings = as.character(warnings),
    fallbacks = as.character(fallbacks),
    excluded_records = excluded_records,
    decision_table = decision_table,
    uncertainty = uncertainty,
    validation_status = if (gate_pass) "passed" else "failed"
  )

  old_class <- class(result)
  old_class <- setdiff(old_class, c("hapblockr_result", "list"))
  class(result) <- unique(c(old_class, "hapblockr_result", "list"))
  result
}

#' Validate a HapBlockR Result Contract
#'
#' Checks the common result schema, immutable identifiers, input hashes,
#' quality-control gates, and decision table. With \code{strict = TRUE}, a
#' failed check stops execution; otherwise the function returns the complete
#' validation report.
#'
#' @param object Object inheriting from \code{"hapblockr_result"}.
#' @param strict Logical. Stop when any validation check fails. Default
#'   \code{TRUE}.
#' @param ... Reserved for future schema versions.
#'
#' @return A data frame containing \code{check}, \code{passed}, and
#'   \code{detail}.
#' @export
validate_hapblockr_result <- function(object, strict = TRUE, ...) {
  if (!inherits(object, "hapblockr_result"))
    stop("object must inherit from 'hapblockr_result'.", call. = FALSE)
  contract <- object$result_contract
  required <- c(
    "schema_version", "method", "call", "parameters", "random_seed",
    "identifiers", "input_hashes", "software", "quality_gates",
    "warnings", "fallbacks", "excluded_records", "decision_table",
    "uncertainty", "validation_status"
  )
  samples <- contract$identifiers$sample_ids
  variants <- contract$identifiers$variant_ids
  checks <- data.frame(
    check = c(
      "required_fields", "schema_version", "method", "sample_ids_unique",
      "variant_ids_unique", "input_hashes_sha256", "quality_gates",
      "validation_status"
    ),
    passed = c(
      is.list(contract) && all(required %in% names(contract)),
      identical(contract$schema_version, "1.0.0"),
      length(contract$method) == 1L && !is.na(contract$method) &&
        nzchar(contract$method),
      !anyNA(samples) && !any(!nzchar(samples)) && !anyDuplicated(samples),
      !anyNA(variants) && !any(!nzchar(variants)) && !anyDuplicated(variants),
      !length(contract$input_hashes) ||
        all(grepl("^[0-9a-f]{64}$", contract$input_hashes)),
      !nrow(contract$quality_gates) ||
        all(!is.na(contract$quality_gates$passed) &
              contract$quality_gates$passed),
      identical(contract$validation_status, "passed")
    ),
    detail = c(
      paste(setdiff(required, names(contract)), collapse = ", "),
      as.character(contract$schema_version),
      as.character(contract$method),
      paste(length(samples), "sample identifier(s)"),
      paste(length(variants), "variant identifier(s)"),
      paste(length(contract$input_hashes), "input hash(es)"),
      paste(nrow(contract$quality_gates), "quality gate(s)"),
      as.character(contract$validation_status)
    ),
    stringsAsFactors = FALSE
  )
  if (isTRUE(strict) && any(!checks$passed))
    stop("hapblockr_result validation failed: ",
         paste(checks$check[!checks$passed], collapse = ", "),
         call. = FALSE)
  checks
}

#' Validate an Object
#'
#' @param object Object to validate.
#' @param ... Additional arguments passed to the class method.
#' @return Class-specific validation output.
#' @export
validate <- function(object, ...) UseMethod("validate")

#' @export
validate.hapblockr_result <- function(object, ...) {
  validate_hapblockr_result(object, ...)
}

#' @export
print.hapblockr_result <- function(x, ...) {
  contract <- x$result_contract
  cat("HapBlockR result\n")
  cat("  Method:", contract$method, "\n")
  cat("  Schema:", contract$schema_version, "\n")
  cat("  Validation:", contract$validation_status, "\n")
  cat("  Decisions:", nrow(contract$decision_table), "\n")
  cat("  Quality gates:", nrow(contract$quality_gates), "\n")
  invisible(x)
}

#' @export
summary.hapblockr_result <- function(object, ...) {
  contract <- object$result_contract
  list(
    method = contract$method,
    schema_version = contract$schema_version,
    validation_status = contract$validation_status,
    quality_gates = contract$quality_gates,
    warnings = contract$warnings,
    fallbacks = contract$fallbacks,
    n_decisions = nrow(contract$decision_table),
    n_uncertainty_records = nrow(contract$uncertainty)
  )
}

#' @export
as.data.frame.hapblockr_result <- function(x, row.names = NULL,
                                           optional = FALSE, ...) {
  out <- x$result_contract$decision_table
  if (!is.data.frame(out))
    stop("The result contract does not contain a decision data frame.",
         call. = FALSE)
  out
}

#' @export
plot.hapblockr_result <- function(x, ...) {
  decision <- x$result_contract$decision_table
  numeric_cols <- names(decision)[vapply(decision, is.numeric, logical(1L))]
  if (!nrow(decision) || !length(numeric_cols))
    stop("No numeric decision column is available to plot.", call. = FALSE)
  values <- decision[[numeric_cols[1L]]]
  labels <- if ("id" %in% names(decision)) decision$id else seq_along(values)
  graphics::barplot(
    values,
    names.arg = labels,
    las = 2L,
    ylab = numeric_cols[1L],
    main = paste("HapBlockR decision:", x$result_contract$method),
    ...
  )
  invisible(x)
}
