# ==============================================================================
# Breeding metadata validation and interoperable recommendation bundles
# ==============================================================================

.hb_nonempty_id <- function(x) {
  is.character(x) && !anyNA(x) && all(nzchar(x))
}

.hb_metadata_issue <- function(table, row, field, code, detail) {
  data.frame(
    table = table,
    row = as.integer(row),
    field = field,
    code = code,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

#' Validate Traceable Breeding Metadata
#'
#' Validates canonical germplasm, trait, environment, and observation tables
#' before modelling or exchange. Referential integrity, ontology identifiers,
#' units, programme, study, trial, environment, and germplasm identities are
#' checked without calling a remote service.
#'
#' @param germplasm Data frame containing \code{germplasm_id} and
#'   \code{breeding_program_id}.
#' @param traits Data frame containing \code{trait_id}, \code{trait_name},
#'   \code{unit}, \code{direction}, and \code{ontology_id}.
#' @param environments Data frame containing \code{environment_id},
#'   \code{study_id}, \code{trial_id}, and \code{location_id}.
#' @param observations Data frame containing \code{observation_id},
#'   \code{germplasm_id}, \code{trait_id}, \code{environment_id},
#'   \code{value}, and \code{unit}.
#' @param strict Logical. Stop if any issue is found. When \code{FALSE}, return
#'   a failed result contract containing the full issue ledger.
#'
#' @return A \code{hapblockr_result} with validated canonical tables,
#'   standards mapping, and an issue ledger.
#' @export
validate_breeding_metadata <- function(
    germplasm,
    traits,
    environments,
    observations,
    strict = TRUE
) {
  result_call <- match.call()
  tables <- list(
    germplasm = as.data.frame(germplasm, stringsAsFactors = FALSE),
    traits = as.data.frame(traits, stringsAsFactors = FALSE),
    environments = as.data.frame(environments, stringsAsFactors = FALSE),
    observations = as.data.frame(observations, stringsAsFactors = FALSE)
  )
  required <- list(
    germplasm = c("germplasm_id", "breeding_program_id"),
    traits = c("trait_id", "trait_name", "unit", "direction", "ontology_id"),
    environments = c(
      "environment_id", "study_id", "trial_id", "location_id"
    ),
    observations = c(
      "observation_id", "germplasm_id", "trait_id", "environment_id",
      "value", "unit"
    )
  )
  missing_columns <- unlist(lapply(names(tables), function(table) {
    missing <- setdiff(required[[table]], names(tables[[table]]))
    if (!length(missing)) return(character())
    paste0(table, ".", missing)
  }))
  if (length(missing_columns))
    stop("Missing canonical metadata column(s): ",
         paste(missing_columns, collapse = ", "), call. = FALSE)

  issues <- list()
  add_issue <- function(x) issues[[length(issues) + 1L]] <<- x
  id_fields <- list(
    germplasm = c("germplasm_id", "breeding_program_id"),
    traits = c("trait_id", "trait_name", "unit", "direction", "ontology_id"),
    environments = c(
      "environment_id", "study_id", "trial_id", "location_id"
    ),
    observations = c(
      "observation_id", "germplasm_id", "trait_id", "environment_id", "unit"
    )
  )
  for (table in names(id_fields)) {
    for (field in id_fields[[table]]) {
      values <- as.character(tables[[table]][[field]])
      bad <- which(is.na(values) | !nzchar(values))
      for (row in bad)
        add_issue(.hb_metadata_issue(
          table, row, field, "missing_or_empty",
          "Required metadata value is missing or empty."
        ))
      tables[[table]][[field]] <- values
    }
  }
  unique_ids <- c(
    germplasm = "germplasm_id",
    traits = "trait_id",
    environments = "environment_id",
    observations = "observation_id"
  )
  for (table in names(unique_ids)) {
    field <- unique_ids[[table]]
    duplicate <- duplicated(tables[[table]][[field]]) |
      duplicated(tables[[table]][[field]], fromLast = TRUE)
    for (row in which(duplicate))
      add_issue(.hb_metadata_issue(
        table, row, field, "duplicate_identifier",
        "Identifier must be unique within its canonical table."
      ))
  }

  invalid_direction <- which(
    !tables$traits$direction %in% c("increase", "decrease")
  )
  for (row in invalid_direction)
    add_issue(.hb_metadata_issue(
      "traits", row, "direction", "invalid_direction",
      "Direction must be 'increase' or 'decrease'."
    ))
  ontology_pattern <- "^[A-Za-z][A-Za-z0-9_.-]*:[A-Za-z0-9_.-]+$"
  invalid_ontology <- which(
    !grepl(ontology_pattern, tables$traits$ontology_id)
  )
  for (row in invalid_ontology)
    add_issue(.hb_metadata_issue(
      "traits", row, "ontology_id", "invalid_ontology_identifier",
      "Use a stable CURIE such as CO_321:0000012."
    ))

  references <- list(
    germplasm_id = tables$germplasm$germplasm_id,
    trait_id = tables$traits$trait_id,
    environment_id = tables$environments$environment_id
  )
  for (field in names(references)) {
    unknown <- which(!tables$observations[[field]] %in% references[[field]])
    for (row in unknown)
      add_issue(.hb_metadata_issue(
        "observations", row, field, "unknown_foreign_key",
        paste(field, "does not occur in its canonical table.")
      ))
  }
  numeric_values <- suppressWarnings(as.numeric(tables$observations$value))
  invalid_value <- which(
    !is.na(tables$observations$value) & !is.finite(numeric_values)
  )
  for (row in invalid_value)
    add_issue(.hb_metadata_issue(
      "observations", row, "value", "non_numeric_value",
      "Observed values must be numeric or explicitly missing."
    ))
  tables$observations$value <- numeric_values

  trait_match <- match(
    tables$observations$trait_id,
    tables$traits$trait_id
  )
  expected_unit <- tables$traits$unit[trait_match]
  unit_mismatch <- which(
    !is.na(trait_match) &
      tables$observations$unit != expected_unit
  )
  for (row in unit_mismatch)
    add_issue(.hb_metadata_issue(
      "observations", row, "unit", "unit_mismatch",
      paste0(
        "Observation unit '", tables$observations$unit[row],
        "' does not match trait unit '", expected_unit[row], "'."
      )
    ))

  issue_table <- if (length(issues)) do.call(rbind, issues) else data.frame(
    table = character(), row = integer(), field = character(),
    code = character(), detail = character(), stringsAsFactors = FALSE
  )
  mapping <- data.frame(
    canonical_table = c(
      "germplasm", "germplasm", "traits", "traits", "traits",
      "environments", "environments", "environments", "observations",
      "observations", "observations"
    ),
    canonical_field = c(
      "germplasm_id", "breeding_program_id", "trait_id", "ontology_id",
      "unit", "environment_id", "study_id", "trial_id", "observation_id",
      "germplasm_id", "value"
    ),
    standard = c(
      rep("BrAPI", 4L), "MIAPPE", rep("BrAPI", 6L)
    ),
    standard_field = c(
      "germplasmDbId", "breedingProgramDbId", "observationVariableDbId",
      "ontologyDbId", "trait/scale unit", "observationUnitDbId",
      "studyDbId", "trialDbId", "observationDbId", "germplasmDbId",
      "value"
    ),
    stringsAsFactors = FALSE
  )
  summary_table <- data.frame(
    table = names(tables),
    rows = vapply(tables, nrow, integer(1L)),
    issues = vapply(names(tables), function(table)
      sum(issue_table$table == table), integer(1L)),
    valid = vapply(names(tables), function(table)
      !any(issue_table$table == table), logical(1L)),
    stringsAsFactors = FALSE
  )
  valid <- !nrow(issue_table)
  result <- list(
    germplasm = tables$germplasm,
    traits = tables$traits,
    environments = tables$environments,
    observations = tables$observations,
    standards_mapping = mapping,
    issues = issue_table,
    summary = summary_table
  )
  result <- .add_hapblockr_contract(
    result = result,
    method = "validate_breeding_metadata",
    call = result_call,
    parameters = list(strict = isTRUE(strict)),
    sample_ids = tables$germplasm$germplasm_id,
    inputs = tables,
    transformations = c(
      "canonical metadata normalisation",
      "referential-integrity validation",
      "unit and ontology validation"
    ),
    quality_gates = c(
      unique_identifiers =
        !any(issue_table$code == "duplicate_identifier"),
      referential_integrity =
        !any(issue_table$code == "unknown_foreign_key"),
      ontology_identifiers =
        !any(issue_table$code == "invalid_ontology_identifier"),
      units_consistent = !any(issue_table$code == "unit_mismatch"),
      required_metadata_complete =
        !any(issue_table$code == "missing_or_empty"),
      numeric_observations =
        !any(issue_table$code == "non_numeric_value")
    ),
    warnings = if (valid) character() else
      paste(nrow(issue_table), "metadata issue(s)"),
    excluded_records = issue_table,
    decision_table = summary_table,
    uncertainty = data.frame()
  )
  if (isTRUE(strict) && !valid)
    stop("Breeding metadata validation failed with ", nrow(issue_table),
         " issue(s). Re-run with strict = FALSE to inspect the ledger.",
         call. = FALSE)
  result
}

.hb_scalar_id <- function(x, name, allow_na = FALSE) {
  if (length(x) != 1L || (!allow_na && (is.na(x) || !nzchar(x))) ||
      (!is.na(x) && !nzchar(x)))
    stop(name, " must be one ", if (allow_na) "non-empty value or NA." else
           "non-empty value.", call. = FALSE)
  as.character(x)
}

.hb_pick_column <- function(data, candidates, default = NA) {
  found <- intersect(candidates, names(data))
  if (!length(found)) rep(default, nrow(data)) else data[[found[1L]]]
}

#' Build a BrAPI/MIAPPE-Oriented Recommendation Bundle
#'
#' Converts a validated HapBlockR decision object into stable tabular
#' recommendations, quality-control, uncertainty, exclusion, provenance, and
#' field-mapping tables. It does not embed remote-service calls.
#'
#' @param result A \code{hapblockr_result}.
#' @param breeding_program_id,study_id,trial_id Stable programme, study, and
#'   trial identifiers.
#' @param environment_id Optional environment identifier.
#' @param trait_id Optional trait or observation-variable identifier.
#' @param unit Optional trait unit.
#' @param ontology_id Optional ontology CURIE.
#' @param require_valid Logical. Refuse a failed decision contract.
#'
#' @return A list of canonical exchange tables with class
#'   \code{"HapBlockR_exchange_bundle"}.
#' @export
build_breeding_exchange <- function(
    result,
    breeding_program_id,
    study_id,
    trial_id,
    environment_id = NA_character_,
    trait_id = NA_character_,
    unit = NA_character_,
    ontology_id = NA_character_,
    require_valid = TRUE
) {
  if (!inherits(result, "hapblockr_result"))
    stop("result must inherit from 'hapblockr_result'.", call. = FALSE)
  breeding_program_id <- .hb_scalar_id(
    breeding_program_id, "breeding_program_id"
  )
  study_id <- .hb_scalar_id(study_id, "study_id")
  trial_id <- .hb_scalar_id(trial_id, "trial_id")
  environment_id <- .hb_scalar_id(
    environment_id, "environment_id", allow_na = TRUE
  )
  trait_id <- .hb_scalar_id(trait_id, "trait_id", allow_na = TRUE)
  unit <- .hb_scalar_id(unit, "unit", allow_na = TRUE)
  ontology_id <- .hb_scalar_id(
    ontology_id, "ontology_id", allow_na = TRUE
  )
  contract <- result$result_contract
  if (isTRUE(require_valid) &&
      !identical(contract$validation_status, "passed"))
    stop("A failed HapBlockR result cannot be exported as a recommendation.",
         call. = FALSE)
  decision <- as.data.frame(contract$decision_table,
                            stringsAsFactors = FALSE)
  if (!nrow(decision))
    stop("The result contract contains no decision rows.", call. = FALSE)

  female <- as.character(.hb_pick_column(
    decision, c("female", "parent1"), NA_character_
  ))
  male <- as.character(.hb_pick_column(
    decision, c("male", "parent2"), NA_character_
  ))
  id <- as.character(.hb_pick_column(
    decision, c("id", "candidate_id", "germplasm_id"), NA_character_
  ))
  is_cross <- !is.na(female) & nzchar(female) &
    !is.na(male) & nzchar(male)
  entity_type <- ifelse(is_cross, "cross", "germplasm")
  entity_id <- ifelse(
    is_cross,
    paste(female, male, sep = " x "),
    ifelse(!is.na(id) & nzchar(id), id,
           paste0("decision_", seq_len(nrow(decision))))
  )
  score_candidates <- names(decision)[
    vapply(decision, is.numeric, logical(1L))
  ]
  score_candidates <- setdiff(
    score_candidates,
    c("rank", "n", "n_crosses", "n_violations",
      "n_binding_constraints")
  )
  score <- if (length(score_candidates))
    as.numeric(decision[[score_candidates[1L]]]) else
      rep(NA_real_, nrow(decision))
  recommendable <- as.logical(.hb_pick_column(
    decision, c("recommendable", "selected", "feasible"),
    identical(contract$validation_status, "passed")
  ))
  recommendable[is.na(recommendable)] <- FALSE
  reason <- as.character(.hb_pick_column(
    decision,
    c("recommendation_reason", "recommendation", "reason", "reasons"),
    ""
  ))
  rank_value <- .hb_pick_column(decision, "rank", NA_integer_)
  if (all(is.na(rank_value)))
    rank_value <- rank(-score, ties.method = "min", na.last = "keep")
  recommendations <- data.frame(
    recommendation_id = sprintf("REC%06d", seq_len(nrow(decision))),
    breeding_program_id = breeding_program_id,
    study_id = study_id,
    trial_id = trial_id,
    environment_id = environment_id,
    trait_id = trait_id,
    ontology_id = ontology_id,
    unit = unit,
    method = contract$method,
    schema_version = contract$schema_version,
    entity_type = entity_type,
    entity_id = entity_id,
    germplasm_id = ifelse(is_cross, NA_character_, entity_id),
    female_germplasm_id = ifelse(is_cross, female, NA_character_),
    male_germplasm_id = ifelse(is_cross, male, NA_character_),
    rank = as.integer(rank_value),
    score = score,
    recommendable = recommendable,
    reason = reason,
    stringsAsFactors = FALSE
  )
  identifiers <- data.frame(
    identifier_type = c(
      rep("sample", length(contract$identifiers$sample_ids)),
      rep("variant", length(contract$identifiers$variant_ids))
    ),
    identifier = c(
      contract$identifiers$sample_ids,
      contract$identifiers$variant_ids
    ),
    stringsAsFactors = FALSE
  )
  provenance <- data.frame(
    field = c(
      "schema_version", "method", "created_at_utc", "call", "random_seed",
      "R_version", "HapBlockR_version"
    ),
    value = c(
      contract$schema_version,
      contract$method,
      contract$created_at_utc,
      contract$call,
      as.character(contract$random_seed),
      contract$software$R,
      contract$software$HapBlockR
    ),
    stringsAsFactors = FALSE
  )
  hashes <- data.frame(
    input = names(contract$input_hashes),
    sha256 = unname(contract$input_hashes),
    stringsAsFactors = FALSE
  )
  mapping <- data.frame(
    canonical_field = c(
      "breeding_program_id", "study_id", "trial_id", "environment_id",
      "trait_id", "germplasm_id", "female_germplasm_id",
      "male_germplasm_id", "score", "unit", "ontology_id"
    ),
    BrAPI_field = c(
      "breedingProgramDbId", "studyDbId", "trialDbId",
      "observationUnitDbId", "observationVariableDbId", "germplasmDbId",
      "parent1DbId/additionalInfo", "parent2DbId/additionalInfo", "value",
      "scale/units", "ontologyDbId"
    ),
    MIAPPE_field = c(
      "investigation identifier", "study identifier", "study identifier",
      "observation unit ID", "observed variable", "biological material ID",
      "crossing metadata", "crossing metadata", "observed value", "unit",
      "trait accession number"
    ),
    stringsAsFactors = FALSE
  )
  uncertainty_table <- as.data.frame(
    contract$uncertainty,
    stringsAsFactors = FALSE
  )
  if (!ncol(uncertainty_table))
    uncertainty_table <- data.frame(
      record_id = character(),
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
  exclusion_table <- as.data.frame(
    contract$excluded_records,
    stringsAsFactors = FALSE
  )
  if (!ncol(exclusion_table))
    exclusion_table <- data.frame(
      record_id = character(),
      reason = character(),
      stringsAsFactors = FALSE
    )
  bundle <- list(
    recommendations = recommendations,
    quality_gates = contract$quality_gates,
    uncertainty = uncertainty_table,
    exclusions = exclusion_table,
    identifiers = identifiers,
    input_hashes = hashes,
    provenance = provenance,
    field_mapping = mapping
  )
  class(bundle) <- c("HapBlockR_exchange_bundle", "list")
  attr(bundle, "bundle_schema_version") <- "1.0.0"
  bundle
}

#' Write a Breeding Exchange Bundle
#'
#' Writes every exchange table as UTF-8 CSV and creates a SHA-256 manifest.
#'
#' @param bundle Object from \code{\link{build_breeding_exchange}}.
#' @param path Target directory.
#' @param overwrite Logical. Replace known bundle files when they exist.
#'
#' @return Invisibly returns the manifest data frame.
#' @export
write_breeding_exchange <- function(bundle, path, overwrite = FALSE) {
  if (!inherits(bundle, "HapBlockR_exchange_bundle"))
    stop("bundle must come from build_breeding_exchange().", call. = FALSE)
  if (length(path) != 1L || is.na(path) || !nzchar(path))
    stop("path must be one non-empty directory path.", call. = FALSE)
  if (!dir.exists(path))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  files <- paste0(names(bundle), ".csv")
  targets <- file.path(path, files)
  manifest_path <- file.path(path, "manifest.csv")
  existing <- c(targets[file.exists(targets)],
                manifest_path[file.exists(manifest_path)])
  if (length(existing) && !isTRUE(overwrite))
    stop("Exchange files already exist. Set overwrite = TRUE to replace ",
         "the known bundle files.", call. = FALSE)
  if (length(existing))
    unlink(existing, force = TRUE)
  for (index in seq_along(bundle)) {
    utils::write.csv(bundle[[index]], targets[index], row.names = FALSE,
                     na = "__HAPBLOCKR_NA__")
  }
  column_type <- function(x) {
    if (is.integer(x)) "integer"
    else if (is.numeric(x)) "numeric"
    else if (is.logical(x)) "logical"
    else "character"
  }
  manifest <- data.frame(
    schema_version = attr(bundle, "bundle_schema_version"),
    file = files,
    rows = vapply(bundle, nrow, integer(1L)),
    column_types = vapply(
      bundle,
      function(x) paste(vapply(x, column_type, character(1L)),
                         collapse = ";"),
      character(1L)
    ),
    sha256 = vapply(
      targets,
      digest::digest,
      character(1L),
      algo = "sha256",
      file = TRUE
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, manifest_path, row.names = FALSE, na = "")
  invisible(manifest)
}

#' Read and Verify a Breeding Exchange Bundle
#'
#' @param path Directory containing \code{manifest.csv}.
#' @param verify Logical. Verify every file's SHA-256 checksum and row count.
#'
#' @return A \code{"HapBlockR_exchange_bundle"}.
#' @export
read_breeding_exchange <- function(path, verify = TRUE) {
  manifest_path <- file.path(path, "manifest.csv")
  if (!file.exists(manifest_path))
    stop("manifest.csv was not found.", call. = FALSE)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  required <- c(
    "schema_version", "file", "rows", "column_types", "sha256"
  )
  if (!all(required %in% names(manifest)) || !nrow(manifest))
    stop("The exchange manifest is malformed.", call. = FALSE)
  if (anyDuplicated(manifest$file) ||
      any(grepl("(^|[\\\\/])\\.\\.([\\\\/]|$)", manifest$file)))
    stop("The exchange manifest contains unsafe or duplicate paths.",
         call. = FALSE)
  targets <- file.path(path, manifest$file)
  if (any(!file.exists(targets)))
    stop("Exchange file(s) listed in the manifest are missing.",
         call. = FALSE)
  if (isTRUE(verify)) {
    observed_hash <- vapply(
      targets,
      digest::digest,
      character(1L),
      algo = "sha256",
      file = TRUE
    )
    if (!identical(
      unname(tolower(observed_hash)),
      unname(tolower(manifest$sha256))
    ))
      stop("Exchange bundle checksum verification failed.", call. = FALSE)
  }
  bundle <- lapply(seq_along(targets), function(index) {
    classes <- strsplit(
      manifest$column_types[index], ";", fixed = TRUE
    )[[1L]]
    utils::read.csv(
      targets[index],
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = "__HAPBLOCKR_NA__",
      colClasses = classes
    )
  })
  names(bundle) <- sub("\\.csv$", "", basename(manifest$file),
                       ignore.case = TRUE)
  if (isTRUE(verify)) {
    rows <- vapply(bundle, nrow, integer(1L))
    if (!identical(as.integer(rows), as.integer(manifest$rows)))
      stop("Exchange bundle row-count verification failed.",
           call. = FALSE)
  }
  class(bundle) <- c("HapBlockR_exchange_bundle", "list")
  attr(bundle, "bundle_schema_version") <- unique(manifest$schema_version)[1L]
  attr(bundle, "manifest") <- manifest
  bundle
}
