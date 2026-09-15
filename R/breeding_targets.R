# ==============================================================================
# Externally analysed breeding targets
# ==============================================================================

.hb_target_definitions <- function() {
  data.frame(
    input_type = c(
      "adjusted_mean", "BLUE", "BLUP_identity", "PBLUP",
      "BV", "GCA", "TGV"
    ),
    estimand = c(
      "model-adjusted entry mean",
      "model-adjusted entry mean",
      "random genotype or entry effect",
      "additive breeding value",
      "additive breeding value",
      "general combining ability",
      "total genetic value"
    ),
    definition = c(
      paste(
        "A model-adjusted entry mean, including a Bayesian posterior adjusted",
        "mean when the analysis does not use the BLUE/BLUP distinction; it is",
        "not a raw arithmetic mean."
      ),
      paste(
        "Best Linear Unbiased Estimate: the genotype or entry was fitted as",
        "a fixed effect, so the estimate is not shrunk towards a population",
        "mean."
      ),
      paste(
        "Best Linear Unbiased Prediction from a random genotype or entry",
        "effect with covariance I multiplied by the genetic variance; I",
        "describes the model covariance and does not assert biological",
        "unrelatedness."
      ),
      paste(
        "Pedigree Best Linear Unbiased Prediction of additive breeding value",
        "with covariance A multiplied by the additive genetic variance, where",
        "A is the numerator relationship matrix."
      ),
      paste(
        "Breeding value: an additive, transmissible genetic estimand rather",
        "than an estimation method. Its estimation basis must be declared;",
        "externally computed genomic breeding values are not accepted."
      ),
      paste(
        "General combining ability in a declared tester or mate population.",
        "The analysis must state whether GCA was fitted as a fixed effect or",
        "as a random effect."
      ),
      paste(
        "Total genetic value containing additive and non-additive components.",
        "Additive and dominance components must be supplied separately so that",
        "HapBlockR can preserve their distinct breeding interpretations."
      )
    ),
    stringsAsFactors = FALSE
  )
}

#' Describe Accepted Breeding-target Inputs
#'
#' Returns the definitions used by HapBlockR for externally analysed
#' genotype-level targets. HapBlockR accepts these summaries instead of raw
#' plot records and does not fit field-trial design factors.
#'
#' @return A data frame defining every accepted input type.
#' @export
breeding_target_types <- function() {
  .hb_target_definitions()
}

#' Prepare Externally Analysed Breeding Targets
#'
#' Validates genotype-level estimates produced outside HapBlockR, computes
#' precision weights, deregresses genuine random-effect predictions where
#' required, and orients every trait so that larger model values are
#' favourable. The returned `model_value` is the response supplied to
#' HapBlockR's marker, haplotype, block, multi-trait, or genotype-by-environment
#' models. External genomic BLUPs and external selection-index values are
#' rejected because genomic effects and selection indices are fitted within
#' HapBlockR.
#'
#' @param data One row per genotype, trait, and optional environment.
#' @param input_type One of the values returned by
#'   \code{breeding_target_types()}. A single type must apply to the call.
#' @param id_col,trait_col,value_col Column names.
#' @param environment_col Optional environment column. Leave `NULL` for an
#'   across-environment estimate.
#' @param unit_col Optional declared-unit column. Units are retained when
#'   supplied but are never mandatory.
#' @param se_col Standard error column for `adjusted_mean` or `BLUE`.
#' @param posterior_sd_col Posterior standard deviation column for a Bayesian
#'   `adjusted_mean`.
#' @param precision_col Externally supplied positive precision column.
#' @param covariance Optional full sampling covariance matrix. Its row and
#'   column names must equal the target record keys returned in `record_key`.
#' @param reliability_col,pev_col Alternative uncertainty inputs for random
#'   predictions. `pev_col` must contain prediction error variances, not
#'   standard errors of differences.
#' @param genetic_variance Additive or genotype variance, either one positive
#'   number, a named value per trait, or a named trait-by-environment value.
#' @param relationship_matrix Required numerator relationship matrix for
#'   `PBLUP`, and for pedigree-based `BV`, `GCA`, or `TGV`.
#' @param parent_average Optional named parent-average vector used in
#'   deregression. The default is zero.
#' @param estimation_basis Required for `BV`; and conditionally for `GCA` and
#'   `TGV`. Allowed values are `"fixed"`, `"identity"`, and `"pedigree"`.
#'   `"genomic"` is explicitly rejected.
#' @param gca_effect For `GCA`, either `"fixed"` or `"random"`.
#' @param tester_population For `GCA`, a non-empty description of the tester or
#'   mate population to which the estimates apply.
#' @param additive_col,dominance_col For `TGV`, columns containing separate
#'   additive and dominance components. Their sum is used as the supplied
#'   total genetic value.
#' @param lower_is_better Trait names for which smaller original values are
#'   favourable. Their signs are reversed in `model_value`.
#' @param heritability Optional provenance only. It is recorded but is not used
#'   to manufacture reliability or precision.
#' @param strict Logical. Stop if any reliability lies outside `(0, 1]`.
#'
#' @return A `hapblockr_result` containing `targets`, a `target_contract`,
#'   precision provenance, optional covariance, and quality gates.
#' @export
prepare_breeding_targets <- function(
    data,
    input_type,
    id_col = "id",
    trait_col = "trait",
    value_col = "value",
    environment_col = NULL,
    unit_col = NULL,
    se_col = NULL,
    posterior_sd_col = NULL,
    precision_col = NULL,
    covariance = NULL,
    reliability_col = NULL,
    pev_col = NULL,
    genetic_variance = NULL,
    relationship_matrix = NULL,
    parent_average = NULL,
    estimation_basis = NULL,
    gca_effect = NULL,
    tester_population = NULL,
    additive_col = NULL,
    dominance_col = NULL,
    lower_is_better = NULL,
    heritability = NULL,
    strict = TRUE
) {
  result_call <- match.call()
  definitions <- .hb_target_definitions()
  rejected <- c(
    "GBLUP", "genomic_BLUP", "GEBV", "selection_index",
    "index", "raw_phenotype", "plot"
  )
  if (length(input_type) != 1L || is.na(input_type)) {
    stop("input_type must be one non-missing value.", call. = FALSE)
  }
  if (input_type %in% rejected) {
    stop(
      input_type,
      " is not an accepted external target. HapBlockR fits genomic effects ",
      "and selection indices internally.",
      call. = FALSE
    )
  }
  if (!input_type %in% definitions$input_type) {
    stop(
      "input_type must be one of: ",
      paste(definitions$input_type, collapse = ", "),
      call. = FALSE
    )
  }
  records <- as.data.frame(data, stringsAsFactors = FALSE)
  required <- c(id_col, trait_col)
  if (input_type != "TGV") required <- c(required, value_col)
  if (!all(required %in% names(records))) {
    stop(
      "data is missing required columns: ",
      paste(setdiff(required, names(records)), collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(environment_col) && !environment_col %in% names(records)) {
    stop("environment_col is absent from data.", call. = FALSE)
  }
  if (!is.null(unit_col) && !unit_col %in% names(records)) {
    stop("unit_col is absent from data.", call. = FALSE)
  }

  targets <- data.frame(
    id = as.character(records[[id_col]]),
    trait = as.character(records[[trait_col]]),
    environment = if (is.null(environment_col)) {
      rep(NA_character_, nrow(records))
    } else {
      as.character(records[[environment_col]])
    },
    unit = if (is.null(unit_col)) {
      rep(NA_character_, nrow(records))
    } else {
      as.character(records[[unit_col]])
    },
    stringsAsFactors = FALSE
  )
  if (anyNA(targets$id) || any(!nzchar(targets$id)) ||
      anyNA(targets$trait) || any(!nzchar(targets$trait))) {
    stop("Every target requires a non-empty genotype ID and trait.",
         call. = FALSE)
  }

  if (input_type == "TGV") {
    if (is.null(additive_col) || is.null(dominance_col) ||
        !all(c(additive_col, dominance_col) %in% names(records))) {
      stop(
        "TGV requires additive_col and dominance_col with separate components.",
        call. = FALSE
      )
    }
    targets$additive_value <- .hb_numeric_column(
      records, additive_col, "additive_col"
    )
    targets$dominance_value <- .hb_numeric_column(
      records, dominance_col, "dominance_col"
    )
    targets$value <- targets$additive_value + targets$dominance_value
  } else {
    targets$value <- .hb_numeric_column(records, value_col, "value_col")
  }
  if (any(!is.finite(targets$value))) {
    stop("Target values must be finite; missing target values are not imputed.",
         call. = FALSE)
  }

  targets$record_key <- paste(
    targets$id,
    targets$trait,
    ifelse(is.na(targets$environment), "ACROSS", targets$environment),
    sep = "::"
  )
  if (anyDuplicated(targets$record_key)) {
    stop(
      "Each genotype-trait-environment target must occur once. Supply ",
      "externally adjusted summaries, not replicate or plot records.",
      call. = FALSE
    )
  }

  random_prediction <- input_type %in% c("BLUP_identity", "PBLUP")
  if (input_type == "BV") {
    estimation_basis <- .hb_validate_basis(estimation_basis, required = TRUE)
    random_prediction <- estimation_basis %in% c("identity", "pedigree")
  }
  if (input_type == "GCA") {
    if (length(tester_population) != 1L || is.na(tester_population) ||
        !nzchar(tester_population)) {
      stop("GCA requires a declared tester_population.", call. = FALSE)
    }
    if (is.null(gca_effect)) {
      stop(
        "GCA requires gca_effect = 'fixed' or 'random'.",
        call. = FALSE
      )
    }
    gca_effect <- match.arg(gca_effect, c("fixed", "random"))
    if (is.null(estimation_basis)) {
      estimation_basis <- if (gca_effect == "fixed") "fixed" else "identity"
    }
    estimation_basis <- .hb_validate_basis(estimation_basis, required = TRUE)
    if (gca_effect == "fixed" && estimation_basis != "fixed") {
      stop("Fixed GCA must use estimation_basis = 'fixed'.", call. = FALSE)
    }
    if (gca_effect == "random" && estimation_basis == "fixed") {
      stop("Random GCA must use identity or pedigree estimation.",
           call. = FALSE)
    }
    random_prediction <- gca_effect == "random"
  }
  if (input_type == "TGV") {
    estimation_basis <- .hb_validate_basis(estimation_basis, required = TRUE)
    random_prediction <- estimation_basis %in% c("identity", "pedigree")
  }
  if (input_type == "BLUP_identity") estimation_basis <- "identity"
  if (input_type == "PBLUP") estimation_basis <- "pedigree"
  if (input_type %in% c("adjusted_mean", "BLUE")) {
    estimation_basis <- if (input_type == "BLUE") "fixed" else
      "model_adjusted"
  }

  relationship_diagonal <- rep(1, nrow(targets))
  if (identical(estimation_basis, "pedigree")) {
    relationship_matrix <- .hb_validate_relationship(
      relationship_matrix, unique(targets$id)
    )
    relationship_diagonal <- diag(relationship_matrix)[
      match(targets$id, rownames(relationship_matrix))
    ]
  } else if (input_type == "PBLUP") {
    stop("PBLUP requires a named numerator relationship matrix A.",
         call. = FALSE)
  }

  uncertainty <- .hb_target_uncertainty(
    records = records,
    targets = targets,
    input_type = input_type,
    random_prediction = random_prediction,
    se_col = se_col,
    posterior_sd_col = posterior_sd_col,
    precision_col = precision_col,
    covariance = covariance,
    reliability_col = reliability_col,
    pev_col = pev_col,
    genetic_variance = genetic_variance,
    relationship_diagonal = relationship_diagonal,
    strict = strict
  )
  targets$reliability <- uncertainty$reliability
  targets$PEV <- uncertainty$PEV
  targets$precision_raw <- uncertainty$precision

  parent_average_value <- rep(0, nrow(targets))
  if (!is.null(parent_average)) {
    if (!is.numeric(parent_average) || is.null(names(parent_average))) {
      stop("parent_average must be a named numeric vector.", call. = FALSE)
    }
    parent_average_value <- as.numeric(parent_average[targets$id])
    if (any(!is.finite(parent_average_value))) {
      stop("parent_average is missing for one or more target IDs.",
           call. = FALSE)
    }
  }
  targets$deregressed_value <- targets$value
  if (random_prediction) {
    targets$deregressed_value <- parent_average_value +
      (targets$value - parent_average_value) / targets$reliability
  }

  precision_group <- .hb_target_precision_group(targets)
  targets$precision_weight <- stats::ave(
    targets$precision_raw,
    precision_group,
    FUN = function(x) x / mean(x)
  )
  direction <- ifelse(targets$trait %in% lower_is_better, -1, 1)
  targets$direction <- direction
  targets$model_value <- direction * targets$deregressed_value
  targets$input_type <- input_type
  targets$estimand <- definitions$estimand[
    match(input_type, definitions$input_type)
  ]
  targets$estimation_basis <- estimation_basis

  covariance_checked <- .hb_validate_target_covariance(
    covariance, targets$record_key
  )
  result <- list(
    targets = targets,
    target_contract = definitions[
      match(input_type, definitions$input_type), , drop = FALSE
    ],
    covariance = covariance_checked,
    precision_provenance = uncertainty$provenance,
    heritability_provenance = heritability,
    tester_population = tester_population,
    relationship_matrix = relationship_matrix
  )
  .add_hapblockr_contract(
    result = result,
    method = "prepare_breeding_targets",
    call = result_call,
    parameters = list(
      input_type = input_type,
      estimation_basis = estimation_basis,
      lower_is_better = lower_is_better,
      units_declared = !is.null(unit_col),
      heritability_used = FALSE
    ),
    sample_ids = unique(targets$id),
    inputs = list(externally_analysed_targets = data),
    transformations = c(
      if (random_prediction) "prediction deregression" else
        "estimate retained without deregression",
      "within-analysis precision normalisation",
      "favourable-direction orientation"
    ),
    quality_gates = c(
      unique_target_records = !anyDuplicated(targets$record_key),
      finite_model_values = all(is.finite(targets$model_value)),
      positive_precision = all(targets$precision_weight > 0),
      reliability_bounded = if (random_prediction)
        all(targets$reliability > 0 & targets$reliability <= 1) else TRUE
    ),
    decision_table = targets,
    uncertainty = targets[
      c("record_key", "PEV", "reliability", "precision_weight")
    ]
  )
}

.hb_numeric_column <- function(records, column, label) {
  value <- suppressWarnings(as.numeric(records[[column]]))
  introduced <- is.na(value) & !is.na(records[[column]])
  if (any(introduced)) {
    stop(label, " contains non-numeric values.", call. = FALSE)
  }
  value
}

.hb_validate_basis <- function(x, required = FALSE) {
  if (is.null(x)) {
    if (required) stop("estimation_basis must be declared.", call. = FALSE)
    return(NULL)
  }
  if (identical(x, "genomic")) {
    stop(
      "estimation_basis = 'genomic' is not accepted for external targets.",
      call. = FALSE
    )
  }
  match.arg(x, c("fixed", "identity", "pedigree"))
}

.hb_validate_relationship <- function(A, ids) {
  if (is.null(A)) {
    stop("A named numerator relationship matrix is required.",
         call. = FALSE)
  }
  A <- as.matrix(A)
  if (!is.numeric(A) || nrow(A) != ncol(A) ||
      is.null(rownames(A)) || is.null(colnames(A)) ||
      !identical(rownames(A), colnames(A)) ||
      !all(ids %in% rownames(A))) {
    stop(
      "relationship_matrix must be a square, consistently named A matrix ",
      "covering every target ID.",
      call. = FALSE
    )
  }
  A <- A[ids, ids, drop = FALSE]
  if (max(abs(A - t(A))) > 1e-8 || any(diag(A) <= 0)) {
    stop("relationship_matrix must be symmetric with a positive diagonal.",
         call. = FALSE)
  }
  A
}

.hb_expand_trait_value <- function(x, targets, name) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0)) {
    stop(name, " must contain positive finite values.", call. = FALSE)
  }
  if (length(x) == 1L) return(rep(as.numeric(x), nrow(targets)))
  group_key <- paste(
    targets$trait,
    ifelse(is.na(targets$environment), "ACROSS", targets$environment),
    sep = "::"
  )
  if (!is.null(names(x)) && all(group_key %in% names(x))) {
    return(as.numeric(x[group_key]))
  }
  if (!is.null(names(x)) && all(targets$trait %in% names(x))) {
    return(as.numeric(x[targets$trait]))
  }
  stop(
    name,
    " must be scalar or named for every trait or trait-environment.",
    call. = FALSE
  )
}

.hb_target_uncertainty <- function(
    records,
    targets,
    input_type,
    random_prediction,
    se_col,
    posterior_sd_col,
    precision_col,
    covariance,
    reliability_col,
    pev_col,
    genetic_variance,
    relationship_diagonal,
    strict
) {
  n <- nrow(targets)
  reliability <- rep(NA_real_, n)
  PEV <- rep(NA_real_, n)
  precision <- rep(NA_real_, n)
  supplied <- Filter(
    Negate(is.null),
    list(
      SE = se_col,
      posterior_SD = posterior_sd_col,
      precision = precision_col,
      reliability = reliability_col,
      PEV = pev_col
    )
  )
  if (length(supplied) > 1L) {
    stop(
      "Supply one primary uncertainty column per call; use covariance ",
      "alongside it only when the full sampling covariance is available.",
      call. = FALSE
    )
  }

  if (!random_prediction) {
    if (!is.null(se_col)) {
      SE <- .hb_positive_uncertainty(records, se_col, "standard errors")
      precision <- 1 / SE^2
      provenance <- "inverse sampling variance from SE"
    } else if (!is.null(posterior_sd_col)) {
      SD <- .hb_positive_uncertainty(
        records, posterior_sd_col, "posterior standard deviations"
      )
      precision <- 1 / SD^2
      provenance <- "inverse posterior variance"
    } else if (!is.null(precision_col)) {
      precision <- .hb_positive_uncertainty(
        records, precision_col, "precision values"
      )
      provenance <- "externally supplied precision"
    } else if (!is.null(covariance)) {
      covariance <- .hb_validate_target_covariance(
        covariance, targets$record_key
      )
      precision <- 1 / diag(covariance)
      provenance <- "inverse diagonal of full sampling covariance"
    } else {
      stop(
        input_type,
        " requires SE, posterior SD, precision, or a full sampling covariance.",
        call. = FALSE
      )
    }
    return(list(
      reliability = reliability,
      PEV = PEV,
      precision = precision,
      provenance = provenance
    ))
  }

  if (!is.null(reliability_col)) {
    reliability <- .hb_numeric_column(
      records, reliability_col, "reliability_col"
    )
    if (!is.null(genetic_variance)) {
      variance <- .hb_expand_trait_value(
        genetic_variance, targets, "genetic_variance"
      )
      marginal_variance <- variance * relationship_diagonal
      PEV <- (1 - reliability) * marginal_variance
    }
    provenance <- "externally supplied reliability"
  } else if (!is.null(pev_col)) {
    variance <- .hb_expand_trait_value(
      genetic_variance, targets, "genetic_variance"
    )
    if (is.null(variance)) {
      stop(
        "genetic_variance is required when reliability is calculated from PEV.",
        call. = FALSE
      )
    }
    marginal_variance <- variance * relationship_diagonal
    PEV <- .hb_positive_uncertainty(records, pev_col, "PEV values")
    reliability <- 1 - PEV / marginal_variance
    provenance <- if (all(relationship_diagonal == 1)) {
      "1 - PEV / genetic variance"
    } else {
      "1 - PEV / (A_ii * additive genetic variance)"
    }
  } else {
    stop(
      input_type,
      " requires reliability or PEV. Standard errors of differences are not ",
      "prediction error variances.",
      call. = FALSE
    )
  }
  invalid <- !is.finite(reliability) | reliability <= 0 | reliability > 1
  if (any(invalid)) {
    message <- paste(
      sum(invalid),
      "reliability value(s) fall outside (0, 1]."
    )
    if (isTRUE(strict)) stop(message, call. = FALSE) else warning(message)
  }
  reliability <- pmin(1, pmax(reliability, .Machine$double.eps))
  precision <- reliability / pmax(1 - reliability, .Machine$double.eps)
  list(
    reliability = reliability,
    PEV = PEV,
    precision = precision,
    provenance = provenance
  )
}

.hb_positive_uncertainty <- function(records, column, label) {
  if (!column %in% names(records)) {
    stop(column, " is absent from data.", call. = FALSE)
  }
  value <- .hb_numeric_column(records, column, label)
  if (any(!is.finite(value)) || any(value <= 0)) {
    stop(label, " must be positive and finite.", call. = FALSE)
  }
  value
}

.hb_validate_target_covariance <- function(covariance, keys) {
  if (is.null(covariance)) return(NULL)
  covariance <- as.matrix(covariance)
  if (!is.numeric(covariance) ||
      any(dim(covariance) != length(keys)) ||
      is.null(rownames(covariance)) ||
      is.null(colnames(covariance)) ||
      !setequal(rownames(covariance), keys) ||
      !setequal(colnames(covariance), keys)) {
    stop(
      "covariance must be a named square matrix covering every record_key.",
      call. = FALSE
    )
  }
  covariance <- covariance[keys, keys, drop = FALSE]
  if (max(abs(covariance - t(covariance))) > 1e-8 ||
      any(diag(covariance) <= 0)) {
    stop("covariance must be symmetric with a positive diagonal.",
         call. = FALSE)
  }
  eigenvalues <- eigen(
    (covariance + t(covariance)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  tolerance <- 1e-8 * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -tolerance) {
    stop("covariance must be positive semidefinite.", call. = FALSE)
  }
  covariance
}

.hb_target_precision_group <- function(targets) {
  required <- c("trait", "environment")
  if (!all(required %in% names(targets))) {
    stop(
      "Prepared targets require trait and environment columns for precision ",
      "normalisation.",
      call. = FALSE
    )
  }
  interaction(
    targets$trait,
    ifelse(is.na(targets$environment), "ACROSS", targets$environment),
    drop = TRUE
  )
}

.hb_unpack_model_targets <- function(x, allow_environment = FALSE) {
  if (!inherits(x, "hapblockr_result") || is.null(x$targets) ||
      is.null(x$target_contract)) {
    return(NULL)
  }
  targets <- x$targets
  required <- c(
    "id", "trait", "environment", "record_key", "model_value",
    "precision_weight"
  )
  if (!all(required %in% names(targets))) {
    stop(
      "The breeding-target object is incomplete; recreate it with ",
      "prepare_breeding_targets().",
      call. = FALSE
    )
  }
  environments <- unique(stats::na.omit(targets$environment))
  if (!isTRUE(allow_environment) && length(environments) > 1L) {
    stop(
      "Targets from multiple environments require fit_gxe_gblup(). To fit ",
      "one environment independently, prepare or subset one environment.",
      call. = FALSE
    )
  }
  if (!isTRUE(allow_environment) && length(environments) == 1L &&
      anyNA(targets$environment)) {
    stop(
      "Do not mix across-environment and environment-specific targets in ",
      "one multivariate model.",
      call. = FALSE
    )
  }
  split_rows <- split(seq_len(nrow(targets)), targets$trait)
  values <- lapply(split_rows, function(index) {
    value <- targets$model_value[index]
    names(value) <- targets$id[index]
    value
  })
  weights <- lapply(split_rows, function(index) {
    value <- targets$precision_weight[index]
    names(value) <- targets$id[index]
    value
  })
  list(
    values = values,
    weights = weights,
    covariance = x$covariance,
    contract = x$target_contract,
    precision_provenance = x$precision_provenance
  )
}
