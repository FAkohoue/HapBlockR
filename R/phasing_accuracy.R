# ==============================================================================
# Truth-set assessment for phased diploid genotypes
# ==============================================================================

.as_phased_truth_set <- function(x, name) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x))
      stop(name, " VCF does not exist: ", x, call. = FALSE)
    x <- read_phased_vcf(x, min_maf = 0, verbose = FALSE)
  }
  if (!is.list(x) || !all(c("hap1", "hap2") %in% names(x)))
    stop(name, " must be a phased VCF path or a list containing hap1 and hap2.",
         call. = FALSE)
  hap1 <- as.matrix(x$hap1)
  hap2 <- as.matrix(x$hap2)
  if (!identical(dim(hap1), dim(hap2)))
    stop(name, " hap1 and hap2 must have identical dimensions.",
         call. = FALSE)
  if (!nrow(hap1) || !ncol(hap1))
    stop(name, " must contain at least one variant and one sample.",
         call. = FALSE)
  if (is.null(rownames(hap1)) || is.null(colnames(hap1)) ||
      is.null(rownames(hap2)) || is.null(colnames(hap2)))
    stop(name, " haplotype matrices must have variant and sample dimnames.",
         call. = FALSE)
  if (!identical(dimnames(hap1), dimnames(hap2)))
    stop(name, " hap1 and hap2 dimnames must be identical.", call. = FALSE)
  if (anyDuplicated(rownames(hap1)) || anyDuplicated(colnames(hap1)))
    stop(name, " variant and sample identifiers must be unique.",
         call. = FALSE)
  values <- c(hap1, hap2)
  if (any(!is.na(values) & !values %in% c(0, 1)))
    stop(name, " must contain only biallelic haploid calls 0, 1, or NA.",
         call. = FALSE)

  snp_info <- x$snp_info
  if (is.null(snp_info)) {
    snp_info <- data.frame(
      SNP = rownames(hap1),
      CHR = "1",
      POS = seq_len(nrow(hap1)),
      REF = NA_character_,
      ALT = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  snp_info <- as.data.frame(snp_info, stringsAsFactors = FALSE)
  required <- c("SNP", "CHR", "POS")
  if (!all(required %in% names(snp_info)) ||
      nrow(snp_info) != nrow(hap1))
    stop(name, " snp_info must contain one SNP/CHR/POS row per variant.",
         call. = FALSE)
  if (!identical(as.character(snp_info$SNP), rownames(hap1)))
    stop(name, " snp_info$SNP must match haplotype-matrix row names.",
         call. = FALSE)
  if (!"REF" %in% names(snp_info)) snp_info$REF <- NA_character_
  if (!"ALT" %in% names(snp_info)) snp_info$ALT <- NA_character_

  list(
    hap1 = matrix(as.numeric(hap1), nrow = nrow(hap1),
                  dimnames = dimnames(hap1)),
    hap2 = matrix(as.numeric(hap2), nrow = nrow(hap2),
                  dimnames = dimnames(hap2)),
    dosage = hap1 + hap2,
    snp_info = snp_info,
    sample_ids = colnames(hap1)
  )
}

.validate_accuracy_threshold <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 0 || x > 1)
    stop(name, " must be one finite number in [0, 1].", call. = FALSE)
  as.numeric(x)
}

.phase_accuracy_stratum <- function(truth_h1, truth_h2, estimate_h1,
                                     estimate_h2) {
  truth_dosage <- truth_h1 + truth_h2
  estimate_dosage <- estimate_h1 + estimate_h2
  truth_observed <- is.finite(truth_dosage)
  called <- truth_observed & is.finite(estimate_dosage)
  dosage_correct <- called & truth_dosage == estimate_dosage

  allele_called <- is.finite(truth_h1) & is.finite(truth_h2) &
    is.finite(estimate_h1) & is.finite(estimate_h2)
  no_swap_matches <- sum(
    (estimate_h1[allele_called] == truth_h1[allele_called]) +
      (estimate_h2[allele_called] == truth_h2[allele_called])
  )
  swap_matches <- sum(
    (estimate_h1[allele_called] == truth_h2[allele_called]) +
      (estimate_h2[allele_called] == truth_h1[allele_called])
  )
  allele_comparisons <- 2L * sum(allele_called)

  informative <- called & truth_dosage == 1 & estimate_dosage == 1
  orientation <- estimate_h1[informative] != truth_h1[informative]
  transitions <- max(length(orientation) - 1L, 0L)
  switches <- if (transitions) sum(orientation[-1L] !=
                                     orientation[-length(orientation)]) else 0L

  list(
    truth_observed = sum(truth_observed),
    called = sum(called),
    dosage_correct = sum(dosage_correct),
    allele_matches = max(no_swap_matches, swap_matches),
    allele_comparisons = allele_comparisons,
    informative_heterozygotes = sum(informative),
    switches = switches,
    transitions = transitions
  )
}

#' Assess Phasing Accuracy against a Truth Set
#'
#' Compares estimated phased, diploid, biallelic genotypes with a phased truth
#' set. Variant, sample, chromosome, position, and available allele identities
#' are checked before calculation. Haplotype labels may be globally swapped
#' within each sample and chromosome. A switch error is counted when the
#' estimated orientation relative to truth changes between consecutive
#' informative heterozygous variants on the same chromosome.
#'
#' @param truth A phased VCF path or a list containing \code{hap1},
#'   \code{hap2}, and preferably \code{snp_info}, as returned by
#'   \code{\link{read_phased_vcf}}.
#' @param estimate Estimated phased data in the same form as \code{truth}.
#' @param max_switch_error_rate Maximum acceptable switch-error rate. Default
#'   \code{1}.
#' @param min_dosage_accuracy Minimum acceptable exact dosage accuracy among
#'   truth genotypes for which an estimated call is available. Default
#'   \code{0}.
#' @param min_allele_concordance Minimum acceptable phased-allele concordance
#'   after the best global haplotype-label orientation is selected within each
#'   sample and chromosome. Default \code{0}.
#' @param min_call_rate Minimum acceptable estimated call rate among observed
#'   truth genotypes. Default \code{0}.
#' @param require_switch_information Logical. Require at least one transition
#'   between informative heterozygous variants. Default \code{TRUE}.
#' @param strict Logical. Stop when a quality gate fails. Default \code{TRUE}.
#'
#' @return A \code{"HapBlockR_phasing_accuracy"} and
#'   \code{"hapblockr_result"} object containing overall metrics, per-sample
#'   metrics, per-chromosome metrics, quality gates, input hashes, and
#'   immutable identifiers.
#'
#' @examples
#' truth <- list(
#'   hap1 = matrix(c(0, 0, 1), ncol = 1,
#'                 dimnames = list(paste0("s", 1:3), "P1")),
#'   hap2 = matrix(c(1, 1, 0), ncol = 1,
#'                 dimnames = list(paste0("s", 1:3), "P1")),
#'   snp_info = data.frame(
#'     SNP = paste0("s", 1:3), CHR = "1", POS = 1:3,
#'     REF = "A", ALT = "G"
#'   )
#' )
#' accuracy <- assess_phasing_accuracy(truth, truth)
#' accuracy$metrics$switch_error_rate
#' validate(accuracy)
#'
#' @export
assess_phasing_accuracy <- function(
    truth,
    estimate,
    max_switch_error_rate = 1,
    min_dosage_accuracy = 0,
    min_allele_concordance = 0,
    min_call_rate = 0,
    require_switch_information = TRUE,
    strict = TRUE
) {
  call <- match.call()
  max_switch_error_rate <- .validate_accuracy_threshold(
    max_switch_error_rate, "max_switch_error_rate"
  )
  min_dosage_accuracy <- .validate_accuracy_threshold(
    min_dosage_accuracy, "min_dosage_accuracy"
  )
  min_allele_concordance <- .validate_accuracy_threshold(
    min_allele_concordance, "min_allele_concordance"
  )
  min_call_rate <- .validate_accuracy_threshold(
    min_call_rate, "min_call_rate"
  )
  if (length(require_switch_information) != 1L ||
      is.na(require_switch_information))
    stop("require_switch_information must be TRUE or FALSE.", call. = FALSE)
  if (length(strict) != 1L || is.na(strict))
    stop("strict must be TRUE or FALSE.", call. = FALSE)

  truth_data <- .as_phased_truth_set(truth, "truth")
  estimate_data <- .as_phased_truth_set(estimate, "estimate")
  truth_samples <- truth_data$sample_ids
  truth_variants <- rownames(truth_data$hap1)
  if (!setequal(truth_samples, estimate_data$sample_ids))
    stop("Truth and estimate sample identities differ.", call. = FALSE)
  if (!setequal(truth_variants, rownames(estimate_data$hap1)))
    stop("Truth and estimate variant identities differ.", call. = FALSE)

  sample_order_identity <- identical(
    truth_samples, estimate_data$sample_ids
  )
  variant_order_identity <- identical(
    truth_variants, rownames(estimate_data$hap1)
  )
  sample_index <- match(truth_samples, estimate_data$sample_ids)
  variant_index <- match(truth_variants, rownames(estimate_data$hap1))
  estimate_data$hap1 <- estimate_data$hap1[
    variant_index, sample_index, drop = FALSE
  ]
  estimate_data$hap2 <- estimate_data$hap2[
    variant_index, sample_index, drop = FALSE
  ]
  estimate_data$dosage <- estimate_data$hap1 + estimate_data$hap2
  estimate_info <- estimate_data$snp_info[variant_index, , drop = FALSE]
  truth_info <- truth_data$snp_info
  coordinate_identity <- identical(
    as.character(truth_info$CHR),
    as.character(estimate_info$CHR)
  ) && identical(as.numeric(truth_info$POS), as.numeric(estimate_info$POS))
  known_alleles <- !is.na(truth_info$REF) & !is.na(truth_info$ALT) &
    !is.na(estimate_info$REF) & !is.na(estimate_info$ALT)
  allele_identity <- !any(known_alleles) || all(
    as.character(truth_info$REF[known_alleles]) ==
      as.character(estimate_info$REF[known_alleles]) &
      as.character(truth_info$ALT[known_alleles]) ==
      as.character(estimate_info$ALT[known_alleles])
  )
  if (!coordinate_identity)
    stop("Truth and estimate chromosome or position identities differ.",
         call. = FALSE)
  if (!allele_identity)
    stop("Truth and estimate REF/ALT allele identities differ.",
         call. = FALSE)

  chromosomes <- unique(as.character(truth_info$CHR))
  rows <- vector("list", length(truth_samples) * length(chromosomes))
  row_index <- 0L
  for (sample_id in truth_samples) {
    sample_column <- match(sample_id, truth_samples)
    for (chromosome in chromosomes) {
      variant_rows <- which(as.character(truth_info$CHR) == chromosome)
      # Switch-error counting in .phase_accuracy_stratum() assumes consecutive
      # array positions are consecutive genomic positions (it compares
      # orientation[-1] against orientation[-length(orientation)] directly).
      # variant_rows is taken in truth_info's stored row order, which is not
      # guaranteed to be position-sorted, so sort explicitly here -- mirrors
      # .block_snp_order()'s explicit sort elsewhere in the package.
      variant_rows <- variant_rows[order(truth_info$POS[variant_rows])]
      stats <- .phase_accuracy_stratum(
        truth_data$hap1[variant_rows, sample_column],
        truth_data$hap2[variant_rows, sample_column],
        estimate_data$hap1[variant_rows, sample_column],
        estimate_data$hap2[variant_rows, sample_column]
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        sample_id = sample_id,
        chromosome = chromosome,
        truth_observed = stats$truth_observed,
        called = stats$called,
        dosage_correct = stats$dosage_correct,
        allele_matches = stats$allele_matches,
        allele_comparisons = stats$allele_comparisons,
        informative_heterozygotes = stats$informative_heterozygotes,
        switches = stats$switches,
        transitions = stats$transitions,
        stringsAsFactors = FALSE
      )
    }
  }
  strata <- do.call(rbind, rows)
  ratio <- function(numerator, denominator) {
    if (sum(denominator) > 0) sum(numerator) / sum(denominator) else NA_real_
  }
  per_sample <- do.call(rbind, lapply(
    split(strata, strata$sample_id),
    function(x) data.frame(
      sample_id = x$sample_id[1L],
      call_rate = ratio(x$called, x$truth_observed),
      dosage_accuracy = ratio(x$dosage_correct, x$called),
      allele_concordance = ratio(x$allele_matches, x$allele_comparisons),
      informative_heterozygotes = sum(x$informative_heterozygotes),
      switches = sum(x$switches),
      transitions = sum(x$transitions),
      switch_error_rate = ratio(x$switches, x$transitions),
      stringsAsFactors = FALSE
    )
  ))
  rownames(per_sample) <- NULL
  per_chromosome <- do.call(rbind, lapply(
    split(strata, strata$chromosome),
    function(x) data.frame(
      chromosome = x$chromosome[1L],
      call_rate = ratio(x$called, x$truth_observed),
      dosage_accuracy = ratio(x$dosage_correct, x$called),
      allele_concordance = ratio(x$allele_matches, x$allele_comparisons),
      informative_heterozygotes = sum(x$informative_heterozygotes),
      switches = sum(x$switches),
      transitions = sum(x$transitions),
      switch_error_rate = ratio(x$switches, x$transitions),
      stringsAsFactors = FALSE
    )
  ))
  rownames(per_chromosome) <- NULL

  metrics <- list(
    sample_order_identity = sample_order_identity,
    variant_order_identity = variant_order_identity,
    n_samples = length(truth_samples),
    n_variants = length(truth_variants),
    call_rate = ratio(strata$called, strata$truth_observed),
    dosage_accuracy = ratio(strata$dosage_correct, strata$called),
    allele_concordance = ratio(
      strata$allele_matches, strata$allele_comparisons
    ),
    informative_heterozygotes = sum(strata$informative_heterozygotes),
    switches = sum(strata$switches),
    transitions = sum(strata$transitions),
    switch_error_rate = ratio(strata$switches, strata$transitions)
  )
  gates <- data.frame(
    gate = c(
      "switch_information", "switch_error_rate", "dosage_accuracy",
      "allele_concordance", "call_rate"
    ),
    passed = c(
      !isTRUE(require_switch_information) || metrics$transitions > 0L,
      if (!metrics$transitions) {
        !isTRUE(require_switch_information)
      } else {
        is.finite(metrics$switch_error_rate) &&
          metrics$switch_error_rate <= max_switch_error_rate
      },
      is.finite(metrics$dosage_accuracy) &&
        metrics$dosage_accuracy >= min_dosage_accuracy,
      is.finite(metrics$allele_concordance) &&
        metrics$allele_concordance >= min_allele_concordance,
      is.finite(metrics$call_rate) && metrics$call_rate >= min_call_rate
    ),
    detail = c(
      paste(metrics$transitions, "informative transition(s)"),
      paste("observed", signif(metrics$switch_error_rate, 6),
            "<=", max_switch_error_rate),
      paste("observed", signif(metrics$dosage_accuracy, 6),
            ">=", min_dosage_accuracy),
      paste("observed", signif(metrics$allele_concordance, 6),
            ">=", min_allele_concordance),
      paste("observed", signif(metrics$call_rate, 6), ">=", min_call_rate)
    ),
    stringsAsFactors = FALSE
  )
  result <- .add_hapblockr_contract(
    result = list(
      metrics = metrics,
      per_sample = per_sample,
      per_chromosome = per_chromosome,
      quality_control = gates
    ),
    method = "assess_phasing_accuracy",
    call = call,
    parameters = list(
      max_switch_error_rate = max_switch_error_rate,
      min_dosage_accuracy = min_dosage_accuracy,
      min_allele_concordance = min_allele_concordance,
      min_call_rate = min_call_rate,
      require_switch_information = require_switch_information
    ),
    sample_ids = truth_samples,
    variant_ids = truth_variants,
    inputs = list(truth = truth_data, estimate = estimate_data),
    transformations = c(
      "aligned estimate to truth by immutable identifiers",
      "optimised only the global haplotype-label orientation within strata"
    ),
    quality_gates = gates,
    decision_table = per_sample
  )
  class(result) <- c(
    "HapBlockR_phasing_accuracy",
    setdiff(class(result), "HapBlockR_phasing_accuracy")
  )
  if (isTRUE(strict) && any(!gates$passed))
    stop(
      "Phasing truth-set quality gate(s) failed: ",
      paste(gates$gate[!gates$passed], collapse = ", "),
      ". Re-run with strict = FALSE to inspect the complete report.",
      call. = FALSE
    )
  result
}
