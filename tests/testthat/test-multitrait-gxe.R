make_test_kernel <- function(n = 16L, p = 30L, seed = 31L) {
  set.seed(seed)
  M <- matrix(rnorm(n * p), n, p)
  K <- tcrossprod(scale(M, center = TRUE, scale = FALSE)) / p
  K <- K + diag(0.25, n)
  ids <- sprintf("G%02d", seq_len(n))
  dimnames(K) <- list(ids, ids)
  K
}

test_that("classical indices accept optional units and check directions", {
  values <- cbind(
    yield = c(5.0, 4.5, 5.4, 4.9),
    disease = c(3.0, 1.5, 2.5, 1.0)
  )
  rownames(values) <- paste0("P", 1:4)
  G <- matrix(c(1.2, -0.2, -0.2, 0.8), 2, 2,
              dimnames = list(colnames(values), colnames(values)))
  P <- matrix(c(1.8, -0.1, -0.1, 1.3), 2, 2,
              dimnames = list(colnames(values), colnames(values)))

  economic <- build_selection_index(
    values, G, P,
    economic_weights = c(yield = 2, disease = 1),
    directions = c(yield = "increase", disease = "decrease"),
    units = c(yield = "t/ha", disease = "score")
  )
  expect_s3_class(economic, "hapblockr_result")
  expect_equal(nrow(economic$scores), 4L)
  expect_true(all(is.finite(economic$scores$selection_index)))
  expect_true(all(validate(economic)$passed))
  expect_true(is.finite(economic$index_accuracy))
  expect_identical(
    names(economic$coefficients),
    c(
      "trait", "unit", "direction", "coefficient", "linear_weight",
      "expected_response", "realised_response_sd_favourable",
      "realised_response_sd_original"
    )
  )

  desired <- build_selection_index(
    values, G, P,
    desired_gains = c(yield = 0.5, disease = 0.2),
    directions = c(yield = "increase", disease = "decrease"),
    units = NULL
  )
  expect_identical(desired$index_type, "pesek_baker")
  expect_true(all(is.na(desired$coefficients$unit)))
  expect_true(all(desired$coefficients$expected_response *
                    c(1, -1) > 0))

  expect_error(
    build_selection_index(
      values, G, P,
      economic_weights = c(yield = 2, disease = 1),
      directions = c(yield = "increase", disease = "sideways"),
      units = c(yield = "t/ha", disease = "score")
    ),
    "increase.*decrease"
  )

  expect_error(
    build_selection_index(
      values, G, P,
      economic_weights = c(yield = 2, disease = -1),
      directions = c(yield = "increase", disease = "decrease")
    ),
    "non-negative favourable-direction magnitudes"
  )
  expect_error(
    build_selection_index(
      values, G, P,
      desired_gains = c(yield = 0.5, disease = -0.2),
      directions = c(yield = "increase", disease = "decrease")
    ),
    "non-negative favourable-direction magnitudes"
  )
  expect_error(
    build_selection_index(
      values, G, P,
      economic_weights = c(yield = 2, disease = 1),
      directions = c(yield = "increase", disease = "decrease"),
      selection_intensity = factor("1")
    ),
    "positive finite real number"
  )
})

test_that("DesiredGainR engines return an internal merit score", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  values <- cbind(
    yield = seq(3.5, 6.5, length.out = 20),
    disease = seq(4, 1, length.out = 20) + sin(seq_len(20))
  )
  rownames(values) <- paste0("P", seq_len(20))
  traits <- colnames(values)
  G <- stats::cov(values)
  P <- G + diag(0.2, 2)
  dimnames(G) <- dimnames(P) <- list(traits, traits)

  dgsi <- build_selection_index(
    values, G, P,
    desired_gains = c(yield = 0.5, disease = 0.3),
    directions = c(yield = "increase", disease = "decrease"),
    method = "dgsi",
    n_select = 5,
    dgsi_control = list(n_iter = 10, n_rep = 3, seed = 11)
  )
  expect_identical(dgsi$index_type, "dgsi")
  expect_equal(dgsi$scores$merit_score, dgsi$scores$selection_index)
  expect_equal(sum(dgsi$engine_result$replicate_diagnostics$Chosen), 1)
  expect_true(all(is.finite(dgsi$coefficients$expected_response)))
  direction_sign <- c(yield = 1, disease = -1)
  D <- diag(direction_sign)
  G_oriented <- D %*% G %*% D
  P_oriented <- D %*% P %*% D
  b <- dgsi$engine_result$coefficients[traits]
  intensity <- stats::dnorm(stats::qnorm(0.75)) / 0.25
  expected_dgsi <- as.numeric(
    intensity * G_oriented %*% b /
      sqrt(as.numeric(crossprod(b, P_oriented %*% b)))
  ) * direction_sign
  expect_equal(
    dgsi$coefficients$expected_response,
    unname(expected_dgsi)
  )
  expect_lt(
    dgsi$coefficients$realised_response_sd_original[
      dgsi$coefficients$trait == "disease"
    ],
    0
  )
  expect_setequal(
    unique(dgsi$response_summary$response_kind),
    c(
      "model_expected_genetic_response",
      "realised_selected_set_differential"
    )
  )

  W <- diag(c(yield = 0.1, disease = 0.05))
  dimnames(W) <- list(traits, traits)
  qgsi <- build_selection_index(
    values, G, P,
    economic_weights = c(yield = 1, disease = 0.5),
    directions = c(yield = "increase", disease = "decrease"),
    method = "qgsi",
    n_select = 5,
    quadratic_weights = W
  )
  expect_identical(qgsi$index_type, "qgsi")
  expect_match(qgsi$contribution_scope, "candidate-specific")
  expect_equal(sum(qgsi$engine_result$ranked_geno$Selected), 5)
  expect_true(all(is.na(qgsi$coefficients$coefficient)))
  expect_equal(qgsi$coefficients$linear_weight, c(1, 0.5))
  expect_true(all(is.finite(qgsi$coefficients$expected_response)))
  expect_lt(
    qgsi$coefficients$expected_response[
      qgsi$coefficients$trait == "disease"
    ],
    0
  )
  engine_gain <- qgsi$engine_result$expected_gain_per_trait
  expected_favourable <- qgsi$coefficients$expected_response * c(1, -1)
  expect_equal(
    engine_gain$Expected_Genetic_Gain[
      match(qgsi$coefficients$trait, engine_gain$Trait)
    ],
    expected_favourable
  )
  expect_true(
    "Expected_Genetic_Gain_LinearSD" %in% names(engine_gain)
  )
  expect_match(
    qgsi$engine_result$expected_gain_basis,
    "total index standard deviation"
  )
  expect_null(qgsi$engine_result$hapblockr_response_reconciliation)
  expect_equal(
    qgsi$selection_intensity,
    qgsi$engine_result$theoretical_parameters$selection_intensity
  )
  expected_coefficient_columns <- c(
    "trait", "unit", "direction", "coefficient", "linear_weight",
    "expected_response", "realised_response_sd_favourable",
    "realised_response_sd_original"
  )
  expect_identical(names(dgsi$coefficients), expected_coefficient_columns)
  expect_identical(names(qgsi$coefficients), expected_coefficient_columns)

  expect_error(
    build_selection_index(
      values, G, P,
      desired_gains = c(yield = 0.5, disease = 0.3),
      directions = c(yield = "increase", disease = "decrease"),
      method = "dgsi",
      n_select = 5,
      selection_intensity = 1
    ),
    "applies only"
  )

  expect_error(
    build_selection_index(
      values, G, P,
      economic_weights = c(yield = 1, disease = 0.5),
      directions = c(yield = "increase", disease = "decrease"),
      method = "smith_hazel",
      n_select = 5
    ),
    "n_select applies only"
  )
})

test_that("DesiredGainR responses retain their units through engine controls", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  set.seed(104)
  values <- cbind(
    yield = rnorm(50, mean = 50, sd = 10),
    disease = rnorm(50, mean = 3, sd = 0.4)
  )
  rownames(values) <- sprintf("P%02d", seq_len(nrow(values)))
  traits <- colnames(values)
  G <- matrix(
    c(80, 1, 1, 0.12), 2, 2,
    dimnames = list(traits, traits)
  )
  P <- matrix(
    c(120, 1.5, 1.5, 0.20), 2, 2,
    dimnames = list(traits, traits)
  )
  directions <- c(yield = "increase", disease = "decrease")

  dgsi <- build_selection_index(
    values, G, P,
    desired_gains = c(yield = 0.5, disease = 0.3),
    directions = directions,
    method = "dgsi",
    n_select = 5,
    dgsi_control = list(
      scale_traits = TRUE, n_iter = 20, n_rep = 2, seed = 104
    )
  )
  expect_equal(
    setNames(dgsi$coefficients$expected_response, traits),
    dgsi$engine_result$theoretical_response$original_units[traits]
  )
  expect_equal(
    dgsi$selection_intensity,
    dgsi$engine_result$theoretical_response$selection_intensity
  )
  W <- diag(c(yield = 0.05, disease = -0.03))
  dimnames(W) <- list(traits, traits)
  qgsi <- build_selection_index(
    values, G, P,
    economic_weights = c(yield = 1, disease = 0.5),
    directions = directions,
    method = "qgsi",
    n_select = 5,
    quadratic_weights = W,
    qgsi_control = list(scale_traits = TRUE)
  )
  engine <- qgsi$engine_result
  gain <- as.data.frame(engine$expected_gain_per_trait)
  expected <- gain$Expected_Genetic_Gain[match(traits, gain$Trait)] *
    engine$transformation$scale[traits] *
    engine$transformation$direction[traits]
  expect_equal(qgsi$coefficients$expected_response, unname(expected))
  expect_equal(
    qgsi$selection_intensity,
    engine$theoretical_parameters$selection_intensity
  )
  expect_error(
    build_selection_index(
      values, G, P,
      economic_weights = c(yield = 1, disease = 0.5),
      directions = directions,
      method = "qgsi",
      n_select = 5,
      quadratic_weights = W,
      qgsi_control = list(n_select = 3)
    ),
    "cannot override.*n_select"
  )
  expect_error(
    build_selection_index(
      values, G, P,
      economic_weights = c(yield = 1, disease = 0.5),
      directions = directions,
      method = "qgsi",
      n_select = 5,
      quadratic_weights = W,
      qgsi_control = list(TRUE)
    ),
    "unique, non-empty argument names"
  )
})

test_that("multivariate GBLUP returns covariance diagnostics, PEV, and reliability", {
  K <- make_test_kernel()
  ids <- rownames(K)
  traits <- c("yield", "disease")
  Sg <- matrix(c(1.0, -0.35, -0.35, 0.7), 2, 2,
               dimnames = list(traits, traits))
  Se <- matrix(c(0.8, 0.10, 0.10, 0.6), 2, 2,
               dimnames = list(traits, traits))
  set.seed(32)
  u <- matrix(
    t(chol(kronecker(Sg, K))) %*% rnorm(length(ids) * 2L),
    nrow = length(ids), ncol = 2L
  )
  e <- matrix(
    t(chol(kronecker(Se, diag(length(ids))))) %*%
      rnorm(length(ids) * 2L),
    nrow = length(ids), ncol = 2L
  )
  Y <- u + e + matrix(rep(c(5, 2), each = length(ids)),
                      nrow = length(ids))
  dimnames(Y) <- list(ids, traits)
  Y[1, 2] <- NA_real_

  fit <- fit_multitrait_gblup(
    Y, K,
    genetic_cov = Sg,
    residual_cov = Se,
    estimate_covariances = FALSE,
    min_reliability = 0.20
  )
  expect_s3_class(fit, "hapblockr_result")
  expect_equal(nrow(fit$predictions), length(ids) * 2L)
  expect_equal(fit$genetic_cov, Sg, tolerance = 1e-8)
  expect_true(all(fit$predictions$PEV >= 0))
  expect_true(all(fit$predictions$reliability >= 0 &
                    fit$predictions$reliability <= 1))
  expect_true(all(validate(fit)$passed))
  expect_equal(
    fit$reml_residual_df,
    fit$n_observed - fit$fixed_effect_rank
  )
})

test_that("multivariate GBLUP estimates positive covariance matrices", {
  K <- make_test_kernel(n = 12L, p = 24L, seed = 33L)
  ids <- rownames(K)
  set.seed(34)
  Y <- cbind(
    trait1 = 0.8 * scale(K %*% rnorm(length(ids)))[, 1] +
      rnorm(length(ids), sd = 0.4),
    trait2 = -0.5 * scale(K %*% rnorm(length(ids)))[, 1] +
      rnorm(length(ids), sd = 0.5)
  )
  rownames(Y) <- ids

  fit <- fit_multitrait_gblup(Y, K, maxit = 300L)
  expect_gt(min(eigen(fit$genetic_cov, symmetric = TRUE)$values), 0)
  expect_gt(min(eigen(fit$residual_cov, symmetric = TRUE)$values), 0)
  expect_true(is.finite(fit$log_likelihood))
})

test_that("environment kernels and reaction-norm GBLUP are auditable", {
  environments <- data.frame(
    environment = c("E1", "E2", "E3"),
    rainfall = c(410, 560, 720),
    temperature = c(25, 22, 19)
  )
  Ke <- build_environment_kernel(environments)
  expect_equal(dim(Ke), c(3L, 3L))
  expect_gte(min(eigen(Ke, symmetric = TRUE)$values), -1e-8)

  K <- make_test_kernel(n = 10L, p = 20L, seed = 35L)
  ids <- rownames(K)
  records <- expand.grid(
    id = ids,
    environment = environments$environment,
    stringsAsFactors = FALSE
  )
  set.seed(36)
  genotype_effect <- setNames(rnorm(length(ids), sd = 0.8), ids)
  sensitivity <- setNames(rnorm(length(ids), sd = 0.4), ids)
  env_score <- setNames(scale(environments$rainfall)[, 1],
                        environments$environment)
  records$phenotype <- 5 +
    genotype_effect[records$id] +
    sensitivity[records$id] * env_score[records$environment] +
    rnorm(nrow(records), sd = 0.25)

  fit <- fit_gxe_gblup(records, K, Ke, min_reliability = 0.10,
                       maxit = 300L)
  expect_s3_class(fit, "hapblockr_result")
  expect_equal(nrow(fit$predictions), nrow(records))
  expect_true(all(fit$variance_components > 0))
  expect_equal(nrow(fit$environment_stability), length(ids))
  expect_true(all(validate(fit)$passed))
  expect_true(is.finite(fit$log_likelihood))
  expect_equal(
    fit$reml_residual_df,
    nrow(records) - fit$fixed_effect_rank
  )
})

test_that("prepared targets propagate precision into multivariate and GxE models", {
  K <- make_test_kernel(n = 10L, p = 20L, seed = 81L)
  ids <- rownames(K)
  set.seed(82)
  wide <- expand.grid(
    id = ids,
    trait = c("yield", "disease"),
    stringsAsFactors = FALSE
  )
  wide$value <- rnorm(nrow(wide))
  wide$SE <- rep(c(0.2, 0.5), each = length(ids))
  targets <- prepare_breeding_targets(
    wide, "BLUE", se_col = "SE"
  )
  Sg <- diag(c(yield = 1, disease = 1))
  Se <- diag(c(yield = 0.5, disease = 0.5))
  dimnames(Sg) <- dimnames(Se) <- list(
    c("yield", "disease"), c("yield", "disease")
  )
  multi <- fit_multitrait_gblup(
    targets, K,
    genetic_cov = Sg,
    residual_cov = Se,
    estimate_covariances = FALSE
  )
  expect_match(
    multi$result_contract$transformations[1],
    "precision-weighted"
  )
  expect_false(is.null(multi$target_provenance))

  records <- expand.grid(
    id = ids,
    environment = c("E1", "E2"),
    stringsAsFactors = FALSE
  )
  records$trait <- "yield"
  records$value <- rnorm(nrow(records))
  records$SE <- rep(c(0.2, 0.6), each = length(ids))
  gxe_targets <- prepare_breeding_targets(
    records,
    "adjusted_mean",
    environment_col = "environment",
    se_col = "SE"
  )
  gxe <- fit_gxe_gblup(gxe_targets, K, maxit = 300L)
  expect_equal(nrow(gxe$across_environment_predictions), length(ids))
  expect_identical(
    anyDuplicated(gxe$across_environment_predictions$id),
    0L
  )
  expect_identical(
    rownames(gxe$across_environment_predictions),
    as.character(seq_len(length(ids)))
  )
  expect_false(is.null(gxe$target_provenance))
})

test_that("full sampling covariance is aligned and has explicit residual semantics", {
  ids <- paste0("G", 1:5)
  traits <- c("yield", "disease")
  records <- expand.grid(
    id = ids,
    trait = traits,
    stringsAsFactors = FALSE
  )
  records <- records[
    !(records$id == "G5" & records$trait == "disease"), , drop = FALSE
  ]
  records$environment <- "E1"
  records$value <- seq_len(nrow(records)) / 10
  records$SE <- ifelse(records$trait == "yield", 0.2, 0.4)
  keys <- paste(
    records$id, records$trait, records$environment, sep = "::"
  )
  S <- diag(records$SE^2)
  dimnames(S) <- list(keys, keys)
  targets <- prepare_breeding_targets(
    records,
    input_type = "BLUE",
    environment_col = "environment",
    se_col = "SE",
    covariance = S
  )
  K <- diag(length(ids))
  dimnames(K) <- list(ids, ids)
  G <- diag(c(yield = 0.5, disease = 0.4))
  dimnames(G) <- list(traits, traits)

  fit <- fit_multitrait_gblup(
    targets,
    K,
    genetic_cov = G,
    estimate_covariances = FALSE
  )
  expect_identical(fit$sampling_covariance_mode, "sampling_only")
  expect_null(fit$residual_cov)
  expect_false(anyNA(rownames(fit$sampling_covariance)))
  observed_covariance <- rownames(fit$sampling_covariance) %in% rownames(S)
  expect_equal(
    diag(fit$sampling_covariance)[observed_covariance],
    diag(S)[match(
      rownames(fit$sampling_covariance)[observed_covariance],
      rownames(S)
    )]
  )
  expect_equal(
    unname(diag(fit$sampling_covariance)[!observed_covariance]),
    0
  )
  sampling_gate <- fit$result_contract$quality_gates
  expect_true(sampling_gate$passed[
    sampling_gate$gate == "sampling_covariance_used"
  ])

  E <- diag(c(yield = 0.2, disease = 0.2))
  dimnames(E) <- list(traits, traits)
  fit_plus <- fit_multitrait_gblup(
    targets,
    K,
    genetic_cov = G,
    residual_cov = E,
    sampling_covariance_mode = "sampling_plus_residual",
    estimate_covariances = FALSE
  )
  expect_equal(fit_plus$residual_cov, E)
  expect_identical(
    fit_plus$sampling_covariance_mode,
    "sampling_plus_residual"
  )
})

test_that("sampling covariance and precision inconsistencies are rejected", {
  ids <- paste0("G", 1:4)
  records <- expand.grid(
    id = ids,
    trait = c("t1", "t2"),
    stringsAsFactors = FALSE
  )
  records$value <- rnorm(nrow(records))
  records$SE <- rep(c(0.2, 0.4), each = length(ids))
  keys <- paste(records$id, records$trait, "ACROSS", sep = "::")
  S <- diag(records$SE^2)
  dimnames(S) <- list(keys, keys)
  targets <- prepare_breeding_targets(
    records, "BLUE", se_col = "SE", covariance = S
  )
  targets$targets$precision_weight[1] <- 99
  K <- diag(length(ids))
  dimnames(K) <- list(ids, ids)

  expect_error(
    fit_multitrait_gblup(targets, K),
    "inconsistent with the target precision weights"
  )
})

test_that("reported multivariate REML likelihood uses residual degrees of freedom", {
  K <- make_test_kernel(n = 6L, p = 12L, seed = 93L)
  ids <- rownames(K)
  traits <- c("t1", "t2")
  Y <- cbind(
    t1 = seq(-1, 1, length.out = 6),
    t2 = c(-0.4, 0.2, 0.8, -0.3, 0.5, 1.1)
  )
  rownames(Y) <- ids
  G <- diag(c(t1 = 0.7, t2 = 0.5))
  E <- diag(c(t1 = 0.3, t2 = 0.4))
  dimnames(G) <- dimnames(E) <- list(traits, traits)

  fit <- fit_multitrait_gblup(
    Y, K,
    genetic_cov = G,
    residual_cov = E,
    estimate_covariances = FALSE
  )
  n <- nrow(Y)
  q <- ncol(Y)
  y <- as.vector(Y)
  X <- kronecker(diag(q), matrix(1, n, 1))
  V <- kronecker(G, K) + kronecker(E, diag(n))
  V_inv <- solve(V)
  XtVinvX <- crossprod(X, V_inv %*% X)
  beta <- solve(XtVinvX, crossprod(X, V_inv %*% y))
  residual <- y - X %*% beta
  expected <- -0.5 * (
    (length(y) - qr(X)$rank) * log(2 * pi) +
      as.numeric(determinant(V, logarithm = TRUE)$modulus) +
      as.numeric(determinant(XtVinvX, logarithm = TRUE)$modulus) +
      as.numeric(crossprod(residual, V_inv %*% residual))
  )
  expect_equal(fit$log_likelihood, expected, tolerance = 1e-8)
})
