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
})

test_that("DesiredGainR engines return an internal merit score", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.2.0")
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

  W <- diag(c(yield = 0.1, disease = 0.05))
  dimnames(W) <- list(traits, traits)
  qgsi <- build_selection_index(
    values, G, P,
    economic_weights = c(yield = 1, disease = 0.5),
    directions = c(yield = "increase", disease = "decrease"),
    method = "qgsi",
    quadratic_weights = W
  )
  expect_identical(qgsi$index_type, "qgsi")
  expect_match(qgsi$contribution_scope, "candidate-specific")
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
  expect_false(is.null(gxe$target_provenance))
})
