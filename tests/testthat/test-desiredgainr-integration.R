# Exercise the dependency's public API, not a mock of its implementation.
desiredgainr_fixture <- function() {
  set.seed(827)
  values <- cbind(yield = rnorm(40, 8, 2), disease = rnorm(40, 3, 1))
  rownames(values) <- sprintf("P%02d", seq_len(nrow(values)))
  traits <- colnames(values)
  G <- matrix(c(2, -0.2, -0.2, 0.5), 2,
              dimnames = list(traits, traits))
  P <- G + diag(c(1, 0.5))
  W <- matrix(c(0.10, -0.03, -0.03, -0.05), 2,
              dimnames = list(traits, traits))
  candidate <- data.frame(id = rownames(values), values, check.names = FALSE)
  list(
    values = values, G = G, P = P, W = W, candidate = candidate,
    traits = traits, directions = c(yield = "increase", disease = "decrease")
  )
}

test_that("DGSI controls and returned results match direct DesiredGainR calls", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  f <- desiredgainr_fixture()
  reference <- f$candidate
  reference$yield <- reference$yield * 1.3 + 2
  reference$disease <- reference$disease * 0.8 - 1
  validation <- reference
  validation$id <- paste0("V", seq_len(nrow(validation)))
  base_control <- list(n_iter = 8L, n_rep = 2L, seed = 827L)
  cases <- list(
    list(replicate_selection = "training"),
    list(scale_traits = TRUE, ref_data = reference, holdout_fraction = 0.25),
    list(scale_traits = TRUE, ref_data = reference,
         validation_data = validation, return_all_reps = FALSE),
    list(replicate_selection = "training", select_mode = "eligible_top_n",
         trait_min = c(yield = max(f$values[, "yield"]),
                       disease = -max(f$values[, "disease"]) - 1))
  )
  for (case in cases) {
    control <- c(base_control, case)
    fit <- build_selection_index(
      f$values, f$G, f$P, desired_gains = c(yield = 0.5, disease = 0.3),
      directions = f$directions, method = "dgsi", n_select = 8,
      dgsi_control = control
    )
    direct <- do.call(DesiredGainR::run_dgsi, c(list(
      init_data = f$candidate["id"], cand_data = f$candidate,
      trait_cols = f$traits, dg = c(yield = 0.5, disease = 0.3),
      P = f$P, G = f$G, id_col = "id", lower_is_better = "disease",
      n_select = 8L
    ), control))
    for (field in setdiff(names(direct), "call")) {
      expect_equal(fit$engine_result[[field]], direct[[field]], info = field)
    }
    scores <- as.data.frame(direct$ranked_geno)
    expect_identical(fit$scores$id, scores$id)
    expect_equal(fit$scores$selection_index, scores$SelectionIndex)
    expect_identical(fit$scores$Selected, scores$Selected)
    expect_equal(fit$coefficients$expected_response,
                 unname(direct$theoretical_response$original_units[f$traits]))
    # These coefficients can be applied to raw trait effects without a
    # second direction reversal or a forgotten reference-scale conversion.
    reconstructed <- as.numeric(
      f$values %*% fit$coefficients_original_units + fit$score_intercept
    )
    expect_equal(
      reconstructed[match(fit$scores$id, rownames(f$values))],
      fit$scores$selection_index, tolerance = 1e-10
    )
    expect_equal(fit$coefficients$coefficient,
                 unname(direct$coefficients[f$traits]))
    expect_equal(fit$result_contract$parameters$n_selected, sum(scores$Selected))
    expect_true(all(validate(fit)$passed))
    expect_true(any(grepl(direct$optimism$selection_rule,
                         fit$result_contract$transformations, fixed = TRUE)))
  }
  expect_equal(sum(fit$scores$Selected), 1L)
})

test_that("DGSI preserves corrected candidate-SD gains across trait units", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  f <- desiredgainr_fixture()
  control <- list(n_iter = 8, n_rep = 2, seed = 827,
                  replicate_selection = "training", ridge_P = 0, ridge_M = 0)
  fit <- build_selection_index(
    f$values, f$G, f$P, desired_gains = c(yield = 0.5, disease = 0.3),
    directions = f$directions, method = "dgsi", n_select = 8,
    dgsi_control = control
  )
  # Change measurement units without changing the scientific objective.
  units_change <- c(yield = 1.2, disease = 0.9)
  S <- diag(units_change)
  G2 <- S %*% f$G %*% S
  P2 <- S %*% f$P %*% S
  dimnames(G2) <- dimnames(P2) <- list(f$traits, f$traits)
  converted <- build_selection_index(
    sweep(f$values, 2L, units_change, "*"), G2, P2,
    desired_gains = c(yield = 0.5, disease = 0.3), directions = f$directions,
    method = "dgsi", n_select = 8, dgsi_control = control
  )
  expect_equal(converted$engine_result$non_iterated$realised_response,
               fit$engine_result$non_iterated$realised_response)
  expect_equal(converted$coefficients_original_units * units_change,
               fit$coefficients_original_units, tolerance = 1e-8)
  expect_identical(converted$scores$id, fit$scores$id)
  expect_equal(converted$scores$selection_index, fit$scores$selection_index)
  expect_equal(converted$coefficients$expected_response,
               fit$coefficients$expected_response * unname(units_change))
})

test_that("QGSI passes covariance, reference and contribution controls exactly", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  f <- desiredgainr_fixture()
  reference <- f$candidate
  reference$yield <- reference$yield + 0.5
  reference$id <- paste0("R", seq_len(nrow(reference)))
  K <- diag(0.8, nrow(reference)) + 0.2
  dimnames(K) <- list(reference$id, reference$id)
  K <- K[rev(reference$id), rev(reference$id)]
  cases <- list(
    list(),
    list(scale_traits = TRUE, Gamma = f$G, true_G = f$P),
    list(scale_traits = TRUE, reference_gebv_data = reference,
         relationship_matrix = K, return_contributions = FALSE)
  )
  for (control in cases) {
    fit <- build_selection_index(
      f$values, f$G, f$P, economic_weights = c(yield = 1, disease = 0.4),
      directions = f$directions, method = "qgsi", n_select = 8,
      quadratic_weights = f$W, qgsi_control = control
    )
    direct <- do.call(DesiredGainR::run_qgsi, c(list(
      init_data = f$candidate["id"], gebv_data = f$candidate,
      trait_cols = f$traits, linear_weights = c(yield = 1, disease = 0.4),
      W = f$W, id_col = "id", lower_is_better = "disease", n_select = 8L
    ), control))
    for (field in setdiff(names(direct), "call")) {
      expect_equal(fit$engine_result[[field]], direct[[field]], info = field)
    }
    ranking <- as.data.frame(direct$ranked_geno)
    expect_identical(fit$scores$id, ranking$id)
    expect_equal(fit$scores$merit_score, ranking$QGSI)
    expect_identical(fit$scores$Selected, ranking$Selected)
    gain <- as.data.frame(direct$expected_gain_per_trait)
    original_gain <- gain$Expected_Genetic_Gain[match(f$traits, gain$Trait)] *
      direct$transformation$scale[f$traits] *
      direct$transformation$direction[f$traits]
    expect_equal(fit$coefficients$expected_response, unname(original_gain))
    expect_null(fit$coefficients_original_units)
    expect_null(fit$score_intercept)
    expect_equal(fit$genetic_cov, f$G)
    expect_true(all(validate(fit)$passed))
  }
})

test_that("delegated controls reject typos, partial names and overrides", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  f <- desiredgainr_fixture()
  for (method in c("dgsi", "qgsi")) {
    args <- list(trait_values = f$values, genetic_cov = f$G,
                 phenotypic_cov = f$P, directions = f$directions,
                 method = method, n_select = 8)
    if (method == "dgsi") {
      args$desired_gains <- c(yield = 0.5, disease = 0.3)
    } else {
      args$economic_weights <- c(yield = 1, disease = 0.4)
      args$quadratic_weights <- f$W
    }
    control_name <- paste0(method, "_control")
    for (bad_name in c("scale_trait", "lower_is_bet", "not_an_argument")) {
      args[[control_name]] <- setNames(list(TRUE), bad_name)
      expect_error(do.call(build_selection_index, args),
                   paste0("unknown argument.*", bad_name))
    }
    args[[control_name]] <- list(lower_is_better = "yield")
    expect_error(do.call(build_selection_index, args), "cannot override")
    args[[control_name]] <- list(scale_traits = TRUE, scale_traits = FALSE)
    expect_error(do.call(build_selection_index, args), "unique, non-empty")
  }
})

test_that("delegated provenance records the objective and engine version", {
  skip_if_not_installed("DesiredGainR", minimum_version = "0.5.0")
  f <- desiredgainr_fixture()
  fit <- build_selection_index(
    f$values, f$G, f$P, economic_weights = c(yield = 1, disease = 0.4),
    directions = f$directions, method = "qgsi", n_select = 8,
    quadratic_weights = f$W
  )
  expect_identical(fit$result_contract$software$external_tools$DesiredGainR,
                   as.character(utils::packageVersion("DesiredGainR")))
  expect_true(all(c("objective", "quadratic_weights", "engine_control") %in%
                    names(fit$result_contract$input_hashes)))
  bad <- f$values
  colnames(bad)[1L] <- "id"
  G <- f$G
  P <- f$P
  dimnames(G) <- dimnames(P) <- list(colnames(bad), colnames(bad))
  expect_error(build_selection_index(
    bad, G, P, desired_gains = c(id = 0.5, disease = 0.3),
    directions = c(id = "increase", disease = "decrease"), method = "dgsi"
  ), "trait name 'id' is reserved")
})
