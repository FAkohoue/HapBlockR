make_index_scenario <- function(values) {
  traits <- colnames(values)
  G <- diag(c(1.0, 0.7))
  P <- diag(c(1.4, 1.1))
  dimnames(G) <- dimnames(P) <- list(traits, traits)
  build_selection_index(
    values, G, P,
    economic_weights = c(yield = 1, disease = 0.5),
    directions = c(yield = "increase", disease = "decrease"),
    units = c(yield = "t/ha", disease = "score")
  )
}

test_that("decision stability compares threshold or holdout scenarios", {
  baseline_values <- cbind(
    yield = c(5.2, 5.0, 4.8, 4.6),
    disease = c(1.5, 1.0, 2.5, 0.8)
  )
  rownames(baseline_values) <- paste0("G", 1:4)
  leave_e1 <- baseline_values
  leave_e1["G1", "yield"] <- 4.7
  leave_e2 <- baseline_values
  leave_e2["G4", "yield"] <- 5.1

  stability <- assess_decision_stability(
    list(
      baseline = make_index_scenario(baseline_values),
      leave_E1_out = make_index_scenario(leave_e1),
      leave_E2_out = make_index_scenario(leave_e2)
    ),
    top_n = 2,
    minimum_selection_frequency = 2 / 3,
    minimum_baseline_jaccard = 0.3
  )
  expect_s3_class(stability, "hapblockr_result")
  expect_equal(nrow(stability$pairwise_overlap), 3L)
  expect_true(all(
    stability$stability$selection_frequency >= 0 &
      stability$stability$selection_frequency <= 1
  ))
  expect_true(any(stability$stability$stable_recommendation))
  expect_true(all(validate(stability)$passed))
})

test_that("decision stability refuses unnamed and failed scenarios", {
  values <- cbind(
    yield = c(5, 4, 3),
    disease = c(1, 2, 3)
  )
  rownames(values) <- paste0("G", 1:3)
  scenario <- make_index_scenario(values)
  expect_error(
    assess_decision_stability(list(scenario, scenario)),
    "uniquely named"
  )

  failed <- scenario
  failed$result_contract$validation_status <- "failed"
  expect_error(
    assess_decision_stability(
      list(baseline = scenario, failed = failed)
    ),
    "Failed result contracts"
  )
})
