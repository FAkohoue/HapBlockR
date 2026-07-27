make_candidate_status <- function() {
  data.frame(
    id = c("A", "B", "C", "D"),
    female_allowed = c(TRUE, TRUE, TRUE, FALSE),
    male_allowed = c(TRUE, TRUE, FALSE, TRUE),
    fertile = TRUE,
    flowering_start = c(5, 6, 7, 7),
    flowering_end = c(7, 8, 9, 9),
    heterotic_group = c("H1", "H2", "H1", "H2"),
    quarantine_group = c("Q1", "Q1", "Q2", "Q2"),
    subpopulation = c("S1", "S1", "S2", "S2"),
    female_capacity = c(30, 20, 20, 0),
    male_capacity = c(20, 30, 0, 20),
    total_capacity = c(30, 30, 20, 20),
    seed_available = c(30, 20, 20, 0),
    pollen_available = c(20, 30, 0, 20),
    stringsAsFactors = FALSE
  )
}

test_that("cross screening reports every operational reason", {
  status <- make_candidate_status()
  quarantine <- matrix(
    c(TRUE, FALSE, FALSE, TRUE),
    2, 2,
    dimnames = list(c("Q1", "Q2"), c("Q1", "Q2"))
  )
  pairs <- rbind(
    c("A", "B"),
    c("A", "C"),
    c("D", "A"),
    c("A", "A")
  )
  screened <- screen_candidate_crosses(
    pairs,
    status,
    require_different_heterotic_groups = TRUE,
    quarantine_compatibility = quarantine,
    forbidden_pairs = data.frame(parent1 = "A", parent2 = "B")
  )

  expect_false(screened$feasible[1])
  expect_match(screened$reasons[1], "forbidden_pair")
  expect_match(screened$reasons[2], "male_role_forbidden")
  expect_match(screened$reasons[2], "same_heterotic_group_forbidden")
  expect_match(screened$reasons[2], "quarantine_incompatible")
  expect_match(screened$reasons[3], "female_role_forbidden")
  expect_match(screened$reasons[4], "selfing_forbidden")
})

test_that("a feasible mating plan receives a complete certificate", {
  status <- make_candidate_status()
  plan <- data.frame(
    female = c("A", "C"),
    male = c("B", "D"),
    family_size = c(20, 20),
    period = c("P1", "P1")
  )
  period_capacity <- expand.grid(
    id = status$id,
    period = "P1",
    stringsAsFactors = FALSE
  )
  period_capacity$total_capacity <- c(20, 20, 20, 20)
  quotas <- data.frame(
    subpopulation = c("S1", "S2"),
    min_contribution = c(40, 40),
    max_contribution = c(50, 50)
  )

  certificate <- certify_mating_plan(
    plan,
    status,
    min_family_size = 20,
    required_pairs = data.frame(parent1 = "A", parent2 = "B"),
    capacity_by_period = period_capacity,
    subpopulation_quotas = quotas,
    require_different_heterotic_groups = TRUE
  )

  expect_s3_class(certificate, "hapblockr_result")
  expect_true(certificate$certificate$feasible)
  expect_length(certificate$violations$constraint, 0L)
  expect_gt(nrow(certificate$binding_constraints), 0L)
  expect_true(all(validate(certificate)$passed))
})

test_that("infeasible plans retain violations and binding information", {
  status <- make_candidate_status()
  plan <- data.frame(
    female = c("A", "B", "A"),
    male = c("A", "A", "B"),
    family_size = c(5, 25, 25),
    period = "P1"
  )
  quotas <- data.frame(
    subpopulation = "S2",
    min_contribution = 20
  )

  certificate <- certify_mating_plan(
    plan,
    status,
    min_family_size = 10,
    required_pairs = data.frame(parent1 = "C", parent2 = "D"),
    subpopulation_quotas = quotas
  )

  expect_false(certificate$certificate$feasible)
  expect_true(all(c(
    "cross_rule", "minimum_family_size", "no_repeated_cross",
    "required_pair", "min_contribution"
  ) %in% certificate$violations$constraint))
  report <- validate(certificate, strict = FALSE)
  expect_true(any(!report$passed))
})
