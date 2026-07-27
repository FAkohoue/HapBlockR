test_that("accepted breeding-target types are defined explicitly", {
  definitions <- breeding_target_types()
  expect_equal(
    definitions$input_type,
    c(
      "adjusted_mean", "BLUE", "BLUP_identity", "PBLUP",
      "BV", "GCA", "TGV"
    )
  )
  expect_true(all(nzchar(definitions$definition)))
})

test_that("BLUEs retain optional units and use inverse-variance precision", {
  values <- data.frame(
    id = paste0("G", 1:4),
    trait = "yield",
    value = c(4.1, 5.0, 4.7, 5.2),
    SE = c(0.2, 0.4, 0.25, 0.5),
    stringsAsFactors = FALSE
  )
  result <- prepare_breeding_targets(
    values,
    input_type = "BLUE",
    se_col = "SE"
  )
  expect_s3_class(result, "hapblockr_result")
  expect_true(all(is.na(result$targets$unit)))
  expect_equal(mean(result$targets$precision_weight), 1)
  expect_equal(
    result$targets$precision_raw,
    1 / values$SE^2
  )
  expect_equal(result$targets$model_value, values$value)
  expect_true(all(validate(result)$passed))
})

test_that("identity BLUPs are deregressed from PEV and genetic variance", {
  values <- data.frame(
    id = paste0("G", 1:3),
    trait = "yield",
    value = c(0.4, 0.8, 1.2),
    PEV = c(0.2, 0.4, 0.1)
  )
  result <- prepare_breeding_targets(
    values,
    input_type = "BLUP_identity",
    pev_col = "PEV",
    genetic_variance = 2
  )
  expected_reliability <- 1 - values$PEV / 2
  expect_equal(result$targets$reliability, expected_reliability)
  expect_equal(
    result$targets$deregressed_value,
    values$value / expected_reliability
  )
})

test_that("supplied reliability does not require genetic variance", {
  values <- data.frame(
    id = paste0("G", 1:3),
    trait = "yield",
    value = c(0.4, 0.8, 1.2),
    reliability = c(0.5, 0.8, 0.9)
  )
  result <- prepare_breeding_targets(
    values,
    input_type = "BLUP_identity",
    reliability_col = "reliability"
  )
  expect_equal(result$targets$reliability, values$reliability)
  expect_true(all(is.na(result$targets$PEV)))
  expect_equal(
    result$targets$deregressed_value,
    values$value / values$reliability
  )
})

test_that("PBLUP reliability uses the A-matrix diagonal", {
  values <- data.frame(
    id = c("G1", "G2"),
    trait = "yield",
    value = c(0.5, 0.9),
    PEV = c(0.2, 0.2)
  )
  A <- matrix(c(1, 0.25, 0.25, 1.5), 2)
  dimnames(A) <- list(values$id, values$id)
  result <- prepare_breeding_targets(
    values,
    input_type = "PBLUP",
    pev_col = "PEV",
    genetic_variance = 2,
    relationship_matrix = A
  )
  expect_equal(
    result$targets$reliability,
    unname(1 - values$PEV / (diag(A) * 2))
  )
})

test_that("lower-is-better traits are oriented to higher model values", {
  values <- data.frame(
    id = paste0("G", 1:3),
    trait = "disease",
    value = c(10, 20, 30),
    SE = rep(1, 3)
  )
  result <- prepare_breeding_targets(
    values,
    input_type = "adjusted_mean",
    se_col = "SE",
    lower_is_better = "disease"
  )
  expect_equal(result$targets$model_value, -values$value)
  expect_equal(result$targets$direction, rep(-1, 3))
})

test_that("external genomic predictions and indices are rejected", {
  values <- data.frame(
    id = "G1", trait = "yield", value = 1, SE = 0.2
  )
  expect_error(
    prepare_breeding_targets(values, "GBLUP", se_col = "SE"),
    "not an accepted external target"
  )
  expect_error(
    prepare_breeding_targets(values, "selection_index", se_col = "SE"),
    "not an accepted external target"
  )
  expect_error(
    prepare_breeding_targets(
      values, "BV", estimation_basis = "genomic", se_col = "SE"
    ),
    "not accepted"
  )
})

test_that("GCA and TGV preserve their breeding interpretation", {
  gca <- data.frame(
    id = paste0("G", 1:2),
    trait = "yield",
    value = c(0.2, 0.4),
    SE = c(0.1, 0.2)
  )
  fixed <- prepare_breeding_targets(
    gca,
    input_type = "GCA",
    se_col = "SE",
    gca_effect = "fixed",
    tester_population = "programme testers in cycle 4"
  )
  expect_equal(unique(fixed$targets$estimation_basis), "fixed")
  expect_error(
    prepare_breeding_targets(
      gca,
      input_type = "GCA",
      se_col = "SE",
      tester_population = "programme testers in cycle 4"
    ),
    "requires gca_effect"
  )

  tgv <- data.frame(
    id = paste0("G", 1:2),
    trait = "yield",
    additive = c(1.0, 1.5),
    dominance = c(0.2, -0.1),
    SE = c(0.2, 0.3)
  )
  total <- prepare_breeding_targets(
    tgv,
    input_type = "TGV",
    additive_col = "additive",
    dominance_col = "dominance",
    estimation_basis = "fixed",
    se_col = "SE"
  )
  expect_equal(
    total$targets$value,
    tgv$additive + tgv$dominance
  )
})

test_that("one environment can be modelled separately but several require GxE", {
  one_environment <- prepare_breeding_targets(
    data.frame(
      id = paste0("G", 1:3),
      trait = "yield",
      environment = "E1",
      value = c(1, 2, 3),
      SE = c(0.2, 0.3, 0.2)
    ),
    input_type = "BLUE",
    environment_col = "environment",
    se_col = "SE"
  )
  expect_silent(HapBlockR:::.hb_unpack_model_targets(one_environment))

  several_environments <- prepare_breeding_targets(
    data.frame(
      id = rep(paste0("G", 1:3), 2),
      trait = "yield",
      environment = rep(c("E1", "E2"), each = 3),
      value = 1:6,
      SE = rep(0.2, 6)
    ),
    input_type = "BLUE",
    environment_col = "environment",
    se_col = "SE"
  )
  expect_error(
    HapBlockR:::.hb_unpack_model_targets(several_environments),
    "multiple environments"
  )
})

test_that("replicate-level records are rejected as duplicate targets", {
  values <- data.frame(
    id = c("G1", "G1"),
    trait = c("yield", "yield"),
    value = c(4, 5),
    SE = c(0.2, 0.2),
    replicate = c("R1", "R2")
  )
  expect_error(
    prepare_breeding_targets(values, "BLUE", se_col = "SE"),
    "externally adjusted summaries"
  )
})
