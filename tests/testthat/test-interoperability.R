make_valid_metadata <- function() {
  germplasm <- data.frame(
    germplasm_id = c("G1", "G2"),
    breeding_program_id = "BP1"
  )
  traits <- data.frame(
    trait_id = "T1",
    trait_name = "Grain yield",
    unit = "t/ha",
    direction = "increase",
    ontology_id = "CO_321:0000012"
  )
  environments <- data.frame(
    environment_id = c("E1", "E2"),
    study_id = "S1",
    trial_id = c("TR1", "TR2"),
    location_id = c("L1", "L2")
  )
  observations <- data.frame(
    observation_id = paste0("O", 1:4),
    germplasm_id = rep(c("G1", "G2"), 2),
    trait_id = "T1",
    environment_id = rep(c("E1", "E2"), each = 2),
    value = c(4.1, 4.5, 3.9, 4.3),
    unit = "t/ha"
  )
  list(
    germplasm = germplasm,
    traits = traits,
    environments = environments,
    observations = observations
  )
}

test_that("canonical breeding metadata validates identity, ontology, and units", {
  x <- make_valid_metadata()
  result <- validate_breeding_metadata(
    x$germplasm, x$traits, x$environments, x$observations
  )
  expect_s3_class(result, "hapblockr_result")
  expect_length(result$issues$code, 0L)
  expect_true(all(result$summary$valid))
  expect_true(all(validate(result)$passed))
  expect_true(all(c("BrAPI", "MIAPPE") %in%
                    result$standards_mapping$standard))
})

test_that("metadata failures return a complete ledger in non-strict mode", {
  x <- make_valid_metadata()
  x$observations$germplasm_id[1] <- "UNKNOWN"
  x$observations$unit[2] <- "kg/plot"
  x$traits$ontology_id <- "not a CURIE"

  result <- validate_breeding_metadata(
    x$germplasm, x$traits, x$environments, x$observations,
    strict = FALSE
  )
  expect_true(all(c(
    "unknown_foreign_key", "unit_mismatch",
    "invalid_ontology_identifier"
  ) %in% result$issues$code))
  report <- validate(result, strict = FALSE)
  expect_true(any(!report$passed))

  expect_error(
    validate_breeding_metadata(
      x$germplasm, x$traits, x$environments, x$observations
    ),
    "validation failed"
  )
})

test_that("recommendation exchange bundles round-trip with checksums", {
  values <- cbind(yield = c(4.0, 5.0, 4.5),
                  disease = c(2.0, 3.0, 1.0))
  rownames(values) <- c("G1", "G2", "G3")
  covariance_names <- colnames(values)
  G <- diag(c(1.0, 0.8))
  P <- diag(c(1.5, 1.2))
  dimnames(G) <- dimnames(P) <- list(
    covariance_names, covariance_names
  )
  decision <- build_selection_index(
    values, G, P,
    economic_weights = c(yield = 1, disease = 0.5),
    directions = c(yield = "increase", disease = "decrease"),
    units = c(yield = "t/ha", disease = "score")
  )
  bundle <- build_breeding_exchange(
    decision,
    breeding_program_id = "BP1",
    study_id = "S1",
    trial_id = "TR1",
    environment_id = "E1",
    trait_id = "INDEX1",
    unit = "index unit",
    ontology_id = "CO_321:0000012"
  )
  expect_s3_class(bundle, "HapBlockR_exchange_bundle")
  expect_equal(nrow(bundle$recommendations), 3L)
  expect_true(all(bundle$recommendations$germplasm_id %in%
                    rownames(values)))

  path <- tempfile("exchange-")
  manifest <- write_breeding_exchange(bundle, path)
  expect_true(all(file.exists(file.path(path, manifest$file))))
  restored <- read_breeding_exchange(path)
  expect_equal(restored$recommendations, bundle$recommendations,
               tolerance = 1e-12)
  expect_identical(
    attr(restored, "bundle_schema_version"),
    "1.0.0"
  )

  writeLines("tampered", file.path(path, "recommendations.csv"))
  expect_error(read_breeding_exchange(path), "checksum")
})

test_that("failed decisions are not exported by default", {
  status <- data.frame(
    id = c("A", "B"),
    female_allowed = TRUE,
    male_allowed = TRUE,
    fertile = TRUE
  )
  failed <- certify_mating_plan(
    data.frame(female = "A", male = "A", family_size = 1),
    status
  )
  expect_error(
    build_breeding_exchange(
      failed,
      breeding_program_id = "BP1",
      study_id = "S1",
      trial_id = "TR1"
    ),
    "failed.*cannot be exported"
  )
})
