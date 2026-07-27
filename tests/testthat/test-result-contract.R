library(testthat)
library(HapBlockR)

test_that("decision outputs carry a valid versioned result contract", {
  ids <- c("A", "B", "C", "D")
  coordinates <- c(0, 1, 4, 8)
  D <- as.matrix(stats::dist(coordinates))
  dimnames(D) <- list(ids, ids)

  result <- select_core_collection(
    G = D, n_core = 3L, type = "distance", verbose = FALSE
  )

  expect_s3_class(result, "hapblockr_result")
  expect_identical(result$result_contract$schema_version, "1.0.0")
  expect_identical(result$result_contract$method, "select_core_collection")
  expect_match(
    unname(result$result_contract$input_hashes),
    "^[0-9a-f]{64}$"
  )
  report <- validate_hapblockr_result(result)
  expect_true(all(report$passed))
  expect_true(all(validate(result)$passed))
})

test_that("result contract methods expose decisions and diagnostics", {
  ids <- c("A", "B", "C")
  D <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3L, 3L,
              dimnames = list(ids, ids))
  result <- select_core_collection(
    G = D, n_core = 2L, type = "distance", verbose = FALSE
  )

  decisions <- as.data.frame(result)
  expect_equal(decisions$id, result$selected)
  expect_equal(summary(result)$validation_status, "passed")
  expect_output(print(result), "select_core_collection")
})

test_that("failed quality gates are machine-readable", {
  result <- HapBlockR:::.add_hapblockr_contract(
    result = list(),
    method = "test_method",
    call = quote(test_method()),
    quality_gates = c(scientific_gate = FALSE)
  )
  expect_identical(result$result_contract$validation_status, "failed")
  expect_error(validate_hapblockr_result(result), "quality_gates")
  report <- validate_hapblockr_result(result, strict = FALSE)
  expect_false(report$passed[report$check == "quality_gates"])
})
