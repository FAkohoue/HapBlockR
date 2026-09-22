## tests/testthat/test-breeder-guide.R
## -----------------------------------------------------------------------------
## Tests for R/breeder_guide.R: open_breeder_guide()
##
## open = FALSE is used throughout so these tests never attempt to launch an
## external application in a non-interactive test/CI environment (the
## function itself also gates opening on interactive(), but open = FALSE is
## used explicitly here to keep the test's intent self-evident).
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

test_that("open_breeder_guide: locates the installed .pdf file", {
  path <- open_breeder_guide(open = FALSE)
  expect_true(is.character(path))
  expect_true(file.exists(path))
  expect_match(path, "HapBlockR_Breeder_Guide\\.pdf$")
})

test_that("open_breeder_guide: locates the editable Word edition", {
  path <- open_breeder_guide(open = FALSE, format = "html")
  expect_true(is.character(path))
  expect_true(file.exists(path))
  expect_match(path, "HapBlockR_Breeder_Guide\\.html$")
})

test_that("open_breeder_guide: returns the path invisibly", {
  expect_invisible(open_breeder_guide(open = FALSE))
})

test_that("open_breeder_guide: validates the requested format", {
  expect_error(
    open_breeder_guide(open = FALSE, format = "txt"),
    "'arg' should be one of"
  )
})

test_that("open_breeder_guide: does not error in a non-interactive session even with open = TRUE", {
  skip_if(interactive(), "This test targets the non-interactive path specifically")
  expect_no_error(open_breeder_guide(open = TRUE))
})
