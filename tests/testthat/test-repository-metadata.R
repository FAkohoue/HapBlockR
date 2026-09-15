.repository_source_root <- function() {
  candidates <- unique(c(
    testthat::test_path("..", ".."),
    getwd(),
    file.path(getwd(), "..")
  ))
  candidates <- normalizePath(candidates, mustWork = FALSE)
  matches <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(matches))
    testthat::skip("Repository-source consistency check.")
  matches[1L]
}

test_that("package, citation, and guide versions are consistent", {
  root <- .repository_source_root()
  description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, ]
  version <- unname(description[["Version"]])
  citation <- utils::readCitationFile(
    file.path(root, "inst", "CITATION"),
    meta = as.list(description)
  )

  cff <- readLines(file.path(root, "CITATION.cff"), warn = FALSE)
  expect_true(any(trimws(cff) == paste("version:", version)))
  expect_s3_class(citation, "bibentry")
  expect_true(any(grepl(version, format(citation), fixed = TRUE)))

  guide <- readLines(
    file.path(root, "inst", "guide", "HapBlockR_Breeder_Guide.md"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  expect_true(any(grepl(
    paste0("Compatible package version: ", version),
    guide,
    fixed = TRUE
  )))
  expect_true(any(grepl("nine decision tools", guide, fixed = TRUE)))
})

test_that("source tree contains no redistributed executables or JAR archives", {
  root <- .repository_source_root()
  extdata <- file.path(root, "inst", "extdata")
  prohibited <- list.files(
    extdata,
    pattern = "\\.(jar|class|exe|dll|so|dylib)$",
    ignore.case = TRUE,
    full.names = TRUE
  )
  if (length(prohibited)) {
    testthat::fail(paste(
      c(
        paste("Repository root:", root),
        paste("Prohibited file:", prohibited)
      ),
      collapse = "\n"
    ))
  }
  expect_length(prohibited, 0L)
})

test_that("breeder guide has reproducible source and a generated PDF", {
  root <- .repository_source_root()
  expect_true(file.exists(file.path(
    root, "inst", "guide", "HapBlockR_Breeder_Guide.md"
  )))
  expect_true(file.exists(file.path(
    root, "tools", "build_breeder_guide_accessible.cjs"
  )))
  guide_pdf <- file.path(
    root, "inst", "extdata", "HapBlockR_Breeder_Guide.pdf"
  )
  expect_true(file.exists(guide_pdf))
  expect_gt(file.info(guide_pdf)$size, 20000)
})

test_that("package declares British English", {
  root <- .repository_source_root()
  description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, ]
  expect_identical(unname(description[["Language"]]), "en-GB")
  expect_true(file.exists(file.path(root, "inst", "WORDLIST")))
})
