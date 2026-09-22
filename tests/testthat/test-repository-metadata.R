.repository_source_root <- function() {
  candidates <- unique(c(
    testthat::test_path("..", ".."),
    getwd(),
    file.path(getwd(), "..")
  ))
  candidates <- normalizePath(candidates, mustWork = FALSE)
  is_checkout <- file.exists(file.path(candidates, "DESCRIPTION")) &
    file.exists(file.path(candidates, "_pkgdown.yml")) &
    file.exists(file.path(candidates, ".git"))
  matches <- candidates[is_checkout]
  if (!length(matches))
    testthat::skip("Repository-layout assertions run only in a Git checkout.")
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

test_that("pkgdown exposes every breeder-facing navigation tab", {
  root <- .repository_source_root()
  config <- paste(
    readLines(file.path(root, "_pkgdown.yml"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(
    config,
    "left:  [reference, tutorial, articles, guide, news]",
    fixed = TRUE
  )
  expect_match(config, "text: Tutorial", fixed = TRUE)
  expect_match(config, "text: Vignettes", fixed = TRUE)
  expect_match(config, "text: Breeder guide", fixed = TRUE)
  expect_match(
    config,
    "href: articles/HapBlockR-breeder-guide.html",
    fixed = TRUE
  )

  article_targets <- c(
    "HapBlockR-intro", "HapBlockR-ld-metrics",
    "HapBlockR-breeding-decisions", "HapBlockR-large-scale",
    "HapBlockR-phasing", "HapBlockR-programme-operations",
    "HapBlockR-workflow", "HapBlockR-breeder-guide"
  )
  expect_true(all(file.exists(file.path(
    root, "vignettes", paste0(article_targets, ".Rmd")
  ))))

  expect_true(file.exists(file.path(
    root, "inst", "guide", "HapBlockR_Breeder_Guide.md"
  )))
  expect_true(file.exists(file.path(
    root, "inst", "extdata", "HapBlockR_Breeder_Guide.pdf"
  )))
  expect_true(file.exists(file.path(
    root, "inst", "extdata", "HapBlockR_Breeder_Guide.docx"
  )))

  workflow <- paste(
    readLines(
      file.path(root, ".github", "workflows", "pkgdown.yaml"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_match(workflow, "docs/breeder-guide.pdf", fixed = TRUE)
  expect_match(workflow, "docs/breeder-guide.docx", fixed = TRUE)
})

test_that("specialist workflows install bounded dependency sets", {
  root <- .repository_source_root()
  workflow_names <- c(
    "beagle-integration.yaml",
    "benchmark-regression.yaml",
    "compiled-valgrind.yaml"
  )
  workflows <- lapply(workflow_names, function(workflow_name) {
    paste(
      readLines(
        file.path(root, ".github", "workflows", workflow_name),
        warn = FALSE
      ),
      collapse = "\n"
    )
  })

  for (workflow in workflows) {
    expect_match(workflow, "dependencies: '\"hard\"'", fixed = TRUE)
    expect_false(grepl("needs: check", workflow, fixed = TRUE))
  }
  expect_match(workflows[[1L]], "SNPRelate", fixed = TRUE)
  expect_match(workflows[[1L]], "gdsfmt", fixed = TRUE)
  expect_match(
    workflows[[3L]],
    '_R_CHECK_FORCE_SUGGESTS_: "false"',
    fixed = TRUE
  )
})

test_that("CI dependency contracts separate core compatibility from integrations", {
  root <- .repository_source_root()
  read_workflow <- function(name) {
    paste(
      readLines(
        file.path(root, ".github", "workflows", name),
        warn = FALSE
      ),
      collapse = "\n"
    )
  }

  check <- read_workflow("R-CMD-check.yaml")
  pkgdown <- read_workflow("pkgdown.yaml")
  integrations <- read_workflow("optional-integrations.yaml")

  expect_match(check, 'r: "release",  dependencies: \'"all"\'', fixed = TRUE)
  expect_match(check, 'r: "devel",    dependencies: \'"hard"\'', fixed = TRUE)
  expect_match(check, 'r: "4.2",      dependencies: \'"hard"\'', fixed = TRUE)
  expect_match(
    check,
    '_R_CHECK_FORCE_SUGGESTS_: ${{ matrix.config.force_suggests }}',
    fixed = TRUE
  )
  expect_match(check, "any::GA", fixed = TRUE)
  expect_false(grepl("Reinstall data.table from source", check, fixed = TRUE))

  for (workflow in c(check, pkgdown)) {
    expect_match(workflow, '"SNPRelate", "gdsfmt", "vsn"', fixed = TRUE)
  }

  expect_match(integrations, "dependencies: '\"hard\"'", fixed = TRUE)
  expect_match(
    integrations,
    "github::FAkohoue/DesiredGainR",
    fixed = TRUE
  )
  expect_false(grepl("needs: check", integrations, fixed = TRUE))
})

test_that("the full-pipeline vignette guards optional exact-solver output", {
  root <- .repository_source_root()
  vignette <- paste(
    readLines(
      file.path(root, "vignettes", "HapBlockR-full-pipeline.Rmd"),
      warn = FALSE
    ),
    collapse = "\n"
  )

  expect_match(
    vignette,
    "if (have_lpsolve) validate(exact_res)",
    fixed = TRUE
  )
  expect_false(grepl("\nvalidate(exact_res)\n", vignette, fixed = TRUE))
})

test_that("workflows pin the archived optiSel source explicitly", {
  root <- .repository_source_root()
  read_workflow <- function(name) {
    paste(
      readLines(
        file.path(root, ".github", "workflows", name),
        warn = FALSE
      ),
      collapse = "\n"
    )
  }

  source <- paste0(
    "github::cran/optiSel@",
    "4317bc1a4468e23f8ceca33b54b308302f984530"
  )
  check <- read_workflow("R-CMD-check.yaml")
  pkgdown <- read_workflow("pkgdown.yaml")
  integrations <- read_workflow("optional-integrations.yaml")

  source_count <- lengths(regmatches(
    check,
    gregexpr(source, check, fixed = TRUE)
  ))
  expect_identical(source_count, 5L)
  expect_match(
    check,
    '${{ matrix.config.optisel_source }}',
    fixed = TRUE
  )
  expect_match(pkgdown, source, fixed = TRUE)
  expect_match(integrations, source, fixed = TRUE)
  expect_false(grepl("any::optiSel", integrations, fixed = TRUE))
})
