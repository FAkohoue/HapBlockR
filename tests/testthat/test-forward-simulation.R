## tests/testthat/test-forward-simulation.R
## -----------------------------------------------------------------------------
## Tests for R/forward_simulation.R (HapSelect gap-analysis, part 3):
##   ga_vs_ts_simulation(), plot_ga_vs_ts_simulation()
##   .gs_normalize_schemes()    -- internal, pure R, no genomicSimulation
##                                  dependency
##
## SECTION 0 below (.gs_normalize_schemes()) runs UNCONDITIONALLY -- it is a
## pure-R helper with no call into genomicSimulation, so it is placed BEFORE
## the skip_if_not_installed("genomicSimulation") line further down and is
## exercised in every environment, including this one (no genomicSimulation
## available). This is the primary hand-verifiable coverage of the OCS/UC-
## informed rapid-cycling `schemes` argument's validation logic.
##
## Everything from the "Shared fixture" section onward wraps the external
## `genomicSimulation` package (not on CRAN -- see ?ga_vs_ts_simulation for
## install instructions). ga_vs_ts_simulation()'s own internal dependency
## guard fires before any of its own input validation, so that section is
## skipped unless genomicSimulation is installed -- there is no way to reach
## its validation-error tests without it. This mirrors test-epistasis.R's
## whole-file skip_if_not_installed("rrBLUP") pattern.
##
## IMPORTANT: the genomicSimulation-gated section is the single most
## important set of tests in this codebase to actually run once
## genomicSimulation is installed -- it is the only executable check on
## the wrapper's real API behaviour, including the newer mating_scheme =
## "ocs"/"uc" code paths added on top of the original GA/TS comparison.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

.fs_ids <- paste0("ind", 1:10)

# ==============================================================================
# 0. .gs_normalize_schemes() (internal, pure R) -- runs unconditionally,
#    without genomicSimulation. Covers the `schemes` argument that
#    generalises ga_vs_ts_simulation() to mating_scheme = "ocs"/"uc".
# ==============================================================================

test_that(".gs_normalize_schemes: NULL schemes builds the legacy GA/TS truncation list", {
  out <- HapBlockR:::.gs_normalize_schemes(
    schemes = NULL, ga_selected = .fs_ids[1:4], ts_selected = .fs_ids[5:8],
    ids = .fs_ids
  )
  expect_named(out, c("GA", "TS"))
  expect_equal(out$GA$founders, .fs_ids[1:4])
  expect_equal(out$GA$mating_scheme, "truncation")
  expect_equal(out$TS$founders, .fs_ids[5:8])
  expect_equal(out$TS$mating_scheme, "truncation")
})

test_that(".gs_normalize_schemes: errors when schemes is NULL and either legacy arg is missing", {
  expect_error(
    HapBlockR:::.gs_normalize_schemes(NULL, ga_selected = NULL,
                                      ts_selected = .fs_ids[5:8], ids = .fs_ids),
    "schemes"
  )
  expect_error(
    HapBlockR:::.gs_normalize_schemes(NULL, ga_selected = .fs_ids[1:4],
                                      ts_selected = NULL, ids = .fs_ids),
    "schemes"
  )
})

test_that(".gs_normalize_schemes: accepts a valid multi-scheme list, preserving extra args and defaulting mating_scheme", {
  schemes <- list(
    GA  = list(founders = .fs_ids[1:4], mating_scheme = "truncation"),
    OCS = list(founders = .fs_ids[3:6], mating_scheme = "ocs", n_crosses = 15L),
    UC  = list(founders = .fs_ids[5:8], mating_scheme = "uc",
              selected_proportion = 0.2),
    NoScheme = list(founders = .fs_ids[7:10])  # mating_scheme omitted
  )
  out <- HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids)
  expect_named(out, c("GA", "OCS", "UC", "NoScheme"))
  expect_equal(out$OCS$mating_scheme, "ocs")
  expect_equal(out$OCS$n_crosses, 15L)
  expect_equal(out$UC$mating_scheme, "uc")
  expect_equal(out$UC$selected_proportion, 0.2)
  expect_equal(out$NoScheme$mating_scheme, "truncation")  # default applied
})

test_that(".gs_normalize_schemes: schemes takes precedence over ga_selected/ts_selected, with a message", {
  schemes <- list(A = list(founders = .fs_ids[1:3], mating_scheme = "truncation"))
  expect_message(
    out <- HapBlockR:::.gs_normalize_schemes(schemes, ga_selected = .fs_ids[1:4],
                                             ts_selected = .fs_ids[5:8], ids = .fs_ids),
    "ignored"
  )
  expect_named(out, "A")
  expect_equal(out$A$founders, .fs_ids[1:3])
})

test_that(".gs_normalize_schemes: errors on an unnamed schemes list", {
  expect_error(
    HapBlockR:::.gs_normalize_schemes(
      list(list(founders = .fs_ids[1:3], mating_scheme = "truncation")),
      NULL, NULL, .fs_ids),
    "named list"
  )
})

test_that(".gs_normalize_schemes: errors on duplicated scheme names", {
  schemes <- list(A = list(founders = .fs_ids[1:3]))
  schemes <- c(schemes, list(A = list(founders = .fs_ids[4:6])))
  expect_error(
    HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids),
    "unique"
  )
})

test_that(".gs_normalize_schemes: errors when a scheme entry has no founders element", {
  schemes <- list(A = list(mating_scheme = "truncation"))
  expect_error(
    HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids),
    "founders"
  )
})

test_that(".gs_normalize_schemes: errors when founders are not a subset of ids", {
  schemes <- list(A = list(founders = c("ind1", "not_an_id")))
  expect_error(
    HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids),
    "not columns"
  )
})

test_that(".gs_normalize_schemes: errors when founders has fewer than 2 individuals", {
  schemes <- list(A = list(founders = "ind1"))
  expect_error(
    HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids),
    ">= 2"
  )
})

test_that(".gs_normalize_schemes: errors on an invalid mating_scheme value", {
  schemes <- list(A = list(founders = .fs_ids[1:3], mating_scheme = "bogus"))
  expect_error(
    HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids),
    "mating_scheme"
  )
})

test_that(".gs_normalize_schemes: all three valid mating_scheme values are accepted", {
  for (ms in c("truncation", "ocs", "uc")) {
    schemes <- list(A = list(founders = .fs_ids[1:3], mating_scheme = ms))
    out <- HapBlockR:::.gs_normalize_schemes(schemes, NULL, NULL, .fs_ids)
    expect_equal(out$A$mating_scheme, ms)
  }
})

skip_if_not_installed("genomicSimulation")

# -- Shared fixture -------------------------------------------------------------

make_sim_inputs <- function(n = 40, p = 20, seed = 11L) {
  phased   <- make_phased(n = n, p = p, seed = seed)
  snp_info <- make_snpinfo(p = p)
  set.seed(seed + 1L)
  snp_effects <- setNames(rnorm(p), snp_info$SNP)
  list(hap1 = phased$hap1, hap2 = phased$hap2, snp_info = snp_info,
      snp_effects = snp_effects)
}

# ==============================================================================
# Smoke test: full pipeline runs end to end
# ==============================================================================

test_that("ga_vs_ts_simulation(): runs end-to-end and returns the documented structure", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)

  sim <- ga_vs_ts_simulation(
    inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
    ga_selected = ids[1:5], ts_selected = ids[6:10],
    n_generations = 2L, pop_size = 10L, selection_intensity = 0.5,
    seed = 1L, verbose = FALSE
  )

  expect_type(sim, "list")
  expect_true(all(c("summary", "ga_final", "ts_final") %in% names(sim)))

  bs <- sim$summary
  expect_true(all(c("scheme", "generation", "mean_gebv", "max_gebv",
                    "sd_gebv", "n_pop") %in% names(bs)))
  # 2 schemes x (n_generations + 1) generations (0, 1, 2) = 6 rows
  expect_equal(nrow(bs), 6L)
  expect_setequal(unique(bs$scheme), c("GA", "TS"))
  expect_equal(sort(unique(bs$generation)), 0:2)

  for (final in list(sim$ga_final, sim$ts_final)) {
    expect_true(all(c("hap1", "hap2", "gebv") %in% names(final)))
    expect_true(is.matrix(final$hap1))
    expect_true(is.matrix(final$hap2))
    expect_equal(dim(final$hap1), dim(final$hap2))
    expect_equal(ncol(final$hap1), length(final$gebv))
    # Phased alleles must stay within {0, 1} -- if genomicSimulation's
    # internal storage/round-trip subtly corrupted phase, this would catch
    # out-of-range or non-integer values.
    expect_true(all(final$hap1 %in% c(0, 1)))
    expect_true(all(final$hap2 %in% c(0, 1)))
  }
})

test_that("ga_vs_ts_simulation(): selection_intensity = NULL mates every offspring generation in full", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)

  sim <- ga_vs_ts_simulation(
    inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
    ga_selected = ids[1:5], ts_selected = ids[6:10],
    n_generations = 2L, pop_size = 10L, selection_intensity = NULL,
    seed = 1L, verbose = FALSE
  )

  # n_pop should be pop_size for every generation > 0 under both schemes
  gens <- sim$summary[sim$summary$generation > 0L, ]
  expect_true(all(gens$n_pop == 10L))
})

# ==============================================================================
# plot_ga_vs_ts_simulation()
# ==============================================================================

test_that("plot_ga_vs_ts_simulation(): runs and returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  sim <- ga_vs_ts_simulation(
    inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
    ga_selected = ids[1:5], ts_selected = ids[6:10],
    n_generations = 2L, pop_size = 10L, seed = 1L, verbose = FALSE
  )
  p <- plot_ga_vs_ts_simulation(sim)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ga_vs_ts_simulation(): errors on a malformed input", {
  expect_error(plot_ga_vs_ts_simulation(list(nope = 1)), "ga_vs_ts_simulation")
})

# ==============================================================================
# Input validation (reached only because genomicSimulation IS installed --
# see file header)
# ==============================================================================

test_that("ga_vs_ts_simulation(): mismatched hap1/hap2 dimensions errors", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  bad_hap2 <- inp$hap2[, -1, drop = FALSE]
  expect_error(
    ga_vs_ts_simulation(inp$hap1, bad_hap2, inp$snp_info, inp$snp_effects,
                       ga_selected = ids[1:3], ts_selected = ids[4:6]),
    "identical dimensions"
  )
})

test_that("ga_vs_ts_simulation(): hap1 rows not matching snp_info errors", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  bad_snp_info <- inp$snp_info[1:5, ]
  expect_error(
    ga_vs_ts_simulation(inp$hap1, inp$hap2, bad_snp_info, inp$snp_effects,
                       ga_selected = ids[1:3], ts_selected = ids[4:6]),
    "nrow"
  )
})

test_that("ga_vs_ts_simulation(): ga_selected IDs not in hap1 columns errors", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  # Validated inside .gs_normalize_schemes() now (schemes[["GA"]] is built
  # internally from ga_selected under the legacy convention), hence "not
  # columns" rather than the literal substring "ga_selected".
  expect_error(
    ga_vs_ts_simulation(inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
                       ga_selected = c("nobody1", "nobody2"),
                       ts_selected = ids[4:6]),
    "not columns"
  )
})

test_that("ga_vs_ts_simulation(): fewer than 2 founders in either scheme errors", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  expect_error(
    ga_vs_ts_simulation(inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
                       ga_selected = ids[1], ts_selected = ids[4:6]),
    ">= 2"
  )
})

test_that("ga_vs_ts_simulation(): selection_intensity outside (0,1] errors", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  expect_error(
    ga_vs_ts_simulation(inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
                       ga_selected = ids[1:3], ts_selected = ids[4:6],
                       selection_intensity = 1.5),
    "selection_intensity"
  )
})

# ==============================================================================
# schemes / mating_scheme = "ocs"/"uc" extension (reached only because
# genomicSimulation IS installed -- see file header)
# ==============================================================================

test_that("ga_vs_ts_simulation(): a 'uc' scheme without blocks errors before any simulation runs", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  expect_error(
    ga_vs_ts_simulation(
      inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
      schemes = list(UC = list(founders = ids[1:4], mating_scheme = "uc")),
      n_generations = 1L, pop_size = 6L, verbose = FALSE
    ),
    "blocks"
  )
})

test_that("ga_vs_ts_simulation(): a 'truncation'-only schemes list runs and labels rows by scheme name", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  sim <- ga_vs_ts_simulation(
    inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
    schemes = list(
      Alpha = list(founders = ids[1:5], mating_scheme = "truncation"),
      Beta  = list(founders = ids[6:10], mating_scheme = "truncation")
    ),
    n_generations = 2L, pop_size = 10L, selection_intensity = 0.5,
    seed = 1L, verbose = FALSE
  )
  expect_setequal(unique(sim$summary$scheme), c("Alpha", "Beta"))
  expect_named(sim$final, c("Alpha", "Beta"))
  # Custom scheme names -> no legacy $ga_final/$ts_final top-level shortcut.
  expect_false("ga_final" %in% names(sim))
})

test_that("ga_vs_ts_simulation(): mating_scheme = 'ocs' runs end-to-end and produces a scoreable trajectory", {
  skip_if_not_installed("optiSel")
  skip_if_not_installed("SimpleMating")
  skip_if_simplemating_too_old("selectCrosses")
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  sim <- ga_vs_ts_simulation(
    inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
    schemes = list(
      TS  = list(founders = ids[1:10], mating_scheme = "truncation"),
      OCS = list(founders = ids[1:10], mating_scheme = "ocs", n_crosses = 8L)
    ),
    n_generations = 1L, pop_size = 8L, selection_intensity = 0.5,
    seed = 1L, verbose = FALSE
  )
  expect_setequal(unique(sim$summary$scheme), c("TS", "OCS"))
  ocs_final <- sim$summary[sim$summary$scheme == "OCS" &
                             sim$summary$generation == 1L, ]
  expect_equal(nrow(ocs_final), 1L)
  expect_true(is.finite(ocs_final$mean_gebv))
})

test_that("ga_vs_ts_simulation(): mating_scheme = 'uc' runs end-to-end given blocks", {
  inp <- make_sim_inputs()
  ids <- colnames(inp$hap1)
  blk <- make_blocks(inp$snp_info, n_blocks = 4L)
  sim <- ga_vs_ts_simulation(
    inp$hap1, inp$hap2, inp$snp_info, inp$snp_effects,
    schemes = list(
      TS = list(founders = ids[1:10], mating_scheme = "truncation"),
      UC = list(founders = ids[1:10], mating_scheme = "uc", n_crosses = 8L)
    ),
    blocks = blk, n_generations = 1L, pop_size = 8L,
    selection_intensity = 0.5, seed = 1L, verbose = FALSE
  )
  expect_setequal(unique(sim$summary$scheme), c("TS", "UC"))
  uc_final <- sim$summary[sim$summary$scheme == "UC" &
                            sim$summary$generation == 1L, ]
  expect_equal(nrow(uc_final), 1L)
  expect_true(is.finite(uc_final$mean_gebv))
})
