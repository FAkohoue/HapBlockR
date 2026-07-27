## tests/testthat/test-ocs.R
## -----------------------------------------------------------------------------
## Tests for R/ocs.R:
##   select_parents_ocs()      -- engine dispatch across THREE engines:
##                                 "alphamate" (true OCS, input-validation-only
##                                 here -- no AlphaMate binary available in
##                                 CI/test env), "optisel" (true OCS via
##                                 optiSel::candes()/opticont()/matings()/
##                                 noffspring() -- NEW, see R/ocs.R), and
##                                 "simplemating" (SimpleMating::planCross()/
##                                 selectCrosses() cross prediction/selection
##                                 -- this is the engine PREVIOUSLY named
##                                 "optisel"; see the breaking-rename NEWS.md
##                                 entry -- do not confuse the two "optisel"
##                                 test sections below with each other)
##   .ocs_validate_inputs()    -- internal, input alignment/validation
##   .rescale_to_nrm()         -- internal, exact rescaling
##   .target_degree_to_n_keep() -- internal, target_degree -> candidate-count
##                                 direction for engine = "simplemating"
##                                 (must match AlphaMate's own TargetDegree
##                                 convention: 0 = max-gain/keep everyone,
##                                 90 = max-diversity/keep fewest -- see
##                                 NEWS.md, this was inverted before)
##
## engine = "alphamate" requires a separately-installed executable. No
## AlphaMate binary is bundled with or auto-detected by HapBlockR (see the
## "Providing the AlphaMate executable" roxygen section on
## select_parents_ocs()) -- alphamate_exe must always be supplied explicitly
## by the caller, so alphamate_exe = NULL (the default) always means
## AlphaMate is unavailable. Tests here therefore fall into four groups:
##   (a) input-validation checks -- always pass an explicit,
##       guaranteed-nonexistent alphamate_exe (or leave it NULL) so they do
##       not depend on any executable being present;
##   (b) engine = "auto" dispatch tests -- prove alphamate_exe resolution
##       gates the engine choice correctly (NULL/nonexistent -> "optisel",
##       the OTHER true-OCS engine, NOT "simplemating" -- "auto" never
##       silently substitutes cross prediction for OCS; a real, user-supplied
##       path -> "alphamate", checked by seeing validation proceed past the
##       alphamate_exe check to the out_dir check, without actually invoking
##       the (platform-specific) executable);
##   (c) engine = "simplemating" (SimpleMating::planCross()/selectCrosses(),
##       NOT true OCS), exercised end-to-end, guarded by
##       skip_if_not_installed() for SimpleMating and optiSel PLUS
##       skip_if_simplemating_too_old() (an installed-but-outdated
##       SimpleMating passes the former but not the latter -- see
##       tests/testthat/helper.R), following this package's existing
##       convention for optional Suggests dependencies (see e.g. the
##       GA-guarded tests in test-parent-selection.R);
##   (d) engine = "optisel" (true OCS via optiSel's own solver), guarded by
##       skip_if_not_installed("optiSel") PLUS skip_if_optisel_missing_fn()
##       for each of candes/opticont/matings/noffspring (this call has not
##       been verified against a real optiSel installation from this
##       development environment -- see R/ocs.R's "Verification status"
##       roxygen section).
## No test here actually launches AlphaMate.exe: doing so is platform-
## specific (Windows only); the wiring/dispatch logic is what's covered.
## -----------------------------------------------------------------------------

library(testthat)
library(HapBlockR)

# -- Shared fixtures ------------------------------------------------------------

.n_ocs <- 10L
.ocs_ids <- paste0("ind", seq_len(.n_ocs))
.ocs_merit <- setNames(rnorm(.n_ocs, 0, 1), .ocs_ids)

set.seed(3L)
.X_ocs <- matrix(rnorm(.n_ocs * 8), .n_ocs, 8)
.ocs_G <- tcrossprod(scale(.X_ocs, scale = FALSE)) / 8
diag(.ocs_G) <- diag(.ocs_G) + 0.5   # ensure positive-definite-ish, away from 0
dimnames(.ocs_G) <- list(.ocs_ids, .ocs_ids)

.ocs_family <- setNames(rep(c("FamA", "FamB"), each = 5L), .ocs_ids)

# ==============================================================================
# 1. .ocs_validate_inputs() (internal)
# ==============================================================================

test_that(".ocs_validate_inputs: errors on unnamed merit", {
  expect_error(HapBlockR:::.ocs_validate_inputs(unname(.ocs_merit), .ocs_G, NULL),
              "merit")
})

test_that(".ocs_validate_inputs: errors on undimnamed G", {
  G_bad <- unname(.ocs_G)
  expect_error(HapBlockR:::.ocs_validate_inputs(.ocs_merit, G_bad, NULL), "G")
})

test_that(".ocs_validate_inputs: errors when G rownames != colnames", {
  G_bad <- .ocs_G
  colnames(G_bad) <- rev(colnames(G_bad))
  expect_error(HapBlockR:::.ocs_validate_inputs(.ocs_merit, G_bad, NULL),
              "row and column names")
})

test_that(".ocs_validate_inputs: errors with fewer than 3 individuals in common", {
  merit_small <- .ocs_merit[1:2]
  expect_error(HapBlockR:::.ocs_validate_inputs(merit_small, .ocs_G, NULL),
              "Fewer than 3")
})

test_that(".ocs_validate_inputs: aligns merit, G, and family to the common ID set", {
  merit_sub <- .ocs_merit[1:8]
  res <- HapBlockR:::.ocs_validate_inputs(merit_sub, .ocs_G, .ocs_family)
  expect_equal(length(res$ids), 8L)
  expect_equal(sort(res$ids), sort(names(merit_sub)))
  expect_equal(dim(res$G), c(8L, 8L))
  expect_equal(names(res$family), res$ids)
})

test_that(".ocs_validate_inputs: messages when candidates are dropped for not being in both merit and G", {
  merit_extra <- c(.ocs_merit, setNames(1, "ind_not_in_G"))
  expect_message(HapBlockR:::.ocs_validate_inputs(merit_extra, .ocs_G, NULL),
                 "excluded")
})

# ==============================================================================
# 2. .rescale_to_nrm() (internal, exact)
# ==============================================================================

test_that(".rescale_to_nrm: rescales to mean diagonal 1", {
  A <- HapBlockR:::.rescale_to_nrm(.ocs_G)
  expect_equal(mean(diag(A)), 1, tolerance = 1e-8)
})

test_that(".rescale_to_nrm: result is exactly symmetric", {
  A <- HapBlockR:::.rescale_to_nrm(.ocs_G)
  expect_equal(A, t(A))
})

test_that(".rescale_to_nrm: scaling is a simple proportional transform", {
  d <- mean(diag(.ocs_G))
  A <- HapBlockR:::.rescale_to_nrm(.ocs_G)
  expect_equal(A, (.ocs_G / d + t(.ocs_G / d)) / 2, tolerance = 1e-10)
})

test_that(".rescale_to_nrm: errors on zero mean diagonal", {
  Z <- matrix(0, 3, 3)
  expect_error(HapBlockR:::.rescale_to_nrm(Z), "zero mean diagonal")
})

# ==============================================================================
# 3. select_parents_ocs(): engine = "alphamate" input validation
# ==============================================================================

# A guaranteed-nonexistent path, used to force the "no valid executable"
# path in tests that want to be explicit about it (equivalent in effect to
# leaving alphamate_exe = NULL, since no binary is ever bundled/detected).
.no_such_exe <- file.path(tempdir(), "definitely_not_a_real_alphamate.exe")

test_that("select_parents_ocs (alphamate): errors when alphamate_exe does not point to an existing file", {
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "alphamate",
                       n_crosses = 5L, out_dir = tempdir(),
                       alphamate_exe = .no_such_exe, verbose = FALSE),
    "alphamate_exe"
  )
})

test_that("select_parents_ocs (alphamate): errors when alphamate_exe is left NULL (no bundled/auto-detected binary)", {
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "alphamate",
                       n_crosses = 5L, out_dir = tempdir(),
                       alphamate_exe = NULL, verbose = FALSE),
    "alphamate_exe"
  )
})

test_that("select_parents_ocs (alphamate): errors when out_dir is missing", {
  fake_exe <- tempfile()
  # Real "MZ" PE header bytes -- the platform warning below is now gated on
  # .is_windows_pe_exe() actually detecting this magic number (see R/ocs.R),
  # not merely on being off Windows, so the fake exe must look like a real
  # Windows binary for these two tests to still exercise that warning path.
  writeBin(as.raw(c(0x4d, 0x5a, 0x90, 0x00)), fake_exe)
  Sys.chmod(fake_exe, "0755")
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "alphamate",
                       n_crosses = 5L, alphamate_exe = fake_exe, verbose = FALSE),
    "out_dir"
  )
})

test_that("select_parents_ocs (auto): falls back to 'optisel' when no valid alphamate_exe resolves", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  res <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "auto",
                            n_crosses = 5L, target_degree = 30,
                            alphamate_exe = .no_such_exe, verbose = FALSE)
  expect_equal(res$engine_used, "optisel")
})

test_that("select_parents_ocs (auto): falls back to 'optisel' when alphamate_exe is left NULL", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  # No binary is ever bundled with or auto-detected by HapBlockR, so the
  # default alphamate_exe = NULL must always resolve to "optisel" (the
  # true-OCS engine, NOT "simplemating") here, regardless of platform.
  res <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "auto",
                            n_crosses = 5L, target_degree = 30,
                            verbose = FALSE)
  expect_equal(res$engine_used, "optisel")
})

# -- engine = "auto" with a user-supplied alphamate_exe --------------------------
#
# These tests prove alphamate_exe resolution correctly gates engine = "auto"
# WITHOUT actually launching the executable (platform-specific and not
# meaningfully testable without a real AlphaMate binary present): validation
# proceeding past the alphamate_exe check to the out_dir check shows "auto"
# picked "alphamate" for a real, user-supplied path.

test_that("select_parents_ocs (auto): selects 'alphamate' when the caller supplies a working alphamate_exe", {
  fake_exe <- tempfile()
  # Real "MZ" PE header bytes so the platform-mismatch warning (if any)
  # fires deterministically; on Windows this would be a genuine PE binary.
  writeBin(as.raw(c(0x4d, 0x5a, 0x90, 0x00)), fake_exe)
  Sys.chmod(fake_exe, "0755")
  expect_error(
    suppressWarnings(
      select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "auto",
                         n_crosses = 5L, alphamate_exe = fake_exe,
                         verbose = FALSE)
    ),
    "out_dir"
  )
})

test_that("select_parents_ocs (auto): a user-supplied alphamate_exe is eligible for auto off Windows", {
  skip_if(.Platform$OS.type == "windows", "this test targets the non-Windows override path")
  fake_exe <- tempfile()
  # Real "MZ" PE header bytes -- the platform warning below is now gated on
  # .is_windows_pe_exe() actually detecting this magic number (see R/ocs.R),
  # not merely on being off Windows, so the fake exe must look like a real
  # Windows binary for these two tests to still exercise that warning path.
  writeBin(as.raw(c(0x4d, 0x5a, 0x90, 0x00)), fake_exe)
  Sys.chmod(fake_exe, "0755")
  # Explicitly supplying alphamate_exe (even off Windows -- e.g. the user's
  # own native build, or a Wine wrapper script) opts back into engine =
  # 'alphamate' for auto-detection; validation proceeds past the platform
  # warning and fails when the (non-functional) fake exe is actually
  # launched, proving auto picked 'alphamate' rather than silently falling
  # back to 'optisel'.
  expect_warning(
    tryCatch(
      select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "auto",
                         n_crosses = 5L, alphamate_exe = fake_exe,
                         out_dir = tempdir(), verbose = FALSE),
      error = function(e) NULL
    ),
    "Windows binary"
  )
})

test_that("select_parents_ocs (alphamate): warns off Windows even when a valid exe path is supplied", {
  skip_if(.Platform$OS.type == "windows", "this test targets the non-Windows warning path")
  fake_exe <- tempfile()
  # Real "MZ" PE header bytes -- the platform warning below is now gated on
  # .is_windows_pe_exe() actually detecting this magic number (see R/ocs.R),
  # not merely on being off Windows, so the fake exe must look like a real
  # Windows binary for these two tests to still exercise that warning path.
  writeBin(as.raw(c(0x4d, 0x5a, 0x90, 0x00)), fake_exe)
  Sys.chmod(fake_exe, "0755")
  # The eventual attempt to actually launch fake_exe will error (it is not a
  # real executable) -- that error is swallowed here via tryCatch so this
  # test isolates and checks only the platform warning, which fires before
  # the launch attempt.
  expect_warning(
    tryCatch(
      select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "alphamate",
                         n_crosses = 5L, alphamate_exe = fake_exe,
                         out_dir = tempdir(), verbose = FALSE),
      error = function(e) NULL
    ),
    "Windows binary"
  )
})

test_that("select_parents_ocs (alphamate): does NOT warn off Windows for a genuine non-PE exe", {
  skip_if(.Platform$OS.type == "windows", "this test targets the non-Windows warning path")
  fake_exe <- tempfile()
  # No "MZ" header -- e.g. a Linux ELF binary would start with 0x7f 'E' 'L'
  # 'F' instead. This stands in for a real native alphamate_exe build (like
  # AlphaGenes' own Linux/macOS AlphaMate release): the platform warning must
  # not fire just because we're off Windows -- only when the file is
  # confirmed to actually be a Windows PE binary. See .is_windows_pe_exe()
  # in R/ocs.R.
  writeBin(as.raw(c(0x7f, 0x45, 0x4c, 0x46)), fake_exe)
  Sys.chmod(fake_exe, "0755")
  # fake_exe is only 4 magic bytes, not a complete/runnable ELF binary, so
  # the system2() call inside .run_alphamate_ocs() still fails to launch it
  # and R itself emits its own unrelated "had status 126" warning -- that is
  # expected noise, not what this test checks. Only the platform-mismatch
  # ("Windows binary") warning from .is_windows_pe_exe() must be absent, so
  # warnings are captured and inspected directly rather than using
  # expect_no_warning(), which would (wrongly) fail on the status-126 one.
  warnings_seen <- character(0)
  withCallingHandlers(
    tryCatch(
      select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "alphamate",
                         n_crosses = 5L, alphamate_exe = fake_exe,
                         out_dir = tempdir(), verbose = FALSE),
      error = function(e) NULL
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("Windows binary", warnings_seen, fixed = TRUE)))
})

# ==============================================================================
# 4. select_parents_ocs(): engine = "simplemating"
#    (SimpleMating::planCross()/selectCrosses() -- NOT true OCS; this is the
#    engine that used to be named "optisel" before the breaking rename. See
#    NEWS.md and R/ocs.R's file header.)
# ==============================================================================

test_that("select_parents_ocs (simplemating): requires SimpleMating and errors clearly when absent", {
  skip_if(requireNamespace("SimpleMating", quietly = TRUE),
         "SimpleMating is installed; this test targets the absent-dependency error path only")
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "simplemating",
                       n_crosses = 5L, verbose = FALSE),
    "SimpleMating"
  )
})

test_that("select_parents_ocs (simplemating): runs and returns the documented structure", {
  skip_if_not_installed("SimpleMating")
  skip_if_not_installed("optiSel")
  skip_if_simplemating_too_old("selectCrosses")
  res <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, family = .ocs_family,
                            engine = "simplemating", n_crosses = 5L,
                            target_degree = 30, allow_selfing = FALSE,
                            verbose = FALSE)
  expect_true(all(c("mating_plan", "contributors", "engine_used", "ok") %in% names(res)))
  expect_equal(res$engine_used, "simplemating")
  expect_true(is.logical(res$ok))
  if (!is.null(res$mating_plan) && nrow(res$mating_plan)) {
    expect_true(all(c("parent1", "parent2", "mean_relationship") %in% names(res$mating_plan)))
    expect_true(all(res$mating_plan$parent1 %in% .ocs_ids))
    expect_true(all(res$mating_plan$parent2 %in% .ocs_ids))
  }
  if (!is.null(res$contributors)) {
    expect_true(all(c("id", "contribution") %in% names(res$contributors)))
  }
})

test_that("select_parents_ocs (simplemating): errors on target_degree outside [0, 90]", {
  skip_if_not_installed("SimpleMating")
  skip_if_not_installed("optiSel")
  skip_if_simplemating_too_old("selectCrosses")
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "simplemating",
                       n_crosses = 5L, target_degree = 120, verbose = FALSE),
    "target_degree"
  )
})

test_that("select_parents_ocs (simplemating): runs at target_degree boundaries 0 and 90 without erroring", {
  skip_if_not_installed("SimpleMating")
  skip_if_not_installed("optiSel")
  skip_if_simplemating_too_old("selectCrosses")
  # Regression test: culling.pairwise.k must be nudged just past the
  # candidate K range's min/max (selectCrosses() keeps rows with
  # data$K < culling.pairwise.k, a strict "<"), or target_degree = 90 sets
  # the cutoff to exactly min(K) -- filtering out every candidate cross and
  # making selectCrosses() error -- and target_degree = 0 would still
  # exclude the maximum-K pair despite being documented as "most
  # permissive, i.e. closest to no culling at all" (direction matches
  # AlphaMate's own TargetDegree convention: 0 = max-gain/no restriction,
  # 90 = max-diversity/most restrictive -- see NEWS.md).
  res0 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "simplemating",
                             n_crosses = 5L, target_degree = 0, verbose = FALSE)
  expect_true(is.logical(res0$ok))
  res90 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "simplemating",
                              n_crosses = 5L, target_degree = 90, verbose = FALSE)
  expect_true(is.logical(res90$ok))
})

# -- .target_degree_to_n_keep(): direction must match AlphaMate's own
# TargetDegree convention -- a pure, deterministic unit test with no
# SimpleMating/optiSel dependency, added specifically because the direction
# was inverted in an earlier version of this file's formula (target_degree
# = 0 was wrongly the max-diversity end instead of the max-gain end -- see
# NEWS.md) and no test at the time actually checked the DIRECTION, only
# that the boundaries ran without erroring.

test_that(".target_degree_to_n_keep(): 0 keeps essentially all candidates (max-gain end)", {
  expect_equal(HapBlockR:::.target_degree_to_n_keep(0, n_cand = 100, min_keep = 20), 100)
})

test_that(".target_degree_to_n_keep(): 90 keeps only the min_keep floor (max-diversity end)", {
  expect_equal(HapBlockR:::.target_degree_to_n_keep(90, n_cand = 100, min_keep = 20), 20)
})

test_that(".target_degree_to_n_keep(): monotonically non-increasing as target_degree rises", {
  n_keep <- vapply(seq(0, 90, by = 5),
                   function(td) HapBlockR:::.target_degree_to_n_keep(td, n_cand = 100, min_keep = 20),
                   integer(1))
  expect_true(all(diff(n_keep) <= 0))
})

test_that(".target_degree_to_n_keep(): a moderate target_degree keeps a middling count, not a boundary value", {
  mid <- HapBlockR:::.target_degree_to_n_keep(45, n_cand = 100, min_keep = 20)
  expect_equal(mid, 50L)   # ceiling((90 - 45) / 90 * 100) = ceiling(50) = 50
  expect_true(mid > 20L && mid < 100L)
})

test_that(".target_degree_to_n_keep(): never returns fewer than min_keep or more than n_cand", {
  for (td in c(0, 30, 60, 90)) {
    n_keep <- HapBlockR:::.target_degree_to_n_keep(td, n_cand = 50, min_keep = 15)
    expect_true(n_keep >= 15L && n_keep <= 50L)
  }
})

test_that("select_parents_ocs (simplemating): honours max_contrib_per_parent (unlike n_parents_max, which is ignored with a message)", {
  skip_if_not_installed("SimpleMating")
  skip_if_not_installed("optiSel")
  skip_if_simplemating_too_old("selectCrosses")
  # n_parents_max has no selectCrosses() equivalent and is always ignored;
  # max_contrib_per_parent, unlike under the old GOCS()-based engine, IS
  # honoured directly (selectCrosses()'s own max.cross argument) -- so no
  # message fires for it specifically.
  expect_message(
    res <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "simplemating",
                              n_crosses = 5L, target_degree = 30,
                              max_contrib_per_parent = 2L, n_parents_max = 6L,
                              verbose = TRUE),
    "not directly enforced"
  )
  if (!is.null(res$contributors) && nrow(res$contributors))
    expect_true(max(res$contributors$contribution) <= 2L)
})

# ==============================================================================
# 5. select_parents_ocs(): engine = "optisel" -- TRUE OCS via optiSel's own
#    solver (optiSel::candes()/opticont()/matings()/noffspring()). NEW as of
#    the engine rename -- see R/ocs.R's "Verification status" roxygen
#    section for exactly what is/isn't independently verified here. Do not
#    confuse this section with section 4 above, which is the DIFFERENT,
#    SimpleMating-based engine now named "simplemating".
# ==============================================================================

test_that("select_parents_ocs (optisel): requires optiSel and errors clearly when absent", {
  skip_if(requireNamespace("optiSel", quietly = TRUE),
         "optiSel is installed; this test targets the absent-dependency error path only")
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                       n_crosses = 5L, verbose = FALSE),
    "optiSel"
  )
})

test_that("select_parents_ocs (optisel): does NOT require SimpleMating", {
  # The whole point of the rename/reimplementation: engine = "optisel" calls
  # optiSel's own solver directly and has no SimpleMating dependency at all,
  # unlike engine = "simplemating". This is checked structurally (the
  # dispatch never touches SimpleMating for this engine) rather than by
  # actually uninstalling SimpleMating in CI.
  expect_false(any(grepl("SimpleMating",
                        deparse(HapBlockR:::.run_optisel_ocs), fixed = TRUE)))
})

test_that("select_parents_ocs (optisel): errors on target_degree outside [0, 90]", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  expect_error(
    select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                       n_crosses = 5L, target_degree = 120, verbose = FALSE),
    "target_degree"
  )
})

test_that("select_parents_ocs (optisel): runs and returns the documented structure", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  res <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, family = .ocs_family,
                            engine = "optisel", n_crosses = 4L,
                            target_degree = 30, allow_selfing = FALSE,
                            verbose = FALSE)
  expect_true(all(c("mating_plan", "contributors", "engine_used", "ok") %in% names(res)))
  expect_equal(res$engine_used, "optisel")
  expect_true(res$ok)
  expect_equal(nrow(res$mating_plan), 4L)
  expect_true(all(unlist(res$constraint_checks)))
  expect_true(all(c("kin_min", "kin_max_gain") %in% names(res$frontier)))
  # kin_min must not exceed kin_max_gain (the max-diversity end of the
  # frontier cannot have a HIGHER mean kinship than the max-gain end) --
  # a basic sanity check on the two opticont()-solved reference points this
  # engine's target_degree interpolates between.
  expect_true(res$frontier$kin_min <= res$frontier$kin_max_gain + 1e-6)
  expect_true(all(c("parent1", "parent2", "mean_relationship") %in% names(res$mating_plan)))
  expect_true(all(res$mating_plan$parent1 %in% .ocs_ids))
  expect_true(all(res$mating_plan$parent2 %in% .ocs_ids))
  expect_true(all(res$mating_plan$parent1 != res$mating_plan$parent2))
  expect_true(all(c("id", "contribution") %in% names(res$contributors)))
  expect_equal(sum(res$contributors$contribution), 1, tolerance = 1e-10)
})

test_that("select_parents_ocs (optisel): runs at target_degree boundaries 0 and 90 without erroring", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  res0 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                             n_crosses = 4L, target_degree = 0, verbose = FALSE)
  expect_true(res0$ok)
  res90 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                              n_crosses = 4L, target_degree = 90, verbose = FALSE)
  expect_true(res90$ok)
})

test_that("select_parents_ocs (optisel): higher target_degree does not increase realised mean kinship", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  # Direction regression test, mirroring .target_degree_to_n_keep()'s own
  # direction tests for the "simplemating" engine: target_degree = 90
  # (max-diversity end) must not produce a HIGHER contribution-weighted mean
  # kinship than target_degree = 0 (max-gain end) for the same candidate
  # set -- this is the exact class of bug fixed in NEWS.md for the other
  # engine, so it is checked here too rather than assumed from the formula
  # alone.
  res0 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                             n_crosses = 4L, target_degree = 0, verbose = FALSE)
  res30 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                              n_crosses = 4L, target_degree = 30, verbose = FALSE)
  res90 <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                              n_crosses = 4L, target_degree = 90, verbose = FALSE)
  # The kinship ceiling actually solved under is non-increasing as
  # target_degree rises, for the same candidate set's frontier (same
  # kin_max_gain/kin_min each time, since they don't depend on
  # target_degree) -- the direct regression check for the direction this
  # engine's target_degree must follow.
  expect_true(res0$frontier$kin_ceiling >= res30$frontier$kin_ceiling - 1e-6)
  expect_true(res30$frontier$kin_ceiling >= res90$frontier$kin_ceiling - 1e-6)
  expect_equal(res0$frontier$kin_ceiling, res0$frontier$kin_max_gain, tolerance = 1e-6)
  expect_equal(res90$frontier$kin_ceiling, res90$frontier$kin_min, tolerance = 1e-6)
})

test_that("select_parents_ocs (optisel): n_parents_max re-solves and is enforced", {
  skip_if_not_installed("optiSel")
  skip_if_optisel_missing_fn("candes")
  skip_if_optisel_missing_fn("opticont")
  skip_if_optisel_missing_fn("matings")
  skip_if_optisel_missing_fn("noffspring")
  # Unlike engine = "simplemating" (where n_parents_max is unsupported),
  # engine = "optisel" retains the leading contributors and re-solves the
  # OCS problem on that restricted candidate set.
  expect_message(
    res <- select_parents_ocs(merit = .ocs_merit, G = .ocs_G, engine = "optisel",
                              n_crosses = 3L, target_degree = 30,
                              n_parents_max = 3L, verbose = TRUE),
    "re-solved"
  )
  expect_true(res$ok)
  expect_lte(nrow(res$contributors), 3L)
  expect_lte(length(unique(c(res$mating_plan$parent1,
                             res$mating_plan$parent2))), 3L)
  expect_true(all(unlist(res$constraint_checks)))
})

test_that("select_parents_ocs (optisel): infeasible hard constraints error", {
  skip_if_not_installed("optiSel")
  skip_if_not_installed("lpSolve")
  expect_error(
    select_parents_ocs(
      merit = .ocs_merit, G = .ocs_G, engine = "optisel",
      n_crosses = 4L, target_degree = 30, n_parents_max = 3L,
      allow_selfing = FALSE, allow_repeated_matings = FALSE,
      verbose = FALSE
    ),
    "No feasible mating plan"
  )
})
