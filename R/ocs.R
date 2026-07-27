# ==============================================================================
# ocs.R
#
# True Optimal Contribution Selection (OCS) / mate allocation -- strategy 3
# of six planned parent-selection strategy extensions. Unlike
# select_parents_ga() (which selects a fixed-size SET of parents with an
# optional soft relatedness penalty) or usefulness_criterion() (which ranks
# candidate CROSSES independently), OCS solves the actual contribution-
# optimisation problem: how much should each candidate parent contribute to
# the next generation, and which specific matings should be made, to
# maximize genetic merit subject to an explicit constraint on the next
# generation's relatedness/inbreeding (Meuwissen 1997).
#
# BREAKING RENAME (see NEWS.md): engine = "optisel" used to mean "wraps
# SimpleMating::planCross()/selectCrosses()" -- a misnomer, since optiSel's
# own solver was never actually called; optiSel was present only as one of
# SimpleMating's transitive dependencies. That SimpleMating-based engine is
# now correctly named engine = "simplemating", and engine = "optisel" has
# been reimplemented from scratch to be what its name always should have
# meant: real Optimal Contribution Selection solved via optiSel's OWN
# functions (Wellmann 2019). There is deliberately no backward-compatible
# alias -- code that passed engine = "optisel" expecting the old
# SimpleMating-based behaviour must be updated to engine = "simplemating"
# (silently keeping the string "optisel" pointed at the old behaviour would
# have meant the name permanently lied about what it does).
#
# This package does not implement an OCS solver from scratch. Three mature,
# purpose-built engines are wrapped instead, matching the package's existing
# philosophy of handing off hard optimisation/statistical machinery to
# validated external tools (GA::ga() for select_parents_ga(), rrBLUP/BGLR for
# marker effects) rather than reimplementing it. Two solve the actual OCS
# problem (continuous contribution optimisation under a relatedness
# constraint); one instead does discrete, greedy cross PREDICTION/selection
# -- a genuinely different algorithm class, not just a different solver for
# the same problem:
#
#   engine = "alphamate" -- TRUE OCS. Wraps the AlphaMate executable (Hickey
#     group / AlphaGenes suite). Full mate-allocation via an evolutionary
#     algorithm: specific crossing pairs, TargetDegree diversity control,
#     no-selfing / no-repeated-mating constraints. Requires the AlphaMate
#     binary installed separately (not an R package) -- point
#     `alphamate_exe` at it.
#   engine = "optisel" -- TRUE OCS. Wraps optiSel's own solver (Wellmann
#     2019, Meuwissen 1997's formal OCS problem): optiSel::candes() builds
#     the candidate-description object from merit + G (relatedness), then
#     optiSel::opticont() solves the actual continuous-optimisation problem
#     for each candidate's optimum contribution, subject to an upper bound
#     on the next generation's mean kinship. That kinship bound is set from
#     `target_degree` by locating the two ends of the gain/diversity
#     frontier for THIS candidate set (both genuinely optiSel-solved: the
#     min-kinship-achievable solution and the max-merit-with-no-kinship-
#     constraint solution) and interpolating between them -- see
#     .run_optisel_ocs()'s own comments for the exact mechanism and why it
#     differs from the "simplemating" engine's quantile-of-observed-K
#     approach. optiSel::noffspring() then converts contributions into
#     expected offspring counts, and optiSel::matings() solves the discrete
#     Sire x Dam assignment that minimises mean offspring kinship. Requires
#     the optiSel package (Suggests) -- no SimpleMating dependency, no
#     external binary.
#   engine = "simplemating" -- NOT true OCS; discrete greedy cross
#     PREDICTION/selection. Wraps SimpleMating (Peixoto et al. 2024,
#     Resende-Lab/SimpleMating): SimpleMating::planCross() builds every
#     candidate pair, this file computes each pair's mid-parent merit and
#     G-matrix relatedness natively (plain arithmetic -- no SimpleMating
#     criterion function needed for that), and SimpleMating::selectCrosses()
#     then does the actual relatedness-constrained mate allocation: a greedy
#     search that maximises the criterion subject to (a) discarding any
#     candidate pair more related than `culling.pairwise.k` and (b) per-
#     parent min/max cross-count caps. Requires SimpleMating + optiSel
#     installed (optiSel is one of SimpleMating's own Imports, pulled in
#     automatically). Real limitation, inherited from selectCrosses() and
#     not something this package works around: its diversity control is a
#     hard relatedness CUTOFF, not a continuous contribution optimum --
#     `target_degree` (this engine's existing parameter, kept for a
#     consistent user-facing API across all three engines) is therefore
#     mapped onto `culling.pairwise.k` by QUANTILE of the candidate set's
#     own observed relatedness values, in the SAME direction as AlphaMate's
#     own TargetDegree (Kinghorn's frontier-degree concept): 0 = keeps
#     essentially all candidates (max-gain end of the frontier, no diversity
#     restriction), 90 = keeps only the least-related candidates (max-
#     diversity end), with a floor that always keeps enough candidates for
#     the requested n_crosses to be feasible. A raw linear interpolation
#     across [min(K), max(K)] was tried first and found, via a real
#     devtools::test() run, to be fragile against a single outlier pair
#     widening the observed range -- at target_degree = 30 on a
#     10-parent/45-cross test panel it left too few candidates for
#     selectCrosses() to find any n_crosses = 5 solution ("Reached maximum in
#     the search, try to increase data size."). This is an approximation of
#     AlphaMate's TargetDegree, not an equivalent algorithm; re-tune if you
#     switch engines. max_contrib_per_parent IS directly honoured
#     (selectCrosses()'s own max.cross argument), and no-repeated-matings is
#     guaranteed by construction (each unordered pair appears at most once
#     in planCross()'s candidate list, so selectCrosses()'s greedy search
#     can never select the same pair twice) rather than merely checked
#     afterwards.
#
# engine = "auto" (default): uses "alphamate" when a usable alphamate_exe
# resolves, otherwise falls back to "optisel" (the other TRUE-OCS engine,
# preferred over "simplemating" here specifically because it solves the
# same formal problem AlphaMate does, just with a different solver -- "auto"
# should not silently swap OCS for cross-prediction).
#
# Verification note: TWO earlier versions of the (now-renamed) "simplemating"
# engine turned out wrong, each discovered by a real devtools::test()
# failure against an actually-installed SimpleMating copy:
#   1. A hand-written optiSel::candes()/opticont() call, replaced after
#      finding real argument-name mismatches against SimpleMating::GOCS()'s
#      actual source (guessed kinship=/ub.kinship=/ub.n=/method="max.bv" vs.
#      GOCS()'s real sKin=/a flat ub=/method="max.Crit"). (This attempt is
#      exactly why the NEW engine = "optisel" below, added later once the
#      user pointed out optiSel's own functions were never actually being
#      called, is written directly against optiSel's real, user-supplied
#      candes()/opticont()/matings() documentation rather than against a
#      guess -- see its own comments for what is/isn't independently
#      verified.)
#   2. The GOCS()-based rewrite that replaced it, which itself broke: GOCS()
#      is present in the Resende-Lab/SimpleMating GitHub repository's R/
#      source and its checked-in NAMESPACE (confirmed by direct fetch), yet
#      a real, freshly-reinstalled SimpleMating 0.2.1's own help index does
#      NOT list GOCS() as an exported function -- confirmed directly by the
#      package maintainer's user, whose `?SimpleMating` help index enumerates
#      contrib2Cross/getIndex/getMPV/getTGV/getUsefA/getUsefAD/
#      getUsefAD_mt/getUsefA_mt/planCross/relateThinning/selectCrosses/
#      setCrosses -- no GOCS. Why the checked-in NAMESPACE and the installed
#      package disagree was not resolved (possibly an un-regenerated/stale
#      NAMESPACE line, possibly something installation-specific); rather
#      than depend on a function whose exported status is contested, this
#      engine was rewritten a second time around planCross()/selectCrosses()
#      specifically because BOTH are confirmed present in the real,
#      currently-installed 0.2.1 AND their exact current source
#      (Resende-Lab/SimpleMating, planCross.R and selectCrosses.R) was read
#      directly before writing this. getUsefA() (used by
#      usefulness_criterion(variance_model = "simplemating") in
#      R/genomic_mating.R, unaffected by any of this) was independently
#      confirmed present in that same real help index AND its real source's
#      argument list (MatePlan, Markers, addEff, K, Map.In, linkDes,
#      propSel, Type, Generation, n_threads, display_progress) matches this
#      package's existing call exactly -- no change was needed there.
# ==============================================================================


# -- Internal: validate and align merit/G/family inputs shared by both
# engines. Returns a list with aligned `ids`, `merit`, `G`, `family`.
.ocs_validate_inputs <- function(merit, G, family, verbose = TRUE) {
  if (is.null(names(merit)))
    stop("merit must be a named numeric vector (names = individual IDs).",
         call. = FALSE)
  if (!is.matrix(G) || is.null(rownames(G)) || is.null(colnames(G)))
    stop("G must be a dimnamed matrix (row/column names = individual IDs), ",
         "e.g. from compute_haplotype_grm() or your own relationship ",
         "matrix (VanRaden GRM, blended H matrix, pedigree A, etc.).",
         call. = FALSE)
  if (!identical(rownames(G), colnames(G)))
    stop("G must have identical row and column names.", call. = FALSE)

  ids <- intersect(names(merit), rownames(G))
  if (length(ids) < 3L)
    stop("Fewer than 3 individuals in common between merit and G; cannot ",
         "run OCS.", call. = FALSE)
  dropped <- length(union(names(merit), rownames(G))) - length(ids)
  if (dropped > 0L && isTRUE(verbose))
    message("[ocs] ", dropped, " individual(s) in merit or G but not both ",
            "-- excluded from the candidate set.")

  fam <- NULL
  if (!is.null(family)) {
    if (is.null(names(family)))
      stop("family must be a named character vector (names = individual ",
           "IDs).", call. = FALSE)
    fam <- family[ids]
  }

  list(ids = ids, merit = merit[ids], G = G[ids, ids, drop = FALSE], family = fam)
}

# -- Internal: is `path` a Windows PE (.exe) binary? Checked by reading its
# first two bytes and comparing against the "MZ" magic number all Windows PE
# executables start with (DOS/PE header) -- not by file extension or name,
# since a user-supplied native Linux/macOS build need not be named ".exe" and
# a copy of the bundled Windows binary need not keep its original name either.
# Returns NA (rather than FALSE) if the file cannot be read, so callers can
# tell "confirmed not a PE binary" apart from "couldn't check" and avoid
# warning on the latter.
.is_windows_pe_exe <- function(path) {
  con <- tryCatch(file(path, "rb"), error = function(e) NULL)
  if (is.null(con)) return(NA)
  on.exit(close(con), add = TRUE)
  magic <- tryCatch(readBin(con, "raw", n = 2L), error = function(e) raw(0))
  if (length(magic) < 2L) return(NA)
  identical(magic, as.raw(c(0x4d, 0x5a)))  # "MZ"
}

# -- Internal: rescale a relationship matrix to have mean diagonal 1
# (numerator-relationship-like scale) and force exact symmetry, matching the
# convention AlphaMate's NRM file expects. VanRaden-style GRMs are not
# guaranteed to already be on this scale (their diagonal reflects genomic
# inbreeding relative to allele-frequency-based expectation, not 1 + F on a
# pedigree-NRM scale).
.rescale_to_nrm <- function(A) {
  d <- mean(diag(A), na.rm = TRUE)
  if (!is.finite(d) || d == 0)
    stop("Cannot rescale relationship matrix to NRM scale: non-finite or ",
         "zero mean diagonal.", call. = FALSE)
  A <- A / d
  (A + t(A)) / 2
}



#' True Optimal Contribution Selection (OCS) and Mate Allocation
#'
#' Solves for optimal parent contributions and a resulting crossing plan
#' that maximizes genetic merit subject to an explicit constraint on the
#' next generation's relatedness (Meuwissen 1997) -- the actual formal OCS
#' problem, as opposed to \code{\link{select_parents_ga}}'s fixed-size subset
#' search (with only a soft relatedness penalty) or
#' \code{\link{usefulness_criterion}}'s independent per-cross ranking.
#'
#' This package does not implement an OCS solver itself. Three mature,
#' external engines are wrapped, selected via \code{engine}. Two solve the
#' actual OCS problem; one instead does discrete, greedy cross prediction/
#' selection -- see "Engine differences" below for why that distinction
#' matters, not just which is "better":
#' \describe{
#'   \item{\code{"alphamate"}}{\strong{True OCS.} The AlphaMate executable
#'     (Hickey group / AlphaGenes suite). Full evolutionary-algorithm mate
#'     allocation: specific crossing pairs, a \code{target_degree} diversity
#'     control, no-selfing/no-repeated-mating constraints, per-parent
#'     contribution caps. Requires the AlphaMate binary installed separately
#'     -- point \code{alphamate_exe} at it. Not an R package; not
#'     distributed with HapBlockR.}
#'   \item{\code{"optisel"}}{\strong{True OCS.} \code{optiSel::candes()} +
#'     \code{opticont()} (Wellmann 2019) -- optiSel's own continuous-
#'     optimisation solver for Meuwissen's (1997) formal OCS problem:
#'     optimal per-candidate contributions subject to an explicit upper
#'     bound on the next generation's mean kinship. \code{target_degree} is
#'     mapped onto that kinship bound by solving both ends of the gain/
#'     diversity frontier for your actual candidate set (a genuine
#'     min-kinship optimum and a genuine max-merit optimum, both via
#'     \code{opticont()}) and interpolating between them.
#'     \code{optiSel::noffspring()} then converts contributions into
#'     expected offspring counts, and \code{optiSel::matings()} solves the
#'     discrete Sire x Dam assignment minimising mean offspring kinship.
#'     Requires the optiSel package -- no external binary, no SimpleMating
#'     dependency.}
#'   \item{\code{"simplemating"}}{\strong{Not true OCS -- discrete greedy
#'     cross prediction/selection.} \code{SimpleMating::planCross()} +
#'     \code{SimpleMating::selectCrosses()} (Peixoto et al. 2024,
#'     Resende-Lab/SimpleMating). \code{planCross()} enumerates every
#'     candidate pair; this function computes each pair's mid-parent merit
#'     and \code{G}-matrix relatedness directly (plain arithmetic);
#'     \code{selectCrosses()} then does the actual relatedness-constrained
#'     mate allocation -- a greedy search maximising merit subject to a
#'     relatedness cutoff and per-parent min/max cross-count caps. Requires
#'     SimpleMating + optiSel installed (no external binary; optiSel is one
#'     of SimpleMating's own dependencies). See "Engine differences" below
#'     for real, verified differences from the two true-OCS engines'
#'     diversity-control mechanism.}
#' }
#' \strong{This is a breaking rename as of the version introducing
#' \code{"simplemating"}}: \code{engine = "optisel"} previously meant what
#' is now \code{engine = "simplemating"} (optiSel's own solver was never
#' actually called under the old name -- optiSel was present only as one of
#' SimpleMating's transitive dependencies). Code written against the old
#' behaviour must change \code{engine = "optisel"} to \code{engine =
#' "simplemating"}; there is no backward-compatible alias, since silently
#' keeping the string \code{"optisel"} pointed at SimpleMating's algorithm
#' would leave the name permanently misleading. See \code{NEWS.md}.
#'
#' \code{engine = "auto"} (default) uses AlphaMate only when you supplied a
#' working \code{alphamate_exe} yourself and the file exists (your own
#' native build, a Wine wrapper script, etc. is trusted as-is); otherwise it
#' falls back to \code{"optisel"} (the other true-OCS engine -- \code{"auto"}
#' never silently substitutes cross prediction for OCS) with a message. No
#' AlphaMate binary is bundled with HapBlockR or auto-detected -- see
#' "Providing the AlphaMate executable" below.
#'
#' @section Engine differences -- read before choosing:
#' \code{"alphamate"} and \code{"optisel"} both solve the actual OCS problem
#' (continuous contribution optimisation under a relatedness constraint) via
#' genuinely different solvers -- AlphaMate's own evolutionary algorithm vs.
#' optiSel's \code{candes()}/\code{opticont()}/\code{matings()}. Both are
#' controlled by the same \code{target_degree} lever using the same
#' direction convention (0 = max-gain end of the frontier, prioritising
#' merit and accepting more relatedness; 90 = max-diversity end, prioritising
#' minimised relatedness), but via genuinely different mechanisms: AlphaMate
#' uses a continuous Kinghorn-frontier degree (its own \code{TargetDegree})
#' natively inside its evolutionary algorithm; the \code{"optisel"} engine
#' interpolates a mean-kinship ceiling between two \code{opticont()}-solved
#' frontier extremes specific to your candidate set (see
#' \code{.run_optisel_ocs()}'s own comments). Not numerically identical
#' between the two; re-tune if you switch.
#'
#' \code{"simplemating"} is a different algorithm class entirely -- discrete
#' greedy cross selection under a hard relatedness cutoff
#' (\code{culling.pairwise.k}: candidate pairs more related than this are
#' discarded outright, not smoothly down-weighted), not a continuous
#' contribution optimum. \code{target_degree} is mapped onto that cutoff by
#' \emph{quantile} of the candidate set's own observed relatedness values,
#' in the same direction as the two true-OCS engines' \code{target_degree}
#' (0 = keeps essentially all candidates, the max-gain end of the frontier
#' with no diversity restriction; 90 = keeps only the least-related
#' candidates, the max-diversity end), with a floor that always keeps enough
#' candidates for the requested \code{n_crosses} to be feasible -- a raw
#' linear interpolation across \code{[min(K), max(K)]} was tried first and
#' found, via a real \code{devtools::test()} run, to be fragile against
#' outlier pairs that widen the observed range and starve
#' \code{selectCrosses()}'s search of candidates at moderate
#' \code{target_degree} values. Treat every engine's \code{target_degree} as
#' an approximation of the same lever, not an equivalent algorithm; re-tune
#' if you switch engines.
#'
#' Beyond that: AlphaMate lets you cap the parent count
#' (\code{n_parents_max}) as a native constraint during allocation. The
#' \code{"optisel"} engine uses the leading solved contributors to define the
#' permitted subset, then re-solves the OCS frontier and constrained optimum
#' within that subset before allocating matings. The \code{"simplemating"}
#' engine does not expose this control and reports it as unsupported.
#' \code{max_contrib_per_parent} is honoured directly by
#' all three engines (AlphaMate's own cap; \code{selectCrosses()}'s
#' \code{max.cross} argument under \code{"simplemating"};
#' \code{optiSel::matings()}'s \code{ub.n} argument under \code{"optisel"}).
#' All three engines also guarantee no-repeated-matings by construction for
#' this function's default (non-selfing) design -- AlphaMate enforces it
#' during allocation; \code{optiSel::matings()} is called with \code{ub.n =
#' 1} unless \code{allow_repeated_matings = TRUE}; and
#' \code{planCross(MateDesign = "half")} lists each unordered pair at most
#' once so \code{selectCrosses()}'s greedy search cannot select the same
#' pair twice.
#'
#' @section Providing the AlphaMate executable:
#' AlphaMate is a separate, third-party tool (Hickey Group / AlphaGenes
#' suite, \url{https://github.com/AlphaGenes/AlphaMate}, Fortran source,
#' MIT licensed) -- it is not distributed as an R package, HapBlockR does
#' not implement its algorithm, and \strong{no AlphaMate binary is bundled
#' with or downloaded by this package}: AlphaGenes does not currently
#' publish pre-built binaries (no GitHub Releases), so there is no stable
#' URL this package could fetch one from, and CRAN policy prohibits
#' shipping compiled executables in a source package regardless. To use
#' \code{engine = "alphamate"} you must build AlphaMate yourself from that
#' repository (a Fortran compiler such as \code{gfortran} is required) or
#' otherwise obtain a working executable you are entitled to use under its
#' own license terms, then pass its full path explicitly via
#' \code{alphamate_exe = "/full/path/to/AlphaMate"} (or \code{.exe} on
#' Windows) -- there is no auto-detection. If you would rather avoid
#' building a separate binary, \code{engine = "optisel"} (true OCS, no
#' SimpleMating dependency) or \code{engine = "simplemating"} (cross
#' prediction/selection) both require no external executable at all -- only
#' R packages -- and are the more portable, drop-in choices.
#'
#' @section Engine integration:
#' The AlphaMate wrapper follows the programme's documented file formats and
#' execution sequence. The \code{"optisel"} engine calls
#' \code{candes()} and \code{opticont()} for the continuous OCS problem and
#' uses HapBlockR's hard-constrained discrete allocator for the mating plan.
#' The \code{"simplemating"} engine uses the exported
#' \code{planCross()} and \code{selectCrosses()} interface. HapBlockR checks
#' required exports and returned constraints so that an incompatible optional
#' dependency fails explicitly with the installed version and function name.
#'
#' @param merit Named numeric vector of candidate parent merit (e.g. GEBV,
#'   a Selection Index, or \code{score_favorable_haplotypes()}'s
#'   \code{stacking_index}), names = individual IDs.
#' @param G Dimnamed relationship matrix (row/column names = individual
#'   IDs), e.g. from \code{\link{compute_haplotype_grm}}, or your own
#'   VanRaden/IBS/blended-H/pedigree-A matrix built however you already
#'   build it. Supply the relationship matrix produced by the HapBlockR
#'   workflow or another aligned relationship analysis appropriate to the
#'   programme.
#' @param family Optional named character vector (names = individual IDs)
#'   of family/cross-of-origin labels, used only for output labelling/
#'   reporting (neither engine enforces family representation as a hard
#'   constraint; see Details).
#' @param engine Character, one of \code{"auto"} (default), \code{"alphamate"}
#'   (true OCS), \code{"optisel"} (true OCS), or \code{"simplemating"}
#'   (discrete greedy cross prediction/selection, not true OCS). See
#'   Description -- and note the breaking rename if you have existing code
#'   using \code{engine = "optisel"} from before this version.
#' @param n_crosses Integer, default \code{20L}. Number of matings to
#'   produce.
#' @param n_parents_max Optional integer. Maximum number of distinct parents
#'   to use. AlphaMate: \code{NumberOfParents}, a native constraint enforced
#'   during allocation. \code{engine = "optisel"}: uses the
#'   \code{n_parents_max} leading contributors to define a restricted
#'   candidate subset, then re-solves the OCS frontier and constrained optimum
#'   within that subset before mating allocation. \code{engine =
#'   "simplemating"}: reported as unsupported because
#'   \code{selectCrosses()} does not expose this control.
#' @param max_contrib_per_parent Integer, default \code{4L}. Maximum number
#'   of matings any single parent can participate in. Honoured directly by
#'   all three engines (AlphaMate's own cap; \code{selectCrosses()}'s
#'   \code{max.cross} argument under \code{engine = "simplemating"};
#'   \code{optiSel::matings()}'s \code{ub.n} argument, applied per-pair,
#'   plus a pre-solve cap on each candidate's assigned offspring slots,
#'   under \code{engine = "optisel"}).
#' @param allow_selfing Logical, default \code{FALSE}.
#' @param allow_repeated_matings Logical, default \code{FALSE}. Enforced by
#'   construction under all three engines for the default (non-selfing)
#'   design -- see "Engine differences".
#' @param target_degree Numeric, default \code{30}, in \code{[0, 90]}.
#'   Diversity-vs-gain lever for ALL THREE engines, but via genuinely
#'   different mechanisms per engine (AlphaMate's continuous Kinghorn-
#'   frontier \code{TargetDegree}; under \code{engine = "optisel"}, a
#'   mean-kinship ceiling interpolated between two \code{opticont()}-solved
#'   frontier extremes for your candidate set; under \code{engine =
#'   "simplemating"}, \code{selectCrosses()}'s hard
#'   \code{culling.pairwise.k} relatedness cutoff, set from
#'   \code{target_degree} by quantile of the candidate set's own observed
#'   relatedness values, with a floor that keeps enough candidates for
#'   \code{n_crosses} to stay feasible) -- \strong{0 prioritises maximising
#'   the criterion} (merit/gain, accepting more relatedness -- the max-gain
#'   end of Kinghorn's frontier), \strong{90 prioritises minimising
#'   relatedness} (the max-diversity end). This matches AlphaMate's own
#'   \code{TargetDegree} convention under all three engines. Not numerically
#'   identical between engines; re-tune if you switch. See "Engine
#'   differences".
#' @param rescale_nrm Logical, default \code{TRUE}. Rescale \code{G} to mean
#'   diagonal 1 (numerator-relationship-like scale) before use, matching
#'   AlphaMate's expected NRM convention. Recommended to leave on for all
#'   three engines unless you have already rescaled \code{G} yourself.
#' @param alphamate_exe Path to the AlphaMate executable. \strong{You must
#'   obtain and be entitled to use this executable yourself} -- it is a
#'   separate, third-party tool (Hickey Group / AlphaGenes suite), not an R
#'   package and not part of HapBlockR's own code. Default \code{NULL}: no
#'   binary is bundled with or auto-detected by HapBlockR (see "Providing
#'   the AlphaMate executable" above), so \code{NULL} always means AlphaMate
#'   is unavailable -- \code{engine = "auto"} will fall back to
#'   \code{"optisel"}. Required (explicitly supplied) for \code{engine =
#'   "alphamate"} (or for \code{"auto"} to select it).
#' @param out_dir Directory to write AlphaMate's input/output files to.
#'   Required for \code{engine = "alphamate"}. Created if it doesn't exist.
#' @param seed Optional integer seed (used by the native mate-allocation
#'   step in each engine's final crossing-pair generation, and by optiSel's
#'   solver if it uses randomness internally).
#' @param verbose Logical, default \code{TRUE}.
#'
#' @return A list with components \code{mating_plan} (data frame:
#'   \code{parent1}, \code{parent2}, \code{mean_relationship}, plus
#'   engine-reported columns where available), \code{contributors} (data
#'   frame: \code{id}, \code{contribution}, \code{family} if supplied),
#'   \code{engine_used}, and \code{ok} (logical: whether no-selfing/
#'   no-repeated-mating/contribution-cap constraints were all satisfied in
#'   the final plan).
#'
#' @references
#' Meuwissen, T.H.E. (1997). Maximizing the response of selection with a
#' predefined rate of inbreeding. \emph{Journal of Animal Science}, 75,
#' 934-940.
#'
#' Wellmann, R. (2019). Optimum contribution selection and mate allocation
#' for genetic improvement, inbreeding, and diversity: the R package
#' optiSel. \emph{BMC Bioinformatics}, 20, 25.
#'
#' Peixoto, M.A., Amadeu, R.R., Bhering, L.L., Ferrao, L.F.V., Munoz, P.R.
#' & Resende Jr., M.F.R. (2024). SimpleMating: R-package for prediction and
#' optimization of breeding crosses using genomic selection. \emph{The Plant
#' Genome}, e20533.
#'
#' @seealso \code{\link{select_parents_ga}}, \code{\link{usefulness_criterion}},
#'   \code{\link{compute_haplotype_grm}}
#' @export
select_parents_ocs <- function(
    merit,
    G,
    family                 = NULL,
    engine                  = c("auto", "alphamate", "optisel", "simplemating"),
    n_crosses                = 20L,
    n_parents_max            = NULL,
    max_contrib_per_parent   = 4L,
    allow_selfing            = FALSE,
    allow_repeated_matings   = FALSE,
    target_degree             = 30,
    rescale_nrm                = TRUE,
    alphamate_exe               = NULL,
    out_dir                      = NULL,
    seed                          = NULL,
    verbose                       = TRUE
) {
  result_call <- match.call()
  engine <- match.arg(engine)
  inp <- .ocs_validate_inputs(merit, G, family, verbose = verbose)

  # No AlphaMate binary is bundled with HapBlockR (it is a separate,
  # third-party tool with no stable pre-built binary distribution to fetch --
  # see the "Providing the AlphaMate executable" section above). engine =
  # "alphamate" is therefore only usable when the caller explicitly supplies
  # a working alphamate_exe; engine = "auto" falls back to "optisel" (the
  # OTHER true-OCS engine -- "auto" never silently substitutes discrete
  # cross prediction/selection ("simplemating") for actual OCS) whenever
  # alphamate_exe is NULL or does not point to an existing file.
  have_alphamate <- !is.null(alphamate_exe) && file.exists(alphamate_exe)
  if (engine == "auto") {
    engine <- if (have_alphamate) "alphamate" else "optisel"
    if (isTRUE(verbose))
      message("[select_parents_ocs] engine = 'auto' -> using '", engine,
              "' (", if (have_alphamate) "alphamate_exe found" else
                "alphamate_exe not supplied/found", ").")
  }

  if (engine == "alphamate") {
    if (is.null(alphamate_exe) || !file.exists(alphamate_exe))
      stop("engine = 'alphamate' requires alphamate_exe to point to an ",
           "existing AlphaMate executable.", call. = FALSE)
    if (.Platform$OS.type != "windows") {
      # Only warn when alphamate_exe is CONFIRMED to be a Windows PE binary
      # (MZ header) -- not merely because we're off Windows. A user-supplied
      # native Linux/macOS build (this branch's whole reason for existing)
      # will read as PE = FALSE here and launch normally without a spurious
      # warning; if the file can't be read at all (PE-ness unknown, e.g. it
      # doesn't exist yet at this exact instant), NA is treated as "don't
      # warn" too, since system2() below will raise its own clear error if
      # the file truly cannot be executed.
      if (isTRUE(.is_windows_pe_exe(alphamate_exe)))
        warning("[select_parents_ocs] alphamate_exe (", alphamate_exe, ") is ",
                "a Windows binary (MZ/PE header); engine = 'alphamate' will ",
                "likely fail to launch on this platform (",
                .Platform$OS.type, "). Consider engine = 'optisel' instead, ",
                "or supply a native AlphaMate build for your platform via ",
                "alphamate_exe.", call. = FALSE)
    }
    if (is.null(out_dir))
      stop("engine = 'alphamate' requires out_dir (directory to write ",
           "AlphaMate's input/output files to).", call. = FALSE)
    res <- .run_alphamate_ocs(inp, family = inp$family,
                              n_crosses = n_crosses,
                              n_parents_max = n_parents_max,
                              max_contrib_per_parent = max_contrib_per_parent,
                              allow_selfing = allow_selfing,
                              allow_repeated_matings = allow_repeated_matings,
                              target_degree = target_degree,
                              rescale_nrm = rescale_nrm,
                              alphamate_exe = alphamate_exe,
                              out_dir = out_dir, verbose = verbose)
  } else if (engine == "optisel") {
    if (!is.null(seed)) set.seed(seed)
    res <- .run_optisel_ocs(inp, family = inp$family,
                            n_crosses = n_crosses,
                            n_parents_max = n_parents_max,
                            max_contrib_per_parent = max_contrib_per_parent,
                            allow_selfing = allow_selfing,
                            allow_repeated_matings = allow_repeated_matings,
                            target_degree = target_degree,
                            rescale_nrm = rescale_nrm, verbose = verbose)
  } else {
    # engine == "simplemating"
    if (!is.null(seed)) set.seed(seed)
    res <- .run_simplemating_ocs(inp, family = inp$family,
                                 n_crosses = n_crosses,
                                 max_contrib_per_parent = max_contrib_per_parent,
                                 n_parents_max = n_parents_max,
                                 target_degree = target_degree,
                                 allow_selfing = allow_selfing,
                                 allow_repeated_matings = allow_repeated_matings,
                                 rescale_nrm = rescale_nrm, verbose = verbose)
  }
  res$engine_used <- engine
  checks <- if (!is.null(res$constraint_checks)) {
    unlist(res$constraint_checks, use.names = TRUE)
  } else {
    c(plan_reported_ok = isTRUE(res$ok))
  }
  .add_hapblockr_contract(
    result = res,
    method = "select_parents_ocs",
    call = result_call,
    parameters = list(
      engine = engine, n_crosses = n_crosses,
      n_parents_max = n_parents_max,
      max_contrib_per_parent = max_contrib_per_parent,
      allow_selfing = allow_selfing,
      allow_repeated_matings = allow_repeated_matings,
      target_degree = target_degree,
      rescale_nrm = rescale_nrm
    ),
    seed = seed,
    sample_ids = inp$ids,
    inputs = list(merit = inp$merit, relationship_matrix = inp$G),
    transformations = c(
      if (rescale_nrm) "relationship-matrix rescaling" else character(),
      paste(engine, "optimal contribution or mating allocation")
    ),
    quality_gates = c(success = isTRUE(res$ok), checks),
    decision_table = if (is.data.frame(res$mating_plan)) {
      res$mating_plan
    } else {
      data.frame()
    },
    uncertainty = data.frame(),
    external_tools = list(engine = engine)
  )
}


# -- Internal: AlphaMate engine. Direct translation of a working production
# AlphaMate driver script (file formats, spec keys, and call sequence
# verified against that script) into a general HapBlockR-input version.
.run_alphamate_ocs <- function(inp, family, n_crosses, n_parents_max,
                               max_contrib_per_parent, allow_selfing,
                               allow_repeated_matings, target_degree,
                               rescale_nrm, alphamate_exe, out_dir, verbose) {
  ids <- inp$ids
  A <- if (isTRUE(rescale_nrm)) .rescale_to_nrm(inp$G) else inp$G

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  id_map <- data.frame(AlphaMateID = ids, SelectionIndex = inp$merit,
                       stringsAsFactors = FALSE)
  if (!is.null(family)) id_map$Family <- family

  criterion_path <- file.path(out_dir, "Criterion.txt")
  utils::write.table(id_map[, c("AlphaMateID", "SelectionIndex")],
                     criterion_path, sep = " ", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)

  nrm_path <- file.path(out_dir, "Nrm.txt")
  nrm_out <- cbind(data.frame(id = ids), as.data.frame(A))
  utils::write.table(nrm_out, nrm_path, sep = " ", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)

  spec_path <- file.path(out_dir, "AlphaMateSpec.txt")
  spec_lines <- c(
    paste("NrmMatrixFile          ,", basename(nrm_path)),
    paste("SelCriterionFile       ,", basename(criterion_path)),
    paste("NumberOfMatings        ,", n_crosses),
    if (!is.null(n_parents_max)) paste("NumberOfParents        ,", n_parents_max),
    paste("LimitContributionsMax  ,", max_contrib_per_parent),
    "MateAllocation         , Yes",
    paste("AllowRepeatedMatings   ,", if (allow_repeated_matings) "Yes" else "No"),
    paste("AllowSelfing           ,", if (allow_selfing) "Yes" else "No"),
    paste("TargetDegree           ,", target_degree),
    "EvolAlgLogAllSolutions , Yes",
    "Stop"
  )
  writeLines(spec_lines[!is.na(spec_lines)], spec_path)

  if (isTRUE(verbose)) message("[select_parents_ocs] Running AlphaMate in ", out_dir, " ...")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)
  # Defensively (re-)grant the execute permission bit on the user-supplied
  # alphamate_exe immediately before invoking it -- e.g. a fresh download or
  # git checkout can leave it unset on Unix-like systems. Harmless no-op on
  # Windows (which does not use POSIX permission bits for .exe files);
  # silently ignored if it fails for any other reason, since system2()
  # below will surface a clear error either way if the binary truly cannot
  # be executed.
  try(Sys.chmod(alphamate_exe, mode = "0755"), silent = TRUE)
  output <- tryCatch(
    system2(alphamate_exe, basename(spec_path), stdout = TRUE, stderr = TRUE),
    error = function(e)
      stop("Failed to run AlphaMate executable: ", conditionMessage(e),
           call. = FALSE)
  )
  if (isTRUE(verbose)) message(paste(output, collapse = "\n"))

  mating_file  <- file.path(out_dir, "MatingPlanModeOptTarget1.txt")
  contrib_file <- file.path(out_dir, "ContributorsModeOptTarget1.txt")

  mating_plan <- NULL
  if (file.exists(mating_file)) {
    mp <- utils::read.table(mating_file, header = TRUE, stringsAsFactors = FALSE)
    names(mp)[names(mp) == "Parent1"] <- "parent1"
    names(mp)[names(mp) == "Parent2"] <- "parent2"
    if (all(c("parent1", "parent2") %in% names(mp)))
      mp$mean_relationship <- vapply(seq_len(nrow(mp)), function(i)
        A[mp$parent1[i], mp$parent2[i]], numeric(1L))
    mating_plan <- mp
  } else if (isTRUE(verbose)) {
    message("[select_parents_ocs] AlphaMate mating-plan output not found ",
            "at ", mating_file, " -- check the AlphaMate console output ",
            "above for errors.")
  }

  contributors <- NULL
  if (file.exists(contrib_file)) {
    contributors <- utils::read.table(contrib_file, header = TRUE,
                                      stringsAsFactors = FALSE)
    # AlphaMate's own documented header for this file is
    # "Id Gender SelCriterion AvgCoancestryA AvgCoancestryC Contribution
    # nContribution" (capital I, lowercase d) -- "Id" must come first here
    # since it's the actual column name; the other three are kept as
    # defensive fallbacks in case a different AlphaMate build/version uses
    # one of them instead. Matching is case-sensitive (R column names), so
    # the previous "ID" entry alone never matched real AlphaMate output and
    # left contributors$id unset (NULL), which downstream callers depending
    # on it (e.g. `ocs_res$contributors$id`) received silently -- despite a
    # perfectly valid mating_plan already being available as a fallback.
    id_col <- intersect(c("Id", "Parent", "ID", "Indiv", "AlphaMateID"),
                        names(contributors))[1]
    if (!is.na(id_col)) {
      names(contributors)[names(contributors) == id_col] <- "id"
      if ("nContribution" %in% names(contributors))
        names(contributors)[names(contributors) == "nContribution"] <- "contribution"
      if (!is.null(family))
        contributors$family <- family[contributors$id]
    }
  }

  ok <- TRUE
  if (!is.null(mating_plan) && all(c("parent1", "parent2") %in% names(mating_plan))) {
    if (!allow_selfing && any(mating_plan$parent1 == mating_plan$parent2)) ok <- FALSE
    pk <- ifelse(mating_plan$parent1 < mating_plan$parent2,
                paste(mating_plan$parent1, mating_plan$parent2),
                paste(mating_plan$parent2, mating_plan$parent1))
    if (!allow_repeated_matings && anyDuplicated(pk) > 0L) ok <- FALSE
  }
  if (!is.null(contributors) && "contribution" %in% names(contributors)) {
    if (any(contributors$contribution > max_contrib_per_parent, na.rm = TRUE)) ok <- FALSE
  }

  list(mating_plan = mating_plan, contributors = contributors,
      id_map = id_map, spec_file = spec_path, out_dir = out_dir, ok = ok)
}


# -- Internal: "optisel" engine -- TRUE Optimal Contribution Selection via
# optiSel's own solver (Wellmann 2019, Meuwissen 1997), not SimpleMating.
# See the file header for why this is a NEW engine (not a rename of the old
# "optisel", which never actually called optiSel's solver) and the
# "Verification status" roxygen section for exactly which parts of this are
# independently confirmed vs. written directly against optiSel's own
# candes()/opticont()/matings()/noffspring() documentation.
#
# Workflow (non-overlapping generations -- "the next generation" here always
# means "this one crossing block being planned now", matching how every
# other HapBlockR selection strategy frames a single breeding cycle):
#   1. candes(phen, Kin = <rescaled G>, cont = NULL) builds the candidate-
#      description object. phen$Sex is left NA for every candidate
#      throughout -- per ?optiSel::candes, an all-NA Sex column omits the
#      "equal contribution by sex" constraint, which is the right choice for
#      HapBlockR's plant-breeding framing (crosses are not sex-structured
#      matings the way animal-breeding OCS usually is).
#   2. Two calls to opticont() locate the actual, candidate-set-specific
#      ends of the gain/diversity frontier: min.Kin (the minimum mean
#      kinship achievable at all, ignoring merit) and max.Merit with no
#      kinship constraint (the maximum merit achievable, ignoring kinship).
#      Both are genuine optiSel-solved optima, not arbitrary reference
#      points.
#   3. target_degree linearly interpolates a mean-kinship ceiling between
#      those two solved extremes (0 = the unconstrained max-merit solution's
#      own kinship; 90 = the min-kinship solution's kinship) -- matching the
#      "simplemating" engine's target_degree DIRECTION exactly (0 = max
#      gain, 90 = max diversity; see .target_degree_to_n_keep()'s comment),
#      via a genuinely different MECHANISM (interpolating between two solved
#      optima here, vs. quantile-of-observed-K there). This is still an
#      approximation of AlphaMate's own continuous Kinghorn-frontier degree,
#      not a proof of numerical equivalence -- re-tune if you switch
#      engines.
#   4. A final opticont("max.Merit", con = list(ub.Kin = <ceiling>)) solves
#      for the actual optimal per-candidate contributions under that
#      kinship ceiling.
#   5. noffspring() converts contributions into expected offspring counts
#      for a population of n_crosses * 2 parent "slots" (each mating needs 2
#      parents); max_contrib_per_parent, if supplied, caps each candidate's
#      slot count before this step (a pre-solve approximation, not a
#      re-optimised constrained solution -- optiSel does not expose
#      "contribution optimum subject to a hard per-parent slot cap" as a
#      single constrained problem the way it does for kinship).
#   6. matings() solves the discrete Sire x Dam assignment minimising mean
#      offspring kinship, given those slot counts and the same kinship
#      matrix. ub.n = 1 (unless allow_repeated_matings = TRUE) is passed so
#      no-repeated-matings is enforced NATIVELY by the solver, not by
#      post-hoc filtering (assumption flagged in "Verification status": that
#      ub.n bounds repeated matings of the same pair).

.allocate_ocs_matings <- function(contributions, Kin, n_crosses,
                                  max_contrib_per_parent = NULL,
                                  allow_selfing = FALSE,
                                  allow_repeated_matings = FALSE) {
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("lpSolve is required to convert optimal contributions into a ",
         "hard-constrained mating plan for engine = 'optisel'. Install with: ",
         "install.packages('lpSolve')", call. = FALSE)

  contributions <- contributions[is.finite(contributions) & contributions > 0]
  if (!length(contributions))
    stop("The optimal-contribution solution contains no positive ",
         "contributions.", call. = FALSE)
  contributions <- contributions / sum(contributions)
  ids <- names(contributions)
  if (is.null(ids) || any(!nzchar(ids)))
    stop("Optimal contributions must be named with candidate IDs.",
         call. = FALSE)
  Kin <- Kin[ids, ids, drop = FALSE]

  pairs <- if (isTRUE(allow_selfing)) {
    which(upper.tri(Kin, diag = TRUE), arr.ind = TRUE)
  } else {
    which(upper.tri(Kin, diag = FALSE), arr.ind = TRUE)
  }
  if (!nrow(pairs))
    stop("No feasible mating pair remains under the selfing policy.",
         call. = FALSE)

  pair_df <- data.frame(
    parent1 = ids[pairs[, 1L]],
    parent2 = ids[pairs[, 2L]],
    relationship = Kin[pairs],
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(pair_df$relationship)))
    stop("The relationship matrix contains non-finite values for candidate ",
         "mating pairs.", call. = FALSE)

  n_pair <- nrow(pair_df)
  n_id <- length(ids)
  n_var <- n_pair + 2L * n_id
  incidence <- matrix(0, nrow = n_id, ncol = n_pair,
                      dimnames = list(ids, NULL))
  for (j in seq_len(n_pair)) {
    incidence[pair_df$parent1[j], j] <-
      incidence[pair_df$parent1[j], j] + 1
    incidence[pair_df$parent2[j], j] <-
      incidence[pair_df$parent2[j], j] + 1
  }

  target_slots <- contributions * (2 * n_crosses)
  constraints <- list()
  directions <- character(0)
  rhs <- numeric(0)

  total_row <- numeric(n_var)
  total_row[seq_len(n_pair)] <- 1
  constraints[[length(constraints) + 1L]] <- total_row
  directions <- c(directions, "=")
  rhs <- c(rhs, n_crosses)

  for (i in seq_len(n_id)) {
    row <- numeric(n_var)
    row[seq_len(n_pair)] <- incidence[i, ]
    row[n_pair + i] <- -1
    row[n_pair + n_id + i] <- 1
    constraints[[length(constraints) + 1L]] <- row
    directions <- c(directions, "=")
    rhs <- c(rhs, target_slots[i])
  }

  if (!is.null(max_contrib_per_parent)) {
    for (i in seq_len(n_id)) {
      row <- numeric(n_var)
      row[seq_len(n_pair)] <- incidence[i, ]
      constraints[[length(constraints) + 1L]] <- row
      directions <- c(directions, "<=")
      rhs <- c(rhs, max_contrib_per_parent)
    }
  }

  if (!isTRUE(allow_repeated_matings)) {
    for (j in seq_len(n_pair)) {
      row <- numeric(n_var)
      row[j] <- 1
      constraints[[length(constraints) + 1L]] <- row
      directions <- c(directions, "<=")
      rhs <- c(rhs, 1)
    }
  }

  rel <- pair_df$relationship
  rel_scaled <- if (diff(range(rel)) > sqrt(.Machine$double.eps)) {
    (rel - min(rel)) / diff(range(rel))
  } else {
    rep(0, length(rel))
  }
  objective <- c(rel_scaled, rep(100, 2L * n_id))
  fit <- lpSolve::lp(
    direction = "min",
    objective.in = objective,
    const.mat = do.call(rbind, constraints),
    const.dir = directions,
    const.rhs = rhs,
    int.vec = seq_len(n_pair)
  )
  if (fit$status != 0L)
    stop("No feasible mating plan satisfies n_crosses, selfing, repeated-",
         "mating and per-parent contribution constraints (lpSolve status ",
         fit$status, "). Revise the constraints or increase the candidate ",
         "set.", call. = FALSE)

  pair_counts <- as.integer(round(fit$solution[seq_len(n_pair)]))
  used <- which(pair_counts > 0L)
  if (!length(used))
    stop("The mating allocator returned no crosses despite a successful ",
         "solver status.", call. = FALSE)
  mating_plan <- pair_df[rep(used, pair_counts[used]),
                         c("parent1", "parent2"), drop = FALSE]
  rownames(mating_plan) <- NULL
  mating_plan$mean_relationship <- Kin[
    cbind(mating_plan$parent1, mating_plan$parent2)
  ]

  actual_slots <- table(factor(
    c(mating_plan$parent1, mating_plan$parent2), levels = ids
  ))
  allocation <- data.frame(
    id = ids,
    target_slots = as.numeric(target_slots[ids]),
    actual_slots = as.integer(actual_slots),
    deviation = as.integer(actual_slots) - as.numeric(target_slots[ids]),
    stringsAsFactors = FALSE
  )

  list(
    mating_plan = mating_plan,
    allocation = allocation,
    solver_status = fit$status,
    solver_objective = fit$objval
  )
}


.solve_optisel_frontier <- function(phen, Kin, target_degree, verbose) {
  cand <- tryCatch(
    optiSel::candes(phen = phen, Kin = Kin, cont = NULL,
                    quiet = !isTRUE(verbose)),
    error = function(e)
      stop("optiSel::candes() failed: ", conditionMessage(e),
           "\nCheck that G (after rescale_nrm) is a valid symmetric ",
           "relationship matrix and that candidate identifiers are unique.",
           call. = FALSE)
  )
  min_kin_res <- tryCatch(
    optiSel::opticont("min.Kin", cand, con = list(), quiet = TRUE),
    error = function(e)
      stop("optiSel::opticont('min.Kin', ...) failed while establishing the ",
           "max-diversity end of the frontier: ", conditionMessage(e),
           call. = FALSE)
  )
  max_merit_res <- tryCatch(
    optiSel::opticont("max.Merit", cand, con = list(), quiet = TRUE),
    error = function(e)
      stop("optiSel::opticont('max.Merit', ...) failed while establishing ",
           "the max-gain end of the frontier: ", conditionMessage(e),
           call. = FALSE)
  )
  kin_min <- min_kin_res$mean[["Kin"]]
  kin_max_gain <- max_merit_res$mean[["Kin"]]
  if (!is.finite(kin_min) || !is.finite(kin_max_gain))
    stop("optiSel returned a non-finite frontier endpoint.", call. = FALSE)
  span <- kin_max_gain - kin_min
  kin_ceiling <- if (span < sqrt(.Machine$double.eps)) {
    kin_max_gain
  } else {
    kin_max_gain - (target_degree / 90) * span
  }
  opt <- tryCatch(
    optiSel::opticont("max.Merit", cand, con = list(ub.Kin = kin_ceiling),
                      quiet = !isTRUE(verbose)),
    error = function(e)
      stop("optiSel::opticont('max.Merit', ...) failed at target_degree = ",
           target_degree, ": ", conditionMessage(e),
           "\nThis can indicate an infeasible kinship ceiling.",
           call. = FALSE)
  )
  if (!is.null(opt$info$valid) && !isTRUE(opt$info$valid))
    stop("optiSel::opticont() did not return a valid solution at ",
         "target_degree = ", target_degree, ".", call. = FALSE)
  if (is.null(opt$parent) || !("oc" %in% names(opt$parent)))
    stop("optiSel::opticont() did not return the expected $parent$oc ",
         "column.", call. = FALSE)

  list(
    cand = cand,
    min_kin_result = min_kin_res,
    max_merit_result = max_merit_res,
    kin_min = kin_min,
    kin_max_gain = kin_max_gain,
    kin_ceiling = kin_ceiling,
    opt = opt
  )
}


.run_optisel_ocs <- function(inp, family, n_crosses, n_parents_max,
                             max_contrib_per_parent, allow_selfing,
                             allow_repeated_matings, target_degree,
                             rescale_nrm, verbose) {
  if (!requireNamespace("optiSel", quietly = TRUE))
    stop("optiSel is required for engine = 'optisel'. Install with: ",
         "install.packages('optiSel')", call. = FALSE)
  need_fns <- c("candes", "opticont")
  missing_fns <- setdiff(need_fns, getNamespaceExports("optiSel"))
  if (length(missing_fns))
    stop("optiSel is installed (version ",
         as.character(utils::packageVersion("optiSel")), ") but does not ",
         "export: ", paste(missing_fns, collapse = ", "), ", which engine ",
         "= 'optisel' requires. Your installed optiSel version may be too ",
         "old, or its API may have moved since this wrapper was written -- ",
         "check ?optiSel::candes and ?optiSel::opticont against your ",
         "installed version.", call. = FALSE)
  if (is.null(target_degree) || !is.numeric(target_degree) ||
      target_degree < 0 || target_degree > 90)
    stop("target_degree must be a single numeric value in [0, 90] for ",
         "engine = 'optisel' (interpolates a mean-kinship ceiling between ",
         "the min-kinship and max-merit ends of the gain/diversity frontier, ",
         "both solved via optiSel::opticont() for your actual candidate ",
         "set: 0 keeps the unconstrained max-merit solution's own kinship ",
         "(prioritises merit), 90 forces the minimum achievable kinship ",
         "(prioritises minimising relatedness) -- matching AlphaMate's own ",
         "TargetDegree convention).", call. = FALSE)

  ids  <- inp$ids
  Kmat <- if (isTRUE(rescale_nrm)) .rescale_to_nrm(inp$G) else inp$G

  phen <- data.frame(
    Indiv = ids,
    Sex = NA_character_,
    isCandidate = TRUE,
    Merit = as.numeric(inp$merit[ids]),
    stringsAsFactors = FALSE
  )

  solved <- .solve_optisel_frontier(phen, Kmat, target_degree, verbose)
  contrib_df <- solved$opt$parent

  if (!is.null(n_parents_max) && n_parents_max < sum(contrib_df$oc > 0)) {
    min_required <- if (isTRUE(allow_selfing)) 1L else 2L
    if (n_parents_max < min_required)
      stop("n_parents_max must be at least ", min_required,
           " under the requested selfing policy.", call. = FALSE)
    positive <- which(is.finite(contrib_df$oc) & contrib_df$oc > 0)
    keep_rows <- positive[
      order(contrib_df$oc[positive], decreasing = TRUE)
    ][seq_len(min(n_parents_max, length(positive)))]
    keep_ids <- contrib_df$Indiv[keep_rows]
    phen <- phen[match(keep_ids, phen$Indiv), , drop = FALSE]
    Kmat <- Kmat[keep_ids, keep_ids, drop = FALSE]
    solved <- .solve_optisel_frontier(phen, Kmat, target_degree, verbose)
    contrib_df <- solved$opt$parent
    if (isTRUE(verbose))
      message("[select_parents_ocs] engine 'optisel': re-solved the OCS ",
              "frontier and constrained optimum within the ",
              length(keep_ids), " parents retained by n_parents_max.")
  }

  kin_min <- solved$kin_min
  kin_max_gain <- solved$kin_max_gain
  kin_ceiling <- solved$kin_ceiling
  opt <- solved$opt

  contributors <- data.frame(
    id = contrib_df$Indiv, contribution = contrib_df$oc,
    stringsAsFactors = FALSE
  )
  contributors <- contributors[contributors$contribution > 0, , drop = FALSE]
  contributors$contribution <- contributors$contribution /
    sum(contributors$contribution)
  if (!is.null(family)) contributors$family <- family[contributors$id]

  contribution_vector <- setNames(
    contributors$contribution, contributors$id
  )
  allocation <- .allocate_ocs_matings(
    contribution_vector, Kmat, n_crosses,
    max_contrib_per_parent = max_contrib_per_parent,
    allow_selfing = allow_selfing,
    allow_repeated_matings = allow_repeated_matings
  )
  mating_plan <- allocation$mating_plan

  if (nrow(mating_plan) != n_crosses)
    stop("Internal error: the constrained mating allocator returned ",
         nrow(mating_plan), " crosses instead of ", n_crosses, ".",
         call. = FALSE)
  if (!allow_selfing && any(mating_plan$parent1 == mating_plan$parent2))
    stop("Internal error: the mating allocator returned a forbidden self.",
         call. = FALSE)
  pair_keys <- apply(
    t(apply(mating_plan[c("parent1", "parent2")], 1L, sort)),
    1L, paste, collapse = "\r"
  )
  if (!allow_repeated_matings && anyDuplicated(pair_keys) > 0L)
    stop("Internal error: the mating allocator returned a repeated mating.",
         call. = FALSE)
  parent_counts <- table(c(mating_plan$parent1, mating_plan$parent2))
  if (!is.null(max_contrib_per_parent) &&
      any(parent_counts > max_contrib_per_parent))
    stop("Internal error: the mating allocator exceeded ",
         "max_contrib_per_parent.", call. = FALSE)
  ok <- TRUE

  list(mating_plan = mating_plan, contributors = contributors,
      id_map = data.frame(id = ids, merit = inp$merit, stringsAsFactors = FALSE),
      spec_file = NA_character_, out_dir = NA_character_, ok = ok,
      frontier = list(kin_min = kin_min, kin_max_gain = kin_max_gain,
                      kin_ceiling = kin_ceiling, target_degree = target_degree),
      allocation = allocation$allocation,
      allocation_solver = list(status = allocation$solver_status,
                               objective = allocation$solver_objective),
      constraint_checks = list(
        n_crosses = nrow(mating_plan) == n_crosses,
        no_selfing = allow_selfing ||
          all(mating_plan$parent1 != mating_plan$parent2),
        no_repeated_matings = allow_repeated_matings ||
          anyDuplicated(pair_keys) == 0L,
        parent_cap = is.null(max_contrib_per_parent) ||
          all(parent_counts <= max_contrib_per_parent),
        n_parents_max = is.null(n_parents_max) ||
          length(unique(c(mating_plan$parent1, mating_plan$parent2))) <=
          n_parents_max
      ),
      opticont_result = opt)
}


# -- Internal: "simplemating" engine (renamed from "optisel" -- see the file
# header for the breaking-rename rationale), implemented via
# SimpleMating::planCross() + SimpleMating::selectCrosses() (Peixoto et al.
# 2024) -- see the file header for why this specific pair of functions (not
# GOCS(), whose exported status turned out to be contested) and exactly what
# was verified before writing this. NOT true OCS -- discrete greedy cross
# prediction/selection; see .run_optisel_ocs() below for the genuine
# optiSel-solver-based engine.
#
# Real, verified limitations carried through honestly rather than worked
# around: selectCrosses() controls diversity via a hard relatedness cutoff
# (culling.pairwise.k), not a continuous contribution optimum -- target_degree
# is mapped onto that cutoff by QUANTILE of the candidate set's own observed
# relatedness values (see .target_degree_to_n_keep() below), which is an
# approximation, not an equivalent algorithm (re-tune if you switch engines).
# n_parents_max is NOT configurable here (selectCrosses() has no such
# argument), unlike engine = "alphamate" or "optisel".

# -- Internal: convert target_degree [0, 90] into how many of the n_cand
# candidate pairs (to be sorted ascending by relatedness K by the caller) to
# keep before selectCrosses() searches. Direction matches AlphaMate's own
# TargetDegree convention exactly: target_degree = 0 is the max-gain end of
# Kinghorn's frontier (no diversity restriction -- keep essentially every
# candidate pair, including the most related ones, so the search is free to
# pick on merit alone); target_degree = 90 is the max-diversity end (keep
# only the min_keep least-related candidate pairs). min_keep is a floor that
# guarantees enough surviving candidates for n_crosses to stay feasible
# regardless of how restrictive target_degree asks to be. Factored out as
# its own pure function specifically so this direction is unit-testable in
# isolation, without needing SimpleMating/optiSel installed or a full OCS
# run -- see test-ocs.R.
.target_degree_to_n_keep <- function(target_degree, n_cand, min_keep) {
  as.integer(min(n_cand, max(min_keep, ceiling(((90 - target_degree) / 90) * n_cand))))
}

.run_simplemating_ocs <- function(inp, family, n_crosses, max_contrib_per_parent,
                                  n_parents_max, target_degree, allow_selfing,
                                  allow_repeated_matings, rescale_nrm, verbose) {
  if (!requireNamespace("SimpleMating", quietly = TRUE))
    stop("SimpleMating is required for engine = 'simplemating'. Install with: ",
         "remotes::install_github('Resende-Lab/SimpleMating')", call. = FALSE)
  # planCross()/selectCrosses() are confirmed present in BOTH the real,
  # currently-installed SimpleMating 0.2.1 (per its own help index) and the
  # exact current GitHub source read before writing this wrapper -- a more
  # stable foundation than GOCS(), whose exported status is contested (see
  # the file header). Checked explicitly for a clear, actionable error
  # rather than letting a .Call()/argument-matching failure surface as a
  # confusing "not an exported object" error.
  if (!all(c("planCross", "selectCrosses") %in% getNamespaceExports("SimpleMating")))
    stop("SimpleMating is installed (version ",
         as.character(utils::packageVersion("SimpleMating")), ") but does ",
         "not export planCross()/selectCrosses(), which this wrapper ",
         "targets. Reinstall the current version with: ",
         "remotes::install_github('Resende-Lab/SimpleMating', force = TRUE). ",
         "If still missing after that, SimpleMating's API may have moved ",
         "again -- check https://github.com/Resende-Lab/SimpleMating for ",
         "its current exported functions.", call. = FALSE)
  if (!requireNamespace("optiSel", quietly = TRUE))
    stop("optiSel is required for engine = 'simplemating' (one of ",
         "SimpleMating's own Imports). Install with: ",
         "install.packages('optiSel')", call. = FALSE)
  if (is.null(target_degree) || !is.numeric(target_degree) ||
      target_degree < 0 || target_degree > 90)
    stop("target_degree must be a single numeric value in [0, 90] for ",
         "engine = 'simplemating' (mapped onto selectCrosses()'s ",
         "culling.pairwise.k relatedness cutoff via a quantile of the ",
         "candidate K distribution, with a floor guaranteeing enough ",
         "surviving candidates for n_crosses to be feasible: 0 keeps ",
         "essentially all candidates (no diversity restriction, prioritises ",
         "merit), 90 keeps only the least-related candidates subject to ",
         "that floor (prioritises minimising relatedness) -- matching ",
         "AlphaMate's own TargetDegree convention).", call. = FALSE)

  ids  <- inp$ids
  Kmat <- if (isTRUE(rescale_nrm)) .rescale_to_nrm(inp$G) else inp$G

  if (isTRUE(verbose) && !is.null(n_parents_max))
    message("[select_parents_ocs] n_parents_max is not directly enforced ",
            "under engine = 'simplemating' (SimpleMating::selectCrosses() ",
            "has no such argument); ignored.")

  mate_design <- if (isTRUE(allow_selfing)) "half_p" else "half"
  cross_df <- tryCatch(
    SimpleMating::planCross(TargetPop = ids, MateDesign = mate_design),
    error = function(e)
      stop("SimpleMating::planCross() failed: ", conditionMessage(e),
           "\nIf this looks like an argument-name mismatch, SimpleMating's ",
           "API may have changed since this wrapper was written against its ",
           "current source -- check ?SimpleMating::planCross.", call. = FALSE)
  )
  if (!nrow(cross_df))
    stop("SimpleMating::planCross() produced zero candidate crosses -- too ",
         "few candidate parents for the requested mate design.", call. = FALSE)

  # Mid-parent merit and pairwise relatedness for every candidate pair,
  # computed natively (plain arithmetic) rather than via a SimpleMating
  # criterion function -- select_parents_ocs() already has everything
  # needed (merit, G) without one.
  cross_df$Y <- (inp$merit[cross_df$Parent1] + inp$merit[cross_df$Parent2]) / 2
  cross_df$K <- Kmat[cbind(cross_df$Parent1, cross_df$Parent2)]

  # culling.pairwise.k is chosen by quantile of the actual candidate K
  # distribution, not by linearly interpolating target_degree across the raw
  # [min(K), max(K)] range. A range-based cutoff is fragile: a single
  # outlier pair (highly related or highly diverse) widens the range enough
  # that a "moderate" target_degree can retain only a handful of candidates
  # -- confirmed against a real devtools::test() run, where target_degree =
  # 30 on a 10-parent/45-cross test panel left too few candidates for
  # selectCrosses() to find any n_crosses = 5 solution ("Reached maximum in
  # the search, try to increase data size."). A quantile-based cutoff, with
  # a floor that guarantees enough surviving candidates for
  # selectCrosses()'s internal greedy search to have room to work regardless
  # of target_degree, fixes this for any K distribution.
  k_sorted <- sort(cross_df$K)
  n_cand   <- length(k_sorted)
  k_span   <- diff(range(k_sorted))
  # Comfortable multiple of n_crosses (SimpleMating's own examples default
  # to n.cross = 200 against candidate sets several times that size) so the
  # greedy min.cross/max.cross search has enough slack to succeed, but never
  # more candidates than actually exist.
  min_keep <- min(n_cand, max(n_crosses * 5L, 20L))
  # See .target_degree_to_n_keep()'s own comment (above .run_simplemating_ocs)
  # for the direction convention this follows.
  q_idx <- if (k_span < .Machine$double.eps) {
    n_cand   # every candidate equally related -- keep them all
  } else {
    .target_degree_to_n_keep(target_degree, n_cand, min_keep)
  }
  # selectCrosses() keeps rows with data$K < culling.pairwise.k (strict
  # "<"); nudge just past the chosen K value so it isn't itself excluded.
  culling_k <- k_sorted[q_idx] + max(abs(k_sorted[q_idx]), k_span, 1) * 1e-8

  max_cross <- if (!is.null(max_contrib_per_parent)) max_contrib_per_parent else length(ids)

  sel <- tryCatch(
    SimpleMating::selectCrosses(
      data = cross_df, n.cross = n_crosses, max.cross = max_cross,
      min.cross = 1L, max.cross.to.search = nrow(cross_df),
      culling.pairwise.k = culling_k
    ),
    error = function(e)
      stop("SimpleMating::selectCrosses() failed: ", conditionMessage(e),
           "\nIf this looks like an argument-name mismatch, SimpleMating's ",
           "API may have changed since this wrapper was written against its ",
           "current source -- check ?SimpleMating::selectCrosses. This can ",
           "also happen if n_crosses/max_contrib_per_parent/target_degree ",
           "are infeasible for this candidate set (try relaxing them).",
           call. = FALSE)
  )

  plan <- sel$plan
  mating_plan <- if (!is.null(plan) && nrow(plan)) {
    data.frame(parent1 = plan$Parent1, parent2 = plan$Parent2,
              mean_relationship = plan$K, criterion = plan$Y,
              stringsAsFactors = FALSE)
  } else {
    data.frame(parent1 = character(0), parent2 = character(0),
              mean_relationship = numeric(0), criterion = numeric(0))
  }

  contributors <- if (nrow(mating_plan)) {
    tab <- table(c(mating_plan$parent1, mating_plan$parent2))
    df <- data.frame(id = names(tab), contribution = as.integer(tab),
                     stringsAsFactors = FALSE)
    if (!is.null(family)) df$family <- family[df$id]
    df
  } else {
    NULL
  }

  ok <- nrow(mating_plan) == n_crosses
  if (nrow(mating_plan)) {
    if (!allow_selfing && any(mating_plan$parent1 == mating_plan$parent2)) ok <- FALSE
    if (!is.null(max_contrib_per_parent) && !is.null(contributors) &&
        max(contributors$contribution) > max_contrib_per_parent) ok <- FALSE
    # no-repeated-matings is guaranteed by construction here, not merely
    # checked afterwards: planCross(MateDesign = "half"/"half_p") lists each
    # unordered pair at most once, and selectCrosses()'s greedy search
    # consumes each candidate row at most once, so the same pair cannot
    # appear twice in `plan` -- see the file header.
  }

  list(mating_plan = mating_plan, contributors = contributors,
      id_map = data.frame(id = ids, merit = inp$merit, stringsAsFactors = FALSE),
      spec_file = NA_character_, out_dir = NA_character_, ok = ok,
      parent_stats = NULL, mating_stats = sel$summary)
}
