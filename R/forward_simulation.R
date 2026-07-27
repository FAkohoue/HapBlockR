# ==============================================================================
# forward_simulation.R
# Forward-in-time recurrent-selection simulation: GA-selected founders vs.
# truncation-selection founders, tracking realised genetic gain over
# generations.
#
# Closes the gap identified against HapSelect's localGEBV_vs_TS_simulation()/
# Haplotype_vs_TS_simulation() (see HapBlockR_vs_HapSelect_comparison.md,
# section 3.2) -- "here is the actual answer, and here is the proof it's
# better than the obvious baseline."
#
# THIS IS A WRAPPER AROUND THE genomicSimulation PACKAGE (the same one
# HapSelect uses: https://github.com/vllrs/genomicSimulation), not a
# self-built meiosis engine. An earlier version of this file implemented its
# own simplified crossover simulator instead, before genomicSimulation's
# real R source (bundled locally under simulation/*.R -- sim-setup.R,
# sim-progression.R, sim-calculators.R, sim-group-utils.R, sim-data-access.R,
# sim-deletors.R, utils.R) was read in full; this version calls the real,
# documented API instead of reimplementing meiosis.
#
# genomicSimulation architecture note: its R API is a thin wrapper around a
# SINGLE package-level mutable pointer (`genomicSimulation:::sim.data$p`) --
# there is no "create an independent session" call. Only one SimData can be
# active at a time. Consequently:
#   - This file's two schemes (GA-founders, TS-founders) are run as two
#     SEPARATE, sequential genomicSimulation sessions
#     (clear.simdata() -> load.data() -> ... -> clear.simdata()), each
#     loading ONLY that scheme's founders. This sidesteps any risk of index
#     collisions if the same individual happens to be selected as a founder
#     by both schemes, and keeps each scheme's simulation fully independent.
#   - ga_vs_ts_simulation() therefore takes exclusive ownership of
#     genomicSimulation's global state for the duration of the call and
#     always leaves it cleared afterwards (via on.exit). Do not call this
#     function from multiple threads/processes sharing one R session, and
#     do not rely on genomicSimulation state surviving after this function
#     returns -- everything needed is copied into this function's return
#     value before the state is cleared.
#
# Biological model:
#   - Requires PHASED haplotypes (hap1/hap2, e.g. from read_phased_vcf()).
#     Unphased 0/1/2 dosage cannot support a haplotype-block-preserving
#     meiosis simulation: without knowing which alleles are in cis at
#     heterozygous sites, simulated recombination would not respect the LD
#     blocks the whole package is built around.
#   - Meiosis/crossing itself is entirely delegated to genomicSimulation's
#     make.random.crosses(), which generates gametes according to the
#     recombination map supplied (via a genetic map file built from
#     snp_info$POS and `recomb_rate`).
#   - Recurrent selection: if `selection_intensity` is supplied, the top
#     fraction of each generation's offspring (by GEBV, via
#     genomicSimulation::break.group.by.GEBV()) becomes the breeding
#     population for the next generation. If NULL, every offspring
#     generation is mated in full (no within-scheme truncation).
#   - GEBVs are calculated by genomicSimulation itself
#     (genomicSimulation::see.GEBVs()), using a marker-effect file built from
#     `snp_effects`, with per-marker centring values set via
#     genomicSimulation::change.eff.set.centres() so that genomicSimulation's
#     raw-dosage GEBV calculation reduces to HapBlockR's own centred-dosage
#     convention: GEBV = sum_t (x_t - 2*p_t) * alpha_t, identical to the
#     formula used throughout compute_local_gebv()/backsolve_snp_effects().
#     This makes genomicSimulation's simulated GEBVs directly comparable in
#     scale to the rest of the package's output.
#
# Two exported functions (original GA-vs-TS two-scheme calling convention
# unchanged, so existing calling code does not need to change):
#   ga_vs_ts_simulation()      -- runs the comparison, returns per-generation
#                                  summary statistics for every scheme.
#   plot_ga_vs_ts_simulation() -- plots realised genetic gain over
#                                  generations for every scheme.
#
# -- OCS/UC-informed rapid-cycling extension -----------------------------------
# ga_selected/ts_selected (both "truncation" mating: random-mate the current
# generation, then keep the top selection_intensity fraction by GEBV) were
# the only two schemes this file could simulate. This closes the "forward
# simulation only compares truncation selection, not OCS-/UC-informed rapid
# cycling" gap: an arbitrary named `schemes` list generalises the comparison
# to N schemes, each with its own `founders` and `mating_scheme` in
# \{"truncation", "ocs", "uc"\}. `ga_selected`/`ts_selected` remain as
# backward-compatible convenience arguments that build the original 2-scheme
# "truncation" list internally when `schemes` is not supplied.
#
#   "truncation" -- unchanged: genomicSimulation::make.random.crosses().
#   "ocs"        -- each generation, the current group's dosage genotypes +
#                    GEBVs are read back from live genomicSimulation state,
#                    a whole-genome additive GRM is built via
#                    compute_haplotype_grm() (generic over any 0/1/2 dosage
#                    matrix, not just haplotype-block features -- see its own
#                    docs), select_parents_ocs(engine = "simplemating",
#                    forced -- a generational loop must not depend on an
#                    OS-specific external executable, and this pins the
#                    exact algorithm/performance profile this scheme was
#                    written and tuned against; engine = "optisel" now also
#                    solves true OCS without an external binary and is a
#                    valid alternative here, but switching would change both
#                    the algorithm and the per-generation solve cost -- not
#                    done automatically by this rename) picks an
#                    optimum-contribution mating_plan, and those EXACT
#                    parent1/parent2 pairs are crossed.
#   "uc"         -- each generation, per-block local GEBVs are computed via
#                    compute_local_gebv() (requires `blocks`, the LD-block
#                    table, in addition to `snp_info`/`snp_effects`),
#                    usefulness_criterion(variance_model = "block_independent")
#                    ranks every candidate pair, and the top-ranked pairs are
#                    crossed.
#
# Both "ocs" and "uc" execute their specific chosen pairs via
# genomicSimulation::make.targeted.crosses(first.parents=, second.parents=),
# NOT the make.random.crosses() used by "truncation".
#
# @section API source: make.random.crosses(), break.group.by.GEBV(),
# see.GEBVs(), see.group.gene.data(), and see.group.data() were already in
# use by this file before this extension. This extension additionally calls
# genomicSimulation::make.targeted.crosses(), whose signature
# (first.parents=, second.parents=, offspring=, give.names=, name.prefix=,
# track.pedigree=, give.ids=, retain=, map=) was read directly from the
# bundled real source (simulation/sim-progression.R,
# make.targeted.crosses(), formerly cross.combinations()) -- not guessed --
# confirming in particular that first.parents/second.parents accept
# individual NAMES (not just indexes), which is what both new mating
# schemes rely on.
# ==============================================================================


# -- Internal: write a genomicSimulation-format genotype matrix file for a
# subset of individuals, preserving hap1/hap2 phase exactly. Uses "0"/"1" as
# the two allele symbols -- genomicSimulation treats allele symbols as
# arbitrary single characters (see simulation/sim-setup.R load.genotypes()
# docs: "Any pair of characters ... each character is an allele"), so there
# is no need to map back to real REF/ALT nucleotides for simulation purposes.
.gs_write_genotype_file <- function(hap1, hap2, ids, founder_ids) {
  idx <- match(founder_ids, ids)
  h1  <- hap1[, idx, drop = FALSE]
  h2  <- hap2[, idx, drop = FALSE]
  pairs <- matrix(paste0(h1, h2), nrow = nrow(h1), ncol = ncol(h1))
  out <- cbind(
    data.frame(name = rownames(hap1), stringsAsFactors = FALSE),
    as.data.frame(pairs, stringsAsFactors = FALSE)
  )
  colnames(out) <- c("name", founder_ids)
  tf <- tempfile("ldxgs_geno_", fileext = ".tsv")
  utils::write.table(out, tf, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE)
  tf
}

# -- Internal: write a genomicSimulation-format genetic map file. Converts
# physical position (bp) to a genetic map distance (cM) using recomb_rate
# (Morgans/bp): pos_cM = POS_bp * recomb_rate * 100. Absolute scale does not
# matter to genomicSimulation (see.genetic.map() re-anchors each chromosome
# to its first marker internally); only relative distances drive the
# Poisson/Haldane crossover process.
.gs_write_map_file <- function(snp_info, recomb_rate) {
  chr    <- .norm_chr_hap(as.character(snp_info$CHR))
  pos_cm <- as.numeric(snp_info$POS) * recomb_rate * 100
  df <- data.frame(marker = as.character(snp_info$SNP), chr = chr, pos = pos_cm,
                   stringsAsFactors = FALSE)
  tf <- tempfile("ldxgs_map_", fileext = ".tsv")
  utils::write.table(df, tf, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE)
  tf
}

# -- Internal: write a genomicSimulation-format marker effect file. Allele
# "1" (our alt-dosage symbol, see .gs_write_genotype_file) gets the
# backsolved SNP effect; allele "0" gets an explicit effect of 0. Per-marker
# centring (to reproduce HapBlockR's 2p-centred convention) is set
# separately afterwards via change.eff.set.centres(), not via this file.
.gs_write_effect_file <- function(common_snp, alpha) {
  df <- data.frame(
    marker = rep(common_snp, 2L),
    allele = rep(c("1", "0"), each = length(common_snp)),
    eff    = c(alpha[common_snp], rep(0, length(common_snp))),
    stringsAsFactors = FALSE
  )
  tf <- tempfile("ldxgs_eff_", fileext = ".tsv")
  utils::write.table(df, tf, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE)
  tf
}

# -- Internal: normalize the `schemes` argument (or the legacy ga_selected/
# ts_selected convenience args) into a single named list, one entry per
# scheme, each with (at least) $founders (character vector, >= 2 members,
# subset of `ids`) and $mating_scheme (one of "truncation"/"ocs"/"uc",
# default "truncation" if not supplied). Pure R, no genomicSimulation
# dependency -- kept as its own function specifically so it can be unit-
# tested without genomicSimulation installed (see test-forward-simulation.R).
.gs_normalize_schemes <- function(schemes, ga_selected, ts_selected, ids) {
  if (is.null(schemes)) {
    if (is.null(ga_selected) || is.null(ts_selected))
      stop("Either `schemes` or both `ga_selected` and `ts_selected` must ",
           "be supplied.", call. = FALSE)
    schemes <- list(
      GA = list(founders = ga_selected, mating_scheme = "truncation"),
      TS = list(founders = ts_selected, mating_scheme = "truncation")
    )
  } else {
    if (!is.null(ga_selected) || !is.null(ts_selected))
      message("[ga_vs_ts_simulation] `schemes` was supplied -- ",
              "`ga_selected`/`ts_selected` are ignored.")
    if (!is.list(schemes) || is.null(names(schemes)) || any(!nzchar(names(schemes))))
      stop("`schemes` must be a named list (one name per scheme, used as ",
           "its label in the returned `summary`).", call. = FALSE)
    if (anyDuplicated(names(schemes)))
      stop("`schemes` names must be unique.", call. = FALSE)
  }

  valid_ms <- c("truncation", "ocs", "uc")
  for (nm in names(schemes)) {
    sc <- schemes[[nm]]
    if (!is.list(sc) || is.null(sc$founders))
      stop("schemes[[\"", nm, "\"]] must be a list with (at least) a ",
           "`founders` element.", call. = FALSE)
    if (!all(sc$founders %in% ids))
      stop("Some founders in schemes[[\"", nm, "\"]] are not columns of ",
           "hap1/hap2.", call. = FALSE)
    if (length(sc$founders) < 2L)
      stop("schemes[[\"", nm, "\"]]$founders must have >= 2 individuals.",
           call. = FALSE)
    if (is.null(sc$mating_scheme)) sc$mating_scheme <- "truncation"
    if (!sc$mating_scheme %in% valid_ms)
      stop("schemes[[\"", nm, "\"]]$mating_scheme must be one of ",
           paste(shQuote(valid_ms), collapse = ", "), ".", call. = FALSE)
    schemes[[nm]] <- sc
  }
  schemes
}

# -- Internal: one generation's OCS-informed crossing for mating_scheme =
# "ocs". Reads the current group's dosage genotypes + GEBVs back from live
# genomicSimulation state, builds a whole-genome additive GRM via
# compute_haplotype_grm() (generic over any individuals x features 0/1/2
# dosage matrix -- see its own docs -- not limited to haplotype-block
# features), calls select_parents_ocs() (engine forced to "simplemating": a
# generational loop must not depend on an OS-specific external executable
# like AlphaMate.exe, and pinning this specific engine -- rather than
# "optisel", which as of the engine rename/re-implementation in R/ocs.R is a
# genuinely different true-OCS algorithm with a different per-generation
# solve cost -- keeps this scheme's algorithm and performance profile
# exactly what it was written and tested against) to get an
# optimum-contribution mating_plan, then
# executes those EXACT parent1/parent2 pairs via
# genomicSimulation::make.targeted.crosses() -- see the file header's
# "Verification status" section for how that function's signature was
# confirmed. `scheme_args` (the scheme list entry minus founders/
# mating_scheme) is passed through to select_parents_ocs(), letting callers
# override n_crosses/max_contrib_per_parent/target_degree/etc per scheme.
.gs_cross_ocs <- function(cur_grp, eff_id, scheme_label, g, pop_size, scheme_args) {
  gm   <- genomicSimulation::see.group.gene.data(cur_grp, count.allele = "1")
  ids  <- colnames(gm)
  dose <- t(gm)
  storage.mode(dose) <- "numeric"
  rownames(dose) <- ids
  gebv <- stats::setNames(
    genomicSimulation::see.GEBVs(cur_grp, effect.set = eff_id), ids)

  G <- compute_haplotype_grm(dose, phased = TRUE)

  ocs_call_args <- utils::modifyList(
    list(merit = gebv, G = G, engine = "simplemating",
        n_crosses = pop_size, verbose = FALSE),
    scheme_args
  )
  ocs_res <- do.call(select_parents_ocs, ocs_call_args)

  mp <- ocs_res$mating_plan
  if (is.null(mp) || !nrow(mp))
    stop("select_parents_ocs() returned an empty mating_plan for scheme '",
         scheme_label, "' at generation ", g, ".", call. = FALSE)

  genomicSimulation::make.targeted.crosses(
    first.parents = mp$parent1, second.parents = mp$parent2,
    offspring = 1L, give.names = TRUE,
    name.prefix = paste0(scheme_label, "_g", g, "_"),
    track.pedigree = TRUE, give.ids = TRUE, retain = TRUE)
}

# -- Internal: one generation's UC-informed crossing for mating_scheme =
# "uc". Reads the current group's dosage genotypes + GEBVs back from live
# genomicSimulation state, computes per-block local GEBV via
# compute_local_gebv() (requires `blocks`, `snp_info`, and the full per-SNP
# additive-effect vector `alpha_full` already built once by
# ga_vs_ts_simulation()), ranks ALL candidate cross pairs by
# usefulness_criterion(variance_model = "block_independent") -- the only
# mode that needs nothing beyond local GEBV, since individuals read back
# from genomicSimulation here are unphased dosage counts, not
# extract_haplotypes()-shaped phased blocks -- and executes the top-ranked
# pairs via make.targeted.crosses(). `scheme_args$n_crosses` controls how
# many top-ranked pairs are actually crossed (default: pop_size, one
# offspring per pair, matching "truncation"'s total offspring count);
# remaining `scheme_args` entries are passed through to
# usefulness_criterion() (e.g. to override selected_proportion).
.gs_cross_uc <- function(cur_grp, eff_id, marker_order, snp_info, blocks,
                         alpha_full, scheme_label, g, pop_size, scheme_args) {
  gm   <- genomicSimulation::see.group.gene.data(cur_grp, count.allele = "1")
  ids  <- colnames(gm)
  dose <- t(gm)
  storage.mode(dose) <- "numeric"
  dimnames(dose) <- list(ids, marker_order)
  gebv <- stats::setNames(
    genomicSimulation::see.GEBVs(cur_grp, effect.set = eff_id), ids)

  loc <- compute_local_gebv(dose, snp_info, blocks, alpha_full,
                            scale = TRUE, complete_decomposition = TRUE)

  n_crosses <- if (!is.null(scheme_args$n_crosses)) scheme_args$n_crosses else pop_size
  uc_call_args <- utils::modifyList(
    list(parent_ids = ids, gebv = gebv, variance_model = "block_independent",
        block_importance = loc$block_importance, local_gebv = loc$local_gebv,
        selected_proportion = 0.1, verbose = FALSE),
    scheme_args[setdiff(names(scheme_args), "n_crosses")]
  )
  uc_res <- do.call(usefulness_criterion, uc_call_args)
  uc_res <- uc_res[!is.na(uc_res$UC), , drop = FALSE]
  if (!nrow(uc_res))
    stop("usefulness_criterion() returned no scorable cross for scheme '",
         scheme_label, "' at generation ", g, ".", call. = FALSE)
  top <- utils::head(uc_res, n_crosses)

  genomicSimulation::make.targeted.crosses(
    first.parents = top$parent1, second.parents = top$parent2,
    offspring = 1L, give.names = TRUE,
    name.prefix = paste0(scheme_label, "_g", g, "_"),
    track.pedigree = TRUE, give.ids = TRUE, retain = TRUE)
}

# -- Internal: run one scheme's full recurrent-selection simulation as its
# own isolated genomicSimulation session (clear -> load -> cross/select loop
# -> extract results -> clear). Returns a per-generation summary data frame
# and the final generation's phased haplotypes (reconstructed from
# genomicSimulation's own internal genotype storage, so marker row order is
# read back from genomicSimulation itself via see.genetic.map() rather than
# assumed to match the input order).
.gs_run_scheme <- function(scheme_label, founder_ids, hap1, hap2, ids,
                           map_file, eff_file, centre_df,
                           n_generations, pop_size, selection_intensity,
                           mating_scheme = "truncation",
                           snp_info = NULL, blocks = NULL, alpha_full = NULL,
                           scheme_args = list(),
                           verbose) {
  genomicSimulation::clear.simdata()
  on.exit(genomicSimulation::clear.simdata(), add = TRUE)

  geno_file <- .gs_write_genotype_file(hap1, hap2, ids, founder_ids)
  on.exit(unlink(geno_file), add = TRUE)

  fmt <- genomicSimulation::define.matrix.format.details(
    has.header = TRUE, markers.as.rows = TRUE, cell.style = "P")

  loaded <- genomicSimulation::load.data(allele.file = geno_file,
                                         map.file = map_file,
                                         effect.file = eff_file,
                                         format = fmt)
  eff_id <- loaded$effectID
  grp    <- loaded$groupNum

  genomicSimulation::change.eff.set.centres(to = centre_df, effect.set = eff_id)

  marker_order <- genomicSimulation::see.genetic.map()$marker

  gebv0 <- genomicSimulation::see.GEBVs(grp, effect.set = eff_id)
  records <- list(data.frame(
    generation = 0L, mean_gebv = mean(gebv0), max_gebv = max(gebv0),
    sd_gebv = stats::sd(gebv0), n_pop = length(gebv0)
  ))

  cur_grp <- grp
  for (g in seq_len(n_generations)) {
    if (isTRUE(verbose))
      message("[ga_vs_ts_simulation] ", scheme_label, ": generation ", g,
              " (mating_scheme = '", mating_scheme, "') ...")

    child_grp <- if (mating_scheme == "truncation") {
      genomicSimulation::make.random.crosses(
        cur_grp, n.crosses = pop_size, cap = 0L, offspring = 1L,
        give.names = TRUE, name.prefix = paste0(scheme_label, "_g", g, "_"),
        track.pedigree = TRUE, give.ids = TRUE, retain = TRUE)
    } else if (mating_scheme == "ocs") {
      .gs_cross_ocs(cur_grp, eff_id, scheme_label, g, pop_size, scheme_args)
    } else if (mating_scheme == "uc") {
      .gs_cross_uc(cur_grp, eff_id, marker_order, snp_info, blocks,
                  alpha_full, scheme_label, g, pop_size, scheme_args)
    } else {
      stop("Unknown mating_scheme '", mating_scheme, "'.", call. = FALSE)
    }

    gebv_g <- genomicSimulation::see.GEBVs(child_grp, effect.set = eff_id)
    records[[length(records) + 1L]] <- data.frame(
      generation = g, mean_gebv = mean(gebv_g), max_gebv = max(gebv_g),
      sd_gebv = stats::sd(gebv_g), n_pop = length(gebv_g)
    )

    if (!is.null(selection_intensity)) {
      sel_grp <- genomicSimulation::break.group.by.GEBV(
        child_grp, low.score.best = FALSE,
        percentage = selection_intensity * 100, effect.set = eff_id)
      # child_grp now refers only to the non-selected remainder.
      genomicSimulation::delete.group(child_grp)
      genomicSimulation::delete.group(cur_grp)
      cur_grp <- sel_grp
    } else {
      genomicSimulation::delete.group(cur_grp)
      cur_grp <- child_grp
    }
  }

  final_ids   <- genomicSimulation::see.group.data(cur_grp, "N")
  final_gebv  <- genomicSimulation::see.GEBVs(cur_grp, effect.set = eff_id)
  final_pairs <- genomicSimulation::see.group.gene.data(cur_grp)

  final_hap1 <- matrix(as.numeric(substr(final_pairs, 1L, 1L)),
                       nrow = nrow(final_pairs), ncol = ncol(final_pairs),
                       dimnames = list(marker_order, final_ids))
  final_hap2 <- matrix(as.numeric(substr(final_pairs, 2L, 2L)),
                       nrow = nrow(final_pairs), ncol = ncol(final_pairs),
                       dimnames = list(marker_order, final_ids))

  list(
    summary    = do.call(rbind, records),
    final_hap1 = final_hap1,
    final_hap2 = final_hap2,
    final_gebv = stats::setNames(final_gebv, final_ids)
  )
}


#' Simulate GA-Selected vs. Truncation-Selected Founders Over Generations
#'
#' @description
#' Runs recurrent-selection forward simulation for two founder sets --
#' typically the output of \code{\link{select_parents_ga}} and
#' \code{\link{truncation_selection}} -- and tracks realised genetic gain
#' (mean/max GEBV of the breeding population) over \code{n_generations}.
#' This is how HapSelect's own \code{localGEBV_vs_TS_simulation()}/
#' \code{Haplotype_vs_TS_simulation()} demonstrate that a GA-optimised
#' founder set outperforms plain truncation selection, and this function is
#' a wrapper around the same underlying engine HapSelect uses for that
#' comparison: the \pkg{genomicSimulation} package
#' (\url{https://github.com/vllrs/genomicSimulation}). See the file header
#' of \code{R/forward_simulation.R} for how the wrapper is structured
#' (each scheme runs as its own isolated genomicSimulation session, since
#' genomicSimulation itself has only one active simulation at a time).
#'
#' @section Installation of genomicSimulation:
#' \pkg{genomicSimulation} is not on CRAN. Install it from a source release:
#' \enumerate{
#'   \item Download the \code{.tar.gz} from
#'     \url{https://github.com/vllrs/genomicSimulation/releases}.
#'   \item \code{install.packages("path/to/genomicSimulation_x.y.z.tar.gz", repos = NULL)}
#' }
#' On a shared HPC system, compile on the same CPU architecture/node type you
#' will run jobs on (see your cluster's documentation for details -- e.g. UQ
#' Bunya requires compiling on an "epyc3" node to avoid "illegal instruction"
#' errors at runtime on other nodes).
#'
#' \strong{Requires phased haplotypes.} Unphased 0/1/2 dosage cannot support
#' a block-preserving meiosis simulation -- without cis/trans information at
#' heterozygous sites, simulated recombination would not respect the LD
#' blocks the founder sets were selected to stack. Use
#' \code{\link{read_phased_vcf}} or another phased source for \code{hap1}/
#' \code{hap2}.
#'
#' @param hap1,hap2 Numeric matrices (SNPs x individuals), values 0/1 --
#'   e.g. \code{read_phased_vcf()$hap1} / \code{$hap2}. Column names are
#'   individual IDs; row names/order must match \code{snp_info}.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS},
#'   same row order as \code{hap1}/\code{hap2}.
#' @param snp_effects Named numeric vector of per-SNP additive effects
#'   (from \code{\link{backsolve_snp_effects}} or
#'   \code{\link{estimate_marker_effects}}), names matching \code{snp_info$SNP}.
#'   GEBVs are calculated by genomicSimulation on the same centred-dosage
#'   scale used elsewhere in HapBlockR (\eqn{(x - 2p) \cdot \alpha}, via
#'   \code{genomicSimulation::change.eff.set.centres()}), so trajectories are
#'   directly comparable to \code{run_haplotype_prediction()} output.
#' @param ga_selected Character vector of founder IDs (from
#'   \code{select_parents_ga()$selected}), must be a subset of
#'   \code{colnames(hap1)}. Ignored (with a message) if \code{schemes} is
#'   supplied. Together with \code{ts_selected}, this is the original
#'   2-scheme calling convention -- kept for backward compatibility; it
#'   builds \code{schemes = list(GA = list(founders = ga_selected,
#'   mating_scheme = "truncation"), TS = list(founders = ts_selected,
#'   mating_scheme = "truncation"))} internally.
#' @param ts_selected Character vector of founder IDs (from
#'   \code{truncation_selection()$selected}), same requirement as
#'   \code{ga_selected}.
#' @param schemes Optional named list, one entry per scheme to simulate
#'   (overrides \code{ga_selected}/\code{ts_selected} when supplied; each
#'   list name is used as the scheme's label in the returned
#'   \code{summary}). Each entry is itself a list with:
#'   \describe{
#'     \item{\code{founders}}{Required. Character vector of founder IDs
#'       (subset of \code{colnames(hap1)}, >= 2 individuals).}
#'     \item{\code{mating_scheme}}{One of \code{"truncation"} (default if
#'       omitted -- random-mate the current generation, keep the top
#'       \code{selection_intensity} fraction by GEBV; the original,
#'       only behaviour of this function), \code{"ocs"} (each generation,
#'       \code{\link{select_parents_ocs}} picks an optimum-contribution
#'       mating plan from the current population, executed via
#'       \code{genomicSimulation::make.targeted.crosses()}; requires
#'       \code{blocks}), or \code{"uc"} (each generation,
#'       \code{\link{usefulness_criterion}} ranks every candidate pair by
#'       predicted progeny usefulness and the top-ranked pairs are crossed;
#'       requires \code{blocks}).}
#'     \item{Any other named element}{Passed through to the underlying
#'       per-generation call: \code{select_parents_ocs()} for
#'       \code{"ocs"} (e.g. \code{n_crosses}, \code{max_contrib_per_parent},
#'       \code{target_degree}), or \code{usefulness_criterion()} for
#'       \code{"uc"} (e.g. \code{selected_proportion}), plus the
#'       \code{"uc"}-only \code{n_crosses} (how many top-ranked pairs are
#'       actually crossed each generation; default \code{pop_size}, one
#'       offspring per pair).}
#'   }
#'   See the "OCS/UC-informed rapid-cycling extension" section of the
#'   \code{R/forward_simulation.R} file header for the full design and its
#'   "Verification status" note.
#' @param blocks LD-block table from \code{\link{run_Big_LD_all_chr}} (or
#'   \code{\link{tune_LD_params}}). Required if any scheme uses
#'   \code{mating_scheme = "uc"} (each generation's per-block local GEBVs are
#'   computed via \code{\link{compute_local_gebv}}, which needs the block
#'   table). Not needed by \code{"ocs"} (which builds a whole-genome GRM
#'   directly from dosage, no block decomposition involved) or
#'   \code{"truncation"}. Validated up front, before any simulation work is
#'   done, whenever at least one scheme requests \code{"uc"}.
#' @param n_generations Integer. Number of generations to simulate. Default
#'   \code{10L}.
#' @param pop_size Integer. Offspring produced per generation, per scheme.
#'   Default \code{100L}.
#' @param recomb_rate Numeric. Used to convert \code{snp_info$POS} (bp) into
#'   the genetic map (cM) that genomicSimulation's own crossing functions
#'   consume: \code{pos_cM = POS_bp * recomb_rate * 100}. Default \code{1e-8}
#'   (~1 cM/Mb, a common genome-wide average approximation). Rescale for a
#'   species with a known genetic map:
#'   \code{recomb_rate = map_length_cM / 100 / genome_length_bp}.
#' @param selection_intensity Numeric in (0,1] or \code{NULL}. Fraction of
#'   each generation's offspring retained as parents for the next generation
#'   (recurrent truncation selection within each scheme, via
#'   \code{genomicSimulation::break.group.by.GEBV()}). Default \code{0.2}
#'   (top 20\%). Set \code{NULL} to mate every offspring generation in full
#'   (no within-scheme truncation).
#' @param seed Integer or \code{NULL}. Random seed for reproducibility
#'   (\code{set.seed()}; genomicSimulation's own crossing functions use R's
#'   random number generator).
#' @param verbose Logical. Print per-generation progress. Default \code{FALSE}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{summary}}{Data frame: \code{scheme} (scheme label -- "GA"/
#'     "TS" under the legacy 2-scheme convention, or the \code{schemes}
#'     list's own names), \code{generation} (0 = founders themselves),
#'     \code{mean_gebv}, \code{max_gebv}, \code{sd_gebv}, \code{n_pop}.}
#'   \item{\code{final}}{Named list (one element per scheme, same names as
#'     \code{summary$scheme}) of lists with \code{hap1}/\code{hap2}
#'     (reconstructed from genomicSimulation's internal storage; markers x
#'     individuals, row names = marker names in genomicSimulation's own
#'     marker order, which may differ from \code{snp_info}'s input row
#'     order) and \code{gebv} (named numeric vector) for the final
#'     generation's breeding population of each scheme, for further analysis
#'     (e.g. diversity metrics).}
#'   \item{\code{ga_final}, \code{ts_final}}{Only present when called via the
#'     legacy \code{ga_selected}/\code{ts_selected} convention (i.e.
#'     \code{schemes} not supplied) -- identical content to
#'     \code{final$GA}/\code{final$TS}, kept as top-level names unchanged
#'     from previous releases so existing calling code does not break.}
#' }
#'
#' @seealso \code{\link{select_parents_ga}}, \code{\link{truncation_selection}},
#'   \code{\link{select_parents_ocs}}, \code{\link{usefulness_criterion}},
#'   \code{\link{plot_ga_vs_ts_simulation}}
#'
#' @examples
#' \dontrun{
#' phased <- read_phased_vcf("mydata_phased.vcf.gz")
#' res    <- run_haplotype_prediction(phased$dosage, phased$snp_info, blocks,
#'                                    blues = blues_vec)
#' top    <- select_top_blocks(res$block_importance, n = 15)
#' vmat   <- res$local_gebv[, top$block_id, drop = FALSE]
#' ga     <- select_parents_ga(vmat, n_founders = 20, seed = 1)
#' ts     <- truncation_selection(res$gebv, n_founders = 20)
#'
#' # Legacy 2-scheme convention (unchanged from previous releases):
#' sim <- ga_vs_ts_simulation(
#'   phased$hap1, phased$hap2, phased$snp_info, res$snp_effects,
#'   ga_selected = ga$selected, ts_selected = ts$selected,
#'   n_generations = 10, seed = 1
#' )
#' plot_ga_vs_ts_simulation(sim)
#'
#' # New: 4-way comparison including OCS- and UC-informed rapid cycling.
#' G <- compute_haplotype_grm(phased$dosage)
#' sim4 <- ga_vs_ts_simulation(
#'   phased$hap1, phased$hap2, phased$snp_info, res$snp_effects,
#'   schemes = list(
#'     GA  = list(founders = ga$selected, mating_scheme = "truncation"),
#'     TS  = list(founders = ts$selected, mating_scheme = "truncation"),
#'     OCS = list(founders = ts$selected, mating_scheme = "ocs",
#'                n_crosses = 20),
#'     UC  = list(founders = ts$selected, mating_scheme = "uc",
#'                selected_proportion = 0.1)
#'   ),
#'   blocks = blocks, n_generations = 10, seed = 1
#' )
#' plot_ga_vs_ts_simulation(sim4)
#' }
#'
#' @export
ga_vs_ts_simulation <- function(hap1, hap2, snp_info, snp_effects,
                                ga_selected = NULL, ts_selected = NULL,
                                schemes = NULL,
                                blocks = NULL,
                                n_generations = 10L,
                                pop_size = 100L,
                                recomb_rate = 1e-8,
                                selection_intensity = 0.2,
                                seed = NULL,
                                verbose = FALSE) {
  if (!requireNamespace("genomicSimulation", quietly = TRUE))
    stop("Package 'genomicSimulation' is required for ga_vs_ts_simulation() ",
         "but is not installed. It is not on CRAN -- see the 'Installation ",
         "of genomicSimulation' section of ?ga_vs_ts_simulation, or ",
         "https://github.com/vllrs/genomicSimulation for source releases ",
         "and installation instructions.", call. = FALSE)

  if (!is.matrix(hap1)) hap1 <- as.matrix(hap1)
  if (!is.matrix(hap2)) hap2 <- as.matrix(hap2)
  if (!identical(dim(hap1), dim(hap2)))
    stop("hap1 and hap2 must have identical dimensions.", call. = FALSE)
  if (nrow(hap1) != nrow(snp_info))
    stop("nrow(hap1)/nrow(hap2) must match nrow(snp_info).", call. = FALSE)
  if (!all(c("SNP", "CHR", "POS") %in% names(snp_info)))
    stop("snp_info must have columns SNP, CHR, POS.", call. = FALSE)

  ids <- colnames(hap1)
  if (is.null(ids))
    stop("hap1/hap2 must have column names (individual IDs).", call. = FALSE)

  schemes_norm <- .gs_normalize_schemes(schemes, ga_selected, ts_selected, ids)
  legacy_2scheme <- is.null(schemes)

  needs_blocks <- vapply(schemes_norm, function(s) identical(s$mating_scheme, "uc"),
                         logical(1))
  if (any(needs_blocks) && is.null(blocks))
    stop("blocks is required when any scheme uses mating_scheme = 'uc' ",
         "(scheme(s): ", paste(names(schemes_norm)[needs_blocks], collapse = ", "),
         "). See ?compute_local_gebv.", call. = FALSE)

  if (!is.null(selection_intensity) &&
      (selection_intensity <= 0 || selection_intensity > 1))
    stop("selection_intensity must be in (0, 1] or NULL.", call. = FALSE)

  snp_names <- as.character(snp_info$SNP)
  rownames(hap1) <- rownames(hap2) <- snp_names

  # Allele frequencies (2p-centring, matching backsolve_snp_effects()/
  # compute_local_gebv()'s diploid convention -- genomicSimulation itself
  # only models diploid crosses, consistent with the rest of the package's
  # phased hap1/hap2 representation being diploid-only).
  dose_all <- t(hap1) + t(hap2)
  p_freq   <- pmax(pmin(colMeans(dose_all, na.rm = TRUE) / 2, 1 - 1e-8), 1e-8)
  names(p_freq) <- snp_names

  alpha_full <- stats::setNames(rep(0, length(snp_names)), snp_names)
  common_snp <- intersect(snp_names, names(snp_effects))
  if (!length(common_snp))
    stop("No SNP names in common between snp_info$SNP and names(snp_effects).",
         call. = FALSE)
  alpha_full[common_snp] <- snp_effects[common_snp]

  # -- Build the shared (map, effect-file, centring) inputs ONCE -------------
  map_file <- .gs_write_map_file(snp_info, recomb_rate)
  eff_file <- .gs_write_effect_file(common_snp, alpha_full)
  on.exit(unlink(c(map_file, eff_file)), add = TRUE)

  centre_df <- data.frame(
    marker = common_snp,
    centre = 2 * p_freq[common_snp] * alpha_full[common_snp],
    stringsAsFactors = FALSE
  )

  if (!is.null(seed)) set.seed(seed)

  scheme_runs <- vector("list", length(schemes_norm))
  names(scheme_runs) <- names(schemes_norm)
  for (nm in names(schemes_norm)) {
    sc <- schemes_norm[[nm]]
    if (isTRUE(verbose))
      message("[ga_vs_ts_simulation] Running '", nm, "' scheme (mating_scheme = '",
              sc$mating_scheme, "') ...")
    extra_args <- sc[setdiff(names(sc), c("founders", "mating_scheme"))]
    scheme_runs[[nm]] <- .gs_run_scheme(
      nm, sc$founders, hap1, hap2, ids,
      map_file, eff_file, centre_df,
      n_generations, pop_size, selection_intensity,
      mating_scheme = sc$mating_scheme,
      snp_info = snp_info, blocks = blocks, alpha_full = alpha_full,
      scheme_args = extra_args, verbose = verbose
    )
  }

  summary_list <- lapply(names(scheme_runs), function(nm) {
    df <- scheme_runs[[nm]]$summary
    df$scheme <- nm
    df
  })
  summary_df <- do.call(rbind, summary_list)
  summary_df <- summary_df[, c("scheme", "generation", "mean_gebv",
                               "max_gebv", "sd_gebv", "n_pop")]

  final_list <- lapply(scheme_runs, function(r)
    list(hap1 = r$final_hap1, hap2 = r$final_hap2, gebv = r$final_gebv))
  names(final_list) <- names(scheme_runs)

  out <- list(summary = summary_df, final = final_list)
  # Backward compatibility: when called the legacy way (schemes = NULL,
  # ga_selected/ts_selected supplied), also expose $ga_final/$ts_final under
  # their original top-level names, unchanged from previous releases.
  if (legacy_2scheme) {
    out$ga_final <- final_list$GA
    out$ts_final <- final_list$TS
  }
  out
}


#' Plot Realised Genetic Gain: GA vs. Truncation Selection (and Beyond)
#'
#' @description
#' Plots mean (solid) and max (dashed) breeding-population GEBV over
#' generations for every scheme from \code{\link{ga_vs_ts_simulation}},
#' mirroring HapSelect's genetic-gain-over-generations plot. Works for the
#' original 2-scheme (GA/TS) comparison as well as an arbitrary-length
#' \code{schemes} comparison (e.g. adding \code{"ocs"}/\code{"uc"}-informed
#' rapid-cycling schemes).
#'
#' @param sim List returned by \code{\link{ga_vs_ts_simulation}}.
#' @param show_max Logical. Overlay the max-GEBV trajectory (dashed) in
#'   addition to the mean (solid). Default \code{TRUE}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{ga_vs_ts_simulation}}
#'
#' @examples
#' \dontrun{
#' plot_ga_vs_ts_simulation(sim)
#' }
#'
#' @export
plot_ga_vs_ts_simulation <- function(sim, show_max = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required: install.packages('ggplot2')", call. = FALSE)
  if (!is.list(sim) || !"summary" %in% names(sim))
    stop("sim must be the list returned by ga_vs_ts_simulation().", call. = FALSE)

  df <- sim$summary
  p <- ggplot2::ggplot(df, ggplot2::aes(x = generation, y = mean_gebv,
                                        colour = scheme)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.8)

  if (isTRUE(show_max)) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(x = generation, y = max_gebv,
                                      colour = scheme),
                        linetype = "dashed", alpha = 0.7) +
      ggplot2::geom_point(ggplot2::aes(x = generation, y = max_gebv,
                                       colour = scheme),
                         shape = 1, alpha = 0.7)
  }

  # Colour-blind-safe (Okabe & Ito 2008) palette, keyed to the ORIGINAL
  # 2-scheme names for exact backward-compatible colours, extended to any
  # further scheme names (e.g. "OCS"/"UC") without erroring the way a fixed
  # 2-entry scale_colour_manual() would on a 3rd+ factor level.
  scheme_names   <- unique(as.character(df$scheme))
  okabe_ito      <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00",
                      "#56B4E9", "#F0E442", "#999999")
  preset         <- c(GA = "#D55E00", TS = "#0072B2")
  unpresetted    <- setdiff(scheme_names, names(preset))
  remaining_pal  <- setdiff(okabe_ito, unname(preset))
  pal <- c(preset[intersect(names(preset), scheme_names)],
          stats::setNames(rep(remaining_pal, length.out = length(unpresetted)),
                          unpresetted))

  p +
    ggplot2::scale_colour_manual(values = pal, name = "Scheme") +
    ggplot2::labs(x = "Generation", y = "GEBV (solid = mean, dashed = max)",
                 title = "Realised genetic gain: GA-selected vs. truncation-selected founders") +
    ggplot2::theme_minimal()
}
