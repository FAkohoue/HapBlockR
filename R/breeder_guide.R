# ==============================================================================
# breeder_guide.R
#
# HapBlockR_Breeder_Guide.pdf is a standalone, non-technical companion to the
# *From Local GEBV to a Crossing Decision* vignette: the same seven parent-
# and cross-selection tools (truncation selection, GA-based parent selection,
# usefulness_criterion(), select_parents_ocs(), select_parents_pareto(),
# validate_crosses_exact(), select_core_collection()), explained in plain
# language -- what each is, when to reach for it, when to be cautious, and
# how it works in practice -- for readers who will make or approve a crossing
# decision but do not themselves run R.
#
# Shipped as PDF (not Word/.docx) deliberately: PDF is not natively editable
# by the reader, which matters for a document meant to be handed to
# breeders/program managers/reviewers as a fixed reference rather than a
# draft they might inadvertently alter. It is built from a source .docx via
# LibreOffice headless conversion (the source is not shipped) and is not a
# vignette (no vignette engine indexes .pdf via vignette()/browseVignettes())
# and not built from source at install time -- it is a static, pre-built
# binary file. The correct place to ship an arbitrary auxiliary file like
# this with an R package is inst/extdata/ (the same convention this package
# already uses for beagle.jar and its example_*.csv/.vcf/.gds/.hmp.txt
# files), which is copied verbatim into the installed package and reachable
# via system.file(). open_breeder_guide() is a thin convenience wrapper
# around that lookup so a user does not need to know the system.file()
# incantation themselves.
# ==============================================================================

#' Locate or Open the HapBlockR Breeder's Guide
#'
#' \code{HapBlockR_Breeder_Guide.pdf} is a standalone, non-technical
#' companion to the \emph{From Local GEBV to a Crossing Decision} vignette --
#' the same seven parent- and cross-selection tools
#' (\code{\link{truncation_selection}}, \code{\link{select_parents_ga}},
#' \code{\link{usefulness_criterion}}, \code{\link{select_parents_ocs}},
#' \code{\link{select_parents_pareto}}, \code{\link{validate_crosses_exact}},
#' \code{\link{select_core_collection}}), explained in plain language -- what
#' each is, when to reach for it, when to be cautious, and how it works in
#' practice -- with breeding-program analogies and a program-shape decision
#' guide, and no R code required to read it. Intended for breeders, program
#' managers, or reviewers who will use or approve a crossing decision but do
#' not necessarily run R themselves.
#'
#' Distributed as a PDF (not an editable Word document) so that the guide
#' reaches readers as a fixed reference rather than something they might
#' unintentionally edit. A \code{.pdf} file has no vignette engine, so it is
#' not indexed by \code{vignette()}/\code{browseVignettes()} the way the
#' package's \code{.Rmd} vignettes are -- it ships as a static file in
#' \code{inst/extdata/} (the same convention this package uses for
#' \code{beagle.jar} and its \code{example_*} datasets) and is located via
#' \code{\link[base]{system.file}}. This function is a thin convenience
#' wrapper around that lookup.
#'
#' @param open Logical, default \code{TRUE}. If \code{TRUE} and the session
#'   is interactive, opens the file with the operating system's default
#'   application for \code{.pdf} files via \code{\link[utils]{browseURL}}.
#'   If \code{FALSE}, or the session is non-interactive, the file is not
#'   opened -- only its path is returned.
#'
#' @return The file path (character, invisibly) to
#'   \code{HapBlockR_Breeder_Guide.pdf} in the installed package. Errors if
#'   the file cannot be found (e.g. the package was installed with
#'   \code{build_vignettes} machinery that stripped \code{inst/extdata} --
#'   not expected under a normal install, but checked explicitly rather than
#'   than silently returning \code{NULL} or an invalid path).
#'
#' @seealso The \emph{From Local GEBV to a Crossing Decision} vignette
#'   (\code{vignette("HapBlockR-breeding-decisions", package = "HapBlockR")})
#'   for the same seven tools with full statistical detail and runnable R
#'   code; \code{\link{select_parents_ga}}, \code{\link{usefulness_criterion}},
#'   \code{\link{select_parents_ocs}}.
#'
#' @examples
#' guide_path <- open_breeder_guide(open = FALSE)
#' file.exists(guide_path)
#'
#' @export
open_breeder_guide <- function(open = TRUE) {
  guide_path <- system.file("extdata", "HapBlockR_Breeder_Guide.pdf",
                            package = "HapBlockR")
  if (!nzchar(guide_path) || !file.exists(guide_path))
    stop("HapBlockR_Breeder_Guide.pdf was not found in the installed ",
         "package (expected under inst/extdata/). Reinstall HapBlockR, or ",
         "download it directly from the package's GitHub repository.",
         call. = FALSE)

  if (isTRUE(open) && interactive()) {
    message("[open_breeder_guide] Opening HapBlockR_Breeder_Guide.pdf ...")
    utils::browseURL(guide_path)
  } else {
    message("[open_breeder_guide] Guide located at: ", guide_path)
  }

  invisible(guide_path)
}
