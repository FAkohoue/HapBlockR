# ==============================================================================
# breeder_guide.R
#
# The PDF and HTML editions of the HapBlockR Breeder's Guide are standalone
# companions to the
# *From Local GEBV to a Crossing Decision* vignette: nine parent- and
# cross-selection tools (truncation selection, GA-based parent selection,
# joint GA+TS parent selection, family- or cluster-quota selection,
# usefulness_criterion(), select_parents_ocs(), select_parents_pareto(),
# validate_crosses_exact(), and select_core_collection()), explained in
# plain language for readers who will make or approve a crossing decision
# but do not themselves run R.
#
# Both installed editions are generated from the single R Markdown source
# inst/guide/HapBlockR_Breeder_Guide.Rmd by
# tools/build_breeder_guide_accessible.cjs (a Node/Playwright pipeline, not
# rmarkdown::render()); open_breeder_guide() locates either edition without
# requiring users to know the system.file() path.
# ==============================================================================

#' Locate or Open the HapBlockR Breeder's Guide
#'
#' The PDF and HTML editions of the HapBlockR Breeder's Guide are standalone
#' companions to the \emph{From Local GEBV to a Crossing Decision} vignette --
#' the same nine parent- and cross-selection tools
#' (\code{\link{truncation_selection}}, \code{\link{select_parents_ga}},
#' \code{\link{select_parents_ga_ts}},
#' \code{\link{select_parents_by_family}},
#' \code{\link{usefulness_criterion}}, \code{\link{select_parents_ocs}},
#' \code{\link{select_parents_pareto}}, \code{\link{validate_crosses_exact}},
#' \code{\link{select_core_collection}}), explained in plain language -- what
#' each is, when to use it, its variants, how it works, and how to interpret
#' its output -- with breeding-programme analogies and a programme-shape
#' decision guide, and no R code required to read it. It is intended for
#' breeders, programme managers, or reviewers who will use or approve a
#' crossing decision but do not necessarily run R themselves.
#'
#' The tagged PDF provides a fixed reference edition and the HTML edition
#' renders the same content for on-screen reading. These standalone files
#' have no vignette engine, so they are not indexed by
#' \code{vignette()}/\code{browseVignettes()} like the package's
#' \code{.Rmd} vignettes. Their version-controlled source is
#' \code{inst/guide/HapBlockR_Breeder_Guide.Rmd}; both generated files ship
#' in \code{inst/extdata/} and are located via
#' \code{\link[base]{system.file}}.
#'
#' @param open Logical, default \code{TRUE}. If \code{TRUE} and the session
#'   is interactive, opens the file with the operating system's default
#'   application for the selected file via \code{\link[utils]{browseURL}}.
#'   If \code{FALSE}, or the session is non-interactive, the file is not
#'   opened -- only its path is returned.
#' @param format Character, one of \code{"pdf"} (default) or \code{"html"}.
#'   Selects the tagged fixed-layout edition or HTML edition.
#'
#' @return The file path (character, invisibly) to the selected guide edition
#'   in the installed package. The function reports an explicit error if the
#'   requested edition is unavailable.
#'
#' @seealso The \emph{From Local GEBV to a Crossing Decision} vignette
#'   (\code{vignette("HapBlockR-breeding-decisions", package = "HapBlockR")})
#'   for the same nine tools with full statistical detail and runnable R
#'   code; \code{\link{select_parents_ga}},
#'   \code{\link{select_parents_ga_ts}},
#'   \code{\link{usefulness_criterion}}, \code{\link{select_parents_ocs}}.
#'
#' @examples
#' pdf_path  <- open_breeder_guide(open = FALSE)
#' html_path <- open_breeder_guide(open = FALSE, format = "html")
#' file.exists(pdf_path)
#' file.exists(html_path)
#'
#' @export
open_breeder_guide <- function(open = TRUE, format = c("pdf", "html")) {
  format <- match.arg(format)
  guide_name <- paste0("HapBlockR_Breeder_Guide.", format)
  guide_path <- system.file("extdata", guide_name,
                            package = "HapBlockR")
  if (!nzchar(guide_path) || !file.exists(guide_path))
    stop(guide_name, " was not found in the installed ",
         "package (expected under inst/extdata/). Reinstall HapBlockR, or ",
         "download it directly from the package's GitHub repository.",
         call. = FALSE)

  if (isTRUE(open) && interactive()) {
    message("[open_breeder_guide] Opening ", guide_name, " ...")
    utils::browseURL(guide_path)
  } else {
    message("[open_breeder_guide] Guide located at: ", guide_path)
  }

  invisible(guide_path)
}
