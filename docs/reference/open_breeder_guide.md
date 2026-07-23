# Locate or Open the HapBlockR Breeder's Guide

`HapBlockR_Breeder_Guide.pdf` is a standalone, non-technical companion
to the *From Local GEBV to a Crossing Decision* vignette – the same
seven parent- and cross-selection tools
([`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`validate_crosses_exact`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
[`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)),
explained in plain language – what each is, when to reach for it, when
to be cautious, and how it works in practice – with breeding-program
analogies and a program-shape decision guide, and no R code required to
read it. Intended for breeders, program managers, or reviewers who will
use or approve a crossing decision but do not necessarily run R
themselves.

## Usage

``` r
open_breeder_guide(open = TRUE)
```

## Arguments

- open:

  Logical, default `TRUE`. If `TRUE` and the session is interactive,
  opens the file with the operating system's default application for
  `.pdf` files via
  [`browseURL`](https://rdrr.io/r/utils/browseURL.html). If `FALSE`, or
  the session is non-interactive, the file is not opened – only its path
  is returned.

## Value

The file path (character, invisibly) to `HapBlockR_Breeder_Guide.pdf` in
the installed package. Errors if the file cannot be found (e.g. the
package was installed with `build_vignettes` machinery that stripped
`inst/extdata` – not expected under a normal install, but checked
explicitly rather than than silently returning `NULL` or an invalid
path).

## Details

Distributed as a PDF (not an editable Word document) so that the guide
reaches readers as a fixed reference rather than something they might
unintentionally edit. A `.pdf` file has no vignette engine, so it is not
indexed by
[`vignette()`](https://rdrr.io/r/utils/vignette.html)/[`browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html)
the way the package's `.Rmd` vignettes are – it ships as a static file
in `inst/extdata/` (the same convention this package uses for
`beagle.jar` and its `example_*` datasets) and is located via
[`system.file`](https://rdrr.io/r/base/system.file.html). This function
is a thin convenience wrapper around that lookup.

## See also

The *From Local GEBV to a Crossing Decision* vignette
([`vignette("HapBlockR-breeding-decisions", package = "HapBlockR")`](https://FAkohoue.github.io/HapBlockR/articles/HapBlockR-breeding-decisions.md))
for the same seven tools with full statistical detail and runnable R
code;
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md).

## Examples

``` r
guide_path <- open_breeder_guide(open = FALSE)
#> Error: HapBlockR_Breeder_Guide.pdf was not found in the installed package (expected under inst/extdata/). Reinstall HapBlockR, or download it directly from the package's GitHub repository.
file.exists(guide_path)
#> Error: object 'guide_path' not found
```
