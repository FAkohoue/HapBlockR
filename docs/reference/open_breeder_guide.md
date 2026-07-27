# Locate or Open the HapBlockR Breeder's Guide

The PDF and HTML editions of the HapBlockR Breeder's Guide are
standalone companions to the *From Local GEBV to a Crossing Decision*
vignette – the same nine parent- and cross-selection tools
([`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md),
[`select_parents_by_family`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`validate_crosses_exact`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
[`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)),
explained in plain language – what each is, when to use it, its
variants, how it works, and how to interpret its output – with
breeding-programme analogies and a programme-shape decision guide, and
no R code required to read it. It is intended for breeders, programme
managers, or reviewers who will use or approve a crossing decision but
do not necessarily run R themselves.

## Usage

``` r
open_breeder_guide(open = TRUE, format = c("pdf", "html"))
```

## Arguments

- open:

  Logical, default `TRUE`. If `TRUE` and the session is interactive,
  opens the file with the operating system's default application for the
  selected file via
  [`browseURL`](https://rdrr.io/r/utils/browseURL.html). If `FALSE`, or
  the session is non-interactive, the file is not opened – only its path
  is returned.

- format:

  Character, one of `"pdf"` (default) or `"html"`. Selects the tagged
  fixed-layout edition or HTML edition.

## Value

The file path (character, invisibly) to the selected guide edition in
the installed package. The function reports an explicit error if the
requested edition is unavailable.

## Details

The tagged PDF provides a fixed reference edition and the HTML edition
renders the same content for on-screen reading. These standalone files
have no vignette engine, so they are not indexed by
[`vignette()`](https://rdrr.io/r/utils/vignette.html)/[`browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html)
like the package's `.Rmd` vignettes. Their version-controlled source is
`inst/guide/HapBlockR_Breeder_Guide.Rmd`; both generated files ship in
`inst/extdata/` and are located via
[`system.file`](https://rdrr.io/r/base/system.file.html).

## See also

The *From Local GEBV to a Crossing Decision* vignette
([`vignette("HapBlockR-breeding-decisions", package = "HapBlockR")`](https://FAkohoue.github.io/HapBlockR/articles/HapBlockR-breeding-decisions.md))
for the same nine tools with full statistical detail and runnable R
code;
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md).

## Examples

``` r
pdf_path  <- open_breeder_guide(open = FALSE)
#> [open_breeder_guide] Guide located at: C:/Users/fakohoue/AppData/Local/R/win-library/4.5/HapBlockR/extdata/HapBlockR_Breeder_Guide.pdf
html_path <- open_breeder_guide(open = FALSE, format = "html")
#> [open_breeder_guide] Guide located at: C:/Users/fakohoue/AppData/Local/R/win-library/4.5/HapBlockR/extdata/HapBlockR_Breeder_Guide.html
file.exists(pdf_path)
#> [1] TRUE
file.exists(html_path)
#> [1] TRUE
```
