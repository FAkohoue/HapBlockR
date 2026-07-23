# Funnel Plot of Block Importance

Visualises every block's local GEBV / haplotype effect values against
its scaled variance, mirroring HapSelect's funnel plot: local
GEBV/haplotype effects on the x-axis, block variance (min-max scaled to
\[0,1\]) on the y-axis. Lets a breeder *see* the ranked-variance elbow
before picking a cutoff for
[`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
or
[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)'s
`importance_threshold`, rather than choosing one blind.

## Usage

``` r
plot_block_funnel(
  local_gebv,
  block_importance,
  highlight_threshold = 0.9,
  max_blocks = 2000L
)
```

## Arguments

- local_gebv:

  Numeric matrix (individuals x blocks) of per-block local
  GEBV/haplotype-effect values, e.g.
  `run_haplotype_prediction()$local_gebv` (or one element of
  `$local_gebv` for a multi-trait result) or
  `compute_local_gebv()$local_gebv`.

- block_importance:

  Data frame with `block_id` and `var_scaled` columns (or
  `var_local_gebv`, rescaled internally), aligned to the columns of
  `local_gebv` by `block_id`, e.g.
  `compute_local_gebv()$block_importance`.

- highlight_threshold:

  Numeric in (0,1\] or `NULL`. If supplied, points from blocks at or
  above this scaled-variance value are drawn in a second colour, so the
  current `importance_threshold` cutoff is visible on the funnel.
  Default `0.9`.

- max_blocks:

  Integer or `NULL`. If the number of blocks in `local_gebv` exceeds
  this, only the `max_blocks` highest-variance blocks are plotted
  (funnel plots with many thousands of points become unreadable and slow
  to render). Default `2000L`.

## Value

A `ggplot2` object.

## See also

[`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md),
[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
plot_block_funnel(res$local_gebv, res$block_importance)
} # }
```
