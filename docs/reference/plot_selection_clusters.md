# Plot Individuals Coloured by Genetic Cluster

Companion plot to
[`cluster_selection_groups`](https://FAkohoue.github.io/HapBlockR/reference/cluster_selection_groups.md):
PC1 vs PC2 of the same (variance-scaled) retained-PC subspace used for
clustering, coloured by cluster membership – distinct from
[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md),
which colours by GA/TS/Both/ Neither selection-group membership instead.
Read them together: this plot shows where the genetic clusters are;
`cluster_selection_groups() $table` shows how well each selection method
represents each cluster.

## Usage

``` r
plot_selection_clusters(cluster_res, save_path = NULL, width = 7, height = 5.5)
```

## Arguments

- cluster_res:

  List returned by
  [`cluster_selection_groups`](https://FAkohoue.github.io/HapBlockR/reference/cluster_selection_groups.md).

- save_path:

  Optional file path ending in `.pdf`. When supplied, the plot is also
  saved via
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  at 300 dpi. `NULL` (default) does not save to disk.

- width, height:

  Numeric, PDF dimensions in inches, used only when `save_path` is
  supplied. Defaults `7`, `5.5`.

## Value

A `ggplot2` object (returned whether or not `save_path` is supplied).

## See also

[`cluster_selection_groups`](https://FAkohoue.github.io/HapBlockR/reference/cluster_selection_groups.md),
[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cl <- cluster_selection_groups(res$G, ga_selected = ga$selected,
                               ts_selected = ts$selected, n_clusters = 3)
plot_selection_clusters(cl, save_path = "clusters.pdf")
} # }
```
