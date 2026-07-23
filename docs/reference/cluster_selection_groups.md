# Cluster Individuals on Retained Principal Components and Cross-Tabulate Selection Coverage

[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)
only ever looks at PC1/PC2, which can be a small fraction of total
genetic variance in a genetically complex/structured population (e.g. 3
distinct clusters captured by only ~26% cumulative variance on PC1+PC2
combined). This function instead retains as many leading principal
components as needed to reach a user-set cumulative-variance threshold
(default 95%), clusters individuals in that higher-dimensional subspace
(hierarchical or k-means), and cross-tabulates an arbitrary number of
named selection strategies (GA, truncation, OCS, core-collection, or any
other ID vector you have) against cluster membership – turning "does
this method represent every genetically distinct sub-group?" into a
table of counts and proportions per cluster, per strategy, instead of a
visual impression from a 2-axis scatter. Individuals may belong to more
than one strategy at once (e.g. selected by both GA and OCS) – each
strategy's representation is counted independently, not as a single
mutually exclusive category.

## Usage

``` r
cluster_selection_groups(
  G = NULL,
  feature_matrix = NULL,
  groups,
  variance_threshold = 0.95,
  method = c("hierarchical", "kmeans"),
  n_clusters,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- G:

  Genomic or haplotype relationship matrix (n x n), dimnames =
  individual IDs. Same duality as
  [`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md):
  leave `NULL` and supply `feature_matrix` instead to cluster in
  target-block feature space rather than genome-wide diversity space.
  **Unlike**
  [`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md),
  the retained PC scores here are eigenvectors scaled by
  `sqrt(eigenvalue)` (the standard PCA-score convention), not raw
  unit-norm eigenvectors – necessary so that Euclidean distance across
  many retained PCs is properly variance-weighted for clustering; a
  2-axis scatter plot does not need this because each axis is
  independently labelled with its own % variance.

- feature_matrix:

  Individuals x target-blocks numeric matrix, e.g. the same
  `value_matrix` passed to
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md).
  Supply exactly one of `G`/`feature_matrix`.

- groups:

  Named list of character vectors, one per selection strategy to check,
  e.g.
  `list(GA = ga$selected, TS = ts$selected, OCS = ocs_parents, Core = core_res$selected)`.
  Names label that strategy's columns in the returned `table`; IDs not
  present in `G`/ `feature_matrix` are ignored (with a message).

- variance_threshold:

  Numeric in (0, 1\]. Retain the smallest number of leading PCs whose
  cumulative proportion of variance is `>=` this threshold (always at
  least 2, so
  [`plot_selection_clusters`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md)
  can still show a PC1/PC2 view of the retained subspace). Default
  `0.95`.

- method:

  `"hierarchical"` (default –
  [`hclust`](https://rdrr.io/r/stats/hclust.html) with Ward's
  minimum-variance linkage (`"ward.D2"`), then
  [`cutree`](https://rdrr.io/r/stats/cutree.html) to `n_clusters`
  groups; deterministic, no seed needed) or `"kmeans"`
  ([`kmeans`](https://rdrr.io/r/stats/kmeans.html), `nstart = 25`; needs
  `seed` for reproducibility).

- n_clusters:

  Integer `>= 2`. Number of clusters to cut the population into.
  **Required** – this package does not silently guess a number of
  clusters for you. If you don't already have a biological reason for a
  specific number (e.g. "this program has 3 founder families"), a
  reasonable starting point is to try a small range (e.g. 2-6) and look
  at `cluster_fit` (the `hclust`/ `kmeans` object returned) with
  standard diagnostics (`plot(cluster_fit)` for a dendrogram; `kmeans`'s
  `tot.withinss` across a range of `k` for an elbow plot) before
  settling on `n_clusters`.

- seed:

  Integer or `NULL`. Only used by `method = "kmeans"`; ignored (with a
  message if supplied) for `"hierarchical"`, which is deterministic.

- verbose:

  Logical. Default `TRUE`.

## Value

Named list:

- `cluster`:

  Named factor (names = individual IDs), cluster membership, levels
  `"Cluster 1"` .. `"Cluster n_clusters"`.

- `table`:

  Data frame, one row per cluster: `n_total` (population size in that
  cluster), `pct_population` (that cluster's % share of the whole
  population), then for every name in `groups` a `n_<name>` (how many of
  that strategy's selections fall in this cluster) and `prop_<name>`
  (`n_<name> / n_total` – what fraction of *this cluster* that strategy
  selected) – this is the actual answer to "which clusters does strategy
  X under-represent?".

- `n_pcs_retained`:

  Integer, how many leading PCs met `variance_threshold`.

- `variance_explained`:

  Numeric, actual cumulative proportion of variance captured by
  `n_pcs_retained` PCs (`>=` `variance_threshold`).

- `variance_threshold`, `method`:

  Echoed back.

- `space`:

  Character, `"population diversity space"` or
  `"target-block feature space"` depending on whether `G` or
  `feature_matrix` was supplied.

- `groups`:

  Character vector, the `names(groups)` echoed back (the strategy labels
  used as `table` column suffixes).

- `pc_scores`:

  The retained, variance-scaled PC score matrix (individuals x
  `n_pcs_retained`) actually clustered on – input to
  [`plot_selection_clusters`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md).

- `cluster_fit`:

  The raw `hclust` or `kmeans` object, for diagnostics (dendrograms,
  within-cluster sum of squares, etc.).

## See also

[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md),
[`plot_selection_clusters`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md),
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
ga  <- select_parents_ga(res$local_gebv, n_founders = 20, seed = 1)
ts  <- truncation_selection(res$gebv, n_founders = 20)
cl  <- cluster_selection_groups(res$G, groups = list(GA = ga$selected,
                                                     TS = ts$selected),
                                variance_threshold = 0.95,
                                method = "hierarchical", n_clusters = 3)
cl$table
} # }
```
