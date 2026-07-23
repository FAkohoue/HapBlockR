# Core-Collection / Diversity-Maximizing Subset Selection

Selects a subset of `n_core` individuals from a candidate panel that
best represents its genetic diversity – the classical core-collection
problem (Schoen & Brown 1993; Gonzalez 1985), useful for genebank
curation, building a diverse training/reference panel, or choosing a
broad founder set for a new diversification programme. Unlike this
package's other parent-selection tools, diversity itself is the primary
objective here, not a penalty on top of a merit criterion – though an
optional merit floor lets you restrict to "the most diverse subset among
candidates that already clear a merit bar" (see `min_sel_value`).

## Usage

``` r
select_core_collection(
  G,
  n_core,
  type = c("relationship", "distance"),
  strategy = c("maximin", "mean_distance"),
  merit = NULL,
  min_sel_value = NULL,
  min_sel_mode = c("value", "percentile", "sd_below_mean"),
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- G:

  Dimnamed relationship or distance matrix (row/column names =
  individual IDs), e.g. from
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md).

- n_core:

  Integer. Size of the core collection to select.

- type:

  Character, one of `"relationship"` (default) or `"distance"`. If
  `"relationship"`, `G` is converted to a genetic distance matrix via
  \\D\_{ij} = G\_{ii} + G\_{jj} - 2G\_{ij}\\ (the exact identity
  relating a Gram/relationship matrix to squared Euclidean distance in
  the space it represents – not an approximation). If `"distance"`, `G`
  is used as a distance matrix directly.

- strategy:

  Character, one of `"maximin"` (default) or `"mean_distance"`. See
  Details above.

- merit:

  Optional named numeric vector (e.g. GEBV or a selection index), names
  = individual IDs. Only used together with `min_sel_value` to
  pre-filter the candidate pool; does not otherwise influence which
  individuals are chosen (this function optimizes diversity, not merit,
  among whichever candidates remain eligible).

- min_sel_value, min_sel_mode:

  Optional merit floor applied to `merit` before diversity selection,
  via the same `.apply_merit_floor()` logic used by
  [`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)/[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  – `min_sel_mode` one of `"value"`, `"percentile"`, `"sd_below_mean"`.
  Both ignored if `merit` is `NULL`.

- seed:

  Optional integer. Currently only relevant for the degenerate
  `n_core = 1` case (no merit supplied), where the single selected
  individual is otherwise chosen at random; included for reproducibility
  and API consistency with this package's other selection functions.

- verbose:

  Logical, default `TRUE`.

## Value

A list with `selected` (character vector of chosen individual IDs, in
selection order), `n_core`, `strategy`, `mean_distance` and
`min_distance` (of the final selected set), and `trace` (data frame, one
row per selection step: `step`, `added`, `criterion` – the maximin
distance or mean-distance value achieved at that step, useful for
plotting how diversity accumulates as the core collection grows).

## Two greedy strategies

- `"maximin"`:

  (Default.) Farthest-point traversal (Gonzalez 1985): start from the
  single most genetically distant pair, then repeatedly add whichever
  remaining candidate has the largest MINIMUM distance to everything
  already selected. Maximizes the smallest pairwise distance in the
  final set (guards against near-duplicate individuals slipping in),
  with a proven 2-approximation guarantee. The standard "M strategy" in
  the core-collection literature.

- `"mean_distance"`:

  At each step, adds whichever remaining candidate most increases the
  selected set's MEAN pairwise distance (the "MD strategy"). Tends to
  spread the selection more evenly across the whole diversity space
  rather than prioritising the single most extreme outliers; more
  expensive to compute (recomputed per candidate per step) but still
  tractable for typical core-collection sizes.

Neither is an exact solver for its objective (both are the standard
greedy heuristics used throughout the core-collection literature, not a
guaranteed-optimal search) – see
[`validate_crosses_exact`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)
if you need a true optimum on a small-enough problem (that function
targets cross selection specifically, not this subset-diversity problem,
but shares the same "exact validation of a heuristic" philosophy).

## References

Schoen, D.J. & Brown, A.H.D. (1993). Conservation of allelic richness in
wild crop relatives is aided by assessment of genetic markers.
*Proceedings of the National Academy of Sciences*, 90, 10623-10627.

Gonzalez, T.F. (1985). Clustering to minimize the maximum intercluster
distance. *Theoretical Computer Science*, 38, 293-306.

## See also

[`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md),
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`validate_crosses_exact`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)
