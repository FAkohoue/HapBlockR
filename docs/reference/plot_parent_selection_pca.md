# PCA of Parent-Selection Groups

Projects every genotyped individual onto the top two principal
components, coloured by which selection group they belong to
(GA-selected, truncation-selected, both, or neither) – mirrors
HapSelect's parent-set PCA plot, which shows where each selection
strategy's founders sit relative to the whole population's diversity.

Two different diversity spaces can be plotted, and they answer different
questions:

- **`G` (default)**: genome-wide (or whole-haplotype-set)
  relationship-matrix PCA, via eigen-decomposition of `G` itself.
  Answers "do these selection groups also look diversity-covering across
  the whole genome?" – a useful cross-check, but `G` is *not* the space
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  actually searched.

- **`feature_matrix`**: PCA (via
  [`prcomp`](https://rdrr.io/r/stats/prcomp.html)) of the individuals x
  target-blocks local-GEBV/haplotype-value matrix – the literal
  `value_matrix` argument
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  optimised coverage over. Answers "does GA's selection actually spread
  out in the space it was asked to spread out in?" directly, at the cost
  of only reflecting the target blocks, not the rest of the genome.

Supply exactly one of the two; `G` remains the default so existing calls
are unaffected.

## Usage

``` r
plot_parent_selection_pca(
  G = NULL,
  ga_selected,
  ts_selected,
  feature_matrix = NULL
)
```

## Arguments

- G:

  Genomic or haplotype relationship matrix (n x n), dimnames =
  individual IDs, e.g. `run_haplotype_prediction()$G` or
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
  output. Leave `NULL` (and supply `feature_matrix` instead) to plot in
  block-feature space instead of genome-wide diversity space.

- ga_selected:

  Character vector of individual IDs selected by
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  (its `$selected`).

- ts_selected:

  Character vector of individual IDs selected by
  [`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
  (its `$selected`).

- feature_matrix:

  Individuals x target-blocks numeric matrix (dimnames rows = individual
  IDs), e.g. the same `value_matrix` passed to
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  (`res$local_gebv[, top_blocks$block_id, drop = FALSE]`). When
  supplied, PCA is computed on this matrix (centred, not scaled –
  local-GEBV columns are already on a common trait-value scale) instead
  of eigen-decomposing `G`. Leave `NULL` (the default) to use `G`
  instead. Supplying both `G` and `feature_matrix` is an error – pick
  one diversity space per plot.

## Value

A `ggplot2` object.

## See also

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
top <- select_top_blocks(res$block_importance, n = 15)
vmat <- res$local_gebv[, top$block_id, drop = FALSE]
ga  <- select_parents_ga(vmat, n_founders = 20, seed = 1)
ts  <- truncation_selection(res$gebv, n_founders = 20)

# Genome-wide diversity cross-check (default):
plot_parent_selection_pca(res$G, ga$selected, ts$selected)

# The actual space GA searched (target-block local GEBV):
plot_parent_selection_pca(ga_selected = ga$selected, ts_selected = ts$selected,
                          feature_matrix = vmat)
} # }
```
