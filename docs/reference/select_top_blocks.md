# Select Top-Ranked Blocks by Count, Percentage, or Cumulative Variance

Filters a block-importance table down to the most informative blocks,
with the same three selection modes as HapSelect's
`select_top_blocks()`:

- `n`:

  Fixed count: the top `n` blocks by rank.

- `perc_total`:

  Fixed percentage of *all* blocks (e.g. `0.5` keeps the top 50% of
  blocks by rank, regardless of how much variance they explain).

- `perc_of_total_var`:

  The smallest top-ranked set whose *cumulative* variance share reaches
  `perc_of_total_var` (e.g. `0.9` keeps just enough top blocks to
  explain 90% of the total scaled variance across all blocks) – lets the
  data decide how many blocks matter, rather than fixing a count or
  percentage up front.

Exactly one of `n`, `perc_total`, `perc_of_total_var` must be supplied.
Unlike a fixed `importance_threshold` cutoff (e.g. in
[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)),
all three modes here always return a usable, non-empty selection (as
long as at least one block has positive variance), which is useful when
a hardcoded threshold would return zero or an unexpectedly large number
of "important" blocks for a given dataset. Use
[`plot_block_funnel`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md)
to see the ranked-variance curve before choosing a cutoff.

## Usage

``` r
select_top_blocks(
  block_importance,
  n = NULL,
  perc_total = NULL,
  perc_of_total_var = NULL,
  var_col = NULL
)
```

## Arguments

- block_importance:

  Data frame with a `var_scaled` or `var_local_gebv` column, e.g. from
  `compute_local_gebv()$block_importance`,
  `run_haplotype_prediction()$block_importance`, or
  `rank_haplotype_blocks()$ranked_blocks` (via its `var_scaled` column,
  when a `pred_result` was supplied).

- n:

  Integer or `NULL`. Fixed count mode.

- perc_total:

  Numeric in (0,1\] or `NULL`. Fixed-percentage-of-blocks mode.

- perc_of_total_var:

  Numeric in (0,1\] or `NULL`. Cumulative-variance mode.

- var_col:

  Character or `NULL`. Column to rank by. Default `NULL`: uses
  `"var_scaled"` if present, else `"var_local_gebv"`.

## Value

`block_importance`, sorted by descending `var_col` and subset to the
selected blocks, with an added integer `rank` column (1 = highest
variance) and a `cum_var_share` column (cumulative fraction of total
variance explained up to and including that block's rank, computed over
the full input before subsetting).

## See also

[`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md),
[`rank_haplotype_blocks`](https://FAkohoue.github.io/HapBlockR/reference/rank_haplotype_blocks.md),
[`plot_block_funnel`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
top15   <- select_top_blocks(res$block_importance, n = 15)
top_50p <- select_top_blocks(res$block_importance, perc_total = 0.5)
top_90v <- select_top_blocks(res$block_importance, perc_of_total_var = 0.9)
} # }
```
