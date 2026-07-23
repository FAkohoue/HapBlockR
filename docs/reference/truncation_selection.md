# Truncation Selection (Top-n by a Single Score)

The simplest possible parent-selection rule and the standard baseline
breeders compare any smarter method against: rank individuals by a
single genome-wide score and keep the top `n_founders`. Used here as the
"TS" comparison arm against
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
(see
[`ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)),
and standalone whenever a plain best-GEBV shortlist is all that's
needed.

## Usage

``` r
truncation_selection(
  score,
  n_founders,
  min_sel_value = NULL,
  min_sel_mode = c("value", "percentile", "sd_below_mean")
)
```

## Arguments

- score:

  Named numeric vector, e.g. whole-genome GEBV
  (`run_haplotype_prediction()$gebv`) or a stacking index
  (`score_favorable_haplotypes()$stacking_index`). Names are individual
  IDs.

- n_founders:

  Integer. Number of individuals to select.

- min_sel_value:

  Numeric or `NULL` (default `NULL` = no floor, every candidate with a
  finite `score` is eligible). Excludes candidates below a merit floor
  *before* ranking, so that a plain top-`n_founders` call never has to
  hit deeper into the population than your program's own quality bar.
  Interpreted according to `min_sel_mode`. Named `min_sel_value` rather
  than `min_selection_index` because `score` need not be a selection
  index – it can be any single genome-wide value.

- min_sel_mode:

  One of `"value"` (default), `"percentile"`, `"sd_below_mean"`. Only
  used when `min_sel_value` is not `NULL`:

  `"value"`

  :   `min_sel_value` is an absolute cutoff on `score` itself (same
      units/scale as `score` – requires knowing that scale in advance).

  `"percentile"`

  :   `min_sel_value` in `(0, 1]` is the fraction of candidates to keep
      from the top, e.g. `0.6` keeps the top 60% by `score`.
      Self-scaling – no need to know `score`'s units.

  `"sd_below_mean"`

  :   `min_sel_value` is the number of standard deviations below
      `mean(score)` the cutoff sits, e.g. `1` keeps everyone
      within/above one SD of the mean. Also self-scaling.

## Value

Named list:

- `selected`:

  Character vector of the top `n_founders` individual IDs (after the
  `min_sel_value` floor, if set), sorted by descending `score`.

- `score`:

  The corresponding scores, same order.

- `cutoff`:

  Numeric. The merit-floor cutoff actually applied (`-Inf` when
  `min_sel_value = NULL`).

## What this does and does not do

This function ranks by a single number and takes the top `n_founders` –
nothing more. It has no notion of which haplotype blocks each individual
carries, no complementarity check between the individuals it selects
(two selected parents could carry identical favourable haplotypes and
none of the ones another candidate has), and no relatedness/coancestry
management. That simplicity is the point: it is the baseline every
smarter method – above all
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
– needs to outperform to justify its added complexity. See
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
*Choosing between this function and truncation_selection()* section for
a full comparison and a decision guide, and the *From Local GEBV to a
Crossing Decision* vignette for a worked example running both side by
side.

## See also

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
ts  <- truncation_selection(res$gebv, n_founders = 20)
ts$selected
} # }
```
