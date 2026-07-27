# Suggest a Starting merit_weight for select_parents_ga_ts()

The block-coverage and whole-genome merit terms used by
[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)
are expressed on dataset-specific scales. This function calibrates their
relative scale directly from the analysed candidate pool. Call it with
`merit_priority` left `NULL` to see the raw diagnostic numbers (the
realistic spread of each term, and the ratio between them), or supply
`merit_priority` (0-100) to also get a literal `merit_weight` value
ready to pass straight into
[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md).
This is the same calculation
[`select_parents_ga_ts()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)'s
own `merit_priority` argument uses internally – calling this function
first just lets you see the numbers before committing to them.

## Usage

``` r
suggest_merit_weight(
  value_matrix,
  merit_score,
  n_founders,
  strategy = c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS"),
  block_weights = NULL,
  merit_priority = NULL
)
```

## Arguments

- value_matrix:

  Numeric matrix (individuals x blocks), identical in shape and meaning
  to
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)'s
  own argument – ideally the exact same, already-filtered matrix you are
  about to pass to that call (after any `min_sel_value`/`top_candidates`
  filtering), so the calibration reflects the real candidate pool the GA
  will search.

- merit_score:

  Named numeric vector, whole-genome merit – same as
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)'s
  argument of the same name. Must cover every individual in
  `value_matrix` and be directionally aligned so that larger values
  always mean greater breeding merit.

- n_founders:

  Integer. Same as
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)'s
  argument of the same name – the founder group size to calibrate for.

- strategy:

  One of `"no_selfing"` (default), `"selfing"`, `"OHS"` (Optimal
  Haplotype Selection), `"OPV"` (Optimal Population Value), or
  `"Haploid_OHS"` – must match the `strategy` you intend to run
  [`select_parents_ga_ts()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)
  with, since it changes how a block's achievable value is computed.

- block_weights:

  Numeric vector, length `ncol(value_matrix)`, or `NULL` (default: equal
  weight 1) – same as
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)'s
  argument of the same name. Must be non-negative.

- merit_priority:

  Numeric in `[0, 100]`, or `NULL` (default). `NULL` returns only the
  diagnostic spread/scale numbers, with `suggested_merit_weight = NULL`.
  A number computes `suggested_merit_weight` too – `0` is always
  equivalent to `merit_weight = 0`; `100` sets merit's spread comparable
  to coverage's spread; values between scale linearly. A zero value is
  useful for scale diagnostics but is not accepted by
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md);
  use
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  for coverage-only selection.

## Value

Named list:

- `merit_span`, `coverage_span`:

  The estimated realistic best-vs-worst achievable spread of each term
  for a group of size `n_founders`, on their own native scales.

- `merit_ceiling`, `merit_floor`, `coverage_ceiling`, `coverage_floor`:

  The four reference values `*_span` is computed from – inspect these
  directly if you want to sanity-check the calibration by hand.

- `scale_factor`:

  `coverage_span / merit_span`. The multiplier that would make merit's
  spread exactly equal to coverage's spread if used as `merit_weight`
  directly (equivalent to `merit_priority = 100`). `NA` if `merit_span`
  was too small to divide by safely (see `ok`).

- `merit_priority`:

  Echoes the argument.

- `suggested_merit_weight`:

  `merit_priority / 100 * scale_factor`, ready to pass to
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)'s
  `merit_weight` argument. `NULL` if `merit_priority` was `NULL`.

- `ok`:

  Logical. `FALSE` if `merit_span` was too close to zero (e.g.
  `merit_score` has almost no spread among these candidates) to safely
  divide by – `suggested_merit_weight` is `NULL` in that case regardless
  of `merit_priority`, with a warning, rather than returning a wild or
  infinite number.

- `trim_n`:

  Integer. How many extreme low points were set aside per term before
  building the "worst achievable" reference group (`0` on small
  candidate pools, where trimming is skipped).

## Details

What "comparable scale" means here, precisely: for a founder group of
size `n_founders`, this function estimates the realistic **best-vs-worst
achievable spread** of the block-coverage sum, and the realistic
best-vs-worst achievable spread of mean `merit_score`, *for this
specific dataset*. `merit_priority = 100` sets `merit_weight` so that
merit's spread becomes comparable in magnitude to coverage's spread;
`merit_priority = 0` is identical to `merit_weight = 0` (no merit term
at all); values in between scale linearly. This gives the breeder a
reproducible definition of relative emphasis based on attainable
contrasts in the supplied data.

The "best achievable" coverage reference is not a naive per-block sum of
each block's own maximum value across all candidates (which is usually
*unreachable* by any single real group, since different blocks' maxima
often belong to different individuals, and would bias the estimate
toward an inflated ceiling). Instead, it is the coverage of an actual,
feasible founder group built by a greedy search: repeatedly add
whichever remaining candidate most improves total coverage. This is not
an arbitrary heuristic – the coverage function is monotone submodular
under a fixed group-size constraint, so greedy is mathematically
guaranteed to land within a known factor (\\1 - 1/e\\, about 63%) of the
true best achievable coverage, using a real, reachable group. The "worst
achievable" references (for both merit and coverage) are built from a
*trimmed* low-ranked group – the bottom `n_founders` candidates after
first setting aside a small, data-size-scaled number of the most extreme
low values – so that a single outlier candidate (e.g. a data-entry
error) cannot single-handedly deflate the floor and distort the
estimated spread. Trimming is skipped entirely on small candidate pools
(under ~20), where it would not be meaningful.

## Interpreting the calibration

The calibration has two interpretation properties:

- Attainable-span scaling:

  `merit_priority` uses the best-vs-worst attainable span of each term.
  The returned diagnostics show the exact contrasts used. A programme
  with an established raw numerical policy may instead supply a positive
  `merit_weight` to
  [`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md).

- A single span number does not capture distribution shape:

  If `merit_score` or the block-coverage values are unusually shaped
  (e.g. strongly bimodal), equal changes in `merit_priority` need not
  yield equal changes in the selected parent set because selection
  depends on candidate combinations, not only marginal distributions.

The calibration scales merit against coverage. Any relationship control
is applied as the separately declared third component of the complete
objective and is reported in `objective_components`.

## See also

[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)
for the joint objective,
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
for coverage-only selection, and
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)
for exploring coverage-relatedness trade-offs.

## Examples

``` r
if (FALSE) { # \dontrun{
res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
top  <- select_top_blocks(res$block_importance, n = 15)
vmat <- res$local_gebv[, top$block_id, drop = FALSE]
cal  <- suggest_merit_weight(vmat, res$gebv, n_founders = 20,
                             merit_priority = 50)
cal$merit_span; cal$coverage_span; cal$suggested_merit_weight
ga_out <- select_parents_ga_ts(
  vmat, n_founders = 20, merit_score = res$gebv,
  merit_weight = cal$suggested_merit_weight
)
} # }
```
