# Suggest a Starting merit_weight for select_parents_ga()

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
`merit_weight` and `coancestry_weight` arguments have no universal
correct value: the block-coverage term and the merit term live on
different, problem-specific scales, so a raw multiplier that works for
one dataset can be meaningless for another. This function estimates a
sensible starting point directly from your own data, in two ways: call
it with `merit_priority` left `NULL` to see the raw diagnostic numbers
(the realistic spread of each term, and the ratio between them), or
supply `merit_priority` (0-100, "how much do you care about merit vs.
coverage") to also get a literal `merit_weight` value ready to pass
straight into
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md).
This is the same calculation
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
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

  Numeric matrix (individuals x blocks), identical in shape/meaning to
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
  own argument – ideally the exact same, already-filtered matrix you are
  about to pass to that call (after any `min_sel_value`/`top_candidates`
  filtering), so the calibration reflects the real candidate pool the GA
  will search.

- merit_score:

  Named numeric vector, whole-genome merit – same as
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
  argument of the same name. Must cover every individual in
  `value_matrix`.

- n_founders:

  Integer. Same as
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
  argument of the same name – the founder group size to calibrate for.

- strategy:

  One of `"no_selfing"` (default), `"selfing"`, `"OHS"`, `"OPV"`,
  `"Haploid_OHS"` – must match the `strategy` you intend to run
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  with, since it changes how a block's achievable value is computed.

- block_weights:

  Numeric vector, length `ncol(value_matrix)`, or `NULL` (default: equal
  weight 1) – same as
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
  argument of the same name. Must be non-negative.

- merit_priority:

  Numeric in `[0, 100]`, or `NULL` (default). `NULL` returns only the
  diagnostic spread/scale numbers, with `suggested_merit_weight = NULL`.
  A number computes `suggested_merit_weight` too – `0` is always
  equivalent to `merit_weight = 0`; `100` sets merit's spread comparable
  to coverage's spread; values between scale linearly.

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
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
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
at all); values in between scale linearly. This is one reasonable,
explicitly-stated definition of "comparable" – not the only possible one
(matching standard deviation instead of spread, for instance, would give
a different number) – see *What this does not solve* below.

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

## What this does not solve

Two limits are inherent to *any* scale-matching approach, not specific
to the method used here, and cannot be resolved by more engineering –
they are documented rather than hidden:

- Matching spread is a choice, not a universal truth:

  A different, equally defensible definition of "comparable" (e.g.
  matching standard deviation across many realistic groups, rather than
  the best-vs-worst achievable span) would produce a different scale
  factor. Treat `merit_priority`'s suggestion as a well-reasoned
  starting point to inspect and adjust, not a uniquely correct answer –
  the literal `merit_weight` argument remains available in
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  for full manual control.

- A single span number does not capture distribution shape:

  If `merit_score` or the block-coverage values are unusually shaped
  (e.g. strongly bimodal), the dial's practical effect may not feel
  perfectly linear across its 0-100 range even though the underlying
  calculation is exact for what it measures.

This calibration also only weighs merit against coverage; if
`coancestry_weight` is also active in your
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
call, its effect is held fixed rather than jointly recalibrated – use
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)'s
sweep to explore that trade-off separately, as already recommended for
tuning `coancestry_weight` on its own.

## See also

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
*Merit-weighted fitness (GA+TS hybrid, optional)* section for the
fitness function this feeds into, and
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)
for exploring the `coancestry_weight` trade-off the same way.

## Examples

``` r
if (FALSE) { # \dontrun{
res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
top  <- select_top_blocks(res$block_importance, n = 15)
vmat <- res$local_gebv[, top$block_id, drop = FALSE]
cal  <- suggest_merit_weight(vmat, res$gebv, n_founders = 20,
                             merit_priority = 50)
cal$merit_span; cal$coverage_span; cal$suggested_merit_weight
ga_out <- select_parents_ga(vmat, n_founders = 20, merit_score = res$gebv,
                            merit_weight = cal$suggested_merit_weight)
} # }
```
