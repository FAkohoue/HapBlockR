# Pareto Frontier of Parent Sets: Merit vs. Relatedness

Sweeps
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
`coancestry_weight` across a grid of values, runs the GA at each one,
and Pareto-filters the resulting (merit, relatedness) points into an
empirical frontier – letting you see the actual gain-vs-diversity
tradeoff curve for your candidate population and pick a point off it
deliberately, rather than guessing a single `coancestry_weight` value
and hoping it was the right one.

## Usage

``` r
select_parents_pareto(
  value_matrix,
  n_founders,
  strategy = c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS"),
  block_weights = NULL,
  top_candidates = NULL,
  G,
  coancestry_weights = c(0, 0.25, 0.5, 1, 2, 4),
  merit = NULL,
  popSize = 100L,
  maxiter = 200L,
  run = 50L,
  pmutation = 0.1,
  pcrossover = 0.8,
  n_reps = 3L,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- value_matrix, n_founders, strategy, block_weights, top_candidates,
  popSize, maxiter, run, pmutation, pcrossover:

  Passed through to
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  at every grid point; see its documentation.

- G:

  Relationship/kinship matrix, required (used both for the coancestry
  penalty during each GA run and to report each resulting set's realised
  `mean_relationship` for the frontier, including at
  `coancestry_weight = 0`).

- coancestry_weights:

  Numeric vector of `coancestry_weight` values to sweep. Default
  `c(0, 0.25, 0.5, 1, 2, 4)` – a broad first pass; problem-specific,
  since block-coverage and relationship scores have no common natural
  scale (same caveat as `coancestry_weight` itself in
  [`?select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)).

- merit:

  Optional named numeric vector (e.g. whole-genome GEBV or a selection
  index) used only to report each resulting set's `mean_merit` on the
  frontier and as the merit axis for Pareto filtering. If `NULL`
  (default), the GA's own block-coverage `fitness` is used as the merit
  axis instead – a valid but less directly interpretable stand-in for
  whole-genome merit.

- n_reps:

  Integer, default `3L`. Passed to
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  at each grid point (lower than that function's own default of 5, to
  keep the sweep's total runtime reasonable).

- seed:

  Optional integer. If supplied, grid point `i` uses `seed + i - 1L` for
  reproducibility across the sweep.

- verbose:

  Logical, default `TRUE`.

## Value

A list with:

- `frontier`:

  Data frame, one row per grid point, with `coancestry_weight`,
  `fitness`, `mean_relationship`, `mean_merit` (`NA` if `merit` not
  supplied), `n_selected`, `converged`, `run_index` (row index into
  `runs`), and the `pareto_optimal`/ `crowding_distance` columns from
  [`pareto_front`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md).

- `runs`:

  List of the full
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  return value at each grid point (in original sweep order, indexed by
  `frontier$run_index`) – use this to get the actual `$selected`
  individual IDs for any frontier point you choose.

## Details

This does not implement a new multi-objective search algorithm (e.g.
NSGA-II) – it reuses
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
(already implemented, already verified) as the underlying solver at each
grid point, and
[`pareto_front`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md)
to filter the resulting points down to the actually non-dominated ones
(a weight sweep is not guaranteed to only produce non-dominated results,
since each run is itself a heuristic GA search – filtering removes sweep
points that turned out to be strictly worse than another point in both
merit and relatedness).

Runtime is the sum of every grid point's GA run – with the default 6
grid points and `n_reps = 3`, that is 18 GA searches. Lower `n_reps` or
narrow `coancestry_weights` for faster, coarser sweeps; the defaults
favour a broad first look over speed.

## See also

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`pareto_front`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
