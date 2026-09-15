# Exact (ILP) Cross Selection: A Validation Check for Heuristic Mating Plans

Solves the cross-selection problem – choose exactly `n_cross` crosses
from a candidate list, maximizing a criterion, subject to a per-parent
maximum contribution and an optional relatedness-based culling threshold
– as a binary integer linear program via
[`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html), guaranteeing
the true optimum (not a heuristic approximation). Use this to check how
close a heuristic mating plan (from
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`SimpleMating::selectCrosses()`](https://rdrr.io/pkg/SimpleMating/man/selectCrosses.html)/
`GOCS()`, or AlphaMate) came to the best achievable plan under the same
constraints, on a small enough candidate set that solving exactly is
practical.

## Usage

``` r
validate_crosses_exact(
  data,
  n_cross,
  max_cross = NULL,
  culling_pairwise_k = NULL,
  parent1_col = "parent1",
  parent2_col = "parent2",
  criterion_col = "UC",
  relatedness_col = NULL,
  G = NULL,
  heuristic_plan = NULL,
  max_vars = 2000L,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame of candidate crosses (e.g.
  [`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
  output, or a
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)/`selectCrosses()`
  table).

- n_cross:

  Integer. Exact number of crosses the plan must contain.

- max_cross:

  Optional integer. Maximum number of crosses any single parent can
  participate in. `NULL` (default) leaves this unconstrained – the true
  optimum may then concentrate heavily on very few parents; supply the
  same value your heuristic plan used for an apples-to-apples
  comparison.

- culling_pairwise_k:

  Optional numeric. Candidate crosses with relatedness above this value
  are excluded before solving (matching
  `SimpleMating`/[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)'s
  culling convention). Requires `relatedness_col` (already in `data`) or
  `G` (to compute relatedness per cross on the fly).

- parent1_col, parent2_col:

  Character, default `"parent1"`/ `"parent2"`. Column names in `data`
  identifying each cross's two parents.

- criterion_col:

  Character, default `"UC"`. Column in `data` to maximize (e.g. `"UC"`
  from
  [`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
  or `"mid_parent_gebv"`).

- relatedness_col:

  Optional character. Column in `data` giving each cross's relatedness
  value, used for `culling_pairwise_k`. If `NULL` and `G` is supplied,
  computed automatically as `G[parent1, parent2]` per row.

- G:

  Optional dimnamed relationship matrix, used to compute relatedness per
  cross when `relatedness_col` is not already in `data`.

- heuristic_plan:

  Optional data frame (same `parent1_col`/ `parent2_col` convention) – a
  heuristic mating plan to compare against the exact optimum. If
  supplied, the return value includes the heuristic plan's total
  criterion and its percentage gap below the exact optimum.

- max_vars:

  Integer, default `2000L`. Safety cap on the number of candidate
  crosses (after culling) the solver will attempt – integer programming
  is exact but can become slow well before this in the worst case; lower
  it if solving is too slow, or pre-filter `data` (e.g. a tighter
  `culling_pairwise_k`, or restrict to your top-ranked crosses by
  `criterion_col`) rather than raising it blindly.

- verbose:

  Logical, default `TRUE`.

## Value

A list inheriting from `HapBlockR_exact_cross_validation` and
`hapblockr_result`, with `exact_plan` (data frame: the optimal cross
selection, a subset of `data`'s rows), `exact_objective` (the true
optimal total criterion), `n_candidates` (candidate crosses considered
after culling), `status` (lpSolve's solver status; `0` = optimal
solution found), and, if `heuristic_plan` was supplied,
`heuristic_objective` and `gap_pct` (the heuristic plan's percentage
shortfall below the exact optimum).

## What this is and is not

This is a validation/sanity-check tool, not a replacement for
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
in normal use. Integer programming does not scale to the size of a real
candidate cross list (all pairwise combinations of a large parent set
can easily run into the tens of thousands of candidate crosses) the way
a GA or greedy heuristic does – see `max_vars`. A per-parent MINIMUM
contribution (SimpleMating's `min.cross`) is deliberately not supported
as a hard constraint; see the source-level comment in
`R/exact_validation.R` for why.

## See also

[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
