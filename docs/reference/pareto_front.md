# Pareto (Non-Dominated) Front of a Set of Candidates

Given a data frame of candidates (parents, crosses, or anything else)
each scored on two or more objectives, flags which candidates are
Pareto-optimal (non-dominated): no other candidate is at least as good
on every objective and strictly better on at least one. This is the
general tool underneath
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
but works on any objective columns you give it – e.g.
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
output scored on `mid_parent_gebv` (maximize) and `predicted_variance`
(context-dependent), or a multi-trait selection index's separate trait
columns.

## Usage

``` r
pareto_front(data, objectives, directions = "max")
```

## Arguments

- data:

  Data frame of candidates.

- objectives:

  Character vector of column names in `data` to treat as objectives.

- directions:

  Character vector, same length as `objectives` (or length 1, recycled),
  each `"max"` or `"min"`. Default `"max"` for every objective.

## Value

`data` with two columns appended: `pareto_optimal` (logical) and
`crowding_distance` (numeric, `NA` for dominated candidates), sorted
with Pareto-optimal candidates first (by descending crowding distance
among them).

## Details

Uses the standard \\O(n^2 \times k)\\ pairwise-dominance algorithm
(\\n\\ = candidates, \\k\\ = objectives) – exact, not a heuristic, but
not intended for huge candidate sets (low hundreds is comfortable;
thousands will be slow). Among the non-dominated set, also computes the
NSGA-II crowding distance (Deb et al. 2002) per objective as a secondary
diagnostic: candidates near the extremes of the front get `Inf`;
candidates in sparsely populated regions of the front get a larger value
than candidates crowded next to near-identical alternatives. This is
purely descriptive (which non-dominated points are most "distinct" from
their neighbours on the front) – it does not change which points are
Pareto-optimal.

## References

Deb, K., Pratap, A., Agarwal, S. & Meyarivan, T. (2002). A fast and
elitist multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on
Evolutionary Computation*, 6, 182-197.

## See also

[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
