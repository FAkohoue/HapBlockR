# Select Parents by Genetic-Algorithm Haplotype Coverage

Selects a fixed-size parent set that jointly covers favourable values
across target haplotype blocks. This is the coverage-only
genetic-algorithm (GA) tool. Use
[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)
when whole-genome merit must also contribute directly to the
optimisation objective, or
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
when parents should be selected only by a whole-genome score.

## Usage

``` r
select_parents_ga(
  value_matrix,
  n_founders,
  strategy = c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS"),
  block_weights = NULL,
  top_candidates = NULL,
  popSize = 100L,
  maxiter = 200L,
  run = 50L,
  pmutation = 0.1,
  pcrossover = 0.8,
  penalty_weight = NULL,
  seed = NULL,
  verbose = FALSE,
  n_reps = 5L,
  G = NULL,
  coancestry_weight = 0,
  target_degree = NULL
)
```

## Arguments

- value_matrix:

  Numeric matrix with candidates in rows and target blocks in columns.
  Row names must contain unique candidate identifiers. Values may be
  local genomic estimated breeding values (local GEBV), directionally
  aligned favourable-haplotype scores or comparable per-block values for
  which larger is better.

- n_founders:

  Positive integer giving the required number of parents.

- strategy:

  One of `"no_selfing"`, `"selfing"`, `"OHS"`, `"OPV"` or
  `"Haploid_OHS"`. See *Crossing-scheme strategies*.

- block_weights:

  Optional non-negative numeric vector of length `ncol(value_matrix)`.
  The default gives every block weight one.

- top_candidates:

  Optional positive integer. Before optimisation, restrict the search to
  this number of candidates ranked by their greatest value across target
  blocks. It must be at least `n_founders`. Leave `NULL` to search all
  candidates.

- popSize, maxiter, run, pmutation, pcrossover:

  Controls passed to
  [`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html). They
  define population size, maximum generations, generations without
  improvement before stopping, mutation probability and crossover
  probability.

- penalty_weight:

  Optional positive coefficient used internally to penalise deviations
  from exactly `n_founders`. The default is calibrated from the block
  weights.

- seed:

  Optional integer seed. With repeated searches, replicate \\i\\ uses
  `seed + i - 1`.

- verbose:

  Logical. If `TRUE`, show the GA monitor for the first replicate.

- n_reps:

  Positive integer number of independent GA searches. Default `5L`. The
  best feasible replicate is selected automatically.

- G:

  Optional named relationship matrix covering all candidates. It is
  required when `coancestry_weight > 0` or `target_degree` is supplied
  and may otherwise be provided to report realised mean relationship.

- coancestry_weight:

  Non-negative numeric weight on mean pairwise relationship. Default
  zero. A positive value produces a weighted coverage-relatedness
  objective. Do not combine it with `target_degree`.

- target_degree:

  Optional number in `[0, 90]`. Zero favours the unconstrained coverage
  end and 90 favours the diversity end. It is converted to a
  relationship ceiling using the supplied `G`. Do not combine it with
  `coancestry_weight`.

## Value

A `hapblockr_result` list. Principal components are `selected`,
`fitness`, `coverage_fitness`, `objective_components`, `per_block`,
`mean_relationship`, `stability`, `feasible`, `converged` and
`result_contract`. `stability$best_rep` identifies the automatically
selected replicate.

## Details

For a candidate set \\S\\, the core objective is
\$\$\mathrm{coverage}(S) = \sum_j w_j\\\mathrm{best}\_j(S),\$\$ where
\\w_j\\ is the weight of target block \\j\\. The definition of
\\\mathrm{best}\_j(S)\\ follows `strategy`. When
`coancestry_weight > 0`, mean pairwise relationship is subtracted from
the objective. When `target_degree` is supplied, the search instead
enforces the corresponding internally calibrated relatedness ceiling.

The function conducts `n_reps` independent searches, checks founder
count and any relatedness ceiling, and automatically returns the
feasible replicate with the greatest complete objective. Breeders do not
choose manually among GA replicates. The `stability` component reports
agreement, convergence and achieved objective values across all
searches.

## Crossing-scheme strategies

- `"no_selfing"`:

  Uses the mean of the two greatest block values among distinct selected
  parents.

- `"OHS"`:

  Optimal Haplotype Selection (OHS); uses the same distinct-parent
  calculation as `"no_selfing"`.

- `"selfing"`:

  Uses the single greatest block value because a selected parent may
  contribute the block without a distinct partner.

- `"OPV"`:

  Optimal Population Value (OPV); uses the same single-best-value
  calculation as `"selfing"`, interpreted as favourable value available
  in the selected population.

- `"Haploid_OHS"`:

  Uses the same single-best-value calculation; one heterozygous parent
  may provide non-homologous gametes.

## See also

[`select_parents_ga_ts`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md),
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)

## Examples

``` r
if (FALSE) { # \dontrun{
prediction <- run_haplotype_prediction(
  geno, snp_info, blocks, blues = phenotype_values
)
target <- select_top_blocks(prediction$block_importance, n = 15)
block_values <- prediction$local_gebv[, target$block_id, drop = FALSE]

selected <- select_parents_ga(
  value_matrix = block_values,
  n_founders = 20,
  strategy = "OHS",
  n_reps = 5,
  seed = 2026
)
selected$selected
selected$stability
} # }
```
