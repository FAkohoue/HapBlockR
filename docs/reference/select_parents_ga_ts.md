# Select Parents by Joint GA and Whole-Genome Merit

Selects parents through a joint GA plus truncation-selection-merit
objective. The search simultaneously rewards complementary favourable
haplotype-block coverage and the selected set's mean whole-genome merit.
It is an explicit alternative to the coverage-only
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
and the merit-only
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
tools.

## Usage

``` r
select_parents_ga_ts(
  value_matrix,
  n_founders,
  merit_score,
  strategy = c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS"),
  block_weights = NULL,
  top_candidates = NULL,
  merit_priority = 50,
  merit_weight = NULL,
  min_sel_value = NULL,
  min_sel_mode = c("value", "percentile", "sd_above_mean", "relaxed_pool",
    "sd_below_mean"),
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

- merit_score:

  Named finite numeric vector covering every candidate in
  `value_matrix`. Larger values must always represent greater
  whole-genome merit.

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

- merit_priority:

  Number greater than zero and no more than 100. Default 50. It is
  converted internally to a dataset-specific `merit_weight` using the
  reachable coverage and merit spans. Omit it when supplying
  `merit_weight`.

- merit_weight:

  Optional advanced positive raw multiplier for mean whole-genome merit.
  When supplied, `merit_priority` must be omitted.

- min_sel_value:

  Optional merit eligibility floor applied before the GA. Its
  interpretation is set by `min_sel_mode`.

- min_sel_mode:

  One of `"value"`, `"percentile"`, `"sd_above_mean"` or
  `"relaxed_pool"`. `"sd_above_mean"` retains scores at least the stated
  number of standard deviations above the mean. The deprecated
  `"sd_below_mean"` alias is accepted temporarily and maps to
  `"relaxed_pool"`.

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

A `hapblockr_result` list with the components documented for
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
plus `mean_merit`, `merit_weight`, `merit_priority` and `cutoff`. The
result contract records `method = "select_parents_ga_ts"`.

## Details

This is a joint optimisation, not a sequential procedure that first
performs truncation selection and then runs a GA. For a candidate parent
set \\S\\, the complete objective is \$\$\mathrm{fitness}(S) =
\mathrm{coverage}(S) + \lambda\\\overline{\mathrm{merit}}(S) -
\gamma\\\overline{G}(S),\$\$ where the final term is included only when
a positive `coancestry_weight` is requested. Because `n_founders` is
fixed, maximising mean merit gives the same ranking as maximising total
merit.

The breeder must provide `merit_score` with larger values consistently
representing better candidates. The score is normally a HapBlockR
prediction or selection index produced after internal genomic,
haplotype, multi-trait or genotype-by-environment modelling.
Lower-is-better objectives must already have been directionally aligned
before this function is called.

Exactly one merit-control form is active:

- `merit_priority` is the recommended breeder-facing control. It is
  greater than zero and no more than 100. The default, 50, asks the
  package to scale merit to half the empirically reachable coverage span
  for the analysed candidate pool.

- `merit_weight` is an advanced positive raw multiplier. If it is
  supplied, omit `merit_priority`.

A value of zero is deliberately unavailable here because it would
silently turn the hybrid tool into coverage-only GA. Use
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
for that analysis.

The optional `min_sel_value` is a separate hard eligibility rule. It
filters the candidate pool using `merit_score` before joint
optimisation; it does not replace the continuous merit term.

As with
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
all `n_reps` searches are assessed automatically and the greatest
complete objective among feasible replicates is returned. The breeder
reviews the combined stability diagnostics rather than choosing an
individual run manually.

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

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`suggest_merit_weight`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
hybrid <- select_parents_ga_ts(
  value_matrix = prediction$local_gebv[, target_blocks, drop = FALSE],
  n_founders = 20,
  merit_score = prediction$gebv,
  strategy = "OHS",
  merit_priority = 50,
  min_sel_value = 0,
  min_sel_mode = "sd_above_mean",
  n_reps = 5,
  seed = 2026
)
hybrid$selected
hybrid$objective_components
hybrid$stability
} # }
```
