# Assess Recommendation Stability Across Scenarios

Compares validated result objects from threshold sweeps, leave-one-
environment-out analyses, leave-one-population-out analyses, or
alternative model assumptions. Scenarios are supplied explicitly so that
the function never refits or silently changes a model.

## Usage

``` r
assess_decision_stability(
  results,
  top_n = 10L,
  score_col = NULL,
  selected_col = NULL,
  minimum_selection_frequency = 0.5,
  minimum_baseline_jaccard = 0.5,
  require_valid = TRUE
)
```

## Arguments

- results:

  Named list of at least two `hapblockr_result` objects.

- top_n:

  Number of top-ranked rows treated as selected when no explicit
  selection column is present. Default `10`.

- score_col:

  Optional common numeric score column.

- selected_col:

  Optional common logical selection column.

- minimum_selection_frequency:

  Minimum scenario frequency required for a stable recommendation.

- minimum_baseline_jaccard:

  Minimum Jaccard overlap required between every scenario and the first
  named baseline scenario.

- require_valid:

  Logical. Refuse failed input result contracts.

## Value

A `hapblockr_result` with item frequencies, pairwise Jaccard overlap,
rank and score summaries, and scenario validation status.
