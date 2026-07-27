# Certify an Operational Mating Plan

Validates a complete plan against cross-level rules, family size, parent
capacities, required pairs, period-specific capacity, and subpopulation
contribution quotas. It returns a feasibility certificate and all
binding constraints.

## Usage

``` r
certify_mating_plan(
  plan,
  candidate_status,
  female_col = "female",
  male_col = "male",
  family_size_col = "family_size",
  period_col = "period",
  min_family_size = 1L,
  required_pairs = NULL,
  capacity_by_period = NULL,
  subpopulation_col = "subpopulation",
  subpopulation_quotas = NULL,
  no_repeated_cross = TRUE,
  ...
)
```

## Arguments

- plan:

  Data frame containing female, male, family size, and optionally
  period.

- candidate_status:

  Candidate status table passed to
  [`screen_candidate_crosses`](https://FAkohoue.github.io/HapBlockR/reference/screen_candidate_crosses.md).
  Optional numeric columns `female_capacity`, `male_capacity`,
  `total_capacity`, `seed_available`, and `pollen_available` are
  enforced.

- female_col, male_col, family_size_col, period_col:

  Plan column names.

- min_family_size:

  Minimum permitted family size.

- required_pairs:

  Optional two-column table of required unordered pairs.

- capacity_by_period:

  Optional table containing `id`, `period`, and one or more of
  `female_capacity`, `male_capacity`, and `total_capacity`.

- subpopulation_col:

  Candidate subpopulation column.

- subpopulation_quotas:

  Optional table with `subpopulation` and optional minimum or maximum
  contribution columns.

- no_repeated_cross:

  Logical.

- ...:

  Additional arguments passed to
  [`screen_candidate_crosses`](https://FAkohoue.github.io/HapBlockR/reference/screen_candidate_crosses.md).

## Value

A `hapblockr_result` with a row-level certificate, violations, capacity
use, quota use, and binding constraints.
