# Screen Crosses for Operational Feasibility

Screens directed crosses, where parent 1 is the female and parent 2 is
the male, against role, fertility, flowering, reciprocal,
heterotic-group, forbidden-pair, and quarantine rules. The function
reports every failed rule rather than stopping at the first failure.

## Usage

``` r
screen_candidate_crosses(
  cross_pairs,
  candidate_status,
  id_col = "id",
  female_allowed_col = "female_allowed",
  male_allowed_col = "male_allowed",
  fertility_col = "fertile",
  flowering_start_col = "flowering_start",
  flowering_end_col = "flowering_end",
  heterotic_group_col = "heterotic_group",
  require_different_heterotic_groups = FALSE,
  quarantine_group_col = "quarantine_group",
  quarantine_compatibility = NULL,
  forbidden_pairs = NULL,
  reciprocal_effects = NULL,
  no_selfing = TRUE
)
```

## Arguments

- cross_pairs:

  Two-column matrix or data frame of directed crosses.

- candidate_status:

  Candidate table with one row per parent.

- id_col:

  Candidate identifier column.

- female_allowed_col, male_allowed_col:

  Optional logical role columns.

- fertility_col:

  Optional logical fertility column.

- flowering_start_col, flowering_end_col:

  Optional numeric or Date flowering-window columns. Both must be
  supplied together.

- heterotic_group_col:

  Optional heterotic-group column.

- require_different_heterotic_groups:

  Logical.

- quarantine_group_col:

  Optional quarantine-group column.

- quarantine_compatibility:

  Optional named logical matrix whose rows are female quarantine groups
  and columns are male quarantine groups.

- forbidden_pairs:

  Optional two-column table. Pairs are treated as unordered unless it
  contains a logical `directional` column.

- reciprocal_effects:

  Optional table with `female`, `male`, and logical `allowed` columns.

- no_selfing:

  Logical.

## Value

A data frame with female, male, feasibility, and a semicolon- separated
reason ledger.
