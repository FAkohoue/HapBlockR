# Prepare Externally Analysed Breeding Targets

Validates genotype-level estimates produced outside HapBlockR, computes
precision weights, deregresses genuine random-effect predictions where
required, and orients every trait so that larger model values are
favourable. The returned \`model_value\` is the response supplied to
HapBlockR's marker, haplotype, block, multi-trait, or
genotype-by-environment models. External genomic BLUPs and external
selection-index values are rejected because genomic effects and
selection indices are fitted within HapBlockR.

## Usage

``` r
prepare_breeding_targets(
  data,
  input_type,
  id_col = "id",
  trait_col = "trait",
  value_col = "value",
  environment_col = NULL,
  unit_col = NULL,
  se_col = NULL,
  posterior_sd_col = NULL,
  precision_col = NULL,
  covariance = NULL,
  reliability_col = NULL,
  pev_col = NULL,
  genetic_variance = NULL,
  relationship_matrix = NULL,
  parent_average = NULL,
  estimation_basis = NULL,
  gca_effect = NULL,
  tester_population = NULL,
  additive_col = NULL,
  dominance_col = NULL,
  lower_is_better = NULL,
  heritability = NULL,
  strict = TRUE
)
```

## Arguments

- data:

  One row per genotype, trait, and optional environment.

- input_type:

  One of the values returned by
  [`breeding_target_types()`](https://FAkohoue.github.io/HapBlockR/reference/breeding_target_types.md).
  A single type must apply to the call.

- id_col, trait_col, value_col:

  Column names.

- environment_col:

  Optional environment column. Leave \`NULL\` for an across-environment
  estimate.

- unit_col:

  Optional declared-unit column. Units are retained when supplied but
  are never mandatory.

- se_col:

  Standard error column for \`adjusted_mean\` or \`BLUE\`.

- posterior_sd_col:

  Posterior standard deviation column for a Bayesian \`adjusted_mean\`.

- precision_col:

  Externally supplied positive precision column.

- covariance:

  Optional full sampling covariance matrix. Its row and column names
  must equal the target record keys returned in \`record_key\`.

- reliability_col, pev_col:

  Alternative uncertainty inputs for random predictions. \`pev_col\`
  must contain prediction error variances, not standard errors of
  differences.

- genetic_variance:

  Additive or genotype variance, either one positive number, a named
  value per trait, or a named trait-by-environment value.

- relationship_matrix:

  Required numerator relationship matrix for \`PBLUP\`, and for
  pedigree-based \`BV\`, \`GCA\`, or \`TGV\`.

- parent_average:

  Optional named parent-average vector used in deregression. The default
  is zero.

- estimation_basis:

  Required for \`BV\`; and conditionally for \`GCA\` and \`TGV\`.
  Allowed values are \`"fixed"\`, \`"identity"\`, and \`"pedigree"\`.
  \`"genomic"\` is explicitly rejected.

- gca_effect:

  For \`GCA\`, either \`"fixed"\` or \`"random"\`.

- tester_population:

  For \`GCA\`, a non-empty description of the tester or mate population
  to which the estimates apply.

- additive_col, dominance_col:

  For \`TGV\`, columns containing separate additive and dominance
  components. Their sum is used as the supplied total genetic value.

- lower_is_better:

  Trait names for which smaller original values are favourable. Their
  signs are reversed in \`model_value\`.

- heritability:

  Optional provenance only. It is recorded but is not used to
  manufacture reliability or precision.

- strict:

  Logical. Stop if any reliability lies outside \`(0, 1\]\`.

## Value

A \`hapblockr_result\` containing \`targets\`, a \`target_contract\`,
precision provenance, optional covariance, and quality gates.
