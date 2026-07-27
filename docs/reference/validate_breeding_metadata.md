# Validate Traceable Breeding Metadata

Validates canonical germplasm, trait, environment, and observation
tables before modelling or exchange. Referential integrity, ontology
identifiers, units, programme, study, trial, environment, and germplasm
identities are checked without calling a remote service.

## Usage

``` r
validate_breeding_metadata(
  germplasm,
  traits,
  environments,
  observations,
  strict = TRUE
)
```

## Arguments

- germplasm:

  Data frame containing `germplasm_id` and `breeding_program_id`.

- traits:

  Data frame containing `trait_id`, `trait_name`, `unit`, `direction`,
  and `ontology_id`.

- environments:

  Data frame containing `environment_id`, `study_id`, `trial_id`, and
  `location_id`.

- observations:

  Data frame containing `observation_id`, `germplasm_id`, `trait_id`,
  `environment_id`, `value`, and `unit`.

- strict:

  Logical. Stop if any issue is found. When `FALSE`, return a failed
  result contract containing the full issue ledger.

## Value

A `hapblockr_result` with validated canonical tables, standards mapping,
and an issue ledger.
