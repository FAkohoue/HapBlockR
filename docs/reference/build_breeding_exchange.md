# Build a BrAPI/MIAPPE-Oriented Recommendation Bundle

Converts a validated HapBlockR decision object into stable tabular
recommendations, quality-control, uncertainty, exclusion, provenance,
and field-mapping tables. It does not embed remote-service calls.

## Usage

``` r
build_breeding_exchange(
  result,
  breeding_program_id,
  study_id,
  trial_id,
  environment_id = NA_character_,
  trait_id = NA_character_,
  unit = NA_character_,
  ontology_id = NA_character_,
  require_valid = TRUE
)
```

## Arguments

- result:

  A `hapblockr_result`.

- breeding_program_id, study_id, trial_id:

  Stable programme, study, and trial identifiers.

- environment_id:

  Optional environment identifier.

- trait_id:

  Optional trait or observation-variable identifier.

- unit:

  Optional trait unit.

- ontology_id:

  Optional ontology CURIE.

- require_valid:

  Logical. Refuse a failed decision contract.

## Value

A list of canonical exchange tables with class
`"HapBlockR_exchange_bundle"`.
