# Validate a HapBlockR Result Contract

Checks the common result schema, immutable identifiers, input hashes,
quality-control gates, and decision table. With `strict = TRUE`, a
failed check stops execution; otherwise the function returns the
complete validation report.

## Usage

``` r
validate_hapblockr_result(object, strict = TRUE, ...)
```

## Arguments

- object:

  Object inheriting from `"hapblockr_result"`.

- strict:

  Logical. Stop when any validation check fails. Default `TRUE`.

- ...:

  Reserved for future schema versions.

## Value

A data frame containing `check`, `passed`, and `detail`.
