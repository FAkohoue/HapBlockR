# Write a Breeding Exchange Bundle

Writes every exchange table as UTF-8 CSV and creates a SHA-256 manifest.

## Usage

``` r
write_breeding_exchange(bundle, path, overwrite = FALSE)
```

## Arguments

- bundle:

  Object from
  [`build_breeding_exchange`](https://FAkohoue.github.io/HapBlockR/reference/build_breeding_exchange.md).

- path:

  Target directory.

- overwrite:

  Logical. Replace known bundle files when they exist.

## Value

Invisibly returns the manifest data frame.
