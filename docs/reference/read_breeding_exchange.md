# Read and Verify a Breeding Exchange Bundle

Read and Verify a Breeding Exchange Bundle

## Usage

``` r
read_breeding_exchange(path, verify = TRUE)
```

## Arguments

- path:

  Directory containing `manifest.csv`.

- verify:

  Logical. Verify every file's SHA-256 checksum and row count.

## Value

A `"HapBlockR_exchange_bundle"`.
