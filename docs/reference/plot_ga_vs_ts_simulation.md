# Plot Realised Genetic Gain: GA vs. Truncation Selection (and Beyond)

Plots mean (solid) and max (dashed) breeding-population GEBV over
generations for every scheme from
[`ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md),
mirroring HapSelect's genetic-gain-over-generations plot. Works for the
original 2-scheme (GA/TS) comparison as well as an arbitrary-length
`schemes` comparison (e.g. adding `"ocs"`/`"uc"`-informed rapid-cycling
schemes).

## Usage

``` r
plot_ga_vs_ts_simulation(sim, show_max = TRUE)
```

## Arguments

- sim:

  List returned by
  [`ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md).

- show_max:

  Logical. Overlay the max-GEBV trajectory (dashed) in addition to the
  mean (solid). Default `TRUE`.

## Value

A `ggplot2` object.

## See also

[`ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_ga_vs_ts_simulation(sim)
} # }
```
