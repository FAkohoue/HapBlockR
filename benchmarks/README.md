# HapBlockR benchmark protocol

`run_ld_benchmark.R` is a deterministic medium-scale correctness and
performance gate. It records the input shape and hash, repeated wall times,
one- and two-thread equality, agreement with a scalar R reference, compiler,
R and package versions, operating system, processor description, and Linux
peak resident memory.

Run from the package root:

```r
Rscript benchmarks/run_ld_benchmark.R
```

The CI workflow retains the resulting CSV artefact and enforces the limits in
`thresholds.csv`. These limits are release gates for gross regressions, not
claims about whole-genome performance. Publishable benchmarks require
representative public data, dedicated hardware, more repetitions, and
independent peak-memory measurement.
