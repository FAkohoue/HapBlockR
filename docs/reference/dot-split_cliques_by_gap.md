# Split cliques at physical gaps, keeping indices and bp positions in lockstep

Splits each clique wherever the gap between consecutive (sorted) member
bp positions exceeds `gapdist`. Takes and returns the SNP column
*indices* alongside their bp positions, sorted/split identically, so
that two SNPs sharing the same bp position (co-located multi-allelic
sites, or an indel and a SNP at the same POS) are never conflated – each
keeps its own distinct column index through the split, unlike an earlier
version that tracked bp values alone and used
[`match()`](https://rdrr.io/r/base/match.html)/[`setdiff()`](https://rdrr.io/r/base/sets.html)
against those values (ambiguous whenever a bp value isn't unique).

## Usage

``` r
.split_cliques_by_gap(cliques_idx, cliques_bp, gapdist)
```

## Arguments

- cliques_idx:

  List of integer vectors: original column indices per clique (into
  `binvector`/`SNPbps` in
  [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md)).

- cliques_bp:

  List of numeric vectors: bp positions, same length and order as the
  corresponding element of `cliques_idx`.

- gapdist:

  Numeric. Gap threshold (bp) above which a clique is split.

## Value

Named list with `idx` and `bp`, each a list of vectors, parallel to each
other, after splitting.
