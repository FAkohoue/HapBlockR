# True Optimal Contribution Selection (OCS) and Mate Allocation

Solves for optimal parent contributions and a resulting crossing plan
that maximizes genetic merit subject to an explicit constraint on the
next generation's relatedness (Meuwissen 1997) – the actual formal OCS
problem, as opposed to
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
fixed-size subset search (with only a soft relatedness penalty) or
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)'s
independent per-cross ranking.

## Usage

``` r
select_parents_ocs(
  merit,
  G,
  family = NULL,
  engine = c("auto", "alphamate", "optisel", "simplemating"),
  n_crosses = 20L,
  n_parents_max = NULL,
  max_contrib_per_parent = 4L,
  allow_selfing = FALSE,
  allow_repeated_matings = FALSE,
  target_degree = 30,
  rescale_nrm = TRUE,
  alphamate_exe = NULL,
  out_dir = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- merit:

  Named numeric vector of candidate parent merit (e.g. GEBV, a Selection
  Index, or
  [`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md)'s
  `stacking_index`), names = individual IDs.

- G:

  Dimnamed relationship matrix (row/column names = individual IDs), e.g.
  from
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md),
  or your own VanRaden/IBS/blended-H/pedigree-A matrix built however you
  already build it. This function does not construct `G` for you – pass
  in whatever relationship matrix you trust.

- family:

  Optional named character vector (names = individual IDs) of
  family/cross-of-origin labels, used only for output labelling/
  reporting (neither engine enforces family representation as a hard
  constraint; see Details).

- engine:

  Character, one of `"auto"` (default), `"alphamate"` (true OCS),
  `"optisel"` (true OCS), or `"simplemating"` (discrete greedy cross
  prediction/selection, not true OCS). See Description – and note the
  breaking rename if you have existing code using `engine = "optisel"`
  from before this version.

- n_crosses:

  Integer, default `20L`. Number of matings to produce.

- n_parents_max:

  Optional integer. Maximum number of distinct parents to use.
  AlphaMate: `NumberOfParents`, a native constraint enforced during
  allocation. `engine = "optisel"`: approximated by zeroing all but the
  `n_parents_max` largest solved contributions after `opticont()` solves
  – not a re-optimised constrained solution (with a message). Ignored
  entirely (with a message) under `engine = "simplemating"` –
  `selectCrosses()` does not expose this control at all.

- max_contrib_per_parent:

  Integer, default `4L`. Maximum number of matings any single parent can
  participate in. Honoured directly by all three engines (AlphaMate's
  own cap; `selectCrosses()`'s `max.cross` argument under
  `engine = "simplemating"`;
  [`optiSel::matings()`](https://rdrr.io/pkg/optiSel/man/matings.html)'s
  `ub.n` argument, applied per-pair, plus a pre-solve cap on each
  candidate's assigned offspring slots, under `engine = "optisel"`).

- allow_selfing:

  Logical, default `FALSE`.

- allow_repeated_matings:

  Logical, default `FALSE`. Enforced by construction under all three
  engines for the default (non-selfing) design – see "Engine
  differences".

- target_degree:

  Numeric, default `30`, in `[0, 90]`. Diversity-vs-gain lever for ALL
  THREE engines, but via genuinely different mechanisms per engine
  (AlphaMate's continuous Kinghorn- frontier `TargetDegree`; under
  `engine = "optisel"`, a mean-kinship ceiling interpolated between two
  `opticont()`-solved frontier extremes for your candidate set; under
  `engine = "simplemating"`, `selectCrosses()`'s hard
  `culling.pairwise.k` relatedness cutoff, set from `target_degree` by
  quantile of the candidate set's own observed relatedness values, with
  a floor that keeps enough candidates for `n_crosses` to stay feasible)
  – **0 prioritises maximising the criterion** (merit/gain, accepting
  more relatedness – the max-gain end of Kinghorn's frontier), **90
  prioritises minimising relatedness** (the max-diversity end). This
  matches AlphaMate's own `TargetDegree` convention under all three
  engines. Not numerically identical between engines; re-tune if you
  switch. See "Engine differences".

- rescale_nrm:

  Logical, default `TRUE`. Rescale `G` to mean diagonal 1
  (numerator-relationship-like scale) before use, matching AlphaMate's
  expected NRM convention. Recommended to leave on for all three engines
  unless you have already rescaled `G` yourself.

- alphamate_exe:

  Path to the AlphaMate executable. **You must obtain and be entitled to
  use this executable yourself** – it is a separate, third-party tool
  (Hickey Group / AlphaGenes suite), not an R package and not part of
  HapBlockR's own code. Default `NULL`: no binary is bundled with or
  auto-detected by HapBlockR (see "Providing the AlphaMate executable"
  above), so `NULL` always means AlphaMate is unavailable –
  `engine = "auto"` will fall back to `"optisel"`. Required (explicitly
  supplied) for `engine = "alphamate"` (or for `"auto"` to select it).

- out_dir:

  Directory to write AlphaMate's input/output files to. Required for
  `engine = "alphamate"`. Created if it doesn't exist.

- seed:

  Optional integer seed (used by the native mate-allocation step in each
  engine's final crossing-pair generation, and by optiSel's solver if it
  uses randomness internally).

- verbose:

  Logical, default `TRUE`.

## Value

A list with components `mating_plan` (data frame: `parent1`, `parent2`,
`mean_relationship`, plus engine-reported columns where available),
`contributors` (data frame: `id`, `contribution`, `family` if supplied),
`engine_used`, and `ok` (logical: whether no-selfing/
no-repeated-mating/contribution-cap constraints were all satisfied in
the final plan).

## Details

This package does not implement an OCS solver itself. Three mature,
external engines are wrapped, selected via `engine`. Two solve the
actual OCS problem; one instead does discrete, greedy cross prediction/
selection – see "Engine differences" below for why that distinction
matters, not just which is "better":

- `"alphamate"`:

  **True OCS.** The AlphaMate executable (Hickey group / AlphaGenes
  suite). Full evolutionary-algorithm mate allocation: specific crossing
  pairs, a `target_degree` diversity control,
  no-selfing/no-repeated-mating constraints, per-parent contribution
  caps. Requires the AlphaMate binary installed separately – point
  `alphamate_exe` at it. Not an R package; not distributed with
  HapBlockR.

- `"optisel"`:

  **True OCS.**
  [`optiSel::candes()`](https://rdrr.io/pkg/optiSel/man/candes.html) +
  `opticont()` (Wellmann 2019) – optiSel's own continuous- optimization
  solver for Meuwissen's (1997) formal OCS problem: optimal
  per-candidate contributions subject to an explicit upper bound on the
  next generation's mean kinship. `target_degree` is mapped onto that
  kinship bound by solving both ends of the gain/ diversity frontier for
  your actual candidate set (a genuine min-kinship optimum and a genuine
  max-merit optimum, both via `opticont()`) and interpolating between
  them.
  [`optiSel::noffspring()`](https://rdrr.io/pkg/optiSel/man/noffspring.html)
  then converts contributions into expected offspring counts, and
  [`optiSel::matings()`](https://rdrr.io/pkg/optiSel/man/matings.html)
  solves the discrete Sire x Dam assignment minimising mean offspring
  kinship. Requires the optiSel package – no external binary, no
  SimpleMating dependency.

- `"simplemating"`:

  **Not true OCS – discrete greedy cross prediction/selection.**
  [`SimpleMating::planCross()`](https://rdrr.io/pkg/SimpleMating/man/planCross.html) +
  [`SimpleMating::selectCrosses()`](https://rdrr.io/pkg/SimpleMating/man/selectCrosses.html)
  (Peixoto et al. 2024, Resende-Lab/SimpleMating). `planCross()`
  enumerates every candidate pair; this function computes each pair's
  mid-parent merit and `G`-matrix relatedness directly (plain
  arithmetic); `selectCrosses()` then does the actual
  relatedness-constrained mate allocation – a greedy search maximising
  merit subject to a relatedness cutoff and per-parent min/max
  cross-count caps. Requires SimpleMating + optiSel installed (no
  external binary; optiSel is one of SimpleMating's own dependencies).
  See "Engine differences" below for real, verified differences from the
  two true-OCS engines' diversity-control mechanism.

**This is a breaking rename as of the version introducing
`"simplemating"`**: `engine = "optisel"` previously meant what is now
`engine = "simplemating"` (optiSel's own solver was never actually
called under the old name – optiSel was present only as one of
SimpleMating's transitive dependencies). Code written against the old
behaviour must change `engine = "optisel"` to `engine = "simplemating"`;
there is no backward-compatible alias, since silently keeping the string
`"optisel"` pointed at SimpleMating's algorithm would leave the name
permanently misleading. See `NEWS.md`.

`engine = "auto"` (default) uses AlphaMate only when you supplied a
working `alphamate_exe` yourself and the file exists (your own native
build, a Wine wrapper script, etc. is trusted as-is); otherwise it falls
back to `"optisel"` (the other true-OCS engine – `"auto"` never silently
substitutes cross prediction for OCS) with a message. No AlphaMate
binary is bundled with HapBlockR or auto-detected – see "Providing the
AlphaMate executable" below.

## Engine differences – read before choosing

`"alphamate"` and `"optisel"` both solve the actual OCS problem
(continuous contribution optimization under a relatedness constraint)
via genuinely different solvers – AlphaMate's own evolutionary algorithm
vs. optiSel's `candes()`/`opticont()`/`matings()`. Both are controlled
by the same `target_degree` lever using the same direction convention (0
= max-gain end of the frontier, prioritising merit and accepting more
relatedness; 90 = max-diversity end, prioritising minimised
relatedness), but via genuinely different mechanisms: AlphaMate uses a
continuous Kinghorn-frontier degree (its own `TargetDegree`) natively
inside its evolutionary algorithm; the `"optisel"` engine interpolates a
mean-kinship ceiling between two `opticont()`-solved frontier extremes
specific to your candidate set (see `.run_optisel_ocs()`'s own
comments). Not numerically identical between the two; re-tune if you
switch.

`"simplemating"` is a different algorithm class entirely – discrete
greedy cross selection under a hard relatedness cutoff
(`culling.pairwise.k`: candidate pairs more related than this are
discarded outright, not smoothly down-weighted), not a continuous
contribution optimum. `target_degree` is mapped onto that cutoff by
*quantile* of the candidate set's own observed relatedness values, in
the same direction as the two true-OCS engines' `target_degree` (0 =
keeps essentially all candidates, the max-gain end of the frontier with
no diversity restriction; 90 = keeps only the least-related candidates,
the max-diversity end), with a floor that always keeps enough candidates
for the requested `n_crosses` to be feasible – a raw linear
interpolation across `[min(K), max(K)]` was tried first and found, via a
real
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
run, to be fragile against outlier pairs that widen the observed range
and starve `selectCrosses()`'s search of candidates at moderate
`target_degree` values. Treat every engine's `target_degree` as an
approximation of the same lever, not an equivalent algorithm; re-tune if
you switch engines.

Beyond that: AlphaMate lets you cap the parent count (`n_parents_max`)
as a native constraint during allocation; the `"optisel"` engine
approximates it by zeroing all but the `n_parents_max` largest solved
contributions *after* solving (not a re-optimised constrained solution –
see `.run_optisel_ocs()`); the `"simplemating"` engine does not support
it at all – ignored (with a message). `max_contrib_per_parent` IS
honoured directly by all three engines (AlphaMate's own cap;
`selectCrosses()`'s `max.cross` argument under `"simplemating"`;
[`optiSel::matings()`](https://rdrr.io/pkg/optiSel/man/matings.html)'s
`ub.n` argument under `"optisel"`). All three engines also guarantee
no-repeated-matings by construction for this function's default
(non-selfing) design – AlphaMate enforces it during allocation;
[`optiSel::matings()`](https://rdrr.io/pkg/optiSel/man/matings.html) is
called with `ub.n = 1` unless `allow_repeated_matings = TRUE`; and
`planCross(MateDesign = "half")` lists each unordered pair at most once
so `selectCrosses()`'s greedy search cannot select the same pair twice.

## Providing the AlphaMate executable

AlphaMate is a separate, third-party tool (Hickey Group / AlphaGenes
suite, <https://github.com/AlphaGenes/AlphaMate>, Fortran source, MIT
licensed) – it is not distributed as an R package, HapBlockR does not
implement its algorithm, and **no AlphaMate binary is bundled with or
downloaded by this package**: AlphaGenes does not currently publish
pre-built binaries (no GitHub Releases), so there is no stable URL this
package could fetch one from, and CRAN policy prohibits shipping
compiled executables in a source package regardless. To use
`engine = "alphamate"` you must build AlphaMate yourself from that
repository (a Fortran compiler such as `gfortran` is required) or
otherwise obtain a working executable you are entitled to use under its
own license terms, then pass its full path explicitly via
`alphamate_exe = "/full/path/to/AlphaMate"` (or `.exe` on Windows) –
there is no auto-detection. If you would rather avoid building a
separate binary, `engine = "optisel"` (true OCS, no SimpleMating
dependency) or `engine = "simplemating"` (cross prediction/selection)
both require no external executable at all – only R packages – and are
the more portable, drop-in choices.

## Verification status (read this)

The AlphaMate wrapper is a direct translation of a working, production
AlphaMate driver script (same file formats, same call sequence) and is
implemented with high confidence.

The `"simplemating"` engine went through two revisions after real
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
failures against an actually-installed SimpleMating: a hand-written
[`optiSel::candes()`](https://rdrr.io/pkg/optiSel/man/candes.html)/`opticont()`
call with real argument-name mismatches, replaced by a wrapper around
`SimpleMating::GOCS()` after reading its actual source – which itself
then failed because a real, freshly-reinstalled SimpleMating 0.2.1 does
not export `GOCS()` despite it being present in the GitHub repository's
checked-in source/NAMESPACE (unresolved discrepancy; possibly a stale
NAMESPACE line). This engine now wraps `planCross()` + `selectCrosses()`
instead, specifically because both are confirmed present in a real,
current SimpleMating 0.2.1 installation's own help index AND their exact
current source was read directly before writing this
(`Resende-Lab/SimpleMating`, `planCross.R` and `selectCrosses.R`). If
`engine = "simplemating"` errors with what looks like an argument-name
mismatch, SimpleMating's own API may have moved again since this wrapper
was written – check
[`?SimpleMating::planCross`](https://rdrr.io/pkg/SimpleMating/man/planCross.html)/
[`?SimpleMating::selectCrosses`](https://rdrr.io/pkg/SimpleMating/man/selectCrosses.html)
against your installed version.

The `"optisel"` engine is written directly against optiSel's own
`candes()`/`opticont()`/`matings()`/`noffspring()` documentation. The
specific assumptions most likely to matter if you hit an argument-name
or return-shape mismatch: (1)
[`optiSel::noffspring()`](https://rdrr.io/pkg/optiSel/man/noffspring.html)'s
exact first positional argument and return shape (assumed to accept the
`opticont()` `$parent` data frame directly and return a list with a
`$parent$n` column); (2)
[`optiSel::matings()`](https://rdrr.io/pkg/optiSel/man/matings.html)'s
`ub.n` argument being an upper bound on repeated matings of the same
Sire x Dam pair (used here to enforce `allow_repeated_matings = FALSE`
natively rather than by post-hoc filtering). Check
[`?optiSel::candes`](https://rdrr.io/pkg/optiSel/man/candes.html),
[`?optiSel::opticont`](https://rdrr.io/pkg/optiSel/man/opticont.html),
[`?optiSel::matings`](https://rdrr.io/pkg/optiSel/man/matings.html), and
[`?optiSel::noffspring`](https://rdrr.io/pkg/optiSel/man/noffspring.html)
against your installed version.

## References

Meuwissen, T.H.E. (1997). Maximizing the response of selection with a
predefined rate of inbreeding. *Journal of Animal Science*, 75, 934-940.

Wellmann, R. (2019). Optimum contribution selection and mate allocation
for genetic improvement, inbreeding, and diversity: the R package
optiSel. *BMC Bioinformatics*, 20, 25.

Peixoto, M.A., Amadeu, R.R., Bhering, L.L., Ferrao, L.F.V., Munoz, P.R.
& Resende Jr., M.F.R. (2024). SimpleMating: R-package for prediction and
optimization of breeding crosses using genomic selection. *The Plant
Genome*, e20533.

## See also

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
