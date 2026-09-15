# Usefulness Criterion (UC) / Genomic Mating: Rank Candidate Crosses

Scores and ranks candidate two-parent crosses by the Usefulness
Criterion (Schnell & Utz 1975; Bernardo 2003; Zhong & Jannink 2007),
i.e. by their predicted potential to produce progeny in the top
`selected_proportion` of the cross – not just by the parents' own
individual merit. This complements
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
and
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md):
those two choose a SET of parents to keep, while
`usefulness_criterion()` instead ranks PAIRS of already- chosen (or
candidate) parents by how good a cross between them is expected to be,
which is the natural next breeding decision once a parent set has been
chosen (which crosses to actually make).

## Usage

``` r
usefulness_criterion(
  parent_ids = NULL,
  cross_pairs = NULL,
  gebv,
  selected_proportion = 0.1,
  n_progeny = NULL,
  n_sim = 20000L,
  seed = NULL,
  variance_model = c("block_independent", "phased", "linked", "simplemating"),
  block_importance = NULL,
  block_ids = NULL,
  local_gebv = NULL,
  segregation_factor = 0.5,
  haplotypes = NULL,
  snp_info = NULL,
  snp_effects = NULL,
  geno_matrix = NULL,
  het_to_na = TRUE,
  G = NULL,
  genetic_map = NULL,
  ld_matrix = NULL,
  type = c("RIL", "DH"),
  generation = 1L,
  n_threads = 1L,
  n_sim_linked = 2000L,
  gebv_se = NULL,
  gebv_reliability = NULL,
  phasing_reliability = NULL,
  downside_quantile = 0.1,
  min_reliability = 0.3,
  verbose = TRUE
)
```

## Arguments

- parent_ids:

  Character vector of candidate parent IDs to generate all pairwise
  crosses from via [`combn`](https://rdrr.io/r/utils/combn.html).
  Ignored if `cross_pairs` is supplied. Either this or `cross_pairs` is
  required.

- cross_pairs:

  Optional data frame or matrix with exactly two columns (parent1,
  parent2) giving a specific, possibly non-exhaustive, list of candidate
  crosses to score. Overrides `parent_ids` if both are given.

- gebv:

  Named numeric vector of whole-genome GEBV (one value per individual,
  names = individual IDs), used for the mid-parent mean term. Required.

- selected_proportion:

  Numeric in (0, 1\]. Proportion of each cross's progeny you intend to
  keep, used to compute the selection intensity \\i\_{sel}\\ (default
  0.1, i.e. keep the top 10%). Smaller values give more weight to
  predicted variance (favour riskier, higher-upside crosses); values
  near 1 make UC converge to the mid-parent value alone.

- n_progeny:

  Optional integer, default `NULL`. If supplied, \\i\_{sel}\\ is
  computed for a FINITE progeny population of this size (by Monte Carlo
  – see "Finite-population selection intensity" below) instead of the
  classical infinite-population asymptotic formula. Use this when your
  realistic cross size is small (tens to a few hundred progeny), where
  the asymptotic formula overstates \\i\_{sel}\\.

- n_sim:

  Integer, default `20000L`. Number of Monte Carlo replicates used for
  the finite-population \\i\_{sel}\\ when `n_progeny` is supplied.
  Ignored otherwise. Larger values reduce Monte Carlo error at the cost
  of runtime.

- seed:

  Optional integer, default `NULL`. Seed for the finite-population Monte
  Carlo (ignored if `n_progeny` is `NULL`). The caller's global RNG
  state is saved and restored, so supplying this does not affect
  random-number generation elsewhere in your script.

- variance_model:

  Character, one of `"block_independent"` (default), `"phased"`,
  `"linked"`, or `"simplemating"`. See "Four variance-prediction modes"
  above. You choose which one to use based on what genotype data you
  have available; this package does not pick one for you.

- block_importance:

  Data frame with (at least) `block_id`, `CHR`, `start_bp`, `end_bp`
  columns identifying which blocks to compute predicted variance over –
  e.g. `run_haplotype_prediction()$block_importance`, optionally
  filtered first via
  [`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  to a manageable, high- importance subset. Required for
  `"block_independent"`/ `"phased"`/`"linked"`; unused (and not
  required) for `"simplemating"`.

- block_ids:

  Optional character vector to further restrict `block_importance` to
  specific `block_id` values before scoring (default `NULL` = use every
  block in `block_importance` as-is). Unused for `"simplemating"`.

- local_gebv:

  Numeric matrix of per-individual, per-block local GEBV (rows =
  individual IDs, columns = `block_id`), e.g.
  `run_haplotype_prediction()$local_gebv`. Required when
  `variance_model = "block_independent"`.

- segregation_factor:

  Numeric scalar, default `0.5`. Scales the single-locus biparental
  segregation-variance formula in `"block_independent"` mode (derived
  from standard F2/RIL additive segregation-variance theory: \\0.5
  \times ((v_i-v_j)/2)^2\\). Exposed so you can adjust it for your
  population's mating design (e.g. a different generation of
  selfing/recombination) rather than being locked to the default
  assumption. Unused outside `"block_independent"`.

- haplotypes:

  Phased haplotype list as returned by
  [`extract_haplotypes`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  (carrying its `block_info` attribute), with `phased = TRUE` blocks.
  Required when `variance_model` is `"phased"` or `"linked"`.

- snp_info:

  Data frame with (at least) `SNP`, `CHR`, `POS` columns, matching what
  was used to build `haplotypes` and fit `snp_effects`. Required when
  `variance_model` is `"phased"` or `"linked"`. Unused for
  `"simplemating"` (use `genetic_map` instead, which carries the same
  SNP identity plus genetic, not physical, position).

- snp_effects:

  Named numeric vector of per-SNP additive marker effects (names = SNP
  IDs), e.g. from
  [`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
  or
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md).
  Required when `variance_model` is `"phased"`, `"linked"`, or
  `"simplemating"`.

- geno_matrix:

  Numeric matrix, individuals x SNPs, dosage-coded 0/1/2/NA. Required
  when `variance_model = "simplemating"`; unused otherwise.
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)
  itself requires strictly 0/2/NA (fully homozygous DH/RIL-style calls)
  – see `het_to_na` for how heterozygous (dosage = 1) calls are handled.

- het_to_na:

  Logical, default `TRUE`. `"simplemating"` only. If `TRUE` (default),
  heterozygous (dosage = 1) calls in `geno_matrix` are automatically
  treated as MISSING (set to `NA`) before calling
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html),
  which requires strictly homozygous 0/2/NA calls – this is a defensible
  statistical choice (equivalent to "no confident homozygous call
  here"), NOT a fabricated one: heterozygous cells are never rounded to
  0 or 2, which would invent a specific allele call the data does not
  support. If `FALSE`, `usefulness_criterion()` errors instead when any
  heterozygous call is found, leaving the decision to you rather than
  converting automatically. Unused otherwise.

- G:

  Dimnamed relationship matrix (row/column names = individual IDs), e.g.
  from
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md).
  Required when `variance_model = "simplemating"` (passed through to
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)'s
  `K` argument – it does not affect the variance calculation, only
  downstream optimisation metadata); unused otherwise.

- genetic_map:

  Data frame with `SNP`, `CHR`, `cM` columns (genetic, not physical,
  position). Preferred (real, validated distances) for `"linked"`;
  required unless `ld_matrix` is supplied instead. For `"simplemating"`,
  likewise required unless `ld_matrix` is supplied instead. Unused for
  `"block_independent"`/`"phased"`.

- ld_matrix:

  Optional square, SNP-dimnamed linkage-disequilibrium (r\\^2\\, in \[0,
  1\]) matrix, e.g. from
  [`compute_r2`](https://FAkohoue.github.io/HapBlockR/reference/compute_r2.md)
  with dimnames set to the SNP IDs used – a recombination-fraction PROXY
  for when no `genetic_map` is available. Used by `"simplemating"` as
  \\1 - LD\\ (passed through to
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)),
  and by `"linked"` as \\0.5 \times (1 - \overline{LD})\\ between
  adjacent target blocks (mean pairwise SNP LD, blocks ordered by
  physical position when no `genetic_map` is given – see "linked" in
  "Four variance-prediction modes" above for the full formula and its
  caveats). **Does not relax the phased-`haplotypes` requirement for
  `"linked"`/`"phased"`** – LD is a population-level statistic and
  cannot substitute for knowing which alleles a specific individual's
  two chromosomes actually carry. Unused for
  `"block_independent"`/`"phased"`.

- type:

  Character, one of `"RIL"` (default) or `"DH"`. `"simplemating"` only:
  the population type derived from each cross (recombinant inbred line
  vs. doubled haploid), which changes the Mendelian-sampling covariance
  formula. Unused otherwise – including by `"linked"`, which models
  single-cross F1-style segregation variance rather than a
  multi-generation RIL/DH population variance (see "Four
  variance-prediction modes" above).

- generation:

  Integer, default `1L`. `"simplemating"` only: the generation at which
  DH lines are generated or RILs extracted (Lehermeier et al. 2017);
  values `>= 10` are treated as an effectively infinite generation.
  Unused otherwise.

- n_threads:

  Integer, default `1L`. `"simplemating"` only: threads used internally
  by
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)'s
  C++ backend. Unused otherwise.

- n_sim_linked:

  Integer, default `2000L`. Number of Monte Carlo progeny simulated PER
  CANDIDATE CROSS for `variance_model = "linked"`. Ignored otherwise.
  Kept separate from `n_sim` (which is for the finite-population
  selection-intensity correction, computed ONCE regardless of cross
  count) because `"linked"` pays this cost once per candidate cross – a
  large value here scales with the number of crosses being scored,
  unlike `n_sim`. `seed`, if supplied, is reused for both Monte Carlo
  procedures.

- gebv_se:

  Optional named numeric vector of GEBV standard errors. When supplied,
  the function propagates parental uncertainty to `mid_parent_SE` and
  the 95 percent UC interval.

- gebv_reliability:

  Optional named numeric vector in \[0, 1\], for example the
  `reliability` column returned by
  [`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md).
  Cross reliability is the conservative minimum of its two parental
  reliabilities.

- phasing_reliability:

  Optional named numeric vector in \[0, 1\]. For
  `variance_model = "phased"`, each cross uses the conservative minimum
  of parental prediction and phasing reliability. Missing phasing
  reliability therefore cannot pass the recommendation gate in phased
  mode.

- downside_quantile:

  Numeric in (0, 0.5). Progeny-distribution quantile reported as
  `downside_value`. Default `0.10`.

- min_reliability:

  Numeric in \[0, 1\]. A cross is marked `recommendation_eligible` only
  when its two-parent reliability meets this threshold. Missing
  reliability never passes the gate. Default `0.30`.

- verbose:

  Logical, default `TRUE`. Print progress and a summary of how many
  candidate crosses could not be scored (e.g. due to missing genotype or
  phase data at target blocks).

## Value

A data frame, one row per candidate cross, sorted by descending UC, with
columns:

- `parent1`, `parent2`:

  The two parent IDs.

- `mid_parent_gebv`:

  Mean of the two parents' whole-genome `gebv`.

- `predicted_variance`:

  Predicted genetic variance of the cross's progeny at the target blocks
  (see Details above for scope). `NA` if the cross could not be scored.

- `selection_intensity`:

  The \\i\_{sel}\\ value used (same for every row, since it depends only
  on `selected_proportion`).

- `UC`:

  `mid_parent_gebv + selection_intensity * sqrt(predicted_variance)`.
  `NA` if the cross could not be scored.

- `rank`:

  Rank by descending UC (`NA` rows sort last).

- `downside_value`:

  The requested lower progeny-distribution quantile, showing downside
  risk alongside expected selected gain.

- `mid_parent_SE`, `UC_SE`, `UC_lower_95`, `UC_upper_95`:

  Propagated uncertainty when `gebv_se` is supplied; the finite-progeny
  contribution is included when `n_progeny` is supplied.

- `prediction_reliability`, `phasing_reliability`, `cross_reliability`,
  `recommendation_eligible`, `eligibility_reason`:

  Explicit reliability gate for promoting a ranked cross to a
  recommendation.

## The Usefulness Criterion

For a cross between parents i and j, UC is defined as \$\$UC =
\mu\_{ij} + i\_{sel} \sqrt{\sigma^2\_{ij}}\$\$ where \\\mu\_{ij}\\ is
the predicted mid-parent value (mean of the progeny distribution),
\\\sigma^2\_{ij}\\ is the predicted genetic variance of the cross's
progeny, and \\i\_{sel}\\ is the standard truncation-selection intensity
for selecting the top `selected_proportion` of that progeny (Falconer &
Mackay 1996). UC therefore favours crosses expected to \*segregate\*
into superior progeny, not just crosses between two already-good parents
– this is what lets it identify transgressive-segregation potential that
whole-genome GEBV alone cannot see.

## Four variance-prediction modes

`variance_model` lets you choose which assumptions to make about the
input data; this package does not force one over the other – pick
whichever matches what you have:

- `"block_independent"`:

  Works on UNPHASED data. Uses each parent's per-block local GEBV
  (`local_gebv`, e.g. from
  [`run_haplotype_prediction`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md))
  and treats each block as a single bi-allelic locus contrasting the two
  parents' block values, via the standard biparental single-locus
  segregation-variance formula \\segregation\\factor \times ((v_i -
  v_j)/2)^2\\, summed across target blocks under an
  independent-assortment assumption. Assumes parents are close to
  homozygous at each block (exact for inbred- line/RIL programs; a
  coarser proxy otherwise). No extra dependency.

- `"phased"`:

  Requires phased haplotypes. Reads each parent's two actual haplotype
  alleles per block via
  [`infer_block_haplotypes`](https://FAkohoue.github.io/HapBlockR/reference/infer_block_haplotypes.md),
  computes each allele's own effect directly from per-SNP effects, and
  exactly enumerates the 4 equally-likely gamete-pair combinations per
  block to get an exact within-block segregation mean/variance. Target
  blocks are still summed independently (no between-block
  linkage/recombination covariance term) – this is a defensible
  block-level extension consistent with the package's "LD block =
  low-recombination unit" premise, but it is NOT a full multi-locus
  LD-aware variance model. No extra dependency.

- `"linked"`:

  Same phased-haplotype input as `"phased"`, plus EITHER a `genetic_map`
  OR an `ld_matrix` – the linkage-aware upgrade of `"phased"` that
  closes its "target blocks summed independently" gap, without the
  external-package requirement of `"simplemating"`. `haplotypes`
  (phased) is required either way – `ld_matrix` only relaxes the
  genetic-distance requirement, it does not substitute for phase
  information (population-level LD tells you nothing about which alleles
  a SPECIFIC individual's two chromosomes carry). Simulates
  `n_sim_linked` progeny by Monte Carlo: for each replicate, one gamete
  per parent is generated by a block-to-block crossover walk, with
  blocks ordered and recombination fractions between adjacent blocks
  obtained either (a) from `genetic_map` via Haldane's (1919) mapping
  function (preferred: a real, validated genetic distance), or (b) when
  no `genetic_map` is supplied, from `ld_matrix` instead: blocks are
  ordered by physical position (CHR, start_bp) and each adjacent pair's
  recombination fraction is PROXIED as \\0.5 \times (1 -
  \overline{LD})\\ (mean pairwise SNP r\\^2\\ between the two blocks'
  member SNPs) – the same \[0, 0.5\] range and boundary behaviour as
  Haldane's mapping function, but a monotonic proxy rather than a
  validated genetic-distance estimator (LD reflects population history,
  not necessarily the true recombination fraction for a specific cross);
  prefer `genetic_map` when you have one. Independent assortment (\\r =
  0.5\\) is used across chromosome boundaries either way. Same scope as
  `"phased"` – single-cross F1-style segregation variance, NOT a
  multi-generation RIL/DH population variance (see `"simplemating"`'s
  `type`/`generation` for that distinct, harder quantity). No extra
  dependency (pure Monte Carlo, no external solver).

- `"simplemating"`:

  The most rigorous option, and the one to prefer when your data fits
  its requirements. Wraps
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)
  (Peixoto et al. 2024, Resende-Lab/SimpleMating – `Suggests`, not
  bundled), which builds a genuine multi-locus Mendelian-sampling
  covariance matrix per chromosome from a genetic map (Haldane-mapped
  recombination fractions) or an LD-matrix proxy when no map is
  available, following Lehermeier et al. (2017) – linkage across ALL
  supplied SNPs, not just within/across independently-treated LD blocks
  like `"phased"`/ `"linked"` above (which only model linkage BETWEEN
  target blocks, not within them). Real constraint:
  [`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)
  requires `geno_matrix` coded strictly 0/2 (fully homozygous calls) –
  built for DH/RIL founder material, not heterozygous outbred parents.
  By default (`het_to_na = TRUE`) heterozygous (dosage = 1) calls are
  automatically treated as missing (set to `NA`) so this mode still runs
  on realistic 0/1/2-coded data rather than erroring outright – see
  `het_to_na` below for the exact behaviour and how to opt back into
  strict error-on-heterozygous validation instead. Does not use
  `block_importance`/`local_gebv`/`haplotypes`/ `n_progeny` at all – see
  the parameter docs below for the arguments this mode actually needs
  (`geno_matrix`, `het_to_na`, `G`, `genetic_map` or `ld_matrix`,
  `type`, `generation`).

## Finite-population selection intensity

The classical UC formula assumes an infinite, normally-distributed
progeny population. Real biparental crosses produce tens to a few
hundred progeny, and the TRUE expected selection intensity for a finite
sample is somewhat smaller than the asymptotic formula gives
(order-statistics theory: the expected mean of the top k of n draws
converges to, but does not equal, the top-p-quantile density ratio as n
grows). Supplying `n_progeny` switches to a Monte Carlo estimate of the
exact finite-sample intensity (see `n_sim`, `seed`) rather than an
approximate closed-form small-sample correction – deliberately, since
several such closed-form corrections exist in the literature and they
are not all mutually consistent, whereas the finite-sample
order-statistic mean itself is estimable directly and exactly (up to
Monte Carlo error) by simulation. Leave `n_progeny = NULL` to use the
standard asymptotic formula (fine for large planned progeny numbers, and
the pre-existing default behaviour).

In both variance modes, `predicted_variance` reflects segregation at the
target blocks in `block_importance` only, not the whole genome – the
same tractability trade-off
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
makes by pre-filtering to a manageable set of top-ranked blocks. If
`block_importance` covers blocks explaining most of local-GEBV variance
(e.g. via
[`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)),
this captures most of the segregating variance that matters in practice,
but is not exactly the whole-genome quantity classical UC formulations
assume. `mid_parent_gebv` (the mean term), by contrast, always uses the
full whole-genome `gebv` you supply, since that is the best available
estimate of the cross's mean.

## What this does and does not do

This ranks pairs of parents for making a cross; it does not choose the
parent set itself (use
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
or
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
for that, optionally first), and it does not allocate differential
numbers of progeny/contributions across crosses or constrain
population-wide inbreeding the way true Optimal Contribution Selection
(OCS) does – that is a separate, not-yet-implemented strategy planned
for this package. Treat `usefulness_criterion()` as a ranking tool for
choosing which of your candidate crosses to prioritise, not as a full
mate-allocation/OCS solver.

## References

Schnell, F.W. & Utz, H.F. (1975). F1-Leistung und Elternwahl in der
Zuchtung von Selbstbefruchtern. *Ber. Arbeitstagung Arbeitsgemeinschaft
Saatzuchtleiter*.

Bernardo, R. (2003). Parental selection, number of breeding populations,
and size of each population in inbred development. *Theoretical and
Applied Genetics*, 107, 1252-1256.

Zhong, S. & Jannink, J.-L. (2007). Using quantitative trait loci results
to discriminate among crosses on the basis of their progeny mean and
variance. *Genetics*, 177, 567-576.

Haldane, J.B.S. (1919). The combination of linkage values and the
calculation of distances between the loci of linked factors. *Journal of
Genetics*, 8, 299-309. (Mapping function used by
`variance_model = "linked"` to convert genetic distance in cM to
recombination fraction.)

Falconer, D.S. & Mackay, T.F.C. (1996). *Introduction to Quantitative
Genetics*, 4th ed. Longman.

Akdemir, D., Beavis, W., Fritsche-Neto, R., Singh, A.K. &
Isidro-Sanchez, J. (2019). Multi-objective optimized genomic breeding
strategies for sustainable food improvement. *Heredity*, 122, 672-683.

Allier, A., Lehermeier, C., Charcosset, A., Moreau, L. & Teyssedre, S.
(2019). Improving short- and long-term genetic gain by accounting for
within-family variance in optimal cross-selection. *Frontiers in
Genetics*, 10, 1006.

## See also

[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
for choosing the parent set itself;
[`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
for building a manageable `block_importance` subset.
