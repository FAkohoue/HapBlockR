# Genetic-Algorithm Founder-Parent Selection

Searches for the set of `n_founders` individuals that jointly maximises
coverage of favourable per-block values (local GEBV, haplotype allele
dosage, or any comparable score) across a set of target blocks, via
`GA::ga(type = "binary")`. This is HapBlockR's answer to HapSelect's
`local_gebv_parent_selection()` / `haplotype_parent_selection()` – the
single biggest capability gap identified against HapSelect (see
`vignette` / gap-analysis notes): turning "here is what's good" (block
importance, stacking scores) into "here is the actual founder set",
instead of leaving a breeder to eyeball a ranked list.

## Usage

``` r
select_parents_ga(
  value_matrix,
  n_founders,
  strategy = c("no_selfing", "selfing", "OHS", "OPV", "Haploid_OHS"),
  block_weights = NULL,
  top_candidates = NULL,
  popSize = 100L,
  maxiter = 200L,
  run = 50L,
  pmutation = 0.1,
  pcrossover = 0.8,
  penalty_weight = NULL,
  seed = NULL,
  verbose = FALSE,
  merit_score = NULL,
  min_sel_value = NULL,
  min_sel_mode = c("value", "percentile", "sd_below_mean"),
  n_reps = 5L,
  G = NULL,
  coancestry_weight = 0,
  merit_weight = 0,
  merit_priority = NULL,
  target_degree = NULL
)
```

## Arguments

- value_matrix:

  Numeric matrix (individuals x blocks): local GEBV
  (`run_haplotype_prediction()$local_gebv`), haplotype allele dosage
  scaled to \[0,1\] (dosage / ploidy), or any comparable per-block
  favourable-value score. Row names = candidate individual IDs.
  **Pre-filter to a manageable set of top-ranked blocks first** (e.g.
  via
  [`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md))
  – exactly as HapSelect's own documented workflow does before calling
  its parent-selection GA – since the fitness function is evaluated
  `popSize * maxiter` times.

- n_founders:

  Integer. Target number of founders to select.

- strategy:

  One of `"no_selfing"` (default), `"selfing"`, `"OHS"`, `"OPV"`,
  `"Haploid_OHS"`. See *Crossing-scheme strategies*.

- block_weights:

  Numeric vector, length `ncol(value_matrix)`, or `NULL` (default: equal
  weight 1 for every block). E.g. pass the `var_scaled` column from
  [`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  to weight higher-variance blocks more heavily.

- top_candidates:

  Integer or `NULL` (default). If supplied, restricts the GA's candidate
  pool to the top `top_candidates` individuals by their own row-max
  value in `value_matrix` before searching – a heuristic prefilter to
  keep the binary chromosome length (and therefore search time)
  manageable for large populations. Set `NULL` to search the full
  population in `value_matrix`.

- popSize, maxiter, run, pmutation, pcrossover:

  Passed to
  [`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html).
  Defaults `100`, `200`, `50`, `0.1`, `0.8` respectively.

- penalty_weight:

  Numeric. Quadratic penalty coefficient for deviating from exactly
  `n_founders` selected individuals. Default `10` times the largest
  single block weight, scaled so the penalty dominates the fitness for
  any deviation while not distorting the optimum for correctly-sized
  subsets. Increase if the GA returns solutions with the wrong count.

- seed:

  Integer or `NULL`. Random seed for reproducibility. When `n_reps > 1`,
  this seeds replicate 1 and replicates `2..n_reps` use `seed + 1`,
  `seed + 2`, ... (deterministic and reproducible as a set, but each
  replicate explores a different starting point).

- verbose:

  Logical. Print
  [`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html)'s
  iteration monitor for the first replicate only (silent for replicates
  `2..n_reps` to avoid flooding the console). Default `FALSE`.

- merit_score:

  Named numeric vector (e.g. `run_haplotype_prediction()$gebv`), or
  `NULL` (default). Whole-genome merit, used for up to two independent
  purposes depending on which of `min_sel_value`/`merit_weight` are
  set: (1) the `min_sel_value` eligibility floor below, a hard
  pre-search exclusion; and (2) when `merit_weight > 0`, a soft,
  additive term *inside* the GA fitness function itself – see
  *Merit-weighted fitness (GA+TS hybrid, optional)*. The two are
  independent and commonly used together (a floor to exclude clearly
  ineligible candidates, plus a weight so merit keeps pulling on the
  search among those who remain), but either can be used alone. Required
  whenever `min_sel_value` is not `NULL`, or `merit_weight > 0`.

- min_sel_value, min_sel_mode:

  Merit floor applied to `value_matrix`'s candidate pool *before* the GA
  searches, using `merit_score` to evaluate it. Same semantics as
  [`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)'s
  arguments of the same name; default `min_sel_value = NULL` applies no
  floor at all (the historical behaviour of this function). Existing
  from an earlier finding: without a floor, a genuinely poor overall
  performer who happens to uniquely carry one target block's favourable
  value will still be selected, purely to cover that block – setting a
  floor here is the fix.

- n_reps:

  Integer \>= 1. Default `5L`. Number of independent GA replicates to
  run (see *GA rigour: replication and convergence* below); the
  replicate with the best fitness is returned as
  `selected`/`fitness`/`per_block`/`ga_fit`. Set to `1L` for the
  fastest, single-run legacy behaviour.

- G:

  Relationship/kinship matrix (n x n, dimnames = individual IDs covering
  every candidate in `value_matrix`), e.g.
  `run_haplotype_prediction()$G` or
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
  output, or `NULL` (default). Required whenever `coancestry_weight > 0`
  or `target_degree` is supplied; see *Coancestry penalty (optional)*
  below. Unused (may be left `NULL`) otherwise.

- coancestry_weight:

  Numeric \>= 0. Default `0` (no coancestry term at all – the historical
  behaviour of this function, block coverage only). When positive,
  subtracts `coancestry_weight` times the mean off-diagonal pairwise
  relationship of the chosen set (from `G`) from the fitness function,
  so the GA trades off block coverage against relatedness directly
  rather than leaving relatedness as an after-the- fact diagnostic.
  There is no universal default scale for this argument – it must be
  tuned against your own `block_weights`/fitness scale (see *Coancestry
  penalty (optional)*). Mutually exclusive with `target_degree` – do not
  supply both.

- merit_weight:

  Numeric \>= 0. Default `0` (no merit term at all – the historical
  behaviour of this function, block coverage only, with `merit_score`
  usable purely as a pre-search eligibility floor via `min_sel_value`).
  When positive, **adds** `merit_weight` times the chosen set's mean
  `merit_score` to the fitness function, so the GA rewards whole-genome
  merit directly inside the search itself, rather than only using it to
  exclude clearly ineligible candidates beforehand. See *Merit-weighted
  fitness (GA+TS hybrid, optional)* below. Requires `merit_score`.
  Mutually exclusive with `merit_priority` – do not supply both.

- merit_priority:

  Numeric in `[0, 100]`, or `NULL` (default – no effect, historical
  behaviour unchanged). The easy alternative to `merit_weight`: instead
  of a raw, unscaled multiplier, state how much you care about merit vs.
  coverage as a plain percentage, and this function calibrates the right
  `merit_weight` from your own data before running the search – see
  [`suggest_merit_weight`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)
  for exactly what "comparable scale" means here and its documented
  limits. `0` is identical to `merit_weight = 0`; `100` scales merit's
  realistic influence on the fitness function to match coverage's;
  values between scale linearly. Requires `merit_score`. Mutually
  exclusive with an explicitly-supplied `merit_weight` (ambiguous
  otherwise) – supplying both is an error.

- target_degree:

  Numeric in `[0, 90]`, or `NULL` (default – no effect, historical
  behaviour unchanged). The easy alternative to `coancestry_weight`:
  instead of a raw, unscaled penalty multiplier, state your
  relatedness-vs-gain preference on the SAME `[0, 90]` scale
  [`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
  already uses (`0` = max gain, prioritising coverage/merit and
  accepting more relatedness; `90` = max diversity, minimising
  relatedness), and this function converts it into an internal
  relatedness-ceiling penalty from your own data before running the
  search – see *Relatedness ceiling (target_degree, optional)* below for
  exactly how, and its documented limits. Requires `G`. Mutually
  exclusive with an explicitly-supplied `coancestry_weight` (ambiguous
  otherwise) – supplying both is an error.

## Value

Named list:

- `selected`:

  Character vector of the selected individual IDs from the best-fitness
  replicate (length `n_founders`, or as close as the GA achieved – check
  `length(selected) == n_founders`).

- `fitness`:

  Numeric. Best fitness value found (raw block-value sum, after
  subtracting the cardinality penalty), best-fitness replicate.

- `per_block`:

  Data frame for the best-fitness replicate: `block_id`, `best_value`,
  `contributor_1`, `contributor_2` (equal to `contributor_1` under a
  self-allowed strategy).

- `strategy`:

  Character, echoes the `strategy` argument.

- `ga_fit`:

  The raw [`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html)
  S4 result object for the best-fitness replicate, for convergence
  diagnostics (e.g. `plot(ga_fit)`).

- `converged`:

  Logical. `TRUE` if the best-fitness replicate's GA stopped because
  fitness plateaued for `run` consecutive generations
  (`ga_fit@iter < maxiter`), rather than being cut off at `maxiter`
  without plateauing. `FALSE` means consider raising `maxiter`.

- `cutoff`:

  Numeric. The `min_sel_value` cutoff actually applied to `merit_score`
  (`-Inf` when `min_sel_value = NULL`).

- `stability`:

  List describing agreement across all `n_reps` replicates – see *GA
  rigour* below. Present even when `n_reps = 1` (trivially: every
  selected individual has `selection_freq = 1`). Also includes
  `$mean_relationship` (realised mean pairwise relationship of each
  replicate's chosen set, `NA` unless `G` was supplied) and
  `$mean_merit` (realised mean `merit_score` of each replicate's chosen
  set, `NA` unless `merit_score` was supplied).

- `mean_relationship`:

  Numeric. The best-fitness replicate's realised mean off-diagonal
  pairwise relationship (from `G`) among the final `selected` set. `NA`
  unless `G` was supplied – present regardless of whether
  `coancestry_weight > 0` or `target_degree` was set, so you can inspect
  relatedness even when the GA wasn't asked to optimise for it (compare
  against an unconstrained run to see the effect).

- `mean_merit`:

  Numeric. The best-fitness replicate's realised mean `merit_score`
  among the final `selected` set. `NA` unless `merit_score` was supplied
  – present regardless of whether `merit_weight > 0`, so you can inspect
  the chosen set's whole-genome merit even when the GA wasn't asked to
  optimise for it directly (compare against a `merit_weight = 0` run to
  see the term's effect).

- `merit_weight`:

  Numeric. The merit weight actually used – either your explicit
  argument, or the value calibrated from `merit_priority` if that was
  supplied instead.

- `merit_priority`:

  Echoes the `merit_priority` argument (`NULL` unless supplied).

- `coancestry_weight`:

  Numeric. Echoes the `coancestry_weight` argument as actually used (`0`
  when `target_degree` was supplied instead – the relatedness constraint
  is enforced via `relatedness_ceiling` in that case, not this
  argument).

- `target_degree`:

  Echoes the `target_degree` argument (`NULL` unless supplied).

- `relatedness_ceiling`:

  Numeric. The interpolated mean- relationship ceiling actually enforced
  when `target_degree` was supplied (`NULL` otherwise) – see
  *Relatedness ceiling (target_degree, optional)*. Compare against
  `mean_relationship` to check whether the search respected it.

## Details

The fitness function, for a candidate founder subset \\S\\ of size
`n_founders`, is: \$\$\text{fitness}(S) = \sum\_{j \in \text{blocks}}
w_j \cdot \text{best}\_j(S)\$\$ where \\\text{best}\_j(S)\\ is the best
value block \\j\\ can achieve from the founders in \\S\\, and \\w_j\\ is
that block's weight (`block_weights`, default 1 for every block). How
\\\text{best}\_j(S)\\ is computed depends on `strategy` – see the
*Crossing-scheme strategies* section. This directly generalises the
pair-average formula shown in HapSelect's own documentation diagram
(\\\sum_j \max((\text{localGEBV}\_{j,1} +
\text{localGEBV}\_{j,2})/2)\\), but is HapBlockR's own implementation –
HapSelect's exact GA fitness function and constraint semantics are not
available beyond its published documentation, so this is not a
byte-for-byte port.

A binary GA chromosome (one bit per candidate individual) is used rather
than a fixed-cardinality representation, with a quadratic penalty for
deviating from `n_founders` selected bits – standard practice for
subset-selection problems with the GA package, since it has no native
hard-cardinality constraint.

## Crossing-scheme strategies

- `"no_selfing"` (default):

  Each block's value is the mean of the **two largest** values among the
  chosen founders – two distinct parents must jointly contribute it.
  Mirrors HapSelect's localGEBV-mode `"no_selfing"`.

- `"selfing"`:

  Each block's value is the **single largest** value among the chosen
  founders – a parent may be selfed to realise its own value alone,
  without needing a complementary partner. Mirrors HapSelect's
  localGEBV-mode `"selfing"`.

- `"OHS"`:

  Same computation as `"no_selfing"` (two distinct parents required),
  named for haplotype-mode use where each parent contributes one
  haplotype copy.

- `"OPV"`:

  Same computation as `"selfing"` (single best value) – "population
  value" framing: any one founder carrying the favourable haplotype is
  enough, poolable across the whole founder set rather than requiring a
  specific complementary pairing.

- `"Haploid_OHS"`:

  Same computation as `"selfing"` – a single (heterozygous) parent may
  donate two non-homologous gametes.

**Why `"selfing"`, `"OPV"`, and `"Haploid_OHS"` compute identically:**
under the max-pair-average formulation, allowing a founder to pair with
itself can never do worse than any distinct pair, because the largest
value in a set is always \\\geq\\ the mean of itself and anything
smaller. This is a direct mathematical consequence of the fitness
formula, not a simplification of convenience – documented here so it
isn't mistaken for the three strategies being unimplemented.

## Coancestry penalty (optional)

By default (`coancestry_weight = 0`) this function optimises block
coverage only, exactly as described in *Details* – relatedness among the
chosen founders plays no role in the search, which is why
[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)
and a manual family-balance check are recommended *after* calling this
function (see *Choosing between this function and
truncation_selection()* below). Setting `coancestry_weight > 0` (with
`G` supplied) moves that check *into* the search itself: the fitness
function becomes \$\$\text{fitness}(S) = \sum_j w_j \cdot
\text{best}\_j(S) - \text{coancestry\\weight} \cdot \overline{G}(S)\$\$
where \\\overline{G}(S)\\ is the mean off-diagonal pairwise relationship
among the chosen set \\S\\ (each individual's own diagonal
self-relationship/inbreeding term is excluded – the penalty reflects
relatedness *between* chosen parents, not their own inbreeding). This is
a direct, minimal extension of the existing block-coverage GA – not a
full optimal-contribution-selection (OCS) formulation (see *Choosing
between this function and truncation_selection()* for what a full OCS/
mate-allocation tool adds beyond this).

**Tuning `coancestry_weight`:** there is no universal correct value –
\\\overline{G}(S)\\ and the block-coverage sum are on different,
problem-specific scales (a VanRaden-style GRM is typically much smaller
in magnitude than a sum of local-GEBV values across many blocks). Start
by running with `coancestry_weight = 0` and inspecting `$fitness`'s
typical magnitude, then choose a weight that makes the coancestry term
large enough to visibly compete with it; increase further if the
resulting `selected` set is still too related for your program, and
compare `$mean_relationship` across a few candidate weights to see the
trade-off curve directly. An easier alternative to tuning this by hand
is `target_degree`, described next.

## Relatedness ceiling (target_degree, optional)

`coancestry_weight` has the same scale problem `merit_weight` does (see
*Merit-weighted fitness* below): there is no universal correct
multiplier, because \\\overline{G}(S)\\ and the block-coverage sum live
on different, dataset-specific scales. Rather than estimating a matching
scale for a weighted trade-off (the approach `merit_priority` takes for
the merit term), `target_degree` sidesteps the scale problem altogether:
it converts your `[0, 90]` preference into a relatedness CEILING – a
bound in `G`'s own real units, not an abstract multiplier – interpolated
between two genuinely constructed reference groups for your actual
candidate pool, exactly mirroring how
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)'s
`engine = "optisel"` interpolates its own kinship bound between two
solved frontier extremes:

- **Gain end (`target_degree = 0`)**: the feasible, greedy-built
  high-coverage group (see *Details*) – extended to also account for
  `merit_weight` when merit is active in this same call, so this
  reference reflects what the search would actually converge to with no
  relatedness constraint, not a coverage-only proxy. Its own mean
  pairwise relationship becomes the ceiling at `target_degree = 0` (i.e.
  no additional constraint beyond what that group naturally has).

- **Diversity end (`target_degree = 90`)**:
  [`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)`(strategy = "maximin")`
  run on `G` for your candidate pool – a proven 2-approximation for
  maximizing the minimum pairwise distance (Gonzalez 1985), built
  ignoring merit entirely. Its own mean pairwise relationship becomes
  the ceiling at `target_degree = 90`.

The ceiling at any `target_degree` in between is a linear interpolation
between these two real, feasible reference points. Inside the search,
exceeding the ceiling is penalised by a squared-violation term scaled to
dominate any possible coverage/merit gain from crossing it (the same
soft-constraint idiom already used for the cardinality penalty) – so no
separate penalty-strength argument is needed.

**Honest limits:** the gain-end reference is a plain greedy hill- climb,
not a formally guaranteed-near-optimal one, once `merit_weight > 0`
(adding a mean-based merit term breaks the strict submodularity argument
that gives the coverage-only greedy its proof) – still a real, feasible
group, just without that guarantee. Separately,
`select_core_collection`'s 2-approximation guarantee covers the MINIMUM
pairwise distance in its chosen group, not directly the MEAN pairwise
relationship used here as the diversity-end reference – a principled,
reused building block, not a proof about the specific number reported.
Neither limitation is unique to this feature; they mirror the same kind
of honestly-stated approximation caveats already documented for
`merit_priority` (see
[`suggest_merit_weight`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md))
and for `select_parents_ocs`'s own `target_degree`.

## Merit-weighted fitness (GA+TS hybrid, optional)

By default (`merit_weight = 0`) this function's search rewards block
coverage only, exactly as described in *Details* – `merit_score`, if
supplied at all, only ever acts as a pre-search eligibility floor via
`min_sel_value` (a hard yes/no cutoff), never as something the search
itself is rewarded for pursuing further. That matters because a floor
alone does not distinguish between a candidate that barely clears it and
one that clears it by a wide margin – both are equally "eligible," and
among the eligible pool the search optimises purely for which target
blocks each candidate happens to carry good value at. On a real panel,
where per-block/local-GEBV estimates carry real estimation noise (they
are themselves statistical estimates, not ground truth), this can let
the search lean on a candidate whose apparent block coverage is partly
or wholly a noisy artefact, provided that candidate clears the floor at
all.

Setting `merit_weight > 0` closes this gap directly: the fitness
function becomes \$\$\text{fitness}(S) = \sum_j w_j \cdot
\text{best}\_j(S) + \text{merit\\weight} \cdot
\overline{\text{merit\\score}}(S)\$\$ (combined with the coancestry
penalty, when both are active, as a third additive term). The search now
has a direct, continuous incentive to prefer higher-merit candidates
*throughout* the eligible pool, not merely to avoid excluding low-merit
ones – pulling the result back toward
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)'s
whole-genome-merit ranking while still searching for joint block
coverage, rather than treating merit and coverage as two entirely
separate stages (a floor, then an unconstrained coverage search). This
is the recommended way to combine this function with a
truncation-selection-style merit signal when you are specifically
concerned that block-coverage-only selection, run on noisy per-block
estimates, could otherwise select on noise: a nonzero `merit_weight`
keeps whole-genome merit – typically a more stable, better-estimated
signal than any single block's local GEBV – pulling on every candidate's
desirability throughout the search, not only at its entry gate.

**Tuning `merit_weight`:** the same scale caveat as `coancestry_weight`
applies, and for the same reason – \\\overline{\text{merit\\score}}(S)\\
and the block-coverage sum are on different, problem-specific scales, so
there is no universal correct value. Start by running with
`merit_weight = 0` and inspecting `$fitness`'s typical magnitude and
`$mean_merit`'s typical magnitude, then choose a weight that makes the
merit term large enough to visibly compete with the block-coverage sum;
compare `$selected` and `$mean_merit` across a few candidate weights
(and against a plain
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
run on the same `merit_score`) to see the effect directly, exactly as
recommended for `coancestry_weight` above. `min_sel_value` and
`merit_weight` address different failure modes and are not substitutes
for each other: keep using `min_sel_value` for a hard floor against
candidates you never want to consider at all, and add `merit_weight` on
top of it when you also want merit to keep influencing the choice among
everyone who clears that floor.

**An easier alternative to tuning `merit_weight` by hand:** set
`merit_priority` instead (a plain 0-100 "how much do I care about merit
vs. coverage" dial) and let this function calibrate the matching
`merit_weight` from your own data automatically – see
[`suggest_merit_weight`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)
for exactly what it computes, why (a greedy, provably near-optimal
reachable high-coverage group as the reference point, not an inflated
unreachable one), and what it honestly does not solve (there is no
single universally "correct" notion of comparable scale –
`merit_priority`'s suggestion is a well-reasoned starting point, not the
only defensible number). `merit_weight` and `merit_priority` are
mutually exclusive – pick whichever style suits you; you can also call
[`suggest_merit_weight`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)
first to inspect the calibration numbers, then pass a literal
`merit_weight` yourself if you want to fine-tune from there.

## GA rigour

replication and convergence: A single GA run, taken at face value, tells
you nothing about whether its answer is a robust optimum or one of
several near-equally-good solutions a stochastic search happened to land
on. This function addresses that directly rather than leaving it to the
caller:

- Convergence:

  Every replicate's `ga_fit@iter` (generations actually run) is checked
  against `maxiter`. Stopping early means GA's own
  `run`-generations-without-improvement criterion triggered – a real
  plateau. Running the full `maxiter` without stopping early means the
  search may not have converged; `$converged` surfaces this instead of
  silently returning a possibly-unconverged result.

- Replication:

  With `n_reps > 1` (default `5`), the search is repeated from `n_reps`
  different starting populations/seeds. `$stability$selection_freq`
  reports, for every candidate selected in at least one replicate, the
  fraction of replicates that selected them – individuals at `1.0` are
  robustly supported regardless of the GA's random starting point;
  individuals selected in only one replicate out of several are
  borderline and worth a second look before committing to them.
  `$stability$fitness_range` shows how much the best achievable fitness
  varied across replicates – a narrow range alongside high selection
  frequencies is the signature of a stable, trustworthy search.

This is deliberately more rigorous by default than a single fixed-seed
GA run: HapSelect's own documented parent-selection GA (see
*Description*) does not report replication stability or a convergence
flag at all. Set `n_reps = 1` only once you have separately confirmed
stability, or for quick iteration during exploratory analysis.

## Choosing between this function and truncation_selection()

[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
ranks by a single whole-genome score and takes the top `n_founders`;
this function instead searches for the *set* of `n_founders` that
jointly covers the most favourable value across a chosen set of target
blocks (typically the output of
[`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)).
They answer genuinely different questions and neither is a strict
improvement on the other:

- [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
  ignores:

  which haplotype blocks each individual carries – its top-N list can
  concentrate on the same favourable blocks while leaving others
  uncovered.

- `select_parents_ga()` ignores:

  whole-genome merit by default, unless `merit_score` is supplied via
  `min_sel_value` (a hard pre-search floor), `merit_weight` (a soft,
  continuous term rewarded throughout the search – see *Merit-weighted
  fitness (GA+TS hybrid, optional)*), or both together. Without either,
  a genuinely poor-performing individual who happens to uniquely carry
  one target block's favourable value will still be selected, purely to
  cover that block – and, even with only a floor and no weight, a
  candidate who barely clears the floor is treated identically to one
  who clears it by a wide margin. Setting `merit_weight > 0` is the
  tool's own built-in "GA+TS hybrid" mode: it keeps searching for joint
  block coverage while also rewarding whole-genome merit directly inside
  the search, which matters most when block-level/local-GEBV estimates
  carry real estimation noise and a coverage-only search risks selecting
  on that noise rather than on genuine signal.

- **Neither function, by default**:

  manages coancestry/ inbreeding risk in the chosen set, assigns
  differential contributions, or decides who mates whom. This function's
  optional `coancestry_weight`/`target_degree`/`G` arguments (see
  *Coancestry penalty (optional)* and *Relatedness ceiling
  (target_degree, optional)*) add a relatedness term to the search
  itself, but still stop at a flat founder set – no differential
  contributions, no mating list. Dedicated
  optimal-contribution-selection (OCS) and mate-allocation tools (e.g.
  AlphaMate, Gorjanc lab) solve all three problems jointly, typically
  maximising gain subject to an explicit coancestry constraint plus
  continuous contribution optimisation. A practical combined workflow:
  use this function (or
  [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md))
  to shortlist candidates worth considering, then hand that shortlist to
  a dedicated OCS tool for the coancestry-managed contribution and
  mating decision, rather than treating either function's output as a
  final mating plan.

In practice, run both and compare with
`intersect(ga_out$selected, ts_out$selected)`: a high overlap means
truncation selection was already close to the GA's answer for this
panel; a low overlap means block coverage and whole-genome ranking
disagree and both founder sets are worth inspecting (including a
family/diversity check via
[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md))
before committing to either. See the *From Local GEBV to a Crossing
Decision* vignette for a full worked comparison.

## References

Scrucca L (2013). GA: A Package for Genetic Algorithms in R. *Journal of
Statistical Software* 53(4):1-37.
[doi:10.18637/jss.v053.i04](https://doi.org/10.18637/jss.v053.i04)

## See also

[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md),
[`plot_parent_selection_pca`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md),
[`ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md),
[`suggest_merit_weight`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
(the `target_degree` convention this function reuses),
[`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
(`target_degree`'s diversity-end reference)

## Examples

``` r
if (FALSE) { # \dontrun{
res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
top  <- select_top_blocks(res$block_importance, n = 15)
vmat <- res$local_gebv[, top$block_id, drop = FALSE]
ga_out <- select_parents_ga(vmat, n_founders = 20, strategy = "no_selfing",
                            seed = 1)
ga_out$selected
} # }
```
