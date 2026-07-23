# Family- or Genetic-Cluster-Quota Parent Selection (Best Groups, Then Best Lines Within Them)

A third parent-shortlist strategy, alongside
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
(rank every candidate by a single score) and
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
(search for the set jointly covering the most target-block value). This
function instead mirrors the shape a breeding program's shortlist
usually already takes in practice: pick the `n_families` best-performing
groups first, then the `n_per_family` best individual lines within each
chosen group. A "group" is either a pedigree/cross-ID label you supply
(`group_by = "family"`, the default) or a genetically data-derived
cluster built directly from a relationship matrix
(`group_by = "genetic_cluster"`) – see *Grouping: pedigree family vs.
genetic cluster* below. Unlike
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
which ranks the whole candidate pool as one list and can, by chance,
draw heavily from a small number of related groups, this function makes
group representation an explicit part of the selection rule itself –
every chosen group contributes exactly `n_per_family` lines (or fewer,
if it has fewer eligible members), regardless of how its individual
lines would have ranked in a single population-wide list.

## Usage

``` r
select_parents_by_family(
  score,
  family = NULL,
  n_families,
  n_per_family = NULL,
  family_select_mode = c("count", "percentage", "sd_threshold", "check_relative"),
  pct_per_family = NULL,
  sd_threshold = NULL,
  check_id = NULL,
  check_value = NULL,
  check_margin_pct = NULL,
  min_sel_value = NULL,
  min_sel_mode = c("value", "percentile", "sd_below_mean"),
  ensure_haplotype_diversity = FALSE,
  value_matrix = NULL,
  haplotypes = NULL,
  block_weights = NULL,
  G = NULL,
  within_group_target_degree = NULL,
  group_by = c("family", "genetic_cluster"),
  rank_k = NULL,
  family_rank_method = c("shrunk_topk_mean", "topk_mean"),
  variance_method = c("reml", "anova"),
  bias_correction = c("order_stats", "none"),
  use_family_relationship = TRUE,
  diversity_method = c("coverage_gain", "dominant_block"),
  n_clusters = NULL,
  cluster_method = "ward.D2",
  verbose = TRUE
)
```

## Arguments

- score:

  Named numeric vector, e.g. whole-genome GEBV
  (`run_haplotype_prediction()$gebv`) or any other selection index.
  Names are individual IDs.

- family:

  Named character or factor vector giving each individual's
  pedigree/family group. Required when `group_by = "family"` (used as
  the actual grouping). Optional when `group_by = "genetic_cluster"`
  (purely echoed in `by_family$family` for cross-referencing; does not
  affect which lines are selected in that mode). Names must cover every
  name in `score` to be considered (individuals in `family` but not
  `score`, or vice versa, are simply not matched).

- n_families:

  Integer. Number of top-ranked groups to select from. If fewer than
  `n_families` distinct groups have at least one eligible member, all of
  them are used and a warning is issued.

- n_per_family:

  Integer, or a NAMED integer vector (names = chosen family/group
  labels) for uneven quotas. Required, and only used, when
  `family_select_mode = "count"` (the default). Number of best lines to
  select from each chosen group. If a chosen group has fewer than its
  quota of eligible members, all of its members are taken and a warning
  is issued. Also the default ranking sample size (`rank_k`) when
  `rank_k` is not supplied separately and `family_select_mode = "count"`
  – see *Details*.

- family_select_mode:

  One of `"count"` (default), `"percentage"`, `"sd_threshold"`, or
  `"check_relative"`. See *family_select_mode: four ways to decide a
  family's take*.

- pct_per_family:

  Numeric in `(0, 100]`. Required, and only used, when
  `family_select_mode = "percentage"`.

- sd_threshold:

  Numeric (signed). Required, and only used, when
  `family_select_mode = "sd_threshold"`.

- check_id:

  Character, a single individual ID present in `score`. Used, and
  mutually exclusive with `check_value`, when
  `family_select_mode = "check_relative"`.

- check_value:

  Numeric, a single fixed benchmark score. Used, and mutually exclusive
  with `check_id`, when `family_select_mode = "check_relative"` (for a
  check that is not itself a scored candidate in this run).

- check_margin_pct:

  Numeric (signed), e.g. `10` for "at least 10% better than the check".
  Required, and only used, when `family_select_mode = "check_relative"`.

- min_sel_value, min_sel_mode:

  Optional merit floor applied to `score` *before* group ranking or
  within-group selection, exactly as in
  [`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)/[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md).
  Default `min_sel_value = NULL` applies no floor.

- ensure_haplotype_diversity:

  Logical, default `FALSE`. See *Haplotype/coverage diversity adjustment
  (optional)*. Requires `value_matrix`.

- value_matrix:

  Numeric matrix (individuals x target blocks), e.g. local GEBV
  (`run_haplotype_prediction()$local_gebv`, typically restricted to
  [`select_top_blocks`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)'s
  output first) or haplotype allele dosage – the same kind of input
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  takes. Row names = individual IDs, must cover every individual under
  consideration once `min_sel_value` and group membership have been
  applied. Required when `ensure_haplotype_diversity = TRUE`; ignored
  (may be left `NULL`) otherwise.

- haplotypes:

  Optional named list, the direct, unmodified return value of
  [`extract_haplotypes`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  (one element per block, each a named character vector of one
  dosage/phased allele string per individual). Only used by
  `diversity_method = "dominant_block"`'s allele-level collision check –
  see *Haplotype/coverage diversity adjustment (optional)*. Ignored
  under `diversity_method = "coverage_gain"`.

- block_weights:

  Numeric vector, length `ncol(value_matrix)`, or `NULL` (default: equal
  weight 1 for every block). Used by
  `diversity_method = "coverage_gain"` exactly as
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
  argument of the same name is.

- G:

  Relationship/kinship matrix (n x n, dimnames = individual IDs), e.g.
  `run_haplotype_prediction()$G` or
  [`compute_haplotype_grm`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
  output, or `NULL` (default). Required when
  `group_by = "genetic_cluster"` or `within_group_target_degree` is
  supplied; optional otherwise (in which case supplying it still unlocks
  the `$mean_relationship` diagnostics – see *Relatedness control
  (optional)* – and, together with `n_clusters`, the optional
  `genetic_group` cross-reference under `group_by = "family"`).

- within_group_target_degree:

  Numeric in `[0, 90]`, or `NULL` (default – no effect). See
  *Relatedness control (optional)*. Requires `G`.

- group_by:

  One of `"family"` (default) or `"genetic_cluster"`. See *Grouping:
  pedigree family vs. genetic cluster*.

- rank_k:

  Integer, or `NULL` (default). The number of each group's own top
  members its ranking criterion (`topk_mean`/ `shrunk_topk_mean`) is
  computed from, independent of the actual take (`n_per_family` or
  otherwise, per `family_select_mode`). Defaults to `n_per_family` (or
  its maximum, if `n_per_family` is a vector) when
  `family_select_mode = "count"` and not supplied separately – the
  historical behaviour; defaults to `1` for the other three modes, which
  have no take-quota to borrow a default from.

- family_rank_method:

  One of `"shrunk_topk_mean"` (default) or `"topk_mean"`. See
  *Shrinkage-corrected family ranking*.

- variance_method:

  One of `"reml"` (default) or `"anova"`. How \\\tau^2\\/\\\sigma^2\\
  are estimated for shrinkage. Ignored under
  `family_rank_method = "topk_mean"`. See *Shrinkage- corrected family
  ranking*.

- bias_correction:

  One of `"order_stats"` (default) or `"none"`. Whether each group's
  `topk_mean` is corrected for finite-population selection-differential
  bias before shrinkage. Ignored under
  `family_rank_method = "topk_mean"`. See *Shrinkage-corrected family
  ranking*.

- use_family_relationship:

  Logical, default `TRUE`. Whether plain i.i.d. shrinkage is replaced by
  genomic-relationship-informed (GBLUP-style) family-effect shrinkage
  when `G` is supplied and covers every eligible group. Ignored under
  `family_rank_method = "topk_mean"` or when `G` is `NULL`. See
  *Shrinkage- corrected family ranking*.

- diversity_method:

  One of `"coverage_gain"` (default) or `"dominant_block"`. See
  *Haplotype/coverage diversity adjustment (optional)*. Ignored unless
  `ensure_haplotype_diversity = TRUE`.

- n_clusters:

  Integer. Number of genetic clusters to create. Required when
  `group_by = "genetic_cluster"`; also enables the optional
  `genetic_group` diagnostic cross-reference under `group_by = "family"`
  when supplied alongside `G`.

- cluster_method:

  Character, default `"ward.D2"`. Linkage method passed to
  [`hclust`](https://rdrr.io/r/stats/hclust.html) when building genetic
  clusters.

- verbose:

  Logical, default `TRUE`. Print informational messages (shrinkage
  fallback, skipped diagnostics).

## Value

Named list:

- `selected`:

  Character vector of all selected individual IDs, ordered by group
  rank, then by selection order within each group.

- `by_family`:

  Data frame, one row per selected individual: `family` (pedigree label,
  `NA` if not available – see *Grouping*), `genetic_group`
  (auto-generated cluster label, `NA` if not computed), `family_rank`,
  `individual`, `score`, `rank_within_family`, and, when
  `ensure_haplotype_diversity = TRUE`, `dominant_block` (when
  `diversity_method = "dominant_block"`) and `collision` (logical – see
  *Haplotype/coverage diversity adjustment*).

- `family_ranking`:

  Data frame, one row per group considered (after the family-size
  eligibility rule has removed undersized groups – see *Family-size
  eligibility rule*): `family` (whichever grouping was ranked – pedigree
  or genetic cluster, per `group_by`), `topk_mean` (raw ranking
  criterion), `topk_mean_bias`/`topk_mean_corrected` (present when
  shrinkage was applied – the finite-population selection-bias estimate
  subtracted, and the corrected value, respectively; bias is 0 and
  corrected equals raw whenever `bias_correction = "none"` or no
  correction could be computed), `shrinkage_weight` (`NA` unless
  `family_rank_method = "shrunk_topk_mean"`, shrinkage was actually
  applied, AND it used the plain scalar i.i.d. formula rather than
  genomic-relationship-informed shrinkage, which has no single scalar
  weight), `rank_score` (the value groups were actually ranked on),
  `n_members`, `rank`, `selected` (logical), and `mean_relationship`
  (that group's own selected picks' realised mean relationship, `NA`
  unless `G` was supplied or the group was not selected).

- `cutoff`:

  Numeric. The `min_sel_value` cutoff actually applied (`-Inf` when
  `min_sel_value = NULL`).

- `n_families`, `n_per_family`, `family_select_mode`, `pct_per_family`,
  `sd_threshold`, `check_id`, `check_value`, `check_margin_pct`,
  `rank_k`, `family_rank_method`, `variance_method`, `bias_correction`,
  `use_family_relationship`, `group_by`, `n_clusters`, `cluster_method`,
  `diversity_method`, `within_group_target_degree`:

  Echo the corresponding arguments (resolved values where applicable;
  `NULL` for whichever of `pct_per_family`/`sd_threshold`/
  `check_id`/`check_value`/`check_margin_pct` were not relevant to the
  `family_select_mode` actually used).

- `variance_method_used`:

  Character, `"reml"` or `"anova"` – which estimator actually produced
  `tau2`/ `sigma2`, which can differ from `variance_method` when
  `"reml"` was requested but fell back to `"anova"` (see
  *Shrinkage-corrected family ranking*). `NA` when shrinkage was not
  applied at all (`family_rank_method = "topk_mean"`, or variance
  components could not be estimated).

- `relationship_informed`:

  Logical. Whether the final ranking actually used
  genomic-relationship-informed (GBLUP-style) shrinkage rather than the
  plain i.i.d. formula (or no shrinkage at all).

- `tau2`, `sigma2`:

  Numeric. The estimated between-group and within-group/residual
  variance components underlying shrinkage (`NULL` when shrinkage was
  not applied).

- `excluded_groups`:

  Character vector of family/group labels removed by the family-size
  eligibility rule before ranking (empty if none were excluded, and
  always empty under
  `family_select_mode %in% c("sd_threshold", "check_relative")`, which
  has no such rule). See *Family-size eligibility rule*.

- `zero_selected_groups`:

  Character vector of ranked family/group labels that were tried but
  contributed zero qualifying members and were passed over via
  auto-backfill (always empty under
  `family_select_mode %in% c("count", "percentage")`, which cannot
  produce a zero-member take). See *family_select_mode: four ways to
  decide a family's take*.

- `ensure_haplotype_diversity`:

  Echoes the argument.

- `mean_relationship`:

  Numeric. Realised mean off-diagonal pairwise relationship among all of
  `selected`. `NA` unless `G` was supplied.

## Details

Groups are ranked by a (by default, shrinkage-corrected) function of
each group's own top `rank_k` members (`rank_k` defaults to
`n_per_family` if not supplied separately) – not by a single best
individual (which would favour one-hit-wonder groups) and not by the
group's whole-membership mean (which would penalise a large group for
having a long tail below its own best lines). See *Shrinkage-corrected
family ranking* below for exactly what "shrinkage-corrected" means and
why. The `n_families` highest-ranked groups by this criterion are kept,
and within each, `n_per_family` members are selected by `score` (subject
to the optional relatedness and haplotype-diversity adjustments below).

## Grouping

pedigree family vs. genetic cluster: `group_by = "family"` (default)
uses your own `family` labels exactly as before.
`group_by = "genetic_cluster"` instead builds groups directly from `G`
(a relationship matrix): hierarchical clustering (Ward's
minimum-variance linkage, `cluster_method = "ward.D2"` by default – the
standard, well-justified choice for compact, genetically coherent
clusters) on the distance matrix implied by `G` (the same exact identity
[`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
uses, not an approximation), cut into `n_clusters` groups. This matters
because pedigree labels are not always a faithful proxy for actual
relatedness – incomplete records, mislabelling, or loosely applied
family IDs can all decouple "same family label" from "actually closely
related," and two DIFFERENT family labels can still correspond to
closely related individuals (e.g. a shared grandparent). Genetic
clustering sidesteps that by building groups directly from the genomic
data itself.

Both labels are reported per selected individual whenever available,
independent of which one actually drove the selection decision (echoed
in `$group_by`): `by_family$family` holds your pedigree label (from the
`family` argument, `NA` if not supplied – supplying it alongside
`group_by = "genetic_cluster"` is purely for cross-referencing, it does
not affect the selection rule), and `by_family$genetic_group` holds the
auto-generated genetic cluster label (`"cluster_1"`, `"cluster_2"`, ...,
relabelled by decreasing mean `score` so `"cluster_1"` is consistently
the strongest-performing genetic group) – populated whenever clustering
was computed, which is always true under `group_by = "genetic_cluster"`
and also true under `group_by = "family"` if you additionally supply `G`
and `n_clusters` purely as a diagnostic cross-check (no effect on which
lines are chosen in that case). `family_ranking`'s own `family` column
always holds whichever grouping was actually ranked/selected on
(pedigree labels or genetic-cluster labels, per `group_by`) – it does
not attempt a group-level crosswalk, since a pedigree family can span
multiple genetic clusters and vice versa; the individual-level crosswalk
in `by_family` is where that correspondence is meaningful.

## Shrinkage-corrected family ranking

A group's raw top-`rank_k` mean score carries two different,
independently-corrected biases rather than being trusted at face value
regardless of group size.

**Reliability shrinkage.** A small group's raw top-`rank_k` mean is an
unreliable estimator of that group's true quality – exactly the kind of
small-sample noise classical quantitative genetics corrects for via
shrinkage / BLUP-style family evaluation.
`family_rank_method = "shrunk_topk_mean"` (the default) implements this:
treating each group's mean score as a random effect (\\\text{score} =
\mu + a_f + e\\, \\a_f \sim N(0, \tau^2)\\, \\e \sim N(0, \sigma^2)\\),
\\\tau^2\\ (true between-group variance) and \\\sigma^2\\
(within-group/residual variance) are estimated from EVERY eligible
member's raw `score` (not from `topk_mean` itself – fitting variance
components on an already-selected/truncated statistic would bias the
estimate) via `variance_method`: either `"reml"` (the default –
`lme4::lmer(score ~ 1 + (1|family))`; more statistically efficient than
method-of-moments under unbalanced group sizes, which this function
always has to deal with, since uneven family size is the entire reason
shrinkage is needed here in the first place) or `"anova"` (the classical
one-way random-effects ANOVA method-of-moments estimator, Searle,
Casella & McCulloch 1992 – also the automatic, non-silent fallback
whenever `"reml"` is requested but lme4 is not installed, there isn't
enough structure to fit the model, or the fit fails). Each group's
(optionally bias-corrected – see below) `topk_mean` is then shrunk
toward the across-group mean by the resulting empirical-Bayes/BLUP
reliability weight \\w_f = n_f \tau^2 / (n_f \tau^2 + \sigma^2)\\ (e.g.
Robinson 1991, "That BLUP Is a Good Thing") – a group of 1 member is
pulled hard toward the pack; a 20-member group's genuine top-k average
is trusted almost fully – UNLESS genomic-relationship-informed shrinkage
(below) is available, in which case it replaces this scalar formula.
`family_rank_method = "topk_mean"` restores the historical, unshrunk
ranking (`variance_method`, `bias_correction` and
`use_family_relationship` are all ignored in that case). Honest limit:
this shrinks the *selection statistic* (`topk_mean`) using reliability
weights estimated from the *full* within-group distribution – a
well-justified, practical empirical-Bayes approach, not a single unified
formal model of `topk_mean` itself. Falls back to unshrunk ranking with
a message, not silently, whenever variance components cannot be
estimated at all (e.g. every eligible group has exactly 1 member).

**Finite-population selection-bias correction.** Separately from
reliability, "mean of a group's own best `rank_k`" is a systematically
*upward*-biased estimator of that group's true mean, purely from having
picked winners out of a finite pool – present even for a group whose
individual scores carry zero estimation error of their own, and
shrinking alone does not remove it. `bias_correction = "order_stats"`
(the default) removes this before shrinkage is applied, by subtracting
each group's expected order-statistics selection differential
(\\\sqrt{\sigma^2} \times\\ the expected mean of the top `rank_k` order
statistics of `n_members` i.i.d. standard normal draws, computed by
direct numerical integration of the order statistic's own density rather
than a tabulated closed-form approximation). `bias_correction = "none"`
restores the uncorrected `topk_mean`. Requires a finite, positive
\\\sigma^2\\ from the variance-component step above; has no effect (bias
= 0) otherwise.

**Genomic-relationship-informed shrinkage.** By default
(`use_family_relationship = TRUE`), whenever `G` is also supplied and
covers every eligible group, the plain scalar \\w_f\\ shrinkage above is
replaced by a full GBLUP-style family-effect BLUP, solved via
Henderson's mixed-model equations from a family-level relationship
matrix (each entry the mean pairwise `G` value between two groups'
members). This lets genetically related groups borrow strength from EACH
OTHER, not only from the grand mean – a strict generalisation of the
scalar \\w_f\\ formula, which it reduces to exactly when groups are
mutually unrelated. Falls back to the plain i.i.d. \\w_f\\ formula
(message, not silently) whenever `G` does not cover every eligible
group, or `use_family_relationship = FALSE`.

## family_select_mode

four ways to decide a family's take: Once the `n_families` best-ranked
groups are chosen, four independent, mutually exclusive modes decide how
many/which of a group's own members are actually taken:

- `"count"` (default):

  `n_per_family` lines per group – a fixed number, identical to the
  original design.

- `"percentage"`:

  `pct_per_family` percent of THAT group's own eligible size (rounded UP
  via [`ceiling()`](https://rdrr.io/r/base/Round.html), so an included
  group never contributes zero – for small groups this means the
  effective share taken can be noticeably higher than the nominal
  percentage, an unavoidable discretisation effect of rounding a
  fraction of a small integer up rather than down).

- `"sd_threshold"`:

  Every group member whose `score` is at least `sd_threshold` standard
  deviations above (or, if negative, below) the MEAN of every eligible
  candidate in this run (population-referenced, not group-referenced) is
  taken – however many that turns out to be.

- `"check_relative"`:

  Every group member whose `score` is at least `check_margin_pct`
  percent above a reference check's score is taken – e.g.
  `check_margin_pct = 10` keeps everyone at least 10% above the check,
  including lines that clear it by 20% or 30%. The check's own score is
  either looked up via `check_id` (a candidate already present in
  `score`) or supplied directly as a fixed `check_value` (for a check
  that is not itself a scored candidate in this run) – exactly one of
  the two is required.

`"count"` and `"percentage"` are QUOTA rules: they always produce a
fixed number of lines (once resolved) and so keep a family-size
eligibility rule (below) that excludes a group outright before ranking
whenever taking that many/that share would mean taking every member
regardless of merit. `"sd_threshold"` and `"check_relative"` are
THRESHOLD rules instead: how many members clear a merit bar is genuine
information, not a degenerate case tied to group size, so there is no
size-based pre-exclusion for them – but a chosen group can then
legitimately contribute zero members if nobody clears the bar. When that
happens, that group is dropped (reported in `$zero_selected_groups`,
with a message unless `verbose = FALSE`) and the next-ranked eligible
group is tried instead, so the final shortlist still ends up with
`n_families` worth of CONTRIBUTING groups (auto-backfill) whenever
enough ranked groups exist.

## Family-size eligibility rule (`"count"`/`"percentage"` only)

A family/group is excluded entirely – before ranking, before it can
count toward `n_families` – whenever its eligible membership does not
exceed the number of lines that would actually be taken from it: under
`"count"`, `n_per_family` itself when it is a single scalar (the common,
flat-quota case), or `rank_k` when `n_per_family` is a named vector for
uneven per-group quotas (since `n_per_family` is then only defined for
whichever groups end up chosen, which is not yet known at exclusion time
– `rank_k`, "the number of lines this group's ranking is based on," is
used as the practical stand-in); under `"percentage"`, that group's OWN
`ceiling(pct_per_family/100 * group_size)`. Below that size there is no
genuine "select the best of" decision for that group at all – every
eligible member would be taken regardless of ranking – so including it
would let a group that was never really competing count toward
`n_families` for free, and would feed a `topk_mean` with zero real
selection content into the shrinkage/bias-correction machinery above
(whose entire premise is that some winnowing happened). Under
`"percentage"`, this is considerably narrower than under `"count"`: the
exclusion condition reduces to `pct_per_family > 100*(n-1)/n` for a
group of size `n`, a threshold that rises toward 100 as `n` grows. At a
realistic rate such as 10%, only single-member groups (`n = 1`) are ever
excluded – any group of 2 or more is untouched. It only starts excluding
larger groups too at high percentages (e.g. 90% excludes every group of
size 9 or smaller). Excluded family/group labels are reported, not
silently dropped, in `$excluded_groups`; an informational message names
them (unless `verbose = FALSE`). An error is raised if nothing remains
once undersized groups are excluded. Does not apply under
`"sd_threshold"`/ `"check_relative"` – see *family_select_mode* above.

## Relatedness control (optional)

Whenever `G` is supplied, the realised mean off-diagonal pairwise
relationship of the final `selected` set (`$mean_relationship`) and of
each selected group's own picks (`family_ranking`'s `mean_relationship`
column) is always reported – present regardless of whether anything was
asked to optimise for it, mirroring
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
own `$mean_relationship`. On top of that, `within_group_target_degree`
(numeric in `[0, 90]`, same convention as
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)/[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md):
`0` = max gain, prioritising `score` and accepting more relatedness;
`90` = max diversity, minimising relatedness) turns this into an actual
dial: filling a group's `n_per_family` slots stops being pure
score-order truncation and instead balances score against a relatedness
ceiling, interpolated between the group's own top-scoring feasible
subset and
[`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)` (strategy = "maximin")`
restricted to that group's own members (a proven 2-approximation,
Gonzalez 1985). Never drops a group's best remaining member purely to
satisfy the relatedness preference – if the quota cannot be filled
without a ceiling violation, the remaining slots are filled from score
order anyway (see
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
for the analogous guarantee in its own `target_degree` documentation).
When `ensure_haplotype_diversity` is ALSO active, relatedness filtering
is layered on top of the coverage/diversity-adjusted preference order,
in that disclosed order – the two constraints are not jointly
co-optimised, consistent with this function's overall greedy,
honestly-scoped design (see *Haplotype/coverage diversity adjustment*
below).

## Haplotype/coverage diversity adjustment (optional)

With `ensure_haplotype_diversity = FALSE` (the default), each chosen
group's `n_per_family` slots are filled by `score` alone. With it set to
`TRUE`, selection proceeds group by group in rank order (best group
first), slot by slot within each (best-scoring remaining member first),
maintaining a registry of what has already been claimed by earlier picks
across *all* groups processed so far. Two methods:

- `diversity_method = "coverage_gain"` (default):

  Reuses
  [`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
  own block-coverage machinery (`.block_best_values()`) to ask whether a
  candidate's addition actually increases total value-weighted coverage
  across every target block in `value_matrix`, given everyone already
  claimed so far – not merely whether one single block collides. The
  first remaining member (in descending `score` order) whose addition
  clears a near-zero gain threshold is taken; if every remaining member
  is already fully redundant with what is claimed, the highest-scoring
  one is taken anyway and flagged `collision = TRUE`.

- `diversity_method = "dominant_block"`:

  The original, coarser heuristic: each individual's *dominant* target
  block – the single target block at which that line's own value in
  `value_matrix` is highest. The first remaining member (score order)
  whose dominant block is not already claimed is taken; if every
  remaining member's dominant block is claimed and `haplotypes` was
  supplied, an allele-level check against the specific claim(s) at that
  block can still accept a candidate whose actual allele differs;
  otherwise the highest-scoring remaining member is taken anyway,
  flagged `collision = TRUE`.

Both are greedy, group-by-group heuristics, not a joint search over the
whole shortlist the way
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)'s
fitness function is – intentionally simpler, at the cost of not
necessarily finding the globally best assignment of lines to slots.

## When to reach for this instead of truncation_selection() or select_parents_ga()

Reach for this function specifically when your program's real shortlist
decision already has the shape "our best few groups, and our best few
lines out of each" – a common, familiar structure in programs organised
around discrete biparental or half-sib families, or when you want a
data-derived genetic-cluster version of that same structure
(`group_by = "genetic_cluster"`) rather than trusting pedigree labels
alone – and you want that structure enforced directly by the selection
rule, rather than checking group balance only as a diagnostic *after*
running
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
or
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md).

Skip it, and use
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
or
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
instead, if your program does not organise around discrete groups in the
first place, or if a strict per-group quota is not actually a constraint
you want enforced.

## References

Searle, S.R., Casella, G. & McCulloch, C.E. (1992). *Variance
Components*. Wiley.

Robinson, G.K. (1991). That BLUP is a Good Thing: The Estimation of
Random Effects. *Statistical Science*, 6(1), 15-32.

Gonzalez, T.F. (1985). Clustering to minimize the maximum intercluster
distance. *Theoretical Computer Science*, 38, 293-306.

Bates, D., Machler, M., Bolker, B. & Walker, S. (2015). Fitting Linear
Mixed-Effects Models Using lme4. *Journal of Statistical Software*,
67(1), 1-48. (`variance_method = "reml"`.)

Falconer, D.S. & Mackay, T.F.C. (1996). *Introduction to Quantitative
Genetics* (4th ed.). Longman. (Within-family selection –
`family_select_mode = "percentage"`/`"sd_threshold"`.)

## See also

[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`select_core_collection`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
(genetic-cluster grouping's distance conversion and
`within_group_target_degree`'s diversity end both reuse it),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)

## Examples

``` r
if (FALSE) { # \dontrun{
haps <- extract_haplotypes(geno, snp_info, blocks)   # same call already
                                                      # used upstream
res  <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues_vec)
top  <- select_top_blocks(res$block_importance, n = 15)
vmat <- res$local_gebv[, top$block_id, drop = FALSE]

# Pedigree-family mode, with relatedness control and coverage-gain
# diversity adjustment.
fam_sel <- select_parents_by_family(
  score         = res$gebv,
  family        = family_id,      # named vector, same names as res$gebv
  n_families    = 6L,
  n_per_family  = 3L,
  ensure_haplotype_diversity = TRUE,
  value_matrix  = vmat,
  G             = res$G,
  within_group_target_degree = 30
)
fam_sel$selected
fam_sel$by_family

# Genetic-cluster mode instead of pedigree family, with the pedigree
# label still echoed for cross-referencing.
clust_sel <- select_parents_by_family(
  score         = res$gebv,
  family        = family_id,          # echoed only, not used for grouping
  group_by      = "genetic_cluster",
  G             = res$G,
  n_clusters    = 8L,
  n_families    = 6L,
  n_per_family  = 3L
)
clust_sel$by_family[, c("individual", "family", "genetic_group")]
} # }
```
