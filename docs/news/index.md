# Changelog

## HapBlockR 0.3.12.9000 (development)

### New: `select_parents_by_family()` gains four opt-in ways to decide a family’s take (`family_select_mode`)

Previously a chosen family’s take was always a fixed number
(`n_per_family`). Four independent, mutually exclusive modes now decide
this, selected via the new `family_select_mode` argument (default
`"count"`, fully backward compatible — nothing changes unless you opt
in):

- `"count"` (default): unchanged, `n_per_family` lines per family.
- `"percentage"`: new `pct_per_family` (0–100\], a fixed share of *that
  family’s own* eligible size rather than a number shared across every
  family, rounded up ([`ceiling()`](https://rdrr.io/r/base/Round.html))
  so an included family never contributes zero. Still a quota, so it
  keeps its own version of the family-size eligibility rule, evaluated
  per family from its own size.
- `"sd_threshold"`: new `sd_threshold` (a signed number of SDs from the
  population mean of every eligible candidate in the run). Every family
  member clearing the bar is taken — however many that is.
- `"check_relative"`: new `check_id` (a candidate already in `score`) or
  `check_value` (a fixed benchmark), plus `check_margin_pct` (e.g. `10`
  for “at least 10% better than the check” — lines clearing it by more,
  e.g. 20% or 30%, are included too, not excluded).

`"sd_threshold"` and `"check_relative"` are merit *bars*, not quotas: a
chosen family can legitimately contribute zero members if nobody clears
the bar, so unlike `"count"`/`"percentage"` they have no size-based
pre-ranking eligibility rule at all. Instead, family selection now
auto-backfills: any ranked family contributing zero members under these
two modes is skipped (reported in the new `$zero_selected_groups` field,
with a message unless `verbose = FALSE`), and the next-ranked eligible
family is tried instead, until `n_families` worth of *contributing*
families are assembled or the ranking is exhausted.

New return fields: `family_select_mode`, `pct_per_family`,
`sd_threshold`, `check_id`, `check_value`, `check_margin_pct`,
`zero_selected_groups`. `n_per_family` is now optional at the R level
(required only when `family_select_mode = "count"`) — no change for
existing callers, who all supply it already.

### BREAKING: `select_parents_by_family()` — more rigorous family-mean correction by default, plus a mandatory family-size eligibility rule

Three statistical refinements to the shrinkage-corrected family ranking
introduced earlier in this same development version, all switched on by
default (a breaking change for anyone relying on the exact numeric
`rank_score`/`family_ranking` values from before), plus one new,
non-optional eligibility rule:

- **REML variance components (new default).** `variance_method` (new
  argument: `"reml"` default, or `"anova"`) controls how the
  between-family variance (`tau2`) and within-family/residual variance
  (`sigma2`) behind shrinkage are estimated. `"reml"` fits
  `lme4::lmer(score ~ 1 + (1|family))`, more statistically efficient
  than the previous ANOVA method-of-moments estimator under unbalanced
  family sizes — the normal case here, since uneven family size is the
  entire reason shrinkage exists in the first place. Falls back to
  `"anova"` automatically (message, not silently) whenever `lme4` isn’t
  installed, there isn’t enough family structure to fit the model, or
  the fit fails. `lme4` is a new `Suggests` dependency.
- **Finite-population selection-bias correction (new default).**
  `bias_correction` (new argument: `"order_stats"` default, or
  `"none"`). “Mean of a family’s own best `rank_k`” is a systematically
  *upward*-biased estimator of that family’s true mean, purely from
  having picked winners out of a finite pool — a different bias than the
  one shrinkage corrects for, present even when a family’s individual
  scores carry zero estimation error. `"order_stats"` removes it before
  shrinkage, via a finite-population order-statistics correction
  computed by direct numerical integration.
- **Genomic-relationship-informed shrinkage (new default).**
  `use_family_relationship` (new argument, default `TRUE`). When `G` is
  supplied and covers every eligible family, plain independent
  (i.i.d.-toward-the-grand-mean) shrinkage is replaced by a GBLUP-style
  family-effect BLUP (Henderson’s mixed-model equations over a
  family-level relationship matrix), letting genetically related
  families borrow strength from each other, not only from the grand
  mean. Falls back to the plain i.i.d. formula (message, not silently)
  whenever `G` doesn’t cover every eligible family, or the argument is
  set to `FALSE`.
- **Family-size eligibility rule (new, always on).** A family/group
  whose eligible membership does not exceed the number of lines that
  would actually be taken from it (`n_per_family`, or `rank_k` for
  uneven, named-vector quotas) is now excluded entirely before ranking,
  since there is no genuine “select the best of” decision possible for
  it — every member would be taken regardless of ranking. Reported in
  the new `$excluded_groups` return field, and messaged (unless
  `verbose = FALSE`).

New return fields: `variance_method`, `variance_method_used`,
`bias_correction`, `use_family_relationship`, `relationship_informed`,
`tau2`, `sigma2`, `excluded_groups`. `family_ranking` gains
`topk_mean_bias`/`topk_mean_corrected` columns. Set
`variance_method = "anova"`, `bias_correction = "none"`,
`use_family_relationship = FALSE` to reproduce the previous default
numeric behaviour exactly (modulo the new, non-optional family-size
exclusion, which has no toggle).

### New: an easy 0-100 dial for `select_parents_ga()`’s `merit_weight`, plus a transparent calibration helper

`merit_weight` (added earlier in this same development version — see
below) has no universal correct value: the block-coverage term and the
merit term live on different, problem-specific scales, so a raw
multiplier that works for one dataset means nothing on another. Two new
pieces close that gap:

- [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  gains an optional `merit_priority` argument (`[0, 100]`, default
  `NULL`, fully backward-compatible, mutually exclusive with an
  explicitly-supplied `merit_weight`). Instead of guessing a raw
  multiplier, state how much you care about merit vs. coverage as a
  plain percentage; the function calibrates the matching `merit_weight`
  from your own data before the GA runs. `0` is identical to
  `merit_weight = 0`; `100` scales merit’s realistic influence on the
  fitness function to match coverage’s; values between scale linearly.
  New output field `$merit_priority` (echoes the argument;
  `$merit_weight` reports whichever value was actually used, whether
  supplied directly or calibrated).
- New exported function
  [`suggest_merit_weight()`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md):
  the same calibration, callable standalone, so you can inspect the
  actual numbers (`merit_span`, `coverage_span`, `scale_factor`, and the
  four reference values they’re built from) before committing to a
  `merit_priority`, or get a literal `merit_weight` back to fine-tune by
  hand.

What “comparable scale” means, precisely: the realistic best-vs-worst
achievable spread of each term, for a founder group of your chosen
`n_founders`. The coverage reference point is not a naive per-block sum
of each block’s own maximum across all candidates, which is usually
*unreachable* by any single real group (different blocks’ maxima often
belong to different individuals) and would bias the calibration toward
an inflated ceiling — instead it is the coverage of an actual, feasible
group built by a greedy search. This is not an arbitrary heuristic: the
coverage function is monotone submodular under a fixed group-size
constraint, so greedy is mathematically guaranteed to land within a
known factor (1 - 1/e, about 63%) of the true best achievable coverage,
using a real, reachable group — not a guess. The “worst achievable”
references (for both merit and coverage) are built from a *trimmed*
low-ranked group (the bottom `n_founders` candidates after setting aside
a small, data-size-scaled number of the most extreme low values first),
so a single outlier candidate can’t single-handedly deflate the floor
and distort the estimate; trimming is skipped on small candidate pools
(under ~20), where it wouldn’t be meaningful.

Two limits are honestly documented rather than engineered around,
because they are inherent to *any* scale-matching approach, not specific
to this method: matching spread is one reasonable, explicitly-stated
definition of “comparable,” not a universal truth (matching standard
deviation instead, for instance, would give a different number); and a
single span number doesn’t capture distribution shape (an
unusually-shaped `merit_score`, e.g. strongly bimodal, may make the
dial’s practical effect feel less than perfectly linear). `merit_weight`
remains available for full manual control, and
[`suggest_merit_weight()`](https://FAkohoue.github.io/HapBlockR/reference/suggest_merit_weight.md)
lets you inspect the calibration before trusting it. The calibration
also only weighs merit against coverage — if `coancestry_weight` is also
active, its effect is held fixed rather than jointly recalibrated; use
[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)’s
sweep to explore that trade-off separately, as already recommended for
tuning `coancestry_weight` on its own.

### New: an easy 0-90 dial for `select_parents_ga()`’s `coancestry_weight`, reusing `select_parents_ocs()`’s `target_degree` convention

`coancestry_weight` has the same scale problem `merit_weight` had, but a
weighted-penalty calibration analogous to `merit_priority` isn’t the
best fit for it: unlike `merit_score` (one independent number per
individual), the coancestry term is a property of *pairs* within the
chosen set, so neither the exact-top-K trick nor a scale-matching
heuristic transfers cleanly. Instead,
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
gains a `target_degree` argument that reuses the SAME lever
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
already has, on the same `[0, 90]` scale and direction (`0` = max gain,
prioritising coverage/merit and accepting more relatedness; `90` = max
diversity, minimising relatedness) — one consistent mental model for the
trade-off across the whole package, rather than a second,
differently-shaped dial.

Mechanically, `target_degree` is not a weighted-sum calibration at all:
it interpolates a mean-relationship CEILING (a bound in the relationship
matrix’s own real units, not an abstract multiplier) between two
genuinely constructed reference groups for the actual candidate pool,
mirroring how `select_parents_ocs(engine = "optisel")` interpolates its
own kinship bound between two solved frontier extremes:

- **Gain end (`target_degree = 0`)**: the feasible, greedy-built
  high-coverage group already used for `merit_priority`’s own
  calibration (`.greedy_coverage_group()`), now extended to also account
  for `merit_weight` when merit is active in the same call — so this
  reference reflects what the search would actually converge to with no
  relatedness constraint, not a coverage-only proxy. Matches the
  `"optisel"` engine’s own “max.Merit, no kinship constraint” extreme.
- **Diversity end (`target_degree = 90`)**:
  `select_core_collection(strategy = "maximin")` run on the relationship
  matrix for the same candidate pool — a proven 2-approximation for
  maximising the minimum pairwise distance (Gonzalez 1985), built
  ignoring merit entirely. Matches the `"optisel"` engine’s own
  “min.Kin, ignoring merit” extreme.

The interpolated ceiling is enforced inside the GA fitness function by a
squared-violation penalty, scaled to dominate any possible
coverage/merit gain from crossing it (the same soft-constraint idiom
already used for the existing cardinality penalty) — so no new
penalty-strength argument is needed. `target_degree` is mutually
exclusive with an explicitly-supplied `coancestry_weight`. New output
fields `$target_degree` (echoes the argument), `$relatedness_ceiling`
(the ceiling actually enforced), and `$coancestry_weight` (echoes that
argument as actually used, `0` when `target_degree` was supplied
instead).

Two limits are honestly documented, consistent with `merit_priority`’s
own: the gain-end reference is a plain greedy hill-climb once
`merit_weight > 0` (adding a mean-based merit term breaks the strict
submodularity argument that gives the coverage-only greedy its proof),
not a formally guaranteed-near-optimal one; and
`select_core_collection`’s 2-approximation guarantee covers the
*minimum* pairwise distance in its chosen group, not directly the *mean*
pairwise relationship used here as the diversity-end reference — a
principled, reused building block, not a proof about the specific number
reported.

### BREAKING: `select_parents_by_family()` overhaul — shrinkage-corrected ranking is now the default, plus relatedness control, genetic-cluster grouping, decoupled rank_k, a second diversity heuristic, and uneven quotas

[`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)’s
family-quota shortlist strategy gains seven related improvements,
addressed together as one coherent rework rather than piecemeal:

- **BREAKING: family ranking is shrinkage-corrected by default.** A
  family’s raw top-`rank_k` mean score is an unreliable estimator of its
  true quality when the family is small — classical small-sample noise
  that quantitative genetics normally corrects for via
  shrinkage/BLUP-style family evaluation rather than trusting a raw
  sample mean regardless of how many observations it’s built from. The
  new default, `family_rank_method = "shrunk_topk_mean"`, estimates
  between-family and within-family/residual variance components from
  every eligible member’s raw `score` via the standard one-way
  random-effects ANOVA method-of-moments estimator (Searle, Casella &
  McCulloch 1992), then shrinks each family’s `topk_mean` toward the
  across-family mean by the resulting empirical-Bayes/BLUP reliability
  weight `w_f = n_f*tau2 / (n_f*tau2 + sigma2)` (Robinson 1991) — a
  family of 1 is pulled hard toward the pack; a large family’s genuine
  top-k average is trusted almost fully. This can change family rank
  order relative to previous releases whenever chosen families differ
  meaningfully in size. Pass `family_rank_method = "topk_mean"` to
  restore the exact historical, unshrunk ranking. Falls back to unshrunk
  ranking (with a message, not silently) whenever variance components
  can’t be estimated at all (e.g. every eligible family has exactly one
  member). New output columns `family_ranking$shrinkage_weight` and
  `$rank_score` (the value families were actually ranked on;
  `$topk_mean` remains the raw, unshrunk figure either way).
- **Relatedness control (`within_group_target_degree`, `G`).** Whenever
  a relationship matrix `G` is supplied, the realised mean pairwise
  relationship of the final shortlist — overall (`$mean_relationship`)
  and per selected family/group (`family_ranking$mean_relationship`) —
  is now always reported, mirroring
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
  own diagnostic, present regardless of whether anything was asked to
  optimise for it. On top of that, the new `within_group_target_degree`
  argument (`[0, 90]`, the SAME convention
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)/[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
  use) turns this into an actual dial: filling a family’s quota stops
  being pure score-order truncation and instead balances score against a
  relatedness ceiling, interpolated between that family’s own
  top-scoring feasible subset and
  `select_core_collection(strategy = "maximin")` restricted to its own
  members (a proven 2-approximation, Gonzalez 1985). Never drops a
  family’s best remaining member purely to satisfy the relatedness
  preference — unfillable slots fall back to score order.
- **`rank_k`, decoupled from `n_per_family`.** The family-ranking
  criterion’s sample size (how many of a family’s own top members its
  ranking mean is computed from) no longer has to equal the actual
  take-quota; defaults to `n_per_family` (its max, if a vector) when not
  supplied, preserving the historical behaviour exactly.
- **A second, more rigorous `diversity_method` for
  `ensure_haplotype_diversity`.** The original `"dominant_block"`
  heuristic (each individual’s single highest-value target block) is now
  one of two choices; the new default, `"coverage_gain"`, reuses
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)‘s
  own block-coverage machinery to ask whether a candidate’s addition
  actually increases total value-weighted coverage across every target
  block, given everyone already claimed so far — not merely whether one
  single block collides. Existing calls that rely on `haplotypes`’
  allele-level dominant-block collision resolution should now pass
  `diversity_method = "dominant_block"` explicitly, since `haplotypes`
  is ignored under the new default.
- **Uneven per-family/per-group quotas.** `n_per_family` may now be a
  single integer (flat quota, unchanged default) or a named integer
  vector (explicit override per chosen family/group).
- **Optional genetic-cluster grouping (`group_by`, `n_clusters`,
  `cluster_method`).** `group_by = "genetic_cluster"` (alongside the
  default, `group_by = "family"`) builds groups directly from `G` via
  hierarchical clustering (Ward’s minimum-variance linkage,
  `cluster_method = "ward.D2"` by default) on the distance matrix
  implied by `G` (the same identity
  [`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
  uses), cut into `n_clusters` groups — for programmes where pedigree
  labels aren’t a fully faithful proxy for actual relatedness. Both a
  pedigree family label and a genetic-cluster label are reported per
  selected individual whenever available (`by_family$family` /
  `by_family$genetic_group`), independent of which one actually drove
  the selection decision (`$group_by` echoes that), so a caller using
  genetic clusters for the quota logic can still see each pick’s
  pedigree family for cross-referencing, and vice versa.
- **Smaller robustness items:** `min_sel_value`/`min_sel_mode` and
  `ensure_haplotype_diversity`’s validation paths are unchanged and
  remain fully backward-compatible; all new arguments are appended after
  the original nine, so existing positional calls keep working
  unmodified.

Two limits are honestly documented, consistent with the rest of the
package’s optional dials: the shrinkage correction is a well-justified,
practical empirical-Bayes approach applied to the selection statistic
(`topk_mean`), not a single unified formal model of that statistic
itself; and `select_core_collection`’s 2-approximation guarantee (used
by both `within_group_target_degree`’s diversity end and genetic-cluster
grouping’s distance conversion) covers the *minimum* pairwise distance
in its chosen group, not directly the *mean* pairwise relationship
reported as the diversity-end reference.

### BREAKING: `select_parents_ocs()` engine names — `"optisel"` is now real optiSel-based OCS; the old `"optisel"` behaviour is now `"simplemating"`

`select_parents_ocs(engine = "optisel")` never actually called optiSel’s
own solver: it wrapped
[`SimpleMating::planCross()`](https://rdrr.io/pkg/SimpleMating/man/planCross.html)/`selectCrosses()`,
with optiSel present only as one of SimpleMating’s own transitive
dependencies. The name was a misnomer from early development (an earlier
hand-written
[`optiSel::candes()`](https://rdrr.io/pkg/optiSel/man/candes.html)/`opticont()`
attempt failed on argument-name mismatches and was replaced by the
SimpleMating-based wrapper, but the engine kept its original name — see
`R/ocs.R`’s file header for the full history). This is now fixed
properly, not just renamed away from:

- **`engine = "optisel"` is a new implementation** that genuinely calls
  optiSel’s own solver (Wellmann 2019):
  [`optiSel::candes()`](https://rdrr.io/pkg/optiSel/man/candes.html)
  builds the candidate-description object from `merit` and `G`;
  [`optiSel::opticont()`](https://rdrr.io/pkg/optiSel/man/opticont.html)
  solves the actual continuous Optimal Contribution Selection problem
  (Meuwissen 1997) — optimal per-candidate contributions subject to an
  explicit upper bound on the next generation’s mean kinship;
  [`optiSel::noffspring()`](https://rdrr.io/pkg/optiSel/man/noffspring.html)
  converts contributions into expected offspring counts;
  [`optiSel::matings()`](https://rdrr.io/pkg/optiSel/man/matings.html)
  solves the discrete Sire x Dam assignment minimising mean offspring
  kinship. `target_degree` (same `[0, 90]` direction convention as
  always: 0 = max-gain end, 90 = max-diversity end) is mapped onto a
  mean-kinship ceiling by solving both ends of the actual gain/diversity
  frontier for your candidate set via `opticont()` and interpolating
  between them. Requires only the `optiSel` package — no SimpleMating
  dependency, no external binary.
- **The previous `"optisel"` behaviour (SimpleMating-based cross
  prediction/selection) is now `engine = "simplemating"`.** Its
  algorithm, `target_degree` mechanism (quantile-of-observed-relatedness
  cutoff via `selectCrosses()`), and all other behaviour are completely
  unchanged – only the engine name changed.
- **There is no backward-compatible alias.** Code that passed
  `engine = "optisel"` expecting the old SimpleMating-based behaviour
  will now get the new optiSel-solver-based engine instead — a different
  algorithm with a different `target_degree` mechanism (though the same
  direction convention). If you have existing code, notebooks, or
  scripts pinning `engine = "optisel"`, update them to
  `engine = "simplemating"` to keep their exact previous behaviour, or
  leave them on `"optisel"` deliberately to switch to real OCS
  (recommended if you were only using `"optisel"` because it required no
  external binary, not because you specifically wanted SimpleMating’s
  greedy cross-selection algorithm).
- **`engine = "auto"`’s fallback (when no usable `alphamate_exe` is
  supplied) now resolves to `"optisel"`** (the new true-OCS engine)
  instead of what is now `"simplemating"` — `"auto"` never silently
  substitutes discrete cross prediction/selection for actual Optimal
  Contribution Selection. This is also a behaviour change for any code
  relying on `engine = "auto"`’s previous fallback target.
- [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
  now offers three engines total: `"alphamate"` and `"optisel"` both
  solve true OCS (via genuinely different solvers); `"simplemating"`
  does discrete greedy cross prediction/selection, a different algorithm
  class entirely. See
  [`?select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)’s
  rewritten “Engine differences” section.
- Internal functions this affects, if you called them directly (not
  recommended, but mentioned for completeness):
  `.run_simplemating_ocs()` is unchanged in behaviour (only reachable
  via `engine = "simplemating"` now); the new `.run_optisel_ocs()` backs
  `engine = "optisel"`.
- `ga_vs_ts_simulation(mating_scheme = "ocs")` and the
  `scripts/CIAT_Peru_parent_selection.R` example script have both been
  pinned explicitly to `engine = "simplemating"` (rather than relying on
  `"auto"`’s new fallback target) so their existing, already-tuned
  behaviour and performance profile are unaffected by this rename.

### Fixed: `select_parents_ocs(engine = "optisel")`’s `target_degree` direction was inverted relative to AlphaMate

*(Note: at the time of the fix described below, `engine = "optisel"`
meant what is now `engine = "simplemating"` — see the breaking-rename
entry above. The bug and fix described here apply to that engine, now
named `"simplemating"`.)*

`target_degree`’s documented convention was, and is meant to be,
identical to AlphaMate’s own `TargetDegree` (Kinghorn’s frontier-degree
concept): `0` is the max-gain end of the frontier (no diversity
restriction, prioritise merit, accept more relatedness), `90` is the
max-diversity end (prioritise minimising relatedness). The
`engine = "optisel"` implementation’s own `target_degree` -\>
`culling.pairwise.k` quantile mapping had this backwards –
`target_degree = 0` was keeping only the *least*-related candidates (the
max-diversity behaviour) and `target_degree = 90` was keeping
*essentially all* candidates (the max-gain behaviour), the opposite of
AlphaMate’s real, published convention (confirmed against the standard
Kinghorn frontier-degree diagram). This was a genuine functional bug,
not just a documentation error: under the previous code, the same
`target_degree` value meant opposite things depending on which engine
you used — `engine = "alphamate"` was unaffected (it passes
`target_degree` straight through to the real AlphaMate executable, which
interprets it correctly), but `engine = "optisel"` was inverted relative
to it. Fixed by flipping the quantile formula from `target_degree / 90`
to `(90 - target_degree) / 90`. If you have existing code or a pinned
`target_degree` value tuned under `engine = "optisel"` against the old
(backwards) behaviour, re-tune it — what used to be a diversity-
prioritising value is now a merit-prioritising one, and vice versa. All
documentation
([`?select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
README, the *From Local GEBV to a Crossing Decision* vignette’s Sections
9-10) has been corrected to match.

### New: two additional parent-selection strategies

- [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  gains an optional `merit_weight` argument (default `0`, fully
  backward-compatible). When positive, it adds a weighted
  whole-genome-merit term directly into the GA’s fitness function,
  alongside the existing block-coverage term (and the coancestry
  penalty, if also active) — the same additive pattern already used for
  `coancestry_weight`. Previously, `merit_score` could only act as a
  hard pre-search eligibility floor (`min_sel_value`), which does not
  distinguish a candidate that barely clears the floor from one that
  clears it by a wide margin. With `merit_weight > 0`, whole-genome
  merit keeps pulling on the search throughout, not only at its entry
  gate — useful specifically when per-block/local-GEBV estimates carry
  real estimation noise and a coverage-only search risks selecting on
  that noise rather than genuine signal. New outputs: `$mean_merit`
  (top-level and per-replicate in `$stability$mean_merit`) and
  `$merit_weight`. See
  [`?select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
  new “Merit-weighted fitness (GA+TS hybrid, optional)” section.
- New function
  [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md):
  a third parent-shortlist strategy alongside
  [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
  and
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
  matching the shape a programme’s shortlist decision often already
  takes – pick the `n_families` best-performing families first (ranked
  by the mean of each family’s own top `n_per_family` members), then the
  `n_per_family` best lines within each chosen family. Unlike
  [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)‘s
  single population-wide ranking, this enforces family balance directly
  in the selection rule rather than checking it as a diagnostic
  afterward. An optional `ensure_haplotype_diversity = TRUE` mode (with
  `value_matrix`, and optionally `haplotypes` — the direct, unmodified
  output of
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md),
  auto-reshaped internally, no manual matrix construction needed — for
  exact allele-level resolution) adjusts within-family picks so that
  representatives of different chosen families are less likely to carry
  their best value at the same target haplotype block — so a cross
  between two chosen families’ representatives is more likely to combine
  two genuinely different favourable haplotypes rather than duplicating
  one. Exported; see
  [`?select_parents_by_family`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md).

### Fixed: two CI/CRAN-check regressions introduced by the rename

- `configure` and `configure.win` lost their executable permission bit
  during the bulk `LDxBlocks` -\> `HapBlockR` text substitution
  (`sed -i` rewrites a file via a temp-file swap that does not always
  preserve the original mode). R requires `configure` to be executable
  on Unix-like systems, so this broke package installation on
  Linux/macOS (Windows is unaffected — it does not enforce POSIX
  permission bits) with
  `ERROR: 'configure' exists but is not executable`. Fixed by restoring
  the executable bit (`git update-index --chmod=+x`, since
  `core.fileMode` is commonly `false` on Windows checkouts and a plain
  `chmod` + `git add` will not pick up the change there).
- `tests/testthat/test-ocs.R`’s “does NOT warn off Windows for a genuine
  non-PE exe” test used `expect_no_warning()` around a call that
  actually attempts to execute a fake AlphaMate binary consisting of
  just 4 ELF magic bytes — not a complete, runnable executable.
  [`system2()`](https://rdrr.io/r/base/system2.html) correctly fails to
  launch it and R itself emits an unrelated “had status 126” warning,
  which `expect_no_warning()` (wrongly) failed on; this test had never
  actually run on Linux/macOS CI before (blocked by the `configure`
  issue above) so the bug was latent. Fixed by capturing warnings
  directly and asserting only that the platform-mismatch (“Windows
  binary”) warning is absent, tolerating the incidental status-126 one.

### BREAKING: No AlphaMate binary is bundled with or auto-detected by HapBlockR anymore

`R CMD check` flagged the bundled `inst/extdata/AlphaMate.exe` and
`inst/extdata/libiomp5md.dll` with a WARNING (“checking for executable
files”) and a NOTE (“checking if this is a source package” — apparent
object files/libraries). Neither is fixable while still shipping the
binaries: CRAN’s own manual (Writing R Extensions) states that source
packages containing binary executables are rejected outright, regardless
of any `BinaryFiles` declaration. AlphaGenes (the AlphaMate maintainers)
also does not currently publish pre-built binaries for
[AlphaGenes/AlphaMate](https://github.com/AlphaGenes/AlphaMate) (Fortran
source, MIT licensed, no GitHub Releases), so there was no stable URL an
`install_alphamate()`-style downloader could have fetched from either.

Both files have been removed from `inst/extdata/`, and
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)’s
auto-detection logic (previously: use the bundled Windows binary
automatically when `alphamate_exe = NULL` and running on Windows) has
been removed along with them. `alphamate_exe` must now always be
supplied explicitly to use `engine = "alphamate"` (build AlphaMate
yourself from its GitHub source, or otherwise obtain/be entitled to use
a working executable); `engine = "auto"` now falls back to `"optisel"`
whenever `alphamate_exe` is `NULL` or does not point to an existing
file, on every platform (the previous Windows-specific bundled-binary
branch of that logic is gone). `engine = "optisel"` itself is unaffected
— it never required an external binary. Updated: `R/ocs.R` (roxygen docs
for
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
in particular “Providing the AlphaMate executable”; the
`@param alphamate_exe` and `engine = "auto"` sections; the auto-detect
code block; the [`Sys.chmod()`](https://rdrr.io/r/base/files2.html)
comment in `.run_alphamate_ocs()`), the README’s “Providing the
AlphaMate executable” section, and `tests/testthat/test-ocs.R` (the
bundled-exe auto-detection tests are replaced with
`alphamate_exe = NULL`/nonexistent-path tests that no longer assume a
bundled file might be present).

### BREAKING: Package renamed from `LDxBlocks` to `HapBlockR`

The package (and GitHub repository) has been renamed from `LDxBlocks` to
`HapBlockR` to better reflect its full scope: LD block detection plus
haplotype analysis, genomic prediction, and breeding-decision support.
The `Title` field is now “Genome-Wide LD Block Detection, Haplotype
Analysis, Genomic Prediction, and Breeding Decision Support”.

This is a mechanical rename with no change in function behaviour, but it
touches every user-facing identifier tied to the old package name:

- [`library(LDxBlocks)`](https://github.com/FAkohoue/LDxBlocks) -\>
  [`library(HapBlockR)`](https://github.com/FAkohoue/HapBlockR);
  `LDxBlocks::`/`LDxBlocks:::` -\> `HapBlockR::`/`HapBlockR:::`;
  `data(..., package = "LDxBlocks")` -\>
  `data(..., package = "HapBlockR")`;
  `system.file(..., package = "LDxBlocks")` -\>
  `system.file(..., package = "HapBlockR")`.
- S3 classes renamed: `LDxBlocks_backend` -\> `HapBlockR_backend`,
  `LDxBlocks_cv` -\> `HapBlockR_cv`, `LDxBlocks_decay` -\>
  `HapBlockR_decay`, `LDxBlocks_diplotype` -\> `HapBlockR_diplotype`,
  `LDxBlocks_effect_concordance` -\> `HapBlockR_effect_concordance`,
  `LDxBlocks_epistasis` -\> `HapBlockR_epistasis`,
  `LDxBlocks_block_epistasis` -\> `HapBlockR_block_epistasis`,
  `LDxBlocks_haplotype_assoc` -\> `HapBlockR_haplotype_assoc`. Code that
  checks `inherits(x, "LDxBlocks_backend")` etc. must be updated.
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  methods dispatch on the new class names. Function names
  ([`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md),
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md),
  etc.) are unchanged.
- Compiled-code registration renamed in lockstep: Rcpp routine symbols
  `_LDxBlocks_*` -\> `_HapBlockR_*` and `R_init_LDxBlocks` -\>
  `R_init_HapBlockR` (`src/RcppExports.cpp`, `R/RcppExports.R`). The
  package must be recompiled (`R CMD INSTALL` /
  [`devtools::install()`](https://devtools.r-lib.org/reference/install.html))
  – a previously compiled `.dll`/`.so` built under the old name will not
  load. `NAMESPACE`’s `useDynLib(LDxBlocks, ...)` -\>
  `useDynLib(HapBlockR, ...)`.
- Internal `bigmemory`-backed cache file prefixes renamed:
  `ldxblocks_bm*` -\> `hapblockr_bm*` (`.bin`/`.desc`/`_snpinfo.rds`/
  `_sampleids.rds`/`_params.rds` sidecar files under `bigmemory_path`,
  used by
  [`read_geno_bigmemory()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno_bigmemory.md),
  `run_ldx_pipeline(use_bigmemory = TRUE)`, and phased-VCF caching). Any
  existing on-disk cache built under the old prefix will not be
  reattached automatically — delete the old `ldxblocks_bm*` files (or
  point `bigmemory_path` at a fresh directory) so the cache rebuilds
  under the new prefix.
- `URL`/`BugReports` in `DESCRIPTION` now point at
  `github.com/FAkohoue/HapBlockR` (assumes the GitHub repository itself
  is also renamed to `HapBlockR`).
- All package source (`R/`, `src/`, `tests/testthat/`, `vignettes/`,
  `scripts/`, `man/*.Rd` content), `README.md`, and this changelog were
  updated for the new name. `NAMESPACE` and `man/*.Rd` were hand-edited
  in this pass; **run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  after pulling this change** to regenerate `NAMESPACE`/`man/`
  authoritatively — it will create correctly-named files such as
  `man/HapBlockR-package.Rd` and `man/print.HapBlockR_backend.Rd`,
  leaving the old `man/LDxBlocks-*.Rd` files behind to be deleted
  manually.
- Files referenced by path from code or docs (`R/breeder_guide.R`’s
  `system.file("extdata", "HapBlockR_Breeder_Guide.pdf", ...)`, and
  README’s `man/figures/HapBlockR_schematic.svg` image) needed matching
  on-disk files, so `HapBlockR_Breeder_Guide.pdf` (top level and under
  `inst/extdata/`) and `man/figures/HapBlockR_schematic.png`/`.svg` were
  added as copies alongside the old `LDxBlocks_*`-named originals (this
  folder does not allow renaming/deleting existing files). Renamed
  copies of `LDxBlocks_formula_audit.md`/`.html`,
  `LDxBlocks_vs_HapSelect_comparison.md`/`.html`, and `LDxBlocks.Rproj`
  were added the same way for consistency, though nothing resolves those
  by path. The old `LDxBlocks_*`-named originals are now redundant and
  can be deleted once you confirm the new copies open correctly.

### Added: `run_ldx_pipeline(phase = TRUE)` now performs genuine Beagle-based imputation, exposed as `beagle_geno_matrix`

Previously, when `phase = TRUE`, the VCF handed to Beagle was built from
`geno_mat` *after* Step 4.5c’s `mean_rounded`/`mode` imputation had
already filled in every missing call. Since Beagle’s HMM only imputes
calls it actually sees as missing, this meant Beagle only ever re-phased
already- decided dosages — its real imputation capability (local
haplotype-cluster model, generally far more accurate than a column
mean/mode) was never exercised, despite
[`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md)
itself being fully capable of it.

[`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)
now snapshots `geno_mat`’s real missingness pattern
(`geno_mat_preimpute`) right before Step 4.5c’s imputation runs, and —
only when `phase = TRUE` — writes that snapshot (not the imputed matrix)
to the cleaned VCF fed to Beagle. Beagle therefore now performs genuine
HMM-based imputation of missing calls in addition to phasing. The
resulting individuals x SNPs dosage matrix is returned as a new field,
`beagle_geno_matrix` (`NULL` unless `phase = TRUE`). The existing
`geno_matrix` field is unchanged and still always holds the
`mean_rounded`/`mode`-imputed matrix used for LD block detection, so
this is purely additive — no backward-compatibility break. A
verbose-gated message now also flags any residual `NA` left in
`beagle_geno_matrix` (e.g. sites Beagle could not resolve without a
reference panel).

`tests/testthat/test-phasing.R`: the existing `phase = TRUE` pipeline
test’s expected-fields check now includes `beagle_geno_matrix`; added a
new test confirming `beagle_geno_matrix` is `NULL` when `phase = FALSE`,
and a new JAR-gated test that injects real missingness into the shared
30-SNP fixture and confirms `beagle_geno_matrix` has no residual `NA`
and differs from `geno_matrix`’s mean/mode fallback at
originally-missing cells (proving genuine Beagle-driven imputation, not
a copy of the fallback).

### Changed: `usefulness_criterion(variance_model = "simplemating")` no longer errors on heterozygous genotypes by default

Previously, any heterozygous (dosage = 1) call anywhere in `geno_matrix`
made this mode error out immediately
([`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)
requires strictly 0/2/NA, fully homozygous DH/RIL-style calls). Real
breeding programmes’ genotypes are usually supplied as ordinary 0/1/2
dosage and realistically carry at least some residual heterozygosity
even in advanced inbred material, so this made the mode impractical to
use without manually pre-processing `geno_matrix` first.

Added `het_to_na` (default `TRUE`): heterozygous calls are now
automatically treated as missing (set to `NA`) before calling
[`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html),
and
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
proceeds normally. This is a defensible statistical choice, not a
fabricated one — a heterozygous cell is marked as “no confident
homozygous call here”, never rounded to 0 or 2 (which would invent a
specific allele call the data does not support), and `NA` was already a
valid value for this argument even before `het_to_na` existed
([`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html)
already receives and handles missing calls the normal way). When
`verbose = TRUE`, a message reports exactly how many calls (count and %)
were converted. Values outside 0/1/2/NA still error unconditionally,
since there’s no defensible automatic interpretation for those. Set
`het_to_na = FALSE` to restore the previous strict behaviour.
`tests/testthat/test-genomic-mating.R` updated: the old “errors loudly
on heterozygous” test now verifies the new default succeeds (with the
expected message), plus a new test confirms `het_to_na = FALSE` still
errors, and a new test confirms genuinely invalid dosage values (e.g. 9)
still error regardless. `scripts/CIAT_Peru_parent_selection.R`’s stage
13b updated to list `het_to_na = TRUE` explicitly and to reflect that
heterozygosity is no longer the likely cause if that cross-check’s
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) safety net ever
fires.

### Updated: `scripts/CIAT_Peru_parent_selection.R` — variance_model cross-checks (stages 13b/13c)

Stage 13’s cross ranking used `variance_model = "block_independent"`
only. Added two cross-checks against the same OCS-selected parents’
candidate crosses: stage 13b runs `variance_model = "simplemating"`
(backsolving `snp_effects` from the whole-genome GEBV via
[`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md),
and supplying `ld_matrix` — via
[`compute_r2()`](https://FAkohoue.github.io/HapBlockR/reference/compute_r2.md)
on the target-block SNPs — since this programme has no real genetic map;
wrapped in [`tryCatch()`](https://rdrr.io/r/base/conditions.html) since
`simplemating` requires strictly homozygous genotypes and this
programme’s data may carry residual heterozygosity, in which case the
cross-check is skipped with a clear message rather than aborting the
script), and reports Spearman rank agreement against stage 13’s ranking.
Stage 13c documents (but does not run) `variance_model = "linked"`: the
new `ld_matrix` fallback above removes the genetic-map barrier, but
“linked” still requires phased haplotypes, which this programme’s
unphased data source does not provide – a ready-to-use call is included
in a comment for once phased data exists.

### Added: `usefulness_criterion(variance_model = "linked")` now accepts `ld_matrix` as an alternative to `genetic_map`

Previously “linked” unconditionally required a real genetic map
(`genetic_map`: SNP, CHR, cM), with no fallback — unlike “simplemating”,
which already accepted an `ld_matrix` (SNP x SNP LD, used as a `1 - LD`
recombination- fraction proxy) when no map was available. Added the same
fallback to “linked”: when `genetic_map` is `NULL`, target blocks are
instead ordered by physical position (CHR, start_bp; new internal
`.block_physical_positions()`) and the recombination fraction between
each pair of adjacent blocks is proxied as `0.5 * (1 - mean(LD))` (mean
pairwise SNP r-squared between the two blocks’ member SNPs; new internal
`.adjacent_r_from_ld()`) — the same `[0, 0.5]` range and boundary
behaviour as Haldane’s mapping function, hand- verified against exact
closed-form values (LD = 1 -\> r = 0; LD = 0 -\> r = 0.5; a known
intermediate LD value) in `tests/testthat/test-genomic-mating.R` section
6c, plus end-to-end
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
tests mirroring the existing `genetic_map`-based ones.

**Important caveat, documented in the roxygen and NEWS here for
visibility:** this is a monotonic PROXY, not a validated
genetic-distance estimator — LD reflects population-level historical
recombination, drift, and selection, not necessarily the true
recombination fraction for a specific cross. Prefer `genetic_map` when
you have one; only fall back to `ld_matrix` when you don’t. Also,
`ld_matrix` only relaxes the genetic-map requirement — “linked” (like
“phased”) still unconditionally requires phased `haplotypes` either way,
since population-level LD cannot substitute for knowing which alleles a
specific individual’s two chromosomes actually carry.

### Fixed: `select_parents_ocs(engine = "alphamate")` — misleading “Windows binary” warning fired for genuine native builds too

The platform warning (“AlphaMate.exe is a Windows binary…”) previously
fired any time `engine = "alphamate"` ran off Windows, regardless of
what `alphamate_exe` actually pointed at — including a real, working
native Linux/macOS AlphaMate build supplied by the user, right after it
had already launched successfully. Fixed by adding an internal
`.is_windows_pe_exe()` check (reads the file’s first two bytes and
compares against the “MZ” DOS/PE header magic number all Windows
executables start with) and gating the warning on that instead of on
`.Platform$OS.type` alone. A confirmed Windows PE binary off Windows
still warns as before; a genuine native build no longer does.
`tests/testthat/test-ocs.R`’s two existing platform-warning tests were
updated to use a fake exe with real “MZ” bytes (previously an empty
file, which happened to still be treated as “presumed Windows” under the
old OS-only check), and a new test confirms a non-PE fake exe (ELF magic
bytes) produces no warning at all.

### Fixed: `select_parents_ocs(engine = "alphamate")` — contributors table’s ID column never matched, silently returning `contributors$id == NULL`

`.run_alphamate_ocs()` parses AlphaMate’s own
`ContributorsModeOptTarget1.txt` output and renames whichever column
holds the individual ID to `id`, by checking for one of
`c("Parent", "ID", "Indiv", "AlphaMateID")` (case-sensitive exact
match). AlphaMate’s actual, documented header for this file is
`Id Gender SelCriterion AvgCoancestryA AvgCoancestryC Contribution nContribution`
– capital I, lowercase d — which never matched any of the four
candidates. The rename was silently skipped, so `contributors` came back
with its original AlphaMate column names and no `id` column at all;
`ocs_res$contributors$id` then silently evaluated to `NULL` even on a
fully successful AlphaMate run with a valid mating plan. Downstream code
that prefers `contributors$id` over deriving parent IDs from
`mating_plan` (as the CIAT Peru script does) would then error out on a
real, valid mating plan. Fixed by adding `"Id"` to the recognised column
names (checked first, since it’s the actual AlphaMate output column; the
other three are kept as defensive fallbacks for other builds/versions).

### Fixed: `scripts/CIAT_Peru_parent_selection.R` — final mating-plan merge dropped most UC/rank values

The final-outputs step merged `ocs_res$mating_plan` (parent1/parent2
direction assigned by optiSel/SimpleMating’s own `selectCrosses()`)
against `uc_res` (parent1/parent2 direction assigned by
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)’s
`combn(unique(ocs_parents), 2)`, which depends only on each ID’s
position in `ocs_parents`) using a direct
`merge(by = c("parent1", "parent2"))`. Because the two sources can
disagree on which parent is listed first for the same unordered pair,
this silently dropped any row where the direction didn’t match — in one
real run, 55 of 60 mating-plan rows ended up with `NA` for
`UC`/`rank`/`predicted_variance` in `mating_plan_next_cycle.csv`, the
programme’s actual crossing-decision deliverable. Fixed by merging on an
order-independent pair key
(`paste(pmin(parent1, parent2), pmax(parent1, parent2))`) instead, while
keeping `ocs_res$mating_plan`’s own parent1/parent2 labels in the
output.

### Updated: `scripts/CIAT_Peru_parent_selection.R` — strategic decision support (stage 16) added

Wires in every function added this cycle: stage 10 now plots both the
genome-wide (`G`) and target-block feature-space (`feature_matrix`)
diversity views side by side, not just the former. New stage 16 (before
final outputs, now stage 17) runs
[`cluster_selection_groups()`](https://FAkohoue.github.io/HapBlockR/reference/cluster_selection_groups.md)
in genome-wide diversity space with BOTH `method = "hierarchical"` and
`method = "kmeans"` reported side by side, cross-tabulating all four
selection strategies used earlier in the script (GA, TS, OCS,
core-collection) against the resulting genetic clusters at once — the
direct, numeric version of the programme’s recurring
population-improvement-vs.-fast-release strategic question (“does this
cycle’s selection spread across every family, or concentrate on a
few?”), rather than inferring it from the stage-10 scatterplot. Each
clustering’s
[`plot_selection_clusters()`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md)
output is saved as a high-resolution PDF. Also fixed the file header’s
pipeline-stage list, which had drifted out of sync with the actual
in-code stage numbers since an earlier session (it still listed only
stages 1-12 under old, since-renumbered labels; now accurately lists all
17 stages).

### Changed: `cluster_selection_groups()` generalised from GA/TS-only to an arbitrary number of named strategies

`ga_selected`/`ts_selected` replaced by a single `groups` argument: a
named list of ID vectors,
e.g. `list(GA = ga$selected, TS = ts$selected, OCS = ocs_parents, Core = core_res$selected)`.
The old `Neither`/`TS-selected`/ `GA-selected`/`Both` mutually-exclusive
category (which cannot generalise past 2 groups without a combinatorial
explosion of categories) is replaced by per-strategy
`n_<name>`/`prop_<name>` columns in `table`, computed independently per
strategy so an individual selected by more than one method (e.g. both GA
and OCS) is correctly counted under each. Not a released API yet (this
function was added earlier in this same development cycle), so no
deprecation path — just an update.

### Fixed: `R CMD check` WARNINGs/NOTEs from `cluster_selection_groups()`’s roxygen docs

Three raw, unescaped `%` characters in the `@param G`/`@description`
prose (“~26% cumulative variance”, “default 95%”, “own % variance”) were
parsed by Rd as LaTeX-style comment markers, silently swallowing the
rest of each line — including closing braces — and cascading into “Lost
braces”, “unexpected section header”, and “unexpected END_OF_INPUT”
across the rest of `cluster_selection_groups.Rd`, which in turn caused
the “checking whether package can be installed” and “checking Rd
sections” WARNINGs (Rd processing runs during install) and the “Rd files
without ” NOTE (the real `@description` content got misplaced by the
cascade). Escaped all three as `\%`, matching the convention already
used everywhere else in the package. Also added `"Cluster"` to this
file’s
[`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
allowlist, fixing the separate “no visible binding for global variable
‘Cluster’” NOTE from
[`plot_selection_clusters()`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md)’s
`ggplot2::aes(colour = Cluster)` (bare column name inside `aes()`, the
same pattern already allowlisted for `PC1`/`PC2`/`group` in
[`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)).
Re-run
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
to regenerate `man/cluster_selection_groups.Rd` from the fixed source
before re-running `R CMD check`.

### Added: `cluster_selection_groups()` and `plot_selection_clusters()` — PC-variance-threshold clustering

[`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)
only ever visualises PC1/PC2, which can be a small fraction of total
variance in a genetically complex population (e.g. ~26% cumulative on
PC1+PC2 alone, with 3 visibly distinct clusters in a real CIAT panel).
`cluster_selection_groups(G, feature_matrix, ga_selected, ts_selected, variance_threshold = 0.95, method = c("hierarchical", "kmeans"), n_clusters)`
instead retains as many leading PCs as needed to reach a
cumulative-variance threshold (properly variance-scaled – eigenvectors x
`sqrt(eigenvalue)` for the `G` path, `prcomp()$x` for the
`feature_matrix` path — so Euclidean distance across many retained PCs
is not dominated equally by near-noise PCs), clusters on that subspace
(hierarchical/Ward or k-means, both available as explicit `method`
choices), and cross-tabulates GA-selected/TS-selected/Both/Neither
membership against cluster membership into a `table` — turning “does
this method represent every genetically distinct sub-group?” into counts
and proportions per cluster rather than a visual read of a 2-axis
scatter.
[`plot_selection_clusters()`](https://FAkohoue.github.io/HapBlockR/reference/plot_selection_clusters.md)
is the companion PC1/PC2 plot of the retained subspace, coloured by
cluster (not selection group); it follows this package’s
`theme_classic()` + explicit readable axis-text-size + high- resolution
(`dpi = 300`) PDF export convention (see `save_path`), matching the
existing style already used for GRM PCA/scree plots in
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)’s
plotting output.

### Added: `plot_parent_selection_pca(feature_matrix = ...)` — plot in the space GA actually searched

Previously `plot_parent_selection_pca(G, ga_selected, ts_selected)` only
supported eigen-decomposing a relationship matrix `G` — typically the
genome-wide GRM, which is a useful diversity cross-check but is *not*
the space
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
actually optimises coverage over (its target- block
local-GEBV/haplotype-value `value_matrix`). Added an alternative
`feature_matrix` argument: when supplied (instead of `G`), PCA is
computed via [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) on that
individuals x target-blocks matrix directly, so the plot answers “does
GA’s selection spread out in the space it was asked to spread out in?”
rather than only the genome-wide cross-check. `G` remains the default;
existing calls are unaffected.

### Fixed: `select_parents_ocs(engine = "auto")` tried to launch the bundled AlphaMate.exe off Windows

`have_alphamate` only checked `file.exists(alphamate_exe)`, not whether
the resolved executable was actually launchable on the current platform.
On Linux/macOS, the auto-detected bundled `AlphaMate.exe` (a Windows PE
binary) satisfies [`file.exists()`](https://rdrr.io/r/base/files.html),
so `engine = "auto"` picked `"alphamate"` anyway and the run failed with
a shell-level `Exec format error` — not a catchable R condition, so
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
could not fall back or give a useful error itself (only the pre-existing
“AlphaMate.exe is a Windows binary”
[`warning()`](https://rdrr.io/r/base/warning.html) hinted at the real
cause). Fixed: `"auto"` now only selects `"alphamate"` for the
*auto-detected bundled* `.exe` when `.Platform$OS.type == "windows"`;
off Windows it falls back to `"optisel"`, matching what
[`?select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
already documented. A user-supplied `alphamate_exe` (a native build, a
Wine wrapper script, etc.) is still always honoured by `"auto"`
regardless of platform.

### Fixed: CI — leftover unresolved git merge-conflict markers broke package parsing

`R CMD build` failed on every platform
(`Error in parse(...) : unexpected input`) because a straggler
`<<<<<<<`/`=======`/`>>>>>>>` conflict block from an earlier merge was
left committed in `NAMESPACE`, `R/forward_simulation.R`,
`R/haplotypes.R`, `R/parent_selection.R`, `tests/testthat/helper.R`,
`tests/testthat/test-forward-simulation.R`, and
`tests/testthat/test-parent-selection.R`. Resolved by keeping the
intended (post-OCS/UC-extension) side of each conflict throughout.

### Fixed: CI — macOS job failed to build SimpleMating (rgl / missing XQuartz)

`rgl`, pulled in transitively via SimpleMating’s dependency tree,
[`dyn.load()`](https://rdrr.io/r/base/dynload.html)s `libGLU.1.dylib` at
package install/lazy-load time. GitHub’s macOS runner images do not ship
XQuartz (which provides `libGLU.1.dylib`) by default, so
`Install R dependencies` failed building SimpleMating from source before
any HapBlockR code was touched — reproducing even on `macos-15`, not
just `macos-latest`/Tahoe. Fixed with two changes: (1) pinned the macOS
job to `macos-15` instead of `macos-latest`, since `macos-latest` now
resolves to macOS 26 (“Tahoe”), which dropped OpenGL support outright
(no XQuartz install can fix that); (2) added a
`brew install --cask xquartz` step before “Install R dependencies” on
macOS, so `libGLU.1.dylib` exists on disk where rgl already looks for
it.

### Updated: `scripts/CIAT_Peru_parent_selection.R` extended to the full mate-allocation layer

Previously stopped at founder-SET selection
([`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)/
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
stages 1-10) with forward simulation explicitly out of scope. Added
stage 5
([`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md),
a k-fold predictive-ability check on the haplotype-GBLUP decomposition,
run before anything is built on top of it) and a new “mate allocation”
block (stages 11-15):
[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)
(merit-vs-diversity frontier, to inform `target_degree` deliberately
rather than guessing it),
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
(true Optimal Contribution Selection — contributions plus an actual
mating plan, `engine = "auto"`, the direct in-R equivalent of this
programme’s existing manual AlphaMate workflow),
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
(`variance_model = "block_independent"`, cross ranking among the
OCS-selected parents),
[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)
(exact ILP optimality-gap check of the OCS plan), and
[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
(a diversity-first cross-check against the merit-aware selections).
Final output is now two CSVs: the founder shortlist (unchanged) and a
new `mating_plan_next_cycle.csv` (OCS mating plan merged with its UC
ranking) — the latter is the actual actionable crossing list, not just a
parent shortlist. Forward simulation
([`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md))
and dominance GBLUP remain explicitly out of scope, each with a stated
reason in the file header (unphased input data; an already-additive
GCA-based index, respectively) rather than a silent omission. Every new
function call follows the script’s existing convention of listing every
argument explicitly, including ones left at their package default.

### Fixed: SimpleMating buried in the optional-dependencies list, unlike genomicSimulation

`SimpleMating` is, like `genomicSimulation`, a GitHub-only remote
resolved via `DESCRIPTION`’s `Remotes:` field and subject to the same
nested- dependency GitHub API timeout risk on a full HapBlockR install —
but unlike genomicSimulation, it only had a terse
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html)/
[`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html)
pair buried inside a code comment in the optional-dependencies block,
with no dedicated section, rationale, or version-check guidance.
README.md’s Installation section now gives it the same treatment as
genomicSimulation: a dedicated “Step 3 (optional)” between
genomicSimulation (Step 2) and the GitHub-token fallback (now Step 4),
covering standalone install, a version/export check
(`getNamespaceExports("SimpleMating")` for `planCross`/`selectCrosses`/
`getUsefA`), the `force = TRUE` reinstall path, and when it’s safe to
skip entirely (`variance_model = "linked"`, `engine = "alphamate"`).
Steps 4-7 renumbered accordingly.

### Fixed: `R CMD check` Rd cross-reference warning in `ga_vs_ts_simulation.Rd`

`@param blocks`’s roxygen text linked to `\link{run_and_tune_Big_LD}`, a
function name that has never existed in this package (the real functions
are
[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
and
[`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md))
— an invented name from an earlier draft of the `blocks` parameter
description that was never caught because no `R CMD check` had been run
against it. Corrected to `\link{tune_LD_params}` in
`R/forward_simulation.R`. A repo-wide scan of every other in-package
`\link{}` target in `R/*.R` turned up no further invented names (the
false positives it did flag — `combn`, `browseURL`, `system.file`, and
the `ldx_*` example-dataset names — all resolve fine; they’re base-R
functions and separately-documented datasets, not missing package
objects).

### Fixed: `select_parents_ocs(engine = "optisel")` errored with “Reached maximum in the search” at ordinary `target_degree` values

A real
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) run
against the `planCross()`/`selectCrosses()` rewrite (below) surfaced
`SimpleMating::selectCrosses() failed: "Reached maximum in the search, try to increase data size."`
at `target_degree = 30` – an ordinary, non-boundary value — on a
10-parent/45-cross test panel requesting only 5 crosses. Root cause: the
`target_degree` -\> `culling.pairwise.k` mapping linearly interpolated
the cutoff across the candidate set’s raw `[min(K), max(K)]` range. That
is fragile against outliers: a single unusually related (or unusually
diverse) candidate pair widens the observed range enough that a
“moderate” `target_degree` can retain only a handful of surviving
candidates — too few for `selectCrosses()`’s internal greedy search to
assemble any `n_crosses` solution. (A first, narrower fix nudged the
cutoff by a small epsilon to stop it landing exactly on
`min(K)`/`max(K)`, which resolved the `target_degree = 0`/`90` boundary
cases specifically, but did not address this broader range-fragility
problem at interior values.) Replaced the range-based interpolation with
a cutoff chosen by *quantile* of the candidate set’s actual observed `K`
values, with a floor (`max(n_crosses * 5L, 20L)`, capped at the
candidate count) that always keeps enough candidates for the requested
`n_crosses` to stay feasible, regardless of how `K` happens to be
distributed. `target_degree = 0` now keeps only the least-related
candidates (subject to that floor); `target_degree = 90` keeps
essentially all of them. New regression test: “select_parents_ocs
(optisel): runs at target_degree boundaries 0 and 90 without erroring”
(`test-ocs.R`); the existing `target_degree = 30` end-to-end tests now
also exercise the fix.

### Fixed: `select_parents_ocs(engine = "optisel")` rewritten around `planCross()`/`selectCrosses()` — `SimpleMating::GOCS()` is not a reliable target

A real
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) run
surfaced
`SimpleMating::GOCS() failed: 'GOCS' is not an exported object from 'namespace:SimpleMating'`.
The first diagnosis — an outdated local install predating a
`GOCS()`-adding 0.2.x rewrite, based on the live
`Resende-Lab/SimpleMating` GitHub `main` branch NAMESPACE/DESCRIPTION
(`export(GOCS)`, `Version: 0.2.1`) — turned out to be wrong: a second
real test run, against a freshly `force = TRUE`-reinstalled
`SimpleMating 0.2.1`, confirmed the *actually installed* package still
does not export `GOCS()`. Its real help index lists `contrib2Cross`,
`getIndex`, `getMPV`, `getTGV`,
`getUsefA`/`getUsefAD`/`getUsefA_mt`/`getUsefAD_mt`, `planCross`,
`relateThinning`, `selectCrosses`, `setCrosses` — no `GOCS`, despite
`R/GOCS.R` and `export(GOCS)` still being present in the repo’s
checked-in source. Rather than keep chasing that unresolved
repo-vs-install discrepancy, `.run_simplemating_ocs()` (`R/ocs.R`) was
rewritten around
[`SimpleMating::planCross()`](https://rdrr.io/pkg/SimpleMating/man/planCross.html)
(builds the candidate cross list) +
[`SimpleMating::selectCrosses()`](https://rdrr.io/pkg/SimpleMating/man/selectCrosses.html)
(selects the final mating plan via a `Parent1`/`Parent2`/`Y`/`K` data
frame and a `culling.pairwise.k` relatedness cutoff) — both
independently confirmed present in the real, installed help index *and*
verified against their actual source fetched directly from
`https://raw.githubusercontent.com/Resende-Lab/SimpleMating/main/R/`.
The `target_degree` `[0, 90]` lever is now approximated by mapping onto
`culling.pairwise.k` over the candidate set rather than AlphaMate’s
continuous optimisation frontier. Two genuine improvements fell out of
the rewrite: `max_contrib_per_parent` is now honoured directly
(previously not configurable under `GOCS()`), and no-repeated-matings is
now guaranteed by construction rather than checked after the fact.
`.run_simplemating_uc()` (`R/genomic_mating.R`, used by
`usefulness_criterion(variance_model = "simplemating")`) was
independently verified unaffected — it calls
[`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html),
confirmed present in the real help index and matching its fetched real
source exactly, so no changes were needed there. Both wrappers still
check their required exports up front via
`getNamespaceExports("SimpleMating")` and stop with an explicit,
actionable message (installed version number + the exact
`remotes::install_github(..., force = TRUE)` reinstall command) instead
of letting a confusing
[`.Call()`](https://rdrr.io/r/base/CallExternal.html)-level “not an
exported object” error surface. Shared test helper
`skip_if_simplemating_too_old()` (`tests/testthat/helper.R`) skips
affected tests cleanly (rather than failing) when a too-old SimpleMating
is installed, alongside the existing
`skip_if_not_installed("SimpleMating")` guards.

### New: additive + dominance GBLUP in `run_haplotype_prediction(include_dominance = TRUE)`

Closes the “no dominance/epistatic relationship matrix — everything is
additive-only” gap. When `include_dominance = TRUE` (default `FALSE`,
requires `marker_effect_method = "gblup"` and the `BGLR` package),
[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
now also computes the Vitezica et al. (2013) genomic dominance
relationship matrix
([`compute_dominance_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md),
from the raw biallelic SNP dosage matrix — dominance coding is undefined
for the multi-allelic haplotype-block feature matrix) and fits a
dual-kernel RKHS model via
`BGLR::BGLR(ETA = list(A = list(K = G, model = "RKHS"), D = list(K = G_dominance, model = "RKHS")))`,
extending this package’s existing, already-used single-kernel BGLR/RKHS
pattern rather than introducing new external-API risk. Returns three new
fields: `dominance_deviation` (per-trait dominance BLUPs),
`total_genetic_value` (additive + dominance), and `G_dominance`.
`snp_effects`/`gebv` continue to reflect the additive component only, so
existing callers see no change in shape or values when
`include_dominance` is left at its default. New tests in
`test-prediction.R` (7 cases: default-off regression guard, argument
validation, full run with field checks, exact `total_genetic_value`
formula check, `G_dominance` symmetry, and an additive-only
shape-unchanged regression check).

### New: `variance_model = "linked"` in `usefulness_criterion()` — native, linkage-aware progeny variance

Closes the “UC’s variance model is block-independent, no linkage between
blocks” gap without adding the GitHub-only `SimpleMating` dependency
that the existing `variance_model = "simplemating"` option requires.
Given the same phased-haplotype/`snp_info`/`snp_effects` input as
`"phased"`, plus a `genetic_map` (SNP/CHR/cM), `"linked"` simulates
`n_sim_linked` (default `2000L`) Monte Carlo progeny: each parent’s
transmitted gamete is generated by a block-to-block crossover walk along
the genetic map, with recombination fraction between adjacent target
blocks from Haldane’s (1919) mapping function (independent assortment, r
= 0.5, across chromosome boundaries). This closes “phased” mode’s own
“target blocks summed independently” simplification — blocks close
together on the genetic map are now correlated the way real linkage
requires. Implemented as first-principles Monte Carlo (Haldane’s mapping
function + a standard Markov-chain crossover walk) rather than a
hand-derived closed-form multi-locus covariance formula, deliberately
mirroring this package’s existing risk-avoidance precedent for
`.run_simplemating_uc()`. Hand- verified in `test-genomic-mating.R`
against two exact closed-form boundary cases: blocks placed far enough
apart that r -\> 0.5 (must recover the independent sum of `"phased"`’s
own per-block variance) and blocks at r = 0 (must recover an exact
4-combo enumeration treating the whole multi-block region as one fused
super-locus), plus a monotonicity check at intermediate r, an
RNG-non-disturbance check, and end-to-end
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
smoke/validation tests. Same scope as `"phased"`: single-cross F1-style
segregation variance, not a multi-generation RIL/DH population variance
(that remains `"simplemating"`’s `type`/`generation` territory).

### New: OCS-/UC-informed rapid-cycling schemes in `ga_vs_ts_simulation()`

Closes the “forward simulation only compares truncation selection, not
OCS-/UC-informed rapid cycling” gap. A new `schemes` argument
generalises the original 2-scheme GA-vs-TS comparison to an arbitrary
named list of schemes, each with a `founders` set and a `mating_scheme`
of `"truncation"` (unchanged: random-mate the current generation, keep
the top `selection_intensity` fraction by GEBV), `"ocs"` (each
generation, `select_parents_ocs(engine = "optisel")` picks an
optimum-contribution mating plan from the live simulated population’s
dosage genotypes + GEBVs, executed via
[`genomicSimulation::make.targeted.crosses()`](https://rdrr.io/pkg/genomicSimulation/man/make.targeted.crosses.html)),
or `"uc"` (each generation, per-block local GEBVs are computed via
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
and `usefulness_criterion(variance_model = "block_independent")` ranks
every candidate pair; requires a new `blocks` argument). `ga_selected`/
`ts_selected` remain as backward-compatible convenience arguments that
build the original 2-scheme `"truncation"`-only list internally when
`schemes` is not supplied — existing calling code, including the
returned `$ga_final`/`$ts_final` fields, is unaffected.
[`plot_ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md)
now uses a colour-blind-safe palette that keeps GA/TS’s original colours
and extends gracefully to any number of additional scheme names.
`make.targeted.crosses()`’s signature (accepting parent NAMES, not just
indexes) was confirmed directly against the bundled real
genomicSimulation source (`simulation/sim-progression.R`) rather than
guessed — see the “API source” note in the `R/forward_simulation.R` file
header for details. New internal `.gs_normalize_schemes()` helper (pure
R, no genomicSimulation dependency) validates/normalises the `schemes`
argument and is unit-tested directly in the new
`test-forward-simulation.R` section 0, which runs unconditionally
(unlike the rest of that file, gated on genomicSimulation being
installed).

### Fixed: 8 failures surfaced by the first real `devtools::test()` run of this cycle

The first
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) run
of this development cycle against a real R installation found genuine
bugs alongside a few test-authoring mistakes:

- **`make_phased()` test fixture (`tests/testthat/helper.R`) was missing
  `dosage`/`sample_ids`.**
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)’s
  phased-list branch requires the full
  [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md)-shaped
  list (`hap1`/`hap2`/`dosage`/`sample_ids`), not just `hap1`/`hap2` —
  omitting `dosage` sent it down a `NULL`-dosage code path that failed
  with `argument is of length zero`. This broke `test-genomic-mating.R`
  and `test-extensions.R` at fixture-construction time, meaning every
  test in `test-genomic-mating.R` had never actually executed (only
  syntax-checked) before this run. Fixed by making `make_phased()`
  return the same shape
  [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md)
  produces.
- **`test-core-collection.R`** had two calls to
  [`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
  passing a distance matrix (`.cc_D`, diagonal = 0) without
  `type = "distance"`, so the default `type = "relationship"`
  misinterpreted it as a Gram matrix and derived all-zero pairwise
  distances. Fixed by adding the explicit `type = "distance"` argument.
- **`test-parent-selection.R`** compared `sort(res$selected)` against a
  numerically-ordered expected vector (`paste0("ind", 6:10)`), but R’s
  default string sort is lexicographic (`"ind10"` sorts before
  `"ind6"`). Fixed by switching to `expect_setequal()`, since selection
  order doesn’t matter for these assertions.
- **`test-association.R`**’s zscore-formula test recomputed mean/sd from
  [`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md)’s
  already-rounded (`normalize = FALSE`) output and rounded again,
  compounding two roundings against the function’s single internal
  rounding — a ~1e-6 mismatch, not a formula error. Fixed by loosening
  the comparison tolerance.
- **`test-genomic-mating.R`**’s `.selection_intensity()` “keeping
  ~everyone” test used `p = 0.9, n_progeny = 5L`, assuming
  `round(0.9 * 5) = round(4.5)` rounds up to `5` (triggering the
  function’s `k >= n_progeny` early return). R’s
  [`round()`](https://rdrr.io/r/base/Round.html) uses round-half-to-even
  (“banker’s rounding”), so `round(4.5) = 4`, not `5` — the test was
  hitting the ordinary Monte Carlo path (selecting the top 4 of 5) and
  correctly getting a nonzero intensity (~0.29), not a bug in
  `.selection_intensity()` itself. This test had also never executed
  before this run (same fixture-construction failure as above). Fixed by
  using `p = 0.99, n_progeny = 10L` (`round(9.9) = 10` unambiguously, no
  `.5`-tie rounding involved) to reliably exercise the intended branch.

### New: `select_parents_ocs(engine = "alphamate")` auto-detects a bundled AlphaMate.exe

A Windows build of the AlphaMate executable (Hickey Group / AlphaGenes
suite), plus the `libiomp5md.dll` it requires at runtime, now ships
under `inst/extdata/` purely as a development/testing convenience. When
`alphamate_exe = NULL` (the default) and that bundled copy is present,
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
finds and uses it automatically via
[`system.file()`](https://rdrr.io/r/base/system.file.html) – a message
is printed when this happens — and `engine = "auto"` now prefers
`"alphamate"` whenever a valid executable resolves. A non-Windows
platform now also triggers an explicit warning when
`engine = "alphamate"` is used, since the bundled binary cannot run
there.

**This auto-detection is a development/testing convenience only.**
AlphaMate is a separate, third-party tool that users must obtain and be
entitled to use themselves — it is not part of HapBlockR and not
guaranteed to be present in any given installation. A new roxygen
section (“Providing the AlphaMate executable”) and a new README section
spell this out explicitly and show how to pass your own `alphamate_exe`
path; `engine = "optisel"` remains the fully portable,
no-external-binary alternative. `test-ocs.R` gained new tests proving
the auto-detection/dispatch logic (guarded so they skip cleanly when the
bundled file is absent) without ever actually launching the executable,
and the two pre-existing tests that assumed no bundled executable would
be present were updated to force the “no valid exe” path explicitly
rather than relying on its absence.

A newer Beagle 5.x release (`beagle.28Jun21.220.jar`) was also added to
`inst/extdata/` alongside the existing `beagle.jar`. No code change was
made for it —
[`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md)/`run_ldx_pipeline(phase = TRUE)`
still require a `beagle.jar` supplied via `beagle_jar` or placed in
`out_dir`, as before; see the new README note under “A note on
`beagle.jar`” for how to use the bundled copy manually.

### New: Breeder’s Guide is now a PDF, not a Word document

`HapBlockR_Breeder_Guide` is now distributed as `.pdf` instead of
`.docx` (built from the same source layout, converted via LibreOffice
headless) so it reaches readers as a fixed, non-editable reference
rather than a document they could inadvertently alter.
[`open_breeder_guide()`](https://FAkohoue.github.io/HapBlockR/reference/open_breeder_guide.md),
`test-breeder-guide.R`, the vignette’s “See also” section, and README’s
Documentation section were all updated accordingly; the old `.docx`
copies (root and `inst/extdata/`) were removed.

### Tests: full coverage added for every function/argument added this cycle

The six new parent/cross-selection strategies,
[`compute_dominance_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md),
[`open_breeder_guide()`](https://FAkohoue.github.io/HapBlockR/reference/open_breeder_guide.md),
and several new arguments on existing functions
(`min_sel_value`/`min_sel_mode` on
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
and
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md);
`coancestry_weight`, `n_reps`/`$stability`/ `$converged` on
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md);
`normalize_method = "zscore"` on
[`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md))
had shipped with no test coverage. This release closes that gap with
nine new/expanded test files: `test-genomic-mating.R`, `test-ocs.R`,
`test-pareto.R`, `test-exact-validation.R`, `test-core-collection.R`,
`test-dominance-grm.R`, `test-breeder-guide.R` (all new), plus
expansions to `test-parent-selection.R`, `test-extensions.R`, and
`test-association.R`. Where possible, tests hand-verify exact expected
values against small fixtures worked out independently of the package’s
own code (e.g. the Vitezica dominance GRM, Pareto dominance, the exact
ILP optimum via exhaustive enumeration, and the farthest-point/maximin
traversal), rather than relying only on structural checks;
structural-only checks are used for the small remainder of functions
that wrap external optimisers (`GA`, `SimpleMating`, `optiSel`,
`lpSolve`), consistent with this package’s existing
`skip_if_not_installed()` convention for optional dependencies.

Writing the
[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
tests surfaced a genuine, previously-untested bug in
`R/core_collection.R`: the merit-floor branch indexed an unnamed
character vector (`common[floor_res$eligible]`) with character names
instead of using the already-filtered ID vector directly, which silently
returned `NA` for every candidate rather than erroring or filtering
correctly. Fixed to `common <- floor_res$eligible`. No prior execution
path had exercised this branch; it is now covered by
`test-core-collection.R`’s merit-floor tests.

Also adds a regression test (`test-extensions.R`) proving the earlier
[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)
phased-haplotype fix (see below) actually works end-to-end on
`make_phased()`-generated data, and hand-verified tests for
`score_favorable_haplotypes(normalize_method = "zscore")`
(`test-association.R`), confirming its output matches the documented
`(S - mean(S)) / sd(S)` formula exactly and preserves the same candidate
ranking as `"minmax"`.

### New: standalone, non-technical Breeder’s Guide + `open_breeder_guide()`

`HapBlockR_Breeder_Guide.pdf` is a new companion document covering the
same seven parent- and cross-selection tools as the *From Local GEBV to
a Crossing Decision* vignette
([`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)),
rewritten for breeders, programme managers, and reviewers who won’t run
the R vignette themselves: plain-language
what-it-is/when-to-use/when-to-be-cautious/how-it- works-in-practice for
each tool, breeding-programme analogies, a programme-shape decision
guide, and a glossary — no R code required to read it. Shipped as PDF
(not an editable Word document): drafted and laid out as a `.docx`
(dot-leader TOC, styled headings/tables), then converted to a fixed,
non-editable PDF via LibreOffice headless conversion for distribution –
the source `.docx` is not shipped.

A `.pdf` has no vignette engine, so it is not built or indexed by
[`vignette()`](https://rdrr.io/r/utils/vignette.html)/[`browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html)
the way the package’s `.Rmd` vignettes are. It ships as a static file
under `inst/extdata/` (the same convention this package already uses for
`beagle.jar` and its `example_*` datasets) and is located via
[`system.file()`](https://rdrr.io/r/base/system.file.html). New exported
function
**[`open_breeder_guide()`](https://FAkohoue.github.io/HapBlockR/reference/open_breeder_guide.md)**
(new file `R/breeder_guide.R`) is a thin convenience wrapper around that
lookup — locates the installed file and opens it with the OS’s default
`.pdf` viewer (`open = TRUE`, interactive sessions only), or just
returns its path. Also linked directly from the vignette’s “See also”
section and from README’s Documentation section (Section 17), and
available as a direct download from the repository root
(`HapBlockR_Breeder_Guide.pdf`) for readers who won’t install the
package at all. New NAMESPACE export; no new dependencies.

### Fixed: `decompose_block_effects()` silently dropped every phased block

Previously flagged as a known issue (not fixed) earlier in this
development cycle — now fixed.
[`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
emits two different per-individual string shapes depending on whether
the input was phased: one fixed-width `0/1/2` dosage string per
individual for unphased blocks, versus a `"gamete1|gamete2"` pair string
(one `0/1` digit per SNP, per gamete) for phased blocks (see
`extract_chr_haplotypes_phased_cpp()` in `src/ld_core.cpp`).
[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)
treated the whole pair string as a single allele label, which made its
length exceed the block’s SNP count by construction, so the digit-count
check silently returned `NA` for every allele and the block was dropped
with no warning — for any dataset built from phased haplotypes, the
function’s output was always empty. Fixed by detecting phased blocks
(via `block_info$phased`, with a `grepl("|", ...)` fallback for older
haplotype objects) and splitting each pair string into its two
constituent gamete alleles before frequency tabulation and effect
decomposition; this also corrects the allele-frequency denominator to 2
x n individuals (gametes), the biologically correct count, rather than
one diplotype label per individual. Unphased behaviour is unchanged.

### Docs: breeding-decisions vignette and README caught up to strategies 2-6

The *From Local GEBV to a Crossing Decision* vignette covered only
strategies 1 (coancestry-penalised GA) and implicitly the
founder-selection baseline through most of this development cycle, even
after
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md),
[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
and
[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)
were all implemented. New Sections 8-13 close that gap with runnable
examples built on the same worked example already in the vignette (no
new example data): UC cross ranking, an optiSel/GOCS()-engine OCS mating
plan (conditionally evaluated on `SimpleMating`+`optiSel` availability),
a narrowed Pareto sweep, an exact-ILP validation of a naive top-N-by-UC
plan against a constrained true optimum (conditionally evaluated on
`lpSolve` availability), a core-collection example, and a closing
decision table mapping each of the seven parent-selection tools to the
breeding question it answers and to a rough
population-improvement-vs-few-elite-families programme-shape guide.
Sections renumbered accordingly (old Section 8 “Translating this into a
breeding decision” -\> 14, old Section 9 “See also” -\> 15). README’s
Section 14.9 function-reference table gained rows for all five new
functions plus
[`pareto_front()`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md);
its now-outdated “What neither one does” paragraph (written when OCS was
still an unimplemented external-tool pointer) was rewritten to reflect
that OCS is now implemented in-package; Sections 1-2
(Motivation/Summary) and Section 4 (installation) gained corresponding
mentions and Suggests-package install commands (`GA`, `glmnet`,
`optiSel`, `SimpleMating`, `lpSolve`). The CIAT Peru production script
(`scripts/CIAT_Peru_parent_selection.R`) was deliberately left unchanged
— it is explicitly scoped to stop at parent selection for that specific
programme, not a general-purpose example, so extending its stages was
judged out of scope for a docs pass rather than a requested feature
change.

### New: strategies 4-6 — Pareto selection, exact ILP validation, core collections

Completes the six-strategy parent-selection roadmap (1:
coancestry-penalised GA; 2: usefulness_criterion(); 3:
select_parents_ocs(); 4-6: this entry).

- **[`pareto_front()`](https://FAkohoue.github.io/HapBlockR/reference/pareto_front.md)**
  (new file `R/pareto_selection.R`) — general, dependency-free
  non-dominated-sort utility (standard pairwise-dominance algorithm plus
  NSGA-II crowding distance as a secondary diagnostic). Works on any
  data frame with objective columns, e.g.
  [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
  output.
- **[`select_parents_pareto()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_pareto.md)**
  (same file) — sweeps
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
  `coancestry_weight` across a grid and Pareto-filters the resulting
  (merit, relatedness) points into an empirical frontier, so you can see
  the actual gain-vs-diversity tradeoff curve instead of guessing a
  single `coancestry_weight` value. Built entirely on top of the
  already-implemented, already-verified GA solver – no new optimisation
  algorithm (e.g. NSGA-II) was hand-built for this.
- **[`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)**
  (new file `R/exact_validation.R`) — a genuine EXACT solver (binary
  integer linear programming via the `lpSolve` package) for the core
  cross-selection problem every other tool in this package approximates
  heuristically (SimpleMating’s `selectCrosses()`/`GOCS()`, AlphaMate,
  the GA). Deliberately scoped as a small-scale validation/sanity-check
  tool (`max_vars` safety cap), not a production workflow replacement.
  Deliberately does NOT support a per-parent minimum-contribution hard
  constraint (SimpleMating’s `min.cross`) — exactly enforcing “min_cross
  applies only to parents used at all” requires a conditional big-M
  formulation, a well-known source of subtle ILP bugs, and this
  function’s entire purpose is to be the trustworthy reference;
  `max_cross` (a simple unconditional upper bound) and the total
  `n_cross` count are implemented as straightforward linear constraints
  instead. New guarded `lpSolve` `Suggests:` dependency.
- **[`select_core_collection()`](https://FAkohoue.github.io/HapBlockR/reference/select_core_collection.md)**
  (new file `R/core_collection.R`) — the classical core-collection /
  diversity-maximising subset problem (Schoen & Brown 1993), where
  earlier strategies treat diversity as a constraint on top of a merit
  objective, this treats diversity as the objective itself (with an
  optional merit floor via the same `.apply_merit_floor()` used by
  [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)/[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)).
  Implemented as the standard farthest-point/maximin greedy heuristic
  (Gonzalez 1985, proven 2-approximation) with a mean-distance (“MD
  strategy”) alternative — both simple, exactly-verifiable greedy
  algorithms, no external dependency, no execution-testing risk (unlike
  this package’s statistical-model-heavy functions). Converts a
  relationship matrix to genetic distance via the exact identity , not
  an approximation.

All four new exports (`pareto_front`, `select_parents_pareto`,
`validate_crosses_exact`, `select_core_collection`) are
balance-verified; `man/*.Rd` still needed via
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
same as every function added this development cycle.

### New: SimpleMating integration — linkage-aware UC and a tested optiSel engine for OCS

Two related additions after reviewing the current (2025) source of
`Resende-Lab/SimpleMating` (Peixoto et al. 2024) directly, rather than
relying on its paper’s description or an older 2021 prototype the
package author had on hand:

**[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
gains `variance_model = "simplemating"`**, wrapping
[`SimpleMating::getUsefA()`](https://rdrr.io/pkg/SimpleMating/man/getUsefA.html).
This is the most rigorous of the three variance modes now available: it
builds a genuine multi-locus Mendelian-sampling covariance matrix per
chromosome from genetic-map recombination fractions (Haldane-mapped) or
an LD-matrix proxy, following Lehermeier et al. (2017) – linkage across
every supplied SNP, not just within/across independently-summed LD
blocks like the existing `"block_independent"`/ `"phased"` modes. Real
constraint, not worked around: requires `geno_matrix` coded strictly 0/2
(fully homozygous DH/RIL-style calls) – not appropriate for heterozygous
outbred parents, for which the two native modes remain the right choice.
New arguments: `geno_matrix`, `G`, `genetic_map`/`ld_matrix`, `type`,
`generation`, `n_threads`.

**[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)’s
`engine = "optisel"` now wraps `SimpleMating::GOCS()`** instead of
calling
[`optiSel::candes()`](https://rdrr.io/pkg/optiSel/man/candes.html)/`opticont()`
directly. The original direct-optiSel implementation was written without
a working R interpreter available and, on comparison against GOCS()’s
actual verified source, turned out to have real argument-name mismatches
(guessed `kinship=`/`ub.kinship=`/`ub.n=`/`method="max.bv"`; GOCS()’s
working code uses `sKin=`/a flat `ub=`/`method="max.Crit"`). Rather than
patch guesses, this engine now wraps GOCS() directly, a tested
implementation. One consequence: GOCS() shares AlphaMate’s
Kinghorn-frontier-degree concept, so `ub_mean_kinship` is gone —
`target_degree` is now the single diversity-vs- gain lever for BOTH
engines (not numerically identical between them, just conceptually the
same dial). Two real, now-verified GOCS() limitations are surfaced
honestly rather than hidden: it fixes each parent’s contribution cap
internally at `2/n_candidates` (`max_contrib_per_parent`/
`n_parents_max` are not configurable for this engine, unlike
`"alphamate"`), and its mating allocator does not guarantee
no-repeated-matings during allocation the way AlphaMate’s evolutionary
algorithm does (checked afterward via `$ok` instead). The package’s own
`.allocate_matings_from_contributions()` native allocator, no longer
used, was removed.

New `SimpleMating` guarded `Suggests:` dependency
(`Remotes: github::Resende-Lab/SimpleMating`, not on CRAN), used by both
features above. `optiSel` remains a separate `Suggests:` entry since
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)’s
optiSel engine still requires it directly (as a dependency of
`SimpleMating::GOCS()`).

Original entry for
[`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
superseded by the above and left here for context: it originally shipped
with a hand-written
[`optiSel::candes()`](https://rdrr.io/pkg/optiSel/man/candes.html)/`opticont()`
call and a native mate-allocator, which turned out to need fixing once
GOCS()’s real source was available.

Not covered by either OCS engine: family/cross-of-origin representation
as a hard constraint (relevant for a population-improvement programme
that wants every family represented each cycle) — `family` is accepted
for output labelling only.

### Scientific-rigor audit of breeding-decision formulas

Full audit of every formula in the package used to inform a breeding
decision (GEBV/GBLUP, GRM, selection index, parent selection, UC/genomic
mating, forward simulation, marker effects, diversity/coancestry
metrics), against current best-established methods. Most held up well
as-is (VanRaden GRM, GBLUP/RR-BLUP/BayesA/B/C marker-effect dispatcher,
Falconer a/d diplotype decomposition, Weir & Cockerham FST). Three gaps
addressed this cycle (a fourth, Smith-Hazel multi-trait selection index,
was explicitly scoped OUT — the user already has a strong external
ASReml-based pipeline for that):

- **New
  [`compute_dominance_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_dominance_grm.md)**
  (`R/haplotypes.R`) — genomic dominance relationship matrix D via the
  Vitezica, Varona & Legarra (2013) natural-and-orthogonal
  parametrisation, diploid/biallelic only, computed from the raw SNP
  genotype matrix (not the multi-allelic haplotype-block feature matrix,
  where a dominance formula isn’t well established). This package does
  not itself fit a dual-kernel (G_A + G_D) model — `D` is meant to be
  paired with
  [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)’s
  `G` and fit externally
  (e.g. [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html),
  ASReml), the same external-handoff pattern
  [`prepare_gblup_inputs()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md)
  already uses for `G` alone. No new dependency.

- **Finite-population selection intensity in
  [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)**
  – new `n_progeny`, `n_sim`, `seed` arguments. The classical UC
  formula’s selection intensity assumes an infinite progeny population,
  which overstates `i_sel` for realistic biparental cross sizes (tens to
  a few hundred progeny). Supplying `n_progeny` switches
  `.selection_intensity()` to a Monte Carlo estimate of the true
  finite-sample expected order-statistic mean, rather than one of
  several mutually-inconsistent closed-form small-sample corrections in
  the literature — exact up to Monte Carlo error, and doesn’t require
  trusting a from-memory constant. Default behaviour
  (`n_progeny = NULL`) is unchanged.

- **`normalize_method` on
  [`score_favorable_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/score_favorable_haplotypes.md)**
  — new `"zscore"` option alongside the existing `"minmax"` (still the
  default, for backward compatibility). Min-max \[0,1\] rescaling has no
  direct genetic interpretation and is sensitive to the panel’s single
  most extreme individual; z-score standardisation (genetic SD units) is
  the more standard quantitative-genetics convention and more robust to
  panel composition changes.

**Deliberately deferred, not done this cycle** (both substantial enough
to warrant their own review pass): a linkage-aware (between-block)
variance mode for
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
— the scientifically strongest version of this would use Monte Carlo
gamete simulation via the already-integrated `genomicSimulation` package
rather than a hand-derived closed-form multi-locus covariance formula,
since the latter is hard to verify without being able to execute R in
this environment; and extending
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)’s
forward-simulation comparison to include a UC-informed (or, once built,
OCS-informed) recurrent-selection strategy alongside truncation
selection, for evaluating diversity-vs-gain tradeoffs over multiple
cycles.

### New: usefulness_criterion() — rank candidate crosses (genomic mating)

Second of six planned parent-selection strategy extensions. New file
`R/genomic_mating.R`, new exported function
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md).
Where
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
and
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
choose a *set* of parents,
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
instead ranks *pairs* of candidate parents (crosses) by the classic
Usefulness Criterion (Schnell & Utz 1975; Bernardo 2003; Zhong & Jannink
2007):

`UC = mid-parent GEBV + selection_intensity * sqrt(predicted cross variance)`

so it can surface crosses with strong transgressive-segregation
potential that whole-genome GEBV ranking alone would miss (a genuinely
good but already-similar pair of parents scores lower than a more
complementary pair predicted to segregate into a wider progeny
distribution).

Two variance-prediction modes, selected by the user via `variance_model`
(both are implemented; the package does not choose one for you — pick
whichever matches the genotype data you have):

- **`"block_independent"`** (default) — unphased-compatible, uses
  `local_gebv` (per-individual, per-block local GEBV) and a standard
  single-locus biparental segregation-variance formula per block, summed
  across target blocks. New `segregation_factor` argument (default
  `0.5`) exposes the formula’s scaling constant rather than hard-coding
  it.
- **`"phased"`** — requires phased haplotypes. Reads each parent’s two
  actual haplotype alleles per block via
  [`infer_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/infer_block_haplotypes.md),
  computes per-allele effects with a new, strictly-validated internal
  helper (deliberately NOT reusing
  [`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md),
  which was found during this work to silently drop every block when
  given phased input – see the “Known issue” note below), and exactly
  enumerates the 4 equally-likely gamete-pair combinations per block for
  an exact within-block segregation mean/variance.

In both modes, `predicted_variance` covers only the blocks in the
`block_importance` argument you supply (the same tractability trade-off
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
already makes) — documented explicitly in
[`?usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
so it isn’t mistaken for a whole-genome quantity. The mean term
(`mid_parent_gebv`) always uses the full whole-genome `gebv` you supply.

Scope, also documented in
[`?usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md):
this ranks candidate crosses, it does not select the parent set itself,
and it does not allocate differential progeny numbers or constrain
population-wide inbreeding the way true Optimal Contribution Selection
(OCS, strategy 3 of 6, not yet implemented) does.

No new package dependencies — uses only `stats` (`qnorm`, `dnorm`,
`var`) and `utils` (`combn`), both already in `Imports`.

**Known issue found while building this feature (now fixed — see the
“Fixed:
[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)…”
entry at the top of this file):**
[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)
silently returned no effects for any block when given phased
(“h1\|h2”-style) haplotype input — its digit-parsing logic assumed
unphased fixed-width dosage strings, so the length check never matched
and the block was dropped with no warning.
[`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)’s
phased mode avoided this entirely with its own helper at the time this
note was written;
[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)
itself has since been fixed directly and is now safe to use on phased
haplotypes.

**Still needed:** `man/usefulness_criterion.Rd` (run
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)),
same as the other functions added this cycle.

### Optional coancestry penalty in select_parents_ga()

First of six planned parent-selection strategy extensions (coancestry-
penalised GA, Usefulness Criterion/genomic mating, true OCS,
multi-objective Pareto selection, ILP exact-validation, core-collection
diversity selection), chosen as the most direct extension of the
existing GA.
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
gains two new arguments, both off by default (fully backward
compatible):

- **`G`** — a relationship/kinship matrix (e.g.
  `run_haplotype_prediction()$G` or
  [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
  output), `NULL` by default.
- **`coancestry_weight`** — numeric, default `0`. When positive (and `G`
  supplied), the fitness function becomes
  `sum(block coverage) - coancestry_weight * mean_pairwise_relationship(chosen set)`,
  so the GA trades off block coverage against relatedness among the
  chosen founders directly, rather than leaving relatedness as an
  after-the-fact
  [`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)
  diagnostic. The relationship term uses only off-diagonal pairwise
  values (relatedness *between* chosen parents, not each one’s own
  inbreeding/self-relationship).

New return fields: `$mean_relationship` (realised mean pairwise
relationship of the final selected set, from `G`; `NA` if `G` not
supplied) and `stability$mean_relationship` (same, per replicate). New
internal helper `.mean_pairwise_relationship()`.

This is explicitly a minimal extension of the existing block-coverage
GA, not a full optimal-contribution-selection (OCS) formulation — no
continuous contribution optimisation, no explicit inbreeding-rate
constraint, no mating list.
[`?select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
new “Coancestry penalty (optional)” section documents the fitness
formula and the (necessarily problem-specific) process for choosing
`coancestry_weight`, since the block-coverage and relationship terms
have no common natural scale.

### Merit-floor filtering and GA replication rigour for parent selection

[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
and
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
gain a shared eligibility floor, and
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
gains built-in replication/convergence diagnostics — both closing gaps
surfaced while comparing HapBlockR’s parent selection against a real
breeding programme’s
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)/
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
output and an external optimal-contribution- selection (OCS) tool.

- **`min_sel_value` / `min_sel_mode`** (new arguments, default
  `min_sel_value = NULL` = no change in behaviour): a whole-genome merit
  floor applied to the candidate pool before ranking/searching, in three
  modes — `"value"` (absolute cutoff), `"percentile"` (keep the top
  fraction by score, self-scaling), `"sd_below_mean"` (cutoff expressed
  in SDs below the mean, self-scaling). Without this,
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  could and did select genuinely poor overall performers (including one
  at a strongly negative selection index) purely because they uniquely
  covered one target haplotype block — its fitness function has no term
  for whole-genome merit.
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  also gains a `merit_score` argument (a named whole-genome value
  vector) since it otherwise has no such input to evaluate the floor
  against; `value_matrix` alone is per-block, not genome-wide.
- **`n_reps`** (new argument on
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
  default `5L`, was implicitly `1`): runs `n_reps` independent GA
  replicates from deterministically-derived seeds and returns the
  best-fitness one, plus a new `$stability` list reporting
  per-individual selection frequency across replicates and the fitness
  range observed — a single fixed-seed GA run gives no information about
  whether its answer is a robust optimum or one of several
  near-equally-good solutions. Set `n_reps = 1` for the previous
  single-run behaviour.
- **`$converged`** (new return field): flags whether the winning
  replicate’s GA stopped because fitness plateaued
  (`ga_fit@iter < maxiter`) rather than being cut off at `maxiter`
  without necessarily converging.
- New internal helpers `.apply_merit_floor()` and `.run_ga_once()` (not
  exported); the latter factors the single-GA-run logic out of
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  so replication can call it repeatedly.
- All additions are backward compatible: every new argument defaults to
  `NULL`/off except `n_reps`, which changes from an implicit `1` to an
  explicit default of `5` — a deliberate behaviour change (slower by
  default, but no longer silently trusting a single GA run) documented
  in
  [`?select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)’s
  new “GA rigour” section.
- Documentation:
  [`?select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  and
  [`?truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
  both gained a full “choosing between these two functions” comparison
  (what each optimises for, what neither does — coancestry management,
  differential contributions, actual mate allocation — and how that maps
  onto dedicated OCS/mate-allocation tools like AlphaMate). The *From
  Local GEBV to a Crossing Decision* vignette’s parent-selection section
  was expanded to match.
- **Not yet done:** `man/select_parents_ga.Rd`,
  `man/truncation_selection.Rd`, and five other exported functions in
  the breeding-decision layer (`plot_parent_selection_pca`,
  `ga_vs_ts_simulation`, `plot_ga_vs_ts_simulation`,
  `select_top_blocks`, `plot_block_funnel`) have no `.Rd` file at all
  yet —
  [`?select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  etc. currently return nothing. Run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  to generate them from the roxygen source (all now complete and
  correct); this was not done by hand here to avoid hand-authoring
  `\usage` blocks that must match the function signatures exactly.

### CI fix: R-CMD-check and pkgdown were both failing at “Install R dependencies”

Every matrix job of `R-CMD-check.yaml` and the `pkgdown.yaml` `build`
job failed identically at the dependency-installation step with
[`pak:: lockfile_create()`](https://pak.r-lib.org/reference/lockfile_create.html)
reporting “Could not solve package dependencies” and listing
`sessioninfo`/`rcmdcheck`/`covr` (or `pkgdown`) as having a “dependency
conflict” — a red herring; those packages were not actually in conflict
with anything. The real cause: `genomicSimulation` was added to
`DESCRIPTION`’s `Suggests:` field in 0.3.9.9000 when
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)
was rewritten as a genuine wrapper around it, but `genomicSimulation` is
not on CRAN and no `Remotes:` field told `pak` where else to find it.
Both CI workflows use `r-lib/actions/setup-r-dependencies@v2`, which
calls
[`pak::lockfile_create()`](https://pak.r-lib.org/reference/lockfile_create.html)
to solve the *entire* dependency graph (all of `Imports`/`Suggests` plus
each workflow’s `extra-packages`) in one pass; an unsolvable entry
anywhere in that graph fails the whole solve, and pak’s error reporting
attributes the failure to whichever packages it happened to be resolving
alongside it rather than to the actual missing package — which is why
the error mentioned `sessioninfo`/`rcmdcheck`/`covr` instead of
`genomicSimulation` and was non-obvious to diagnose from the log alone.

Fixed properly by keeping `genomicSimulation` in `Suggests:` (it is a
key dependency of
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md),
not an incidental one, and removing it from the formal dependency graph
was the wrong fix even though it made the immediate CI error disappear)
and adding:

``` R
Remotes:
    github::vllrs/genomicSimulation
```

[`pak::lockfile_create()`](https://pak.r-lib.org/reference/lockfile_create.html)
reads `Remotes:` to learn where to fetch packages that are declared in
`Imports`/`Suggests` but not resolvable from CRAN/PPM.
`remotes::install_github('vllrs/genomicSimulation')` is
genomicSimulation’s own documented supported install path for its
development version (per its README’s “The Very Latest Features”
section) — the Bunya HPC guide’s tarball-download procedure is a
workaround for offline/restricted compute nodes, not evidence that the
live GitHub source fails to build normally; on a standard,
internet-connected build machine (including GitHub Actions runners)
`install_github()` is expected to work directly. `NAMESPACE` still has
no `import(genomicSimulation)` directive — every call site remains a
guarded, fully-qualified `genomicSimulation::fn()` inside code that
checks
[`requireNamespace("genomicSimulation", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html)
first — so this change only affects what CI (and end users running
`remotes::install_github("FAkohoue/HapBlockR", dependencies = TRUE)`)
can resolve automatically; it does not change runtime behaviour when the
package is absent. This is also expected to let
`test-forward-simulation.R` actually run in CI for the first time,
rather than being whole-file-skipped via
`skip_if_not_installed("genomicSimulation")` as it always has been until
now.

Both `genomicSimulation` and HapBlockR itself require a C/C++ compiler
to build from source; all five R-CMD-check matrix runners and the
pkgdown Ubuntu runner already compile HapBlockR’s own Rcpp/RcppArmadillo
code successfully, so the toolchain is confirmed present. What is not
yet confirmed from within this environment is whether
genomicSimulation’s build succeeds cleanly on every one of those runners
without additional system libraries — if a CI run surfaces a *new*,
different failure specifically inside genomicSimulation’s own
compilation step (as opposed to the “dependency conflict” solve failure
this fixes), that would be the next thing to diagnose, separately, from
that run’s actual log.

## HapBlockR 0.3.11.9000 (development)

### `run_haplotype_stability()`: the `complete_decomposition = FALSE` fix (0.3.10.9000, item 1 below) was necessary but not sufficient — second, deeper bug found and fixed

The second
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) run
(2546 passing) showed the *same* four
[`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
failures at the *same* [`merge()`](https://rdrr.io/r/base/merge.html)
error (`'by' must specify a uniquely valid column'`), unchanged by the
`complete_decomposition = FALSE` fix. That fix was correct (it removed a
real column-count mismatch) but it was not the reason every block’s
regression was degenerate.

The actual cause: the function’s Finlay-Wilkinson environmental index
`I` was computed as the grand mean of each environment’s entire
`local_gebv` matrix
(`env_means <- vapply(local_by_env, function(m) mean(as.numeric(m)), ...)`).
But every entry of `local_gebv` is a `(x - ploidy*p) * alpha` term,
which is mean-zero by construction across individuals within a single
fit — this is a general property of BLUP/backsolved random-effect
predictions, not specific to any one block. So `env_means` was always
approximately `(0, 0, ...)` regardless of real between-environment
differences, `stats::var(I)` collapsed below the `1e-10` degeneracy
guard for *every* block, `rows` ended up empty,
`out <- do.call(rbind, rows)` became `NULL`, and the final
`merge(bk, out, by = "block_id", all.y = TRUE)` failed because `NULL`
coerces to a 0-row/0-column data frame with no `block_id` column to
merge on. Substituting the whole-genome
[`kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) GEBV
(`fit$g`) instead would not have helped either — BLUP predictions are
mean-zero by construction for the same reason.

Fixed by computing the environmental index the way Finlay & Wilkinson
(1963) originally defined it: the mean *observed* (raw phenotype/BLUE)
performance across all genotypes in each environment, not a derived
GEBV/local-GEBV quantity.
[`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
now tracks `pheno_env_means[env] <- mean(y[common], na.rm = TRUE)` (the
same `y` already pulled from `blues_list[[env]]`) inside its existing
per-environment loop, and `I` is now built from `pheno_env_means`,
restricted to whichever environments actually made it into
`local_by_env`. Nothing downstream of `I`’s assignment (the regression,
the `stable`/`p_b1` test, or the final
[`merge()`](https://rdrr.io/r/base/merge.html)) needed to change — only
the index’s source. With the shipped `ldx_blues_list`-style
two-environment test fixture (environments differing by a genuine +0.4
mean offset), `var(I)` is now comfortably above the degeneracy threshold
and every block produces a real row.

## HapBlockR 0.3.10.9000 (development)

### First real `devtools::test()` run: three bugs found and fixed

This is the first entry in this development cycle written after the test
suite actually ran
([`devtools::test()`](https://devtools.r-lib.org/reference/test.html),
456s, 2543 passing, 6 failing). All six failures traced back to three
distinct issues, all fixed:

1.  **[`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
    regression from
    [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)’s
    new `complete_decomposition = TRUE` default (4 test failures).**
    [`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
    (`R/haplotype_analysis.R`) calls
    [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
    internally and then [`merge()`](https://rdrr.io/r/base/merge.html)s
    the result against `blocks` by a `block_id` it reconstructs directly
    from `blocks`’ own `CHR`/`start.bp`/`end.bp` columns — a design that
    assumes `colnames(local_gebv)` corresponds 1:1 with rows of
    `blocks`. Once `complete_decomposition` defaulted to `TRUE`
    (0.3.4.9000), every call to
    [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
    that didn’t explicitly override it started adding extra
    singleton-pseudo-block columns for SNPs outside every block window —
    columns with no matching row in `blocks`, and (with the shipped
    `ldx_geno`/`ldx_blocks` example data, which has documented “isolated
    singleton” SNPs between blocks) enough of them to change the grand
    mean used in the function’s environmental-index calculation, which
    in turn made every block’s Finlay-Wilkinson regression degenerate
    and produced an empty results table — surfacing downstream as
    [`merge()`](https://rdrr.io/r/base/merge.html)’s
    `'by' must specify a uniquely valid column` (merging against a
    0-row/0-column empty data frame). This is exactly the kind of “the
    default changed underneath an existing caller” regression that
    inspection alone did not catch: I checked whether any *test*
    hardcoded an exact block count, but did not separately audit every
    other *internal caller* of
    [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
    for a same-default assumption. Fixed by passing
    `complete_decomposition = FALSE` explicitly in
    [`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)’s
    call, restoring its pre-0.3.4.9000 behaviour (it was never designed
    to handle singleton pseudo-blocks, and doesn’t need to — its whole
    point is comparing the same set of real LD blocks across
    environments).

2.  **`estimate_marker_effects(method = "bayesr")` was never valid (2
    test failures/errors).**
    [`BGLR::BGLR()`](https://rdrr.io/pkg/BGLR/man/BGLR.html) does not
    implement a model called “BayesR” — it supports
    `FIXED`/`BRR`/`BL`/`BayesA`/`BayesB`/`BayesC`/ `RKHS`. The “BayesR”
    method name was added under the mistaken assumption that BGLR’s
    model menu included it (BayesR, as in Erbe et al. 2012, is a
    four-component-mixture method not implemented by BGLR). This was
    invisible to code review because
    [`estimate_marker_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
    correctly *dispatches* to `BGLR::BGLR(..., model = "BayesR")` — the
    bug is one level down, inside BGLR itself refusing that model string
    at runtime. Fixed by replacing `"bayesr"` with `"bayesa"` everywhere
    (function signatures,
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)’s
    `marker_effect_method`, `solver_used` labelling, all roxygen docs,
    DESCRIPTION, and the HapSelect gap-analysis document). BayesA
    (marker-specific variances, scaled-inverse-chi-squared prior, no
    point mass at zero) is a real BGLR model and a genuinely different
    flavour from BayesB/BayesC (which add a point-mass-at-zero
    variable-selection component), so the “three Bayesian options”
    design intent is preserved with real, working models.

3.  **A test’s own expectation was mathematically wrong, not the source
    code (1 test failure).** `test-marker-effects.R` asserted that
    `estimate_marker_effects(method = "rrblup")`’s reported `gebv`
    should differ between `ploidy = 2` and `ploidy = 4`. It does not,
    and that turns out to be *correct* behaviour: the centring term used
    throughout this cycle’s ploidy generalisation is `ploidy * p_hat`,
    where `p_hat = colMeans(dosage) / ploidy` — and
    `ploidy * (colMeans(dosage) / ploidy)` is algebraically identical to
    `colMeans(dosage)` for *any* ploidy value (the ploidy cancels out of
    the centring step itself, as long as the `[1e-8, 1-1e-8]` safety
    clamp on `p_hat` doesn’t engage, which it does not for realistic
    dosage ranges). Ploidy’s real, discriminating effect is on the
    *scaling/denominator* terms (`ploidy * sum(p*(1-p))`), which is why
    the
    [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)
    ploidy tests correctly showed a difference across ploidy values
    while this one, wrongly, expected one for a step where there isn’t
    one. Fixed by correcting the test’s expectation (now documents and
    checks the ploidy-invariance of the centring step directly) rather
    than touching the source, which was already correct. Added a note to
    this effect in
    [`estimate_marker_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)’s
    documentation so this doesn’t need rediscovering.

None of these fixes change any test that was previously passing; the
[`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)/[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)/[`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
ploidy tests,
[`select_top_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)/[`plot_block_funnel()`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md)
tests, and
[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)/[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
tests (`GA` installed) all passed on the first run with no changes
needed.
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)
could not be exercised in this run (`genomicSimulation` not installed in
the test environment) and remains the top priority to verify next.

------------------------------------------------------------------------

## HapBlockR 0.3.9.9000 (development)

### ga_vs_ts_simulation() now wraps genomicSimulation, as HapSelect does

Supersedes the previous entry’s design decision. The 0.3.7.9000 release
shipped a self-built meiosis engine for
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)
instead of wrapping the `genomicSimulation` package, since that
package’s exact R API wasn’t yet confirmed against real documentation.
Real source of `genomicSimulation`’s R layer (`sim-setup.R`,
`sim-progression.R`, `sim-calculators.R`, `sim-group-utils.R`,
`sim-data-access.R`, `sim-deletors.R`, `utils.R`) has since been made
available locally and was read in full. `R/forward_simulation.R` has
been rewritten accordingly:
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)
and
[`plot_ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md)
keep their exact previous signatures (no change needed in calling code),
but the crossing/selection/GEBV engine underneath is now
`genomicSimulation` itself — the same package HapSelect uses for its own
`localGEBV_vs_TS_simulation()`/`Haplotype_vs_TS_simulation()`
comparison.

**How the wrapper works:** Each scheme (GA-founders, TS-founders) is run
as its own isolated `genomicSimulation` session: `clear.simdata()`, then
`load.data()` with a genotype-matrix file (founders’ phased alleles,
written with “0”/“1” as allele symbols to preserve `hap1`/`hap2`
cis/trans phase exactly), a genetic map file (`snp_info$POS` converted
to cM via `recomb_rate`), and a marker-effect file (`snp_effects`,
allele “1” gets the effect, allele “0” gets 0). Per-marker centring is
then set via
[`genomicSimulation::change.eff.set.centres()`](https://rdrr.io/pkg/genomicSimulation/man/change.eff.set.centres.html)
so that genomicSimulation’s raw-dosage GEBV calculation reduces exactly
to HapBlockR’s own `(x - 2p) * alpha` convention used throughout
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)/
[`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
— genomicSimulation was clearly designed with this exact use case in
mind (its own docs for `change.eff.set.centres.of.allele()` describe
“marker effects … calculated using GBLUP etc. … a count matrix centred
by allele frequency” as the intended input). Each generation:
`make.random.crosses()` produces `pop_size` offspring from the current
parent group; `see.GEBVs()` records that generation’s mean/max/sd; if
`selection_intensity` is set, `break.group.by.GEBV()` selects the top
fraction as next generation’s parents (matching HapSelect’s own
recurrent truncation-selection comparison scheme). `genomicSimulation`’s
single global `SimData` pointer means the two schemes cannot run
concurrently in one session — they run sequentially, each starting from
and ending with a fully cleared state (`on.exit(clear.simdata())`), so a
failure or interruption in one scheme cannot leak state into the other
or into unrelated `genomicSimulation` use elsewhere in the same R
session.

**Why each scheme gets its own session instead of sharing one population
split into two groups:**
[`genomicSimulation::make.group()`](https://rdrr.io/pkg/genomicSimulation/man/make.group.html)
reassigns group membership by index; if the same individual happened to
be selected as a founder by both the GA and TS methods, splitting one
shared population into two groups would silently move that individual
into whichever group’s `make.group()` call ran last, quietly shrinking
the other founder set by one. Loading each scheme’s founders into its
own freshly-cleared session avoids this edge case entirely rather than
working around it.

**What changed in the return value:** `ga_final`/`ts_final` still
contain `hap1`/`hap2` (now reconstructed from `genomicSimulation`’s own
internal genotype storage via `see.group.gene.data()`, with row order
read back from `see.genetic.map()` rather than assumed to match the
input `snp_info` order — documented explicitly, since
`genomicSimulation` does not guarantee it preserves input file row order
internally) plus a new `gebv` element (named numeric vector for the
final generation).

**Dependency:** `genomicSimulation` is not on CRAN and has been added to
`Suggests` without a version constraint. Install instructions (GitHub
source release `.tar.gz`, `install.packages(path, repos = NULL)`) are in
[`?ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)’s
“Installation of genomicSimulation” section, including a note that HPC
users must compile on the same CPU architecture they’ll run jobs on
(e.g. UQ Bunya’s “epyc3” node requirement) to avoid “illegal
instruction” crashes at runtime.

**Note:** this function depends on a second external package
(`genomicSimulation`) alongside its own. The exact row order
`see.genetic.map()`/`see.group.gene.data()` return, the exact rounding
`break.group.by.GEBV(percentage=)` uses, and how `load.data()`’s format
auto-detection or the explicit
`define.matrix.format.details(cell.style = "P")` format used here parses
single-digit “0”/“1” allele-pair strings as phase-preserving pairs are
all taken from `genomicSimulation`’s own roxygen documentation. Before
relying on this function: install `genomicSimulation`, run a small
example (e.g. the shipped `ldx_geno`/`ldx_snp_info` example data with
`n_generations = 2`, `pop_size = 10`), and confirm `sim$summary` looks
sane (non-degenerate GEBV spread, monotonic-ish trend under selection)
before scaling up.

------------------------------------------------------------------------

## HapBlockR 0.3.8.9000 (development)

### Ploidy-aware dosage centering (HapSelect gap-analysis, part 4, scoped)

The last of the four gap-analysis phases. Original scoping question was
whether to add full polyploid support; the answer implemented here is
deliberately narrower, for reasons explained below.

New `ploidy` parameter (default `2L`, fully backward-compatible) added
to:
[`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md),
[`estimate_marker_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md),
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md),
[`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md),
and
[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
(which passes it through to the first three). Generalises the VanRaden
(2008) centering/scaling used throughout the marker-effect and
local-GEBV pipeline from the diploid-specific `2p`/ `2*sum(p(1-p))`
convention to `ploidy*p`/`ploidy*sum(p(1-p))`, following the standard
dosage-scaling generalisation for polyploid genomic prediction (Endelman
et al. 2018). With `ploidy = 2L` every formula reduces exactly to its
previous form — this is a pure generalisation, not a behaviour change,
confirmed by inspection (each new formula collapses to the old hardcoded
`2`/`2p` expression when `ploidy` is 2).

**What this does NOT cover, and why:** HapBlockR’s phased-data
representation (`hap1`/`hap2` from
[`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md),
the “\|”-delimited 2-gamete haplotype strings built throughout
[`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
and
[`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md),
and the crossover meiosis model in
[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md))
is diploid by construction — it stores exactly two gametes per
individual throughout. Generalising that to arbitrary ploidy would mean
redesigning the core phased-data representation (hap1..hapK instead of
hap1/hap2) across every function that touches phased strings, not adding
a parameter to a formula — a much larger, higher-risk change to the
package’s central data model. It is also of debatable practical value:
reliable haplotype phasing is itself largely unsolved for polyploid
species, so most polyploid genomic selection in practice is dosage-based
(0..ploidy counts from array/GBS calls), not phase-based — which is
exactly the workflow this release’s `ploidy` parameter targets. The
compiled r2/rV2 LD kernels (`compute_r2_cpp`/`compute_rV2_cpp` in
`src/ld_core.cpp`) remain diploid/biallelic-only as originally scoped,
since they cannot be safely modified without the ability to recompile
and test in this environment.

[`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md)’s
He/Shannon/n_eff_alleles formulas were left unchanged: they operate on
categorical haplotype-string frequencies, not numeric dosage, so
“ploidy” doesn’t enter their arithmetic the same way — and since they
consume the same 2-gamete phased strings noted above, a `ploidy`
argument there would be unexercisable dead code given the current data
representation.

**Caveat:** as with every entry in this development cycle, none of this
has been executed against real or example data (no R available in the
environment it was written in). The claim that `ploidy = 2L` exactly
reproduces prior output was verified by direct algebraic inspection of
each changed formula (substituting `ploidy = 2` recovers the original
expression character-for-character), not by running the code.
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
should be run before trusting this, with particular attention to
[`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md)’s
`denom` calculation, which changed from an unconditional
`2 * sum(p*(1-p))` to `ploidy * sum(p*(1-p))` — identical at
`ploidy = 2` but worth a direct numeric check.

------------------------------------------------------------------------

## HapBlockR 0.3.7.9000 (development)

### New module: forward-in-time GA-vs-truncation-selection simulation (HapSelect gap-analysis, part 3)

Closes the remaining piece of HapSelect’s core differentiator: not just
selecting a founder set (previous entry), but demonstrating whether it
actually outperforms the obvious baseline over generations. New file
`R/forward_simulation.R`, two new exported functions:

**[`ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/ga_vs_ts_simulation.md)**
— runs recurrent-selection forward simulation for two founder sets
(typically
[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
vs. [`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
output) and tracks realised mean/max breeding-population GEBV over
`n_generations`. Mirrors what HapSelect’s own
`localGEBV_vs_TS_simulation()`/`Haplotype_vs_TS_simulation()`
demonstrate, but is HapBlockR’s own, self-contained simulation engine
rather than a `genomicSimulation` wrapper — at the time, that package’s
exact API hadn’t yet been confirmed against documentation, and wrapping
an unfamiliar external API before being able to test it seemed worse
than building a small, fully-owned, clearly-documented simulator. The
meiosis model is a standard simplified crossover simulation:
Poisson-distributed crossover count (rate = physical chromosome span x
`recomb_rate`, default ~1 cM/Mb), uniformly placed breakpoints,
alternating strand transmission, no crossover interference. Requires
**phased** haplotypes (`hap1`/`hap2`, e.g. from
[`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md))
— unphased 0/1/2 dosage cannot support a block-preserving meiosis
simulation, since without cis/trans information at heterozygous sites,
simulated recombination would not respect the LD blocks the founder sets
were selected to stack in the first place. Recurrent truncation
selection within each scheme (`selection_intensity`, default top 20%) is
applied each generation by default, matching HapSelect’s own documented
simulation scheme; set to `NULL` to instead re-mate the original founder
set every generation.

**[`plot_ga_vs_ts_simulation()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md)**
— plots mean (solid) and max (dashed) breeding-population GEBV over
generations for both schemes.

**Caveat:** as with all entries in this development cycle, this has
never been executed (no R available in the environment it was written
in). The crossover simulation, recurrent-selection loop, and GEBV
bookkeeping should all be validated against a real phased dataset before
being relied on – particularly runtime (a full generation loop with
`pop_size` individuals x `n_generations` x two schemes, each requiring
per-chromosome crossover simulation, is the most computationally heavy
new code added this cycle).

------------------------------------------------------------------------

## HapBlockR 0.3.6.9000 (development)

### New module: genetic-algorithm founder-parent selection (HapSelect gap-analysis, part 2)

Closes the single biggest capability gap identified against HapSelect:
turning block-importance/local-GEBV output into an actual founder-parent
*combination*, rather than stopping at a ranked list a breeder has to
manually pick from. New file `R/parent_selection.R`, three new exported
functions:

**[`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)**
— searches for the set of `n_founders` individuals jointly maximising
coverage of favourable per-block values (local GEBV, haplotype allele
dosage, or any comparable score) across a set of target blocks, via
`GA::ga(type = "binary")` (new `GA` Suggests dependency). Mirrors
HapSelect’s
`local_gebv_parent_selection()`/`haplotype_parent_selection()` in spirit
and crossing-scheme vocabulary (`strategy = "no_selfing"` (default),
`"selfing"`, `"OHS"`, `"OPV"`, `"Haploid_OHS"`), but is HapBlockR’s own
implementation: HapSelect’s exact GA fitness function and constraint
semantics aren’t available beyond its published documentation. The
fitness function is directly inspired by the pair-average formula shown
in HapSelect’s own documentation diagram
(`sum_j max((localGEBV_j1 + localGEBV_j2) / 2)`), generalised so a
block’s achievable value is either the mean of the two largest founder
values (`"no_selfing"`/`"OHS"` — two distinct parents required) or the
single largest founder value (`"selfing"`/`"OPV"`/`"Haploid_OHS"` — a
parent may supply both copies itself). The help page documents why the
latter three strategies compute identically (a direct consequence of the
max-pair-average formulation: self-pairing can never do worse than any
distinct pair), so this isn’t mistaken for the strategies being
unimplemented. A binary GA chromosome with a quadratic cardinality
penalty is used (standard practice for GA subset selection, since
[`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html) has no
native hard-cardinality constraint), seeded with an initial population
of exactly-`n_founders` suggestions so the search starts in the feasible
region. An optional `top_candidates` prefilter keeps the chromosome
length manageable for large populations.

**[`truncation_selection()`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)**
— trivial top-`n` baseline by a single genome-wide score (whole-genome
GEBV, stacking index), used standalone and as the comparison arm for the
forward-simulation module (next entry).

**[`plot_parent_selection_pca()`](https://FAkohoue.github.io/HapBlockR/reference/plot_parent_selection_pca.md)**
— PCA of a genomic/haplotype relationship matrix, coloured by selection
group (GA-only, TS-only, both, neither) — mirrors HapSelect’s parent-set
PCA plot showing where each strategy’s founders sit relative to overall
population diversity.

**Caveat:** as with all entries in this development cycle, none of this
could be executed/tested in the environment it was written in (no R
available this session). The GA-based search in particular has never
been run — please validate convergence behaviour, timing, and output
sanity via
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) and
a real dataset before relying on it, and inspect `ga_fit`
(e.g. `plot(select_parents_ga(...)$ga_fit)`) to confirm convergence.

------------------------------------------------------------------------

## HapBlockR 0.3.5.9000 (development)

### Cleanup and block-selection ergonomics (HapSelect gap-analysis, part 1)

Small, mechanical fixes plus block-selection tooling identified in the
HapBlockR-vs-HapSelect gap analysis. None change existing default
behaviour except where noted.

- **Stale “eleven C++ functions” claim fixed.** `README.md` said
  `src/ld_core.cpp` exports eleven compiled functions; it exports
  thirteen (`compute_r2_cpp`, `compute_rV2_cpp`, `maf_filter_cpp`,
  `build_adj_matrix_cpp`, `col_r2_cpp`, `compute_r2_sparse_cpp`,
  `boundary_scan_cpp`, `build_hap_strings_cpp`, `resolve_overlap_cpp`,
  `block_snp_ranges_cpp`, `extract_chr_haplotypes_cpp`,
  `extract_chr_haplotypes_phased_cpp`, `impute_and_filter_cpp` —
  confirmed against `R/RcppExports.R`). Both mentions in `README.md` §2
  and §16.1 updated and the full function list spelled out.
- **`.Rbuildignore` deduplicated.** `^cran-comments\.md$` and
  `^pkgdown$` were each listed twice. Also added `^.*\.bak$` so stray
  backup files never ship in the source tarball regardless of whether
  they’re removed from the working tree.
- **`R/ld_decay.R.bak` flagged for removal.** This stray backup
  (predating several fixes now in `R/ld_decay.R`) could not be deleted
  from this environment (no shell access this session) — `^.*\.bak$` in
  `.Rbuildignore` neutralises the packaging risk, but the file itself
  still needs manual removal from the repository.
- **Fixed a duplicate-bp-position bug in
  [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md)**
  (Density/Maximal Bron-Kerbosch path). Clique assignment previously
  round-tripped through bp *values* (`match(best, SNPbps)`,
  `setdiff(x, best)`), which silently misassigns SNPs whenever two
  columns share the same bp position (co-located multi-allelic sites, or
  an indel and a SNP recorded at the same POS — not rare in real VCFs):
  [`match()`](https://rdrr.io/r/base/match.html) returns only the first
  matching column, and [`setdiff()`](https://rdrr.io/r/base/sets.html)
  on bp values could strip a *different* SNP that merely shares a value
  with one actually selected. Cliques now carry the original column
  index alongside its bp position in lockstep throughout (including
  through
  [`.split_cliques_by_gap()`](https://FAkohoue.github.io/HapBlockR/reference/dot-split_cliques_by_gap.md),
  rewritten to split index/bp pairs together), so assignment no longer
  depends on bp values being unique.
- **Importance threshold parameterised.** The `>= 0.9` scaled-variance
  cutoff for the `important` flag was hardcoded in
  [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
  and in
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)’s
  multi-trait `importance_rule = "mean"` path. Both gain an
  `importance_threshold` parameter (default `0.9`, unchanged behaviour).
- **[`select_top_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/select_top_blocks.md)
  (new, exported)** — HapSelect-style block selection with three modes:
  fixed count (`n`), fixed percentage of all blocks (`perc_total`), or
  the smallest top-ranked set explaining a given percentage of
  cumulative variance (`perc_of_total_var`). Unlike a fixed
  `importance_threshold`, all three modes always return a usable
  selection.
- **[`plot_block_funnel()`](https://FAkohoue.github.io/HapBlockR/reference/plot_block_funnel.md)
  (new, exported)** — HapSelect-style funnel plot: local GEBV/haplotype
  effect values (x-axis) vs. scaled block variance (y-axis), with the
  current `importance_threshold` cutoff overlaid, so a breeder can see
  the ranked-variance elbow before choosing where to cut.

**Caveat:** as with prior entries, none of this could be executed/tested
in the environment these changes were written in (no R, no shell
available this session). Run
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
before relying on these changes, and manually remove `R/ld_decay.R.bak`
(and consider the stray bare `README` file alongside `README.md`, which
is a stale duplicate not otherwise flagged in `.Rbuildignore`).

------------------------------------------------------------------------

## HapBlockR 0.3.4.9000 (development)

### New: pluggable marker-effect models, and a complete GEBV decomposition matching HapSelect

Two changes closing the gap identified against HapSelect’s
`localGEBV`/marker- effect workflow.

**1.
[`estimate_marker_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
(new, exported)** — unified “Step B” marker- effect estimator, matching
the flexibility of HapSelect’s own `create_marker_effects_file()` (which
accepts rrBLUP / BGLR / Sommer / ASReml-R). `method =`:

- `"gblup"` (default) — the existing
  [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
  against a supplied relationship matrix, then backsolved to per-SNP
  effects via
  [`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  (VanRaden 2008; Tong et al. 2025). Identical to what
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
  has always done internally.
- `"rrblup"` —
  [`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)
  fit directly as a marker-effect (SNP-BLUP / ridge-regression BLUP)
  model (Meuwissen et al. 2001; Endelman 2011). No relationship matrix
  or backsolving step.
- `"bayesb"` / `"bayesc"` / `"bayesa"` — Bayesian variable-selection
  regression via
  [`BGLR::BGLR()`](https://rdrr.io/pkg/BGLR/man/BGLR.html) (Perez & de
  los Campos 2014). Requires the new `BGLR` Suggests dependency.
  (Originally shipped as `"bayesr"`; BGLR does not actually implement a
  “BayesR” model — this failed at runtime with BGLR’s own “model BayesR
  not implemented” error the first time the new test suite was actually
  run. Fixed by replacing it with the real BGLR model `"bayesa"` — see
  the 0.3.9.9000 entry below for the full story.)

All four methods return `gebv` on the same scale (`M %*% alpha`, the
centred genotype matrix times the estimated marker effects), so
block-importance rankings are comparable across method choices.
[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
gains a `marker_effect_method` parameter (default `"gblup"`, fully
backward-compatible: `solver_used` remains exactly `"rrBLUP"` for the
default) plus `n_iter`/`burn_in`/`seed` for the Bayesian methods.

**2.
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
— complete decomposition (`complete_decomposition = TRUE`, new
default).** Previously, any SNP not covered by a row of `blocks` —
isolated markers dropped at LD-detection time
(`singleton_as_block = FALSE`, the
[`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md)/[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
default), or blocks below `min_snps` — silently vanished from the
decomposition: it had a real, non-zero backsolved effect, but that
effect was never attributed to any block, so `rowSums(local_gebv)` fell
short of an individual’s true genome-wide GEBV. HapSelect does not have
this gap: its `def_blocks()` always keeps an isolated marker as its own
one-marker haploblock rather than dropping it, so its `localGEBV` values
are complete by construction.
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
now matches that guarantee: every SNP with a known effect that isn’t
covered by an existing block is appended as its own “singleton”
pseudo-block (`block_id` prefixed `"singleton_"`, `n_snps = 1`, flagged
via a new `singleton` column in `block_importance`), so
`rowSums(local_gebv)` reconstructs the full backsolved GEBV regardless
of how `blocks` was generated. Ranking (`var_scaled`/`important`) is
recomputed over the full combined set so singleton and LD blocks compete
on the same scale. Set `complete_decomposition = FALSE` to restore the
previous (possibly incomplete) behaviour.
[`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)’s
`n_blocks` now reflects the actual post-decomposition block count (LD
blocks + any singletons added), which may exceed `nrow(blocks)`.

The completeness guarantee depends on `geno_matrix` being the full,
unfiltered, genome-wide genotype matrix (the same one used to derive
`snp_effects`) rather than one already subset to block-member SNPs
(e.g. a haplotype feature matrix) — a SNP dropped upstream of
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
was never a column there to begin with and can’t be recovered inside it.
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
now warns when this is detectable: if `snp_info`’s own coordinates place
SNPs outside every block window, but none of those SNPs are present in
`geno_matrix` at all, `geno_matrix` is almost certainly pre-filtered and
the decomposition will look complete without actually being so.

**Caveat:** as with the formula-audit fixes above, none of this could be
executed/tested in the environment these changes were written in (no R
available). Run
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
before relying on
[`estimate_marker_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)
(especially the BGLR-backed methods) or the new
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
default in production.

------------------------------------------------------------------------

## HapBlockR 0.3.3.9000 (development)

### Correctness: formula audit fixes across LD, diversity, prediction, and association modules

A full audit of every estimator implemented in the package (LD metrics,
diversity, GEBV backsolving/local prediction, association testing,
cross-population concordance) against its cited source (VanRaden 2008,
Nei 1973, Weir & Cockerham 1984, Gao et al. simpleM 2008/2010/2011, Tong
et al. 2024/2025) turned up six formula-level defects. All six are fixed
in this release; none change the public function signatures.

**1. `prepare_geno(method = "rV2")` — missing-genotype NA propagation
(critical).** `V_inv_sqrt %*% geno_centered` is a full matrix product
across all `n` individuals. A single unimputed `NA` in one genotype
column silently contaminated that column’s whitened values for *every*
individual (not just the one with the missing call), which downstream
made `compute_r2_cpp()` treat the corrupted column as monomorphic and
zero it out — silently dropping SNPs from clique detection with no
warning. Fixed by mean-imputing (to 0 on the centred scale) immediately
before whitening. The `r2` path was never affected (its C++ kernel
already mean-imputes per column with no cross-individual mixing).

**2.
[`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md)
and
[`scan_diversity_windows()`](https://FAkohoue.github.io/HapBlockR/reference/scan_diversity_windows.md)
— Shannon entropy computed in nats instead of bits.** Both used natural
log ([`log()`](https://rdrr.io/r/base/Log.html)) while the documented
formula is `H' = -sum(p_i * log2(p_i))`. Every reported Shannon entropy
value was too small by a constant factor of `ln(2) ≈ 0.693`. Fixed by
switching both call sites to
[`log2()`](https://rdrr.io/r/base/Log.html).

**3.
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
— centering convention inconsistent with
[`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md).**
The backsolved per-SNP effects (`alpha`) are defined on mean-centred
dosages (`M = geno - 2*p`), but
[`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
was combining them with genotypes scaled to `[0, 1]` (`h = G/2`) via
`h %*% alpha + (1 - h) %*% (-alpha)`, an inconsistent parameterisation
that distorts local GEBV whenever allele frequencies deviate from 0.5.
Fixed to use the same centred convention as the backsolving step:
`local_GEBV = (G - 2*p) %*% alpha`.

**4.
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
— approximated EMMAX instead of implementing it.** The function fit the
null REML model, then subtracted the estimated fixed effects and
polygenic BLUP (`y - X*beta - u_hat`) and ran plain OLS on the
remainder, without ever transforming the marker columns. `u_hat` is an
estimate, not the true random effect, so the leftover residual is not
fully whitened, and OLS standard errors on untransformed markers do not
account for the correlation structure surviving in the residual —
weakest in exactly the related/structured populations this machinery
exists for. Replaced with a genuine GLS/EMMAX transform: both the
phenotype and every haplotype-allele column are pre-multiplied by
`V_hat^(-1/2)` (via the null model’s GRM eigendecomposition, computed
once per trait and reused for every marker), the transformed null design
(intercept + PCs) is projected out via the Frisch-Waugh-Lovell theorem,
and both the per-allele Wald scan and the per-block omnibus F-test now
run on the fully GLS-transformed quantities with correspondingly
corrected degrees of freedom (`n_ind - ncol(X_null) - 1` per allele;
`n_ind - ncol(X_null) - df_lrt` per block). A residual-scan fallback
(the previous approximation) is retained and used, with a warning, only
if the null model fit, its GRM eigendecomposition, or the transformed
null-design inversion fails. Also fixed a latent bug in the same
fallback path and in the per-block omnibus test: the per-allele degrees
of freedom were hardcoded assuming a single intercept (`n_ind - 2`,
`n_ind - df_lrt - 1`) rather than discounting the full null design size,
which under-counted degrees of freedom whenever `n_pcs > 0`
(anti-conservative p-values in the Q+K setting).

**5.
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)
— Cochran’s Q degrees of freedom off by one.** `Q_df` was computed as
`n_shared_alleles - 1`, but Cochran’s Q for `k` independent per-allele
effect-size comparisons has `k` degrees of freedom, not `k - 1` (there
is no shared common-mean estimation step consuming a degree of freedom
here, unlike a single-outcome meta-analysis). Fixed to
`Q_df = n_shared_alleles`.

**6.
[`compare_haplotype_populations()`](https://FAkohoue.github.io/HapBlockR/reference/compare_haplotype_populations.md)
— FST reduced to a simple allele-frequency-variance ratio instead of the
documented Weir & Cockerham (1984) estimator.** The previous formula
omitted the finite-sample bias correction and the within-population
heterozygosity component entirely, overstating FST for small sample
sizes. Replaced with the full three-variance- component WC1984 estimator
(`a`, `b`, `c`, with `n_bar`, `n_C`, `p_bar`, `s_squared`, `h_bar`),
verified against a well-established reference implementation
(scikit-allel). This also fixed a previously unflagged latent bug:
phased `"gamete1|gamete2"` diplotype strings were never split into
individual alleles before frequency/FST computation, so phased-data FST
was silently computed over diplotype combinations rather than true
allele frequencies. Unphased data (no heterozygosity component
available) falls back to `h_bar = 0`, equivalent to treating the data as
haploid, as before.

**Caveat:** these are formula-level corrections verified by manual
derivation and cross-referencing against the cited literature and, for
FST, an external reference implementation. Run
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) /
`R CMD check` before relying on results computed with this version,
particularly for
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md).

------------------------------------------------------------------------

## HapBlockR 0.3.2.9000 (development)

### New module: epistasis detection

Three new exported functions implement within-block and between-block
epistasis detection. All operate on GRM-corrected REML residuals from
the same null model as
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md),
ensuring population-structure-corrected tests throughout.

**[`scan_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_epistasis.md)**
— Within-block pairwise SNP epistasis scan.

Tests all C(p, 2) SNP pairs within significant blocks for the
interaction term `aa_ij` in
`y = mu + ai*xi + aj*xj + aa_ij*(xi*xj) + e`. Restricted to significant
blocks (controlled by `sig_blocks`) to avoid the genome-wide
combinatorial explosion. Multiple-testing correction: Bonferroni and
simpleM Sidak within each block; `Meff` estimated from the eigenspectrum
of the pairwise interaction column matrix. The `sig_metric` parameter
(default `"p_simplem_sidak"`) controls which correction drives the
`significant` flag. Output columns: `p_wald`, `p_bonf`, `p_simplem`,
`p_simplem_sidak`, `Meff`, `significant`, `significant_bonf`,
`significant_simplem`, `significant_simplem_sidak`. Returns
`HapBlockR_epistasis`.

**[`scan_block_by_block_epistasis()`](https://FAkohoue.github.io/HapBlockR/reference/scan_block_by_block_epistasis.md)**
— Trans-haplotype between-block epistasis scan.

Tests significant haplotype alleles against every allele at all other
blocks: O(n_sig × n_total_alleles) tests, Bonferroni corrected.
Identifies genetic background dependence where a resistance haplotype at
one locus only functions in the presence of a specific background at
another locus — a form of epistasis that single-block and single-SNP
analyses cannot detect. With 25 significant alleles × 17,943 total
alleles this scan involves ~450,000 tests. Returns
`HapBlockR_block_epistasis`.

**[`fine_map_epistasis_block()`](https://FAkohoue.github.io/HapBlockR/reference/fine_map_epistasis_block.md)**
— Single-block epistasis fine-mapping.

Identifies the specific interacting SNP pairs within one block.
Dispatches to exhaustive pairwise scan for blocks with p ≤ 200 SNPs
(`method = "pairwise"`), or LASSO with pairwise interaction terms via
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
at `lambda.1se` for larger blocks (`method = "lasso"`).
`method = "auto"` (default) selects automatically. Requires pre-computed
REML residuals (`y_resid`).

**Shared internal helpers added:** `.fit_null_reml()`,
`.pairwise_interaction_scan()`.

**`glmnet` added to `Suggests`** (required only for
`fine_map_epistasis_block(method = "lasso")`).

**43 new tests** in `tests/testthat/test-epistasis.R`.

------------------------------------------------------------------------

## HapBlockR 0.3.1.9000 (development)

### Enhancement: PC model selection and expanded plot outputs in test_block_haplotypes()

[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
now supports automatic selection of the number of population structure
covariates (n_pcs) via a new BIC/lambda/hybrid optimisation framework,
and produces three new diagnostic plots in PDF format.

**New parameters:**

- `optimize_pcs = FALSE` — when `TRUE`, fits REML null models for
  `n_pcs = 0..optimize_pcs_max` and selects the value minimising the
  score criterion chosen by `optimize_method`. Uses only the first trait
  (GRM is shared across traits). When `FALSE` (default), `n_pcs` is used
  directly.
- `optimize_pcs_max = 10L` — upper bound on the number of PCs evaluated
  during optimisation.
- `optimize_method = c("bic_lambda", "bic", "lambda")` — criterion for
  PC model selection (only used when `optimize_pcs = TRUE`):
  - `"bic"` — minimise BIC of the null REML model:
    `−2·logLik + k·log(n)`, where k = intercept + n_pcs + σ²_g + σ²_e.
  - `"lambda"` — minimise \|λ_GC − 1\|, where λ is estimated from a fast
    500-allele scan on GRM-corrected residuals. Most directly targets
    genomic control calibration.
  - `"bic_lambda"` (**default**, recommended for GWAS) — hybrid score:
    \|λ − 1\| + 0.01 × scaled_BIC. Minimises inflation/deflation while
    BIC breaks ties toward the simpler (fewer PCs) model.

**New return element:**

- `$pc_model_selection` — `data.frame` with one row per k tested:
  `n_pcs`, `BIC`, `lambda_gc`, `score`, `selected`. Printed by
  [`print()`](https://rdrr.io/r/base/print.html) and visible in the run
  log. `NULL` when `optimize_pcs = FALSE`.

**Plot changes (breaking):**

All plots now saved as **PDF** instead of PNG. Three plots are produced
per run:

- `manhattan_<trait>.pdf` — unchanged content, format changed to PDF.
- `qq_<trait>.pdf` — unchanged content, format changed to PDF.
- `pca_grm.pdf` — **new**: individuals in PC1 × PC2 space (GRM
  eigenvectors), coloured by the first trait’s phenotype using a
  blue-white-red diverging palette. Subtitle reports the number of PCs
  used in the model.
- `grm_scree.pdf` — **new**: scree plot of GRM eigenvalues (up to PC30),
  with the selected PC cutoff highlighted in red and the cumulative
  variance curve overlaid as a dashed line. When `optimize_pcs = TRUE`,
  the subtitle reports the selection criterion, selected k, and
  lambda_GC.

**Breaking change note:** any downstream code that expected
`manhattan_*.png` or `qq_*.png` file names must be updated to use the
`.pdf` extension.

### Enhancement: full multiple-testing correction set in estimate_diplotype_effects()

[`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md)
now computes the same four correction columns as
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md),
making the two modules symmetric. The `significant` flag is now driven
by a user-selected criterion rather than always using plain Bonferroni.

**New parameters:**

- `sig_threshold = 0.05` — significance cutoff applied to the p-value
  selected by `sig_metric`.
- `sig_metric = c("p_omnibus_adj", "p_omnibus_fdr", "p_omnibus_simplem", "p_omnibus_simplem_sidak")`
  — which correction drives the `significant` flag:
  - `"p_omnibus_adj"` — plain Bonferroni × n_blocks_per_trait (old
    default, retained for backward compatibility).
  - `"p_omnibus_fdr"` — Benjamini-Hochberg FDR.
  - `"p_omnibus_simplem"` — simpleM Bonferroni-style: min(p × Meff, 1),
    where Meff is estimated from the block-summary PC1 eigenspectrum.
  - `"p_omnibus_simplem_sidak"` (**recommended**) — simpleM Šidák-style:
    1 − (1 − p)^Meff. Consistent with
    [`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md).
- `meff_percent_cut = 0.995` — variance threshold for simpleM Meff
  estimation.
- `meff_max_cols = 1000L` — chunk size for large eigendecompositions.

**New output columns in `$omnibus_tests`** (all always present):

| Column | Description |
|----|----|
| `p_omnibus_adj` | Plain Bonferroni × n_blocks (unchanged) |
| `p_omnibus_fdr` | BH-FDR per trait |
| `p_omnibus_simplem` | simpleM Bonferroni-style: min(p × Meff, 1) |
| `p_omnibus_simplem_sidak` | simpleM Šidák-style: 1 − (1−p)^Meff |
| `Meff` | Effective number of independent tests (block-summary PC1 eigenspectrum) |

## HapBlockR 0.3.1.9000 (development)

### New feature: block_match = “position” in compare_block_effects() and compare_gwas_effects()

Both cross-population comparison functions now support matching LD
blocks between populations by **genomic interval overlap** rather than
by `block_id` string equality.

**Why this matters.** LD block boundaries are population-specific: the
same causal QTL region may be carved into a 100 kb block in Population A
and a 130 kb block in Population B, producing different `block_id`
strings (`"block_1_10000_85000"` vs `"block_1_10000_91000"`). With the
default `block_match = "id"`, these would not be compared — the block
appears as Pop1-only and the replication signal is lost.

**New parameters:**

- `block_match = c("id", "position")` — `"id"` (default) preserves
  backward compatibility; `"position"` matches by
  Intersection-over-Union (IoU) in base pairs.
- `overlap_min = 0.50` — minimum IoU for two blocks to be considered the
  same region. Blocks below this threshold are labelled `"pop1_only"`.

**New output column `match_type`** in `$concordance`:

- `"exact"` — same `block_id` string (boundaries identical)
- `"position"` — matched by genomic overlap (boundaries differ but IoU ≥
  `overlap_min`)
- `"pop1_only"` — no Pop2 block overlaps this Pop1 block at the
  threshold
- `NA` — no block tables were supplied

**New internal function `.match_blocks_by_position()`** performs the
interval join using a CHR-filtered IoU search. For each Pop1 block, it
finds the best-matching Pop2 block by IoU and records the match type.
Accessible via `HapBlockR:::.match_blocks_by_position()`.

8 new tests in `test-association.R` (122 → 130 total).

## HapBlockR 0.3.1.9000 (development)

### New function: compare_gwas_effects()

Cross-population effect concordance from **external GWAS results**.
Complements
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)
for users who ran association analysis outside HapBlockR (GAPIT, TASSEL,
FarmCPU, PLINK, or any other tool).

**Two input paths:**

- **Pre-mapped (recommended):** supply
  [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
  output for each population. Block assignment is explicit and
  auditable.
- **Raw GWAS + blocks (convenience):** supply raw GWAS data frames and
  block tables;
  [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
  is called internally.

**SE derivation.** When the SE column is absent (common in GAPIT/FarmCPU
output), SE is derived from the z-score: `SE = |BETA| / |Φ⁻¹(P/2)|`. The
`se_derived_pop1` / `se_derived_pop2` output columns flag which
population required this step.

**Key differences from
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md):**
External GWAS produces one lead SNP per block rather than multiple
haplotype allele effects. Consequently: - `effect_correlation` is always
`NA` (needs ≥ 3 alleles). - `direction_agreement` is 0 or 1 only. -
`Q_stat`, `Q_p`, `I2` are always `NA` (Cochran Q undefined with df =
0). - `replicated` uses `meta_p ≤ 0.05` instead of `Q_p > 0.05`.

**GWAS-specific output columns** added to `$concordance`:
`lead_snp_pop1`, `lead_snp_pop2`, `lead_p_pop1`, `lead_p_pop2`,
`se_derived_pop1`, `se_derived_pop2`, `both_pleiotropic`.

**Flexible column naming.** The `beta_col`, `se_col`, and `p_col`
arguments accept any column names (e.g. `"effect"`, `"std_err"`,
`"pvalue"`). The `Marker` column is accepted as an alias for `SNP`.

**Output class.** Returns `HapBlockR_effect_concordance` — the same
class as
[`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md).
The existing [`print()`](https://rdrr.io/r/base/print.html) method works
immediately.

18 new tests in `test-association.R` (104 → 122 total).

## HapBlockR 0.3.1.9000 (development)

### New function: compare_block_effects()

Cross-population haplotype effect concordance. Takes two
[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
result objects and returns per-block statistics for systematic GWAS
replication:

- **IVW meta-analysis**: inverse-variance weighted combined effect and
  SE per block, identical framework to two-sample Mendelian
  randomisation.
- **Cochran Q heterogeneity**: tests whether effect sizes differ
  significantly between populations. Significant Q (low Q_p) flags G×E
  interaction or population-specific LD structure differences.
- **I² inconsistency**: 0–100% measure of between-population
  heterogeneity. Values \> 50% indicate that the two populations are
  telling a different biological story at that block.
- **Direction agreement**: fraction of shared alleles with the same
  effect sign. Controlled by `direction_threshold` (default 0.75).
- **`replicated` flag**: composite criterion —
  `enough_shared AND directionally_concordant AND Q_p > 0.05`.
- **Block boundary diagnostics**: when `blocks_pop1` and `blocks_pop2`
  are supplied, `boundary_overlap_ratio` is **automatically computed**
  (not user-set) as bp(intersection) / bp(union) for every block. This
  quantifies how similarly the two populations carved the region into LD
  blocks. `boundary_overlap_warn` (input parameter, default `0.80`) is
  the threshold below which `boundary_warning = TRUE` is set in the
  output.
- **`$shared_alleles`**: per-allele detail table with `ivw_effect`,
  `ivw_SE`, `direction_agree`, and raw effects/SEs from both
  populations.
- 15 new tests in `test-association.R` (89 → 104 total).

### Enhancement: simpleM multiple-testing correction in test_block_haplotypes()

[`test_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/test_block_haplotypes.md)
now implements the simpleM procedure (Gao et al. 2008, 2010, 2011) for
LD-aware multiple-testing correction of correlated haplotype allele
tests. simpleM estimates the effective number of independent tests
(Meff) from the eigenspectrum of the haplotype allele dosage correlation
matrix, replacing raw test counting with an LD-aware effective count
that is less conservative than Bonferroni while providing family-wise
error control.

**New parameters:**

- `sig_metric`: which p-value drives `significant` /
  `significant_omnibus`. One of `"p_wald"`, `"p_fdr"`, `"p_simplem"`
  (Bonferroni-style), or `"p_simplem_sidak"` (Šidák-style, recommended).
  Default `"p_wald"`.
- `meff_scope`: scope for Meff estimation — `"chromosome"`
  (recommended), `"global"`, or `"block"`. Default `"chromosome"`.
- `meff_percent_cut`: variance threshold for simpleM eigendecomposition.
  Default `0.995` (99.5%), following the original simpleM
  recommendation.
- `meff_max_cols`: chunk size for large eigendecompositions. Default
  `1000L`.

**New output columns — always present regardless of `sig_metric`:**

- `allele_tests`: `Meff`, `alpha_simplem`, `alpha_simplem_sidak`,
  `p_simplem`, `p_simplem_sidak`, `p_fdr`
- `block_tests`: `Meff`, `alpha_simplem`, `alpha_simplem_sidak`,
  `p_omnibus_fdr`, `p_omnibus_simplem`, `p_omnibus_simplem_sidak` (plus
  `p_omnibus_adj` retained for backward compatibility)

**New return list elements:** `meff_scope`, `meff_percent_cut`, `meff`
(nested list of Meff summaries per trait: `$allele$global`,
`$allele$chromosome`, `$allele$block`, `$block$global`,
`$block$chromosome`).

13 new tests in `test-association.R` (76 → 89).

### Bug fix: read_phased_vcf() dot-ID synthesis

When a phased VCF has `.` or empty string in the ID column (column 3),
[`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md)
now synthesises `CHR_POS` identifiers for those rows, matching the
behaviour of
[`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md).
Without this fix, all-dot VCFs produced duplicate `rownames(hap1)` that
caused
[`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
to fail when slicing by SNP name.

- Handles `NA`, `"."`, and `""` (empty string) row-by-row.
- Real rsIDs in mixed VCFs are preserved; only missing IDs are replaced.
- Verbose message reports the count of synthesised IDs.
- 3 new regression tests in `test-phasing.R` (30 → 33).

Note: this gap only affects direct calls to
[`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md)
on externally phased VCFs with dot IDs. The
[`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)
path was already safe because
[`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
synthesises `CHR_POS` IDs before writing the cleaned VCF to Beagle.

## HapBlockR 0.3.1.9000 (development patch)

### Breaking change: SNP ID separator changed from `:` to `_`

When a VCF or GDS file has no rsID (i.e. the ID field is `.` or empty),
[`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
previously generated fallback SNP identifiers in `CHR:POS` format
(e.g. `"1:4106"`). These are now generated as `CHR_POS` format
(e.g. `"1_4106"`).

**Reason:** The colon `:` is a reserved character in many genomic file
formats (BED, VCF INFO field, R data frames used as rownames) and causes
silent parsing errors downstream. The underscore `_` is safe in all
contexts.

**Impact:** Any cached bigmemory backing files (`hapblockr_bm*.rds`,
`hapblockr_bm*.bin`, `hapblockr_bm*.desc`) built with the old format
must be deleted before rerunning. The pipeline will rebuild them
automatically.

**Files to delete before rerunning:**

``` r
results_dir <- "/your/results/dir"
for (f in list.files(results_dir, pattern = "^hapblockr_bm", full.names = TRUE))
  file.remove(f)
```

## HapBlockR 0.3.1.9000 (development patch)

### Bug fixes

- **`extract_chr_haplotypes_cpp()` — `retained_idx` field added
  (critical)** The C++ extractor now returns a `retained_idx` integer
  vector (1-based) giving the row index in the input
  `block_sb`/`block_eb` arrays for each retained block. Previously the
  compact `hap_strings` and `n_snps` arrays (indexed over retained
  blocks only) were paired in R with `chr_blk_srt[b, ]` (indexed over
  all blocks), causing a systematic mismatch: singleton coordinates were
  paired with multi-SNP haplotype strings and inflated `n_snps` values.
  Both the backend streaming path and the matrix path in
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  now use `ret_idx <- cpp_res$retained_idx` to fetch block metadata,
  eliminating the mismatch entirely. Sanity checks stop immediately with
  an informative message if `retained_idx` is missing or malformed
  (indicating the old compiled binary is still in use). Both extraction
  paths also pre-filter `chr_blk_srt` with `findInterval` before calling
  C++, providing a first-pass guard independent of compilation.

- **Pipeline QC singleton check** —
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)
  now asserts that no `block_info` row has `start_bp == end_bp` with
  `n_snps > 1` immediately after
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md),
  stopping with a clear message if the stale binary is still loaded.

- **`.validate_hap_output()` added** — called at end of pipeline;
  reports NA count in feature matrix, duplicated `hap_id` values, n_snps
  range, and hap_id / column alignment.

### Improvements

- **`alleles` column in writer** —
  [`write_haplotype_numeric()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_numeric.md)
  now decodes nucleotide sequences directly from `hap_info$hap_string` +
  `snp_info` per row, without calling
  [`decode_haplotype_strings()`](https://FAkohoue.github.io/HapBlockR/reference/decode_haplotype_strings.md).
  Single source of truth: the dosage string stored in `hap_info` is
  always the one used to build the matrix column.

- **Backend extraction list accumulation** — the backend streaming path
  in
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  now accumulates `block_info` rows in a pre-allocated list (`bi_rows`)
  and binds once at the end with
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html),
  matching the matrix path and eliminating O(n²) rbind growth.

- **`.prefilter_blocks_by_span()` helper** — shared `findInterval`-based
  pre-filter extracted to a single internal function used by both
  extraction paths; eliminates duplicated logic and future drift risk.

- **`extract_chr_haplotypes_phased_cpp()` — new C++ phased extractor
  (Item 5)** A dedicated phased haplotype extractor is now compiled into
  `ld_core.cpp`. It accepts `hap1_chr` and `hap2_chr` (0/1 gamete
  matrices) and builds `"gamete1|gamete2"` strings in one OpenMP pass,
  counting gamete frequencies correctly (each individual contributes two
  observations). Returns the same contract as
  `extract_chr_haplotypes_cpp()`: `retained_idx`, `retained_start_bp`,
  `retained_end_bp`, `hap_freq`, `freq_dominant`, etc. The matrix path
  in
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  now dispatches to this function when `isp = TRUE`, eliminating the
  previous R loop that called `build_hap_strings_cpp()` twice per block
  and concatenated results with
  [`paste0()`](https://rdrr.io/r/base/paste.html). The R `for(b)` loop
  is now identical for phased and unphased: both read
  `cpp_res$hap_strings[[b]]` from their respective C++ extractor.

- **Imputed backend parameter fingerprinting (Item 6)** —
  `.make_imputed_backend()` now accepts `maf_cut`, `min_callrate`, and
  `impute_method` parameters and saves them as
  `hapblockr_bm_imputed_params.rds` alongside the other four backing
  files. On reattach the saved fingerprint is compared to the current
  parameters; if they differ (or the fingerprint file is absent,
  indicating an old-format cache), all five files are cleaned and the
  backend is rebuilt. Without this, a run with `maf_cut = 0.10` would
  silently reuse an imputed backend built with `maf_cut = 0.05`.
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)
  passes `maf_cut`, `min_callrate`, and the `impute` argument through to
  `.make_imputed_backend()`.

- **C++ as bp coordinate truth — `retained_start_bp` and
  `retained_end_bp` (Item 8)** Both `extract_chr_haplotypes_cpp()` and
  `extract_chr_haplotypes_phased_cpp()` now return `retained_start_bp`
  and `retained_end_bp` integer vectors: for each retained block, the
  values are direct copies of `block_sb[b]` and `block_eb[b]` from the
  C++ arrays — no R-side index lookup into `chr_blk_srt`. Both R
  extraction loops check `!is.null(cpp_res$retained_start_bp)` and read
  coordinates directly from C++ when available. The R-side fallback
  (`chr_blk_srt[ret_idx[b], ]`) remains for binary-compatibility with
  old compiled objects. This removes the last avenue for a coordinate
  mismatch between the C++ compact array index and the R
  full-block-table row index.

------------------------------------------------------------------------

## HapBlockR 0.3.1 (development)

### Performance: C++ overlap resolution and WGS stall fixes

The
[`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md)
post-segment pipeline now routes through C++ for all overlap resolution,
eliminating several O(n²)–O(n·p) bottlenecks that caused chromosome
processing to stall after segment detection completed.

**`resolve_overlap_cpp()` — new exported C++ function**

Replaces the two sequential R `.resolve_overlap()` calls with a single
C++ pass implementing the identical cumulative-score split rule. Four
improvements over the R version:

- **BLAS DGEMM scoring** — for each overlapping block pair, disputed SNP
  scores are computed as `rowMeans(C_L²) − rowMeans(C_R²)` where
  `C_L = (Z_overlap.t() × Z_left_reps) / (n−1)` via Armadillo matrix
  multiply. This replaces a
  [`vapply()`](https://rdrr.io/r/base/lapply.html) loop calling
  `col_r2_cpp()` per SNP against all p columns (314k for chr1). Cost per
  disputed SNP drops from O(n × p) to O(n × 2k_rep) — a 15,700×
  reduction for chr1.
- **Lazy column cache** — each column of `adj_mat` is standardised at
  most once. Columns never accessed as representatives or overlap SNPs
  are never touched. For ~1,000 overlapping pairs using 20 reps each,
  fewer than 1% of chr1’s 314k columns are ever standardised.
- **Single pass** — one call replaces two. The second R call was needed
  because the first pass could create new adjacent overlaps; the C++
  while loop handles this naturally within the same pass.
- **OpenMP over pairs** — overlapping pairs are collected first in an
  O(n_blocks) scan, then resolved in `#pragma omp parallel for` (scores
  computed in parallel; boundaries applied serially in reverse to
  preserve index stability).

**R-level running counters**

Four additional O(n) scan operations replaced with O(1) running
counters:

- `sum(!is.na(LDblocks[,1L]))` called 419 times per chromosome →
  `ld_count`
- `max(which(!is.na(newLDblocks[,1L])))` in re-merge loop →
  `remerge_count`
- `min(which(is.na(newLDblocks[,1L])))` in re-merge loop →
  `remerge_count`
- `done <- LDblocks[!is.na(LDblocks[,1L]),]` →
  `LDblocks[seq_len(ld_count),]`

**Re-merge skip when modeNum=2**

The optional re-merge loop across forced cut-points is now skipped when
all cut-points are forced (`length(atfcut) >= length(cutpoints) - 2`),
which occurs when `cutsequence.modi` switches to `modeNum=2`. In this
mode every cut is forced by definition, so re-running CLQD across all
boundaries produces no benefit. With 419 forced segments on chr1 this
eliminated potentially hundreds of CLQD calls in a serial R loop.

**[`intersect()`](https://rdrr.io/r/base/sets.html) →
[`findInterval()`](https://rdrr.io/r/base/findInterval.html) in re-merge
check**

The per-pair forced-cut check
`length(intersect(eb[2L]:nb[1L], atfcut)) > 0L` — O(n_blocks × n_atfcut)
total — replaced with `findInterval(eb[2L], atfcut)` binary search:
O(log n_atfcut) per call.

**`score_overlap_cpp()` and `resolve_seam_cpp()` (internal static
helpers)**

Two additional C++ helpers in `src/ld_core.cpp` implement the BLAS
scoring kernel and a seam-local resolver for future use at truly massive
scale (pangenome panels where the global O(n_blocks) scan itself becomes
a bottleneck). These are `static` functions not exported to R.

### Tests

- `test-resolve-overlap.R` (307 lines, new file): parallel property
  tests for both `.resolve_overlap()` (R reference) and
  `resolve_overlap_cpp()` (C++), plus cross-validation tests verifying
  identical boundaries on all six biological cases: non-overlapping
  passthrough, union merge (B inside A), all-left path, all-right path,
  tie-breaking (zero-variance overlap), mixed cumulative-score split.
  Output invariant tests: no new blocks created, output is
  non-overlapping, global span preserved.
- `test-cpp.R`: five additional tests for `resolve_overlap_cpp()`
  covering dimensionality, non-overlapping passthrough, no remaining
  overlaps, start ≤ end invariant, and R vs C++ agreement on random
  overlapping input.

### Haplotype association testing and breeding decision functions

Four new exported functions complete the statistical inference and
breeding decision layer.

**Haplotype association testing**

- **`test_block_haplotypes(haplotypes, blues, blocks, n_pcs, top_n, min_freq, id_col, blue_col, blue_cols, alpha, verbose)`**
  — Block-level haplotype association tests via a unified Q+K mixed
  linear model (EMMAX/GAPIT3 formulation): y = mu + alpha·x_hap +
  sum(beta_k·PC_k) + g + e, where PC_k are GRM-derived eigenvectors
  (fixed effects for population structure) and g ~ MVN(0, sigma_g^2 G)
  is the polygenic kinship random effect. GRM is inverted once per trait
  via
  [`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)
  (O(n^3)); per-allele scan across all blocks is fully vectorised in a
  single
  [`crossprod()`](https://rdrr.io/pkg/Matrix/man/matmult-methods.html)
  call (O(n\*p), same BLAS trick as marginal SNP screening). Returns a
  `HapBlockR_haplotype_assoc` object with `$allele_tests` (per-allele
  Wald tests: effect, SE, t, p_wald, p_wald_adj, significant) and
  `$block_tests` (omnibus F-test per block: F_stat, p_omnibus,
  p_omnibus_adj, var_explained, significant_omnibus). `n_pcs = 0L`
  (default): EMMAX pure GRM; `n_pcs > 0`: Q+K model; `n_pcs = NULL`:
  auto-select via scree elbow.

- **`estimate_diplotype_effects(haplotypes, blues, blocks, min_freq, min_n_diplotype, id_col, blue_col, blue_cols, verbose)`**
  — Estimates additive and dominance effects from diplotype class means
  at each LD block, after GRM kinship correction. For each allele pair
  (A, B): additive effect a = (mean_BB - mean_AA) / 2; dominance
  deviation d = mean_AB - midpoint; dominance ratio d/a (0 = additive,
  +/-1 = complete dominance, \|d/a\| \> 1 = overdominance). Returns a
  `HapBlockR_diplotype` object with `$diplotype_means`,
  `$dominance_table` (a, d, d_over_a, overdominance per allele pair per
  block per trait), and `$omnibus_tests` (F-test per block).

**Breeding decision tools**

- **`score_favorable_haplotypes(haplotypes, allele_effects, min_freq, missing_string, normalize)`**
  — Scores each individual’s genome-wide haplotype portfolio against a
  table of known per-allele effects. Block score = sum(allele_effect ×
  dosage) per block. Genome-wide stacking index = sum across all scored
  blocks, normalised to \[0,1\] when `normalize = TRUE`. Returns a
  ranked data frame (one row per individual) with `stacking_index`,
  `n_blocks_scored`, `mean_block_score`, `rank`, and one
  `score_<block_id>` column per scored block for detailed inspection.

- **`summarize_parent_haplotypes(haplotypes, candidate_ids, allele_effects, blocks, min_freq, missing_string)`**
  — Produces a tidy long-format allele inventory (one row per individual
  × block × allele) for candidate parents. Reports allele dosage (0/1
  unphased; 0/1/2 phased), population allele frequency, optional allele
  effect, and a `is_rare` flag (freq \< 0.10). Includes rows with dosage
  = 0 so all candidates can be compared on the same rows. Primary tool
  for identifying complementary rare alleles across candidates and
  designing haplotype stacking crosses.

### New example dataset

- **`ldx_blues_list`** — named list of two environments (`env1`, `env2`)
  of named numeric BLUEs (120 individuals each), generated from the same
  polygenic architecture as `ldx_blues` with environment-specific
  offsets. Used in examples for
  [`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)
  and
  [`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md).
  Flat-file copy at `inst/extdata/example_blues_env.csv`.

### Analysis extension functions (v0.3.1 additions)

Ten new exported functions extend the pipeline without modifying any
existing function. All are in two new files: `R/analysis_extensions.R`
and `R/haplotype_inference.R`.

**Cross-validation and model evaluation**

- **[`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md)**
  — k-fold cross-validation for the haplotype GBLUP model. Masks
  phenotypes fold-by-fold, predicts via
  [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
  using the shared haplotype GRM, and returns predictive ability
  (Pearson r) and RMSE per trait per fold. Supports multiple
  replications and multiple traits. Returns an `HapBlockR_cv` object
  with `pa_summary`, `pa_mean`, `k`, `n_rep`.

**Population comparison**

- **[`compare_haplotype_populations()`](https://FAkohoue.github.io/HapBlockR/reference/compare_haplotype_populations.md)**
  — computes Weir-Cockerham (1984) FST and allele frequency differences
  per block between two named sample groups. Returns `FST`,
  `max_freq_diff`, dominant allele per group, chi-squared p-value (Monte
  Carlo, B=2000), and a `divergent` flag (FST \> 0.1 AND p \< 0.05).
  Suitable for breeding cycle monitoring and wild/elite panel
  comparisons.

**Visualisation**

- **[`plot_haplotype_network()`](https://FAkohoue.github.io/HapBlockR/reference/plot_haplotype_network.md)**
  — draws a minimum-spanning network of haplotype alleles within one LD
  block using [`igraph::mst()`](https://r.igraph.org/reference/mst.html)
  with Hamming-distance edge weights. Node size proportional to
  frequency. Optional group colouring via a `groups` named vector.
  Returns the `igraph` MST object invisibly.

**Multi-environment stability**

- **[`run_haplotype_stability()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_stability.md)**
  — Finlay-Wilkinson (1963) regression of per-block local GEBV
  contributions against the environmental index across environments.
  Returns slope b (stability coefficient), SE, R², deviation mean square
  s²d, and a `stable` flag (H0: b=1 not rejected at alpha=0.05).
  Requires at least 2 environments.

**Annotation export**

- **[`export_candidate_regions()`](https://FAkohoue.github.io/HapBlockR/reference/export_candidate_regions.md)**
  — converts
  [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
  output to BED (0-based, UCSC/BEDtools-compatible), CSV, or a named
  list ready for `biomaRt::getBM()`. Supports `chr_prefix`, LD-extended
  windows via `use_lead_snp`, and `padding_bp`.

**Effect decomposition**

- **[`decompose_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/decompose_block_effects.md)**
  — aggregates per-SNP additive effects (from
  [`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md))
  into a per-haplotype-allele effect table. Effect = sum(SNP_effect ×
  allele_dosage) per allele position. Returns `allele_effect`,
  `effect_rank`, and `frequency` per allele per block. Directly links
  prediction model output to selection index construction.

**Genome-wide diversity scanning**

- **[`scan_diversity_windows()`](https://FAkohoue.github.io/HapBlockR/reference/scan_diversity_windows.md)**
  — sliding-window He / Shannon / n_eff_alleles scan across the genome
  independent of LD block boundaries. Window size and step controlled by
  `window_bp` and `step_bp`. Returns a data frame with one row per
  window including `sweep_flag` (freq_dominant \>= 0.90).

**True haplotype inference and harmonisation**

- **[`infer_block_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/infer_block_haplotypes.md)**
  — converts raw haplotype strings to a structured per-individual,
  per-block diplotype table with explicit `hap1`, `hap2`, `diplotype`
  (canonical sorted string), `heterozygous`, `phase_ambiguous`, and
  `missing` columns. Handles both phased input (from
  [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md),
  where `phase_ambiguous` is always `FALSE`) and unphased input (where
  heterozygous genotypes set `phase_ambiguous = TRUE` unless
  `resolve_unphased = TRUE` triggers a maximum-parsimony heuristic).

- **[`collapse_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/collapse_haplotypes.md)**
  — merges rare haplotype alleles (below `min_freq`) into biologically
  meaningful groups rather than dropping them. Three strategies:
  `"rare_to_other"` (pool into `<other>`), `"nearest"` (merge with most
  similar common allele by Hamming distance), `"tree_based"` (UPGMA
  dendrogram; merges rare alleles at the coarsest cut that avoids
  merging common alleles with each other). Preserves a `label_map`
  attribute recording every original→collapsed mapping for use by
  [`harmonize_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/harmonize_haplotypes.md).

- **[`harmonize_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/harmonize_haplotypes.md)**
  — makes haplotype allele labels transferable across
  training/validation splits, populations, or environments. Builds a
  reference dictionary from alleles above `min_freq_ref` in the
  reference panel; matches target alleles by exact string first, then
  nearest Hamming neighbour (up to `max_hamming`), then labels unmatched
  alleles `"<novel>"`. Attaches a `harmonization_report` attribute
  reporting `n_exact`, `n_nearest`, `n_novel`, and `mean_hamming_dist`
  per block.

### New functions

- **[`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)**
  — LD decay analysis per chromosome. Estimates the distance at which r²
  (or rV²) drops below a critical threshold using a memory-efficient
  position-first random sampling strategy: pair indices are sampled from
  SNP positions only, then
  [`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md)
  loads only the unique SNP columns involved (~2% of a WGS chromosome
  for 50k random pairs). `compute_r2_sparse_cpp()` (Armadillo BLAS +
  OpenMP) is called once per chromosome in place of an R pair-by-pair
  loop. Two threshold approaches: fixed numeric (e.g. 0.1, standard GWAS
  practice) and parametric (95th percentile of r² between unlinked
  markers on different chromosomes, measuring the background
  kinship-induced LD level). Optional LOESS and Hill-Weir (1988)
  nonlinear decay model fitting. Chromosome-specific decay distances can
  be passed directly to `define_qtl_regions(ld_decay = decay)` to define
  biologically justified candidate gene windows. Censored distances
  (threshold never crossed within `max_dist`) are flagged with
  `censored = TRUE` in output and emitted as warnings. Requires no
  change to existing code.

- **[`plot_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ld_decay.md)**
  — ggplot2 visualisation of r² vs physical distance per chromosome,
  with optional raw points, threshold line, per-chromosome decay
  distance markers, and facet option.

### Algorithm improvements

- **LD-informed overlap resolution** replacing blind union merge in
  [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md).
  Overlapping blocks at sub-segment seams are now resolved by a
  cumulative-score boundary rule: for each disputed SNP,
  `score = mean_r2(with left core) - mean_r2(with right core)`. The
  cumulative sum is tracked across the overlap zone and the split
  boundary is placed at the last position where the cumulative score is
  \>= 0. Representatives are selected from boundary-adjacent SNPs of
  each core (nearest to the overlap zone, not first-k). Three clean
  cases: all-left (block A keeps overlap, block B shrinks), all-right
  (block A shrinks, block B keeps overlap), mixed (split at cumulative
  boundary). Falls back to union merge only when one core is empty (one
  block fully inside the other).

### C++ updates

- **`compute_r2_sparse_cpp()`** gains an `n_threads` parameter (OpenMP).
  Thread-local vector accumulation prevents contention; results are
  merged after the parallel loop. The R wrapper in `RcppExports.R` is
  updated to expose the new argument with default `n_threads = 1L`.

### Bug fixes and robustness

- **[`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  pair sampling** — `snp_info` is now sorted by POS within each
  chromosome before any pair index sampling. Previously, unsorted input
  could produce incorrect pair distances silently.
- **[`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  column precomputation** — prepared columns (mean-imputed, optionally
  whitened) are now computed once per chunk, not per pair. Monomorphic
  columns (sd \< 1e-6 after preparation) are skipped before any r²
  computation.
- **[`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  `both` mode** — sliding-window and random pairs are now kept separate:
  sliding window feeds the decay curve, random pairs feed the parametric
  threshold only. Previously they were merged into one pool, biasing the
  curve shape estimate.
- **[`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  fitting failures** — LOESS and Hill-Weir nonlinear failures now emit
  [`warning()`](https://rdrr.io/r/base/warning.html) with chromosome
  name and error message instead of being silently swallowed by
  `tryCatch(..., error = function(e) NULL)`.
- **[`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  threshold crossing** — linear interpolation between the last-above and
  first-below bins for sub-bin precision (was: first-below bin
  midpoint). Non-monotone curves fall back to first-below.
- **`r2_threshold` validation** — numeric thresholds outside \[0,1\], NA
  values, and unrecognised strings now produce informative errors
  instead of silent wrong results.

### Performance improvements

- **[`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)
  C++ optimisation** — haplotype strings now built via
  `build_hap_strings_cpp()`, replacing an R
  [`vapply()`](https://rdrr.io/r/base/lapply.html) loop that incurred
  one function call per individual per block. For a 3M-SNP panel with
  17,078 blocks and 204 individuals: ~3.5 million R calls -\> one C++
  call per block. Expected speedup: 20-50x (3.5 hours -\> ~5-10 minutes
  for that panel).
- **[`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md)
  bigmemory NA fix** — char-type big.matrix stores NA as -128; these are
  now correctly restored to `NA_integer_`. Also removed a redundant
  integer-\>double type conversion.

### WGS-scale acceleration — three new approaches

- **Louvain/Leiden community detection** (`CLQmode = "Louvain"` or
  `"Leiden"`): Polynomial-time O(n log n) community detection replaces
  Bron-Kerbosch
  ([`igraph::max_cliques()`](https://r.igraph.org/reference/cliques.html))
  which has exponential worst-case complexity. The original Big-LD run
  on a 3M-SNP WGS panel found 4.26 million maximal cliques in a single
  1500-SNP window and ran for \> 1 hour without completing.
  Louvain/Leiden finish the same window in \< 1 second. Block boundaries
  are equivalent to or better than the Density mode for WGS panels where
  LD structure is dense. Available in
  [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md),
  [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md),
  [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md),
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md).
  Edge weights are inversely proportional to base-pair distance so local
  LD structure is respected.

- **Sparse LD computation** (`max_bp_distance` parameter): When
  `max_bp_distance > 0`, only SNP pairs within that physical distance
  have their r² computed via `compute_r2_sparse_cpp()`. Pairs beyond the
  threshold are assumed to be in negligible LD and set to zero in the
  adjacency matrix. This reduces O(p²) LD computation to near-O(p) for
  WGS panels — at 500 kb, approximately 70–90% of pairs are skipped.
  Available in
  [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md),
  propagated through all callers. Default `0L` (disabled) preserves
  original behaviour.

- **Memory-mapped genotype store**
  ([`read_geno_bigmemory()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno_bigmemory.md)):
  Wraps any genotype source in a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  backed by a binary file on disk.
  [`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md)
  retrieves columns via OS page faults — only the requested bytes are
  loaded into RAM. Peak RAM is proportional to
  `n_samples × subSegmSize × 8 bytes`, not the full genome matrix.
  Backing files persist across R sessions: supply the `.desc` path to
  reattach without re-loading source data. Storage type `"char"` (1 byte
  per cell) saves 8× RAM vs double for 0/1/2 dosage values. Requires
  `bigmemory` (added to Suggests).

### Recommended configuration for 3M-SNP WGS panels

``` r
blocks <- run_Big_LD_all_chr(
  be,
  CLQmode         = "Leiden",    # polynomial — guaranteed connected communities
  CLQcut          = 0.70,        # sparser LD graph
  max_bp_distance = 500000L,     # skip pairs > 500 kb
  subSegmSize     = 500L,        # smaller windows for safety
  leng            = 50L,         # narrow boundary scan for WGS density
  checkLargest    = TRUE,        # belt-and-suspenders for Density mode
  n_threads       = n_threads
)
```

------------------------------------------------------------------------

## HapBlockR 0.3.0

### Breaking changes

- [`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md):
  removed `pheno_ids` parameter. Sample subsetting is the caller’s
  responsibility — use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype
  handling is not part of LD block detection and these functions had no
  algorithmic role.
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md):
  parameter `rV2method` renamed to `kin_method`.
- `digits` default changed from `6` to `-1` (no rounding) throughout.
- SeqArray backend replaced by SNPRelate throughout. GDS files are now
  created via
  [`SNPRelate::snpgdsVCF2GDS()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsVCF2GDS.html)
  and read via
  [`SNPRelate::snpgdsGetGeno()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsGetGeno.html).
  Existing `.gds` files created by SeqArray are not compatible; delete
  and re-convert with `read_geno(..., verbose = TRUE)`.

### New features — haplotype module (comprehensive rewrite)

- **Phasing functions**:
  - [`read_phased_vcf()`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md):
    read pre-phased VCF with `0|1` GT fields; returns `hap1`/`hap2`
    gamete matrices plus combined dosage.
  - [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md):
    call Beagle 5.x for statistical phasing of WGS data; returns path to
    phased VCF.gz.
  - [`unphase_to_dosage()`](https://FAkohoue.github.io/HapBlockR/reference/unphase_to_dosage.md):
    collapse phased gamete matrices back to 0/1/2 (internal helper; not
    exported).
- **Haplotype extraction**
  ([`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md)):
  - Auto-detects phased vs unphased input from list structure.
  - Unphased mode: diploid string `"012201"` per individual per block.
  - Phased mode: two-gamete string `"011|100"` per individual per block.
  - Blocks are strictly per-chromosome; cross-chromosome blocks are
    architecturally impossible.
- **Enhanced diversity metrics**
  ([`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md)):
  - `n_eff_alleles`: effective number of alleles = 1/Σpᵢ² (Hill 1973),
    ranging from 1 (monomorphic) to k (equal frequencies).
  - `sweep_flag`: TRUE when `freq_dominant ≥ 0.90`, flagging possible
    selective sweeps or strong founder effects (Difabachew et al. 2023).
  - `He` is now sample-size corrected following Nei (1973): He = n/(n-1)
    × (1 − Σpᵢ²).
  - Output now includes `CHR`, `start_bp`, `end_bp`, `n_snps` columns
    for direct use in genomic region analyses.
- **Haplotype string decoder**
  ([`decode_haplotype_strings()`](https://FAkohoue.github.io/HapBlockR/reference/decode_haplotype_strings.md)):
  - Converts raw dosage strings (e.g. `"02110"`) to nucleotide sequences
    (e.g. `"AGTTA"`) using REF/ALT from `snp_info`.
  - Returns a data frame: block_id, CHR, start_bp, end_bp, hap_rank,
    hap_id, dosage_string, nucleotide_sequence, frequency, n_carriers,
    snp_positions, snp_alleles.
- **Feature matrix**
  ([`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md)):
  - New `encoding =` parameter: `"additive_012"` (default) or
    `"presence_01"`.
  - Phased + `"additive_012"`: true 0/1/2 allele counts per gamete
    (0=absent, 1=one gamete, 2=both gametes).
  - Unphased + `"additive_012"`: 0/1/NA — 1=present, 0=absent. The value
    2 is not used because homozygosity cannot be inferred from unphased
    strings.
  - `"presence_01"`: 0/1 presence/absence for kernel methods or random
    forests (formerly `"presence_02"`; `"presence_02"` accepted as a
    backward-compat alias).
  - New `min_freq =` parameter to drop rare haplotype allele columns.
  - `top_n` default changed from `5L` to `NULL` (retain all alleles
    above `min_freq`).
- **QTL region definition**
  ([`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)):
  - Maps significant GWAS markers onto LD blocks post-GWAS.
  - Flags pleiotropic blocks (significant hits from multiple traits).
  - Implements the haploblock-based QTL cataloguing approach of Tong et
    al. (2024) *Theoretical and Applied Genetics* 137:274.
  - Now accepts optional `BETA` column in `gwas_results`. When supplied,
    output includes: `lead_beta` (effect of lead SNP), `sig_snps` (all
    significant SNP IDs, semicolon-separated), `sig_betas` (their
    marginal effects). See
    [`?define_qtl_regions`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
    for block effect estimation approaches.
- **Output writers** (all produce rows = haplotype alleles, cols =
  individuals, with metadata columns hap_id, CHR, start_bp, end_bp,
  n_snps before individual columns):
  - [`write_haplotype_numeric()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_numeric.md):
    tab-delimited dosage matrix (0/1/2/NA). Metadata column `alleles`
    shows the nucleotide sequence of this allele. Compatible with
    rrBLUP, BGLR, ASReml-R.
  - [`write_haplotype_character()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_character.md):
    tab-delimited nucleotide matrix. Each cell shows the nucleotide
    sequence if the individual carries that allele, `"-"` if absent,
    `"."` if missing. Heterozygous SNP positions are encoded with IUPAC
    ambiguity codes (R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C), keeping
    the sequence the same length as `n_snps`. `Alleles` column shows all
    REF/ALT at every block SNP (semicolon-separated).
  - [`write_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/write_haplotype_diversity.md):
    CSV of per-block diversity metrics with optional genome-wide mean
    summary row.
  - `write_haplotype_hapmap()` is **removed**. The diploid AA/AT/TT
    HapMap encoding was not meaningful for multi-SNP haplotype alleles.
- **Haplotype string decoder**
  ([`decode_haplotype_strings()`](https://FAkohoue.github.io/HapBlockR/reference/decode_haplotype_strings.md)):
  exported function mapping dosage strings to nucleotide sequences with
  full SNP-level metadata.
- **Multi-trait haplotype prediction** (extended
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)):
  - Extends
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
    to k traits simultaneously.
  - One trait-agnostic GRM (computed once) shared across all traits.
  - GBLUP solver: attempts
    [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html)
    multi-trait model first; automatically falls back to
    [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
    per-trait loop if sommer is not installed or the model fails to
    converge.
  - Block importance aggregated across traits: `var_scaled_mean`,
    `n_traits_important`, `important_any`, `important_all` per block.
  - `importance_rule` argument controls combined flag: `'any'` (≥ 1
    trait), `'all'` (all traits), `'mean'` (mean ≥ 0.9).
  - `blues` accepts a wide data frame (one column per trait) or a named
    list of named numeric vectors (different individuals per trait).
  - `sommer` added to `Suggests`; `rrBLUP` is the required fallback.
  - `run_haplotype_prediction_mt()` is removed; its functionality is
    fully absorbed into the unified
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
    API.
  - Cites Covarrubias-Pazaran (2016) for sommer and Endelman (2011) for
    rrBLUP.
- **Haplotype-based genomic prediction** (Tong et al. 2024, 2025):
  - [`compute_haplotype_grm()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_grm.md):
    VanRaden (2008) GRM from haplotype feature matrix. Mean-imputes NA,
    clamps frequencies, returns symmetric n×n matrix.
  - [`backsolve_snp_effects()`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md):
    derives per-SNP additive effects from GEBV without refitting the
    marker model: α̂ = M′G⁻¹ĝ / 2Σp(1−p) (Tong et al. 2025 *Theor Appl
    Genet* 138:267).
  - [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md):
    local haplotype GEBV per block = sum of SNP effects within block.
    Ranks blocks by Var(local GEBV); `important` flag when scaled
    variance ≥ 0.90 (Tong et al. 2024 *Theor Appl Genet* 137:274).
  - [`prepare_gblup_inputs()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_gblup_inputs.md):
    aligns haplotype feature matrix with a phenotype data frame;
    computes and returns a bended GRM ready for rrBLUP, sommer,
    ASReml-R, or BGLR.
  - [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md):
    end-to-end Tong et al. (2024/2025) pipeline from pre-adjusted
    phenotype values (BLUEs/BLUPs) to block importance. Accepts `blues`
    as a named numeric vector or a data frame with `id_col` and
    `blue_col` arguments. Uses
    [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
    for REML-based GBLUP. Returns all standard pipeline outputs plus
    GEBV, per-SNP effects, local haplotype GEBV matrix, and block
    importance table.
  - [`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/integrate_gwas_haplotypes.md):
    combines GWAS evidence (`has_gwas_hit`), variance evidence
    (`is_important`, scaled Var(local GEBV) ≥ 0.9), and diversity
    evidence (`is_diverse`, He ≥ threshold) into a `priority_score`
    (0–3) with plain-language `recommendation` per block.
  - [`rank_haplotype_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/rank_haplotype_blocks.md):
    unified block ranking across three use cases:
    1.  diversity-only — rank by He; (2) GWAS-only — binary GWAS-hit
        flag then He as tiebreaker, p-value not used for ranking; (3)
        phenotype — rank by scaled Var(local GEBV). Returns all standard
        pipeline outputs plus `ranked_blocks` data frame.
- **End-to-end pipeline**
  ([`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md)):
  - `hap_format = "hapmap"` renamed to `hap_format = "character"`.
  - All Big_LD arguments now exposed: `clstgap`, `split`, `appendrare`,
    `singleton_as_block`, `checkLargest`, `digits`, `kin_method`,
    `CLQmode`. Previously 8 of these were silently using Big_LD
    defaults.
  - `min_freq` parameter exposed (was previously hardcoded inside
    `build_haplotype_feature_matrix`).
  - Return list now includes `geno_matrix` (individuals × SNPs,
    MAF-filtered), enabling direct use with
    [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md)
    and
    [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
    without reloading the genotype file.
  - [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md)
    default grid now jointly optimises `CLQcut` (4 values) and
    `min_freq` (2 values), giving 8 combinations. Weber et al. (2023)
    show both are hyperparameters requiring dataset-specific tuning.

### New features — I/O and backend

- **Multi-format I/O**:
  [`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md)
  supports numeric dosage CSV, HapMap, VCF / bgzipped VCF, SNPRelate
  GDS, PLINK BED/BIM/FAM, and plain R matrices.
- **SNPRelate GDS streaming backend**: peak RAM proportional to one
  window. Replaces the previous SeqArray backend entirely.
- **PLINK BED backend**: `BEDMatrix`-backed memory-mapped access.
- **[`read_chunk()`](https://FAkohoue.github.io/HapBlockR/reference/read_chunk.md)**,
  **[`close_backend()`](https://FAkohoue.github.io/HapBlockR/reference/close_backend.md)**:
  unified accessor interface.
- **[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)**
  accepts `HapBlockR_backend` directly.
- **`min_snps_chr`** parameter: skip scaffolds with fewer SNPs (default
  10).
- **`clean_malformed`** parameter in
  [`read_geno()`](https://FAkohoue.github.io/HapBlockR/reference/read_geno.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md),
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md):
  streams input files to remove malformed lines before GDS conversion
  (for NGSEP and some GATK-origin VCFs).
- **`chr`** parameter in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/HapBlockR/reference/run_ldx_pipeline.md):
  restrict processing to specific chromosomes.
- **`ldx_gwas`** example dataset now includes a `trait` column
  (`"TraitA"` / `"TraitB"`) to demonstrate multi-trait pleiotropic block
  detection via `define_qtl_regions(trait_col = "trait")`.

### Performance improvements — C++ core

- `boundary_scan_cpp()`: C++ boundary scan replaces R inner loop (~20x
  faster).
- `maf_filter_cpp()`: combined MAF + monomorphic filter, single O(np)
  pass (~10x faster for panels \> 100 k SNPs).
- `build_adj_matrix_cpp()`: in-place adjacency construction, eliminates
  intermediate allocation.
- `compute_r2_sparse_cpp()`: sparse r² within bp window, avoids O(p²)
  cost for large sub-segments.

### Performance improvements — never-full-genome memory model

- **Chunked pre-allocated numeric CSV reader**: two-pass strategy —
  header scan only in Pass 1, then fixed 50,000-row chunks fill a single
  pre-allocated matrix in Pass 2. Peak RAM = one chunk, not 2x the file.
- **VCF and HapMap mandatory SNPRelate GDS auto-conversion**:
  transparent streaming GDS backend on first call; cache reused on
  subsequent calls.
- **Explicit per-chromosome gc()** in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
  and
  [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md):
  prevents heap fragmentation across 20-30 chromosome passes.

### Tests

- `test-ld-decay.R`: 33 tests for
  [`compute_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld_decay.md)
  and
  [`plot_ld_decay()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ld_decay.md)
  covering: input validation, `HapBlockR_decay` object structure, all
  sampling modes, unsorted `snp_info` handling, all threshold types
  (fixed, parametric, both, NULL), `fit_model` options, decay distance
  correctness, censored flag, all three backend types (matrix,
  `HapBlockR_backend`, bigmemory), print and plot methods, and
  [`define_qtl_regions()`](https://FAkohoue.github.io/HapBlockR/reference/define_qtl_regions.md)
  integration with `ld_decay=`.
- `test-basic.R`: smoke tests covering all major functions via the
  230-SNP example dataset (3 chromosomes, 9 LD blocks).
- `test-algorithm.R`: property-based tests for the Big-LD segmentation
  algorithm: `CLQD`, `Big_LD`, `run_Big_LD_all_chr`, `summarise_blocks`,
  `plot_ld_blocks`.
- `test-haplotypes.R`: comprehensive tests for haplotype extraction
  (phased and unphased), diversity metrics (`n_eff_alleles`,
  `sweep_flag`), QTL region definition (`lead_beta`, `sig_snps`,
  `sig_betas`), feature matrix encoding (`additive_012` /
  `presence_01`), output writers, `rank_haplotype_blocks`, and
  `integrate_gwas_haplotypes`.

------------------------------------------------------------------------

## HapBlockR 0.2.0 (2025-04-07)

### New features

- C++/Armadillo core: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()`.
- OpenMP parallelism via `n_threads` parameter.
- Dual LD metric: `method = "r2"` or `"rV2"`.
- [`prepare_geno()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_geno.md),
  [`compute_ld()`](https://FAkohoue.github.io/HapBlockR/reference/compute_ld.md).

### Bug fixes

- [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md):
  fixed `re_idx` tracking after dense-core pre-pass.
- [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md):
  [`vapply()`](https://rdrr.io/r/base/lapply.html) used throughout for
  type-stability.

------------------------------------------------------------------------

## HapBlockR 0.1.0 (2025-04-07)

Initial release.

- [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md),
  [`CLQD()`](https://FAkohoue.github.io/HapBlockR/reference/CLQD.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md),
  [`tune_LD_params()`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md).
- [`extract_haplotypes()`](https://FAkohoue.github.io/HapBlockR/reference/extract_haplotypes.md),
  [`compute_haplotype_diversity()`](https://FAkohoue.github.io/HapBlockR/reference/compute_haplotype_diversity.md),
  [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/HapBlockR/reference/build_haplotype_feature_matrix.md).
- [`summarise_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/summarise_blocks.md),
  [`plot_ld_blocks()`](https://FAkohoue.github.io/HapBlockR/reference/plot_ld_blocks.md).
