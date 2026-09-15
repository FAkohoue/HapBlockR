# DesiredGainR 0.5.0 (in development)

- Extended `compare_selection_methods()` to accept classical, restricted,
  general, iteratively optimised desired-gain and quadratic genomic fits.
  `comparison_objective()` fixes one target, utility, response scale and
  genetic covariance before methods are inspected. Internal adapters harmonise
  candidates and response units. Common validation values support per-trait
  response and linear or quadratic utility comparisons. QGSI approximations
  remain labelled separately from exact linear-index responses.

- Added the original common-condition comparison of expected response, target
  attainment, classical criteria, Mahalanobis alignment, rank agreement and
  selected-set agreement. The multiple-trait
  selection vignette now defines the major strategy classes, explains every
  implemented method, gives a method-choice framework and develops a complete
  worked comparison. Independent-culling and Elston limits now use original
  trait units. Economic indices now accept valid negative economic weights.

## Critical remediation of 1 August 2026

- Added a reproducible empirical validation suite using six real breeding-
  programme sources. It estimates common-check and year-adjusted recurrent-
  cycle trends, freezes a desired-gain index before later INIA rice cohorts,
  compares it with yield-only selection under a hierarchical temporal
  bootstrap, tests Genomes-to-Fields transport across non-overlapping year
  windows, and records the CIMMYT wheat predicted-versus-realized aggregate
  anchor. Every endpoint ships with its evidence tier and interpretation,
  making the package's mathematical, genomic and transport evidence directly
  auditable. A prospective comparison of competing vectors is identified as a
  future opportunity to estimate incremental field response among strategies.
  The sixth source is the Zhang et al. CIMMYT maize RCGS programme: raw
  two-location cycle gains reproduce the published trend, a nested
  leave-one-environment-out analysis tests desired-gain directions, and a
  marker-subsampled repeated five-fold GBLUP reproduction includes a permuted
  negative control. Cross-file checks expose the 43-versus-44 C4 discrepancy
  and the byte-identical C3/C4 selected genotype files. Unphased bulk-family
  calls are explicitly rejected as actual AlphaSimR founders.

- The population-driven recommender now uses `G_target` for breeder-facing
  genetic-SD units and treats `G_realised` only as a finite-QTL calibration
  diagnostic. It confirms multiple finalists, applies exact simultaneous lower
  and upper probability bounds, reports supported/uncertain/not-supported
  decisions, runs independent multi-start and half-budget stability checks,
  and validates finalists across newly sampled genetic architectures. Optional
  covariance draws are integrated when genetic/residual degrees of freedom are
  supplied; otherwise support is explicitly conditional on the point
  covariance.
- Added leakage-free cross-fitted RR-BLUP selection to
  `simulate_selection_cycles()` through its compact `prediction` list. The
  simulator reports held-out GEBV accuracy and never uses hidden simulated
  breeding values to construct a selection index.
- The exact box-QP active-set solver is now independently validated in CI
  against Clarabel (primary interior-point oracle) and polished/unpolished OSQP
  across analytical, random, ill-conditioned, degenerate and invariance cases.
  Full KKT, feasibility, objective, status and discrepancy results ship under
  `inst/validation/`.
- Added an independent multivariate-normal Monte Carlo test of the desired-gain
  response direction and achievable-response ellipsoid. It calls neither the
  package coefficient helper nor AlphaSimR, providing a separate check of the
  one-cycle selection-response mathematics.
- Renamed misleading `probability_best` outputs to
  `bootstrap_selection_frequency`; this is a resampling stability diagnostic,
  not a posterior probability that a vector is truly best.

- Added `suggest_desired_gains()` for breeders who cannot defend a complete
  desired-gain vector. It accepts trait-specific minimum gains in estimated
  `G_target` genetic-standard-deviation units and returns separate
  minimum-attainment and maximum-balanced recommendations. An exact convex
  active-set solution verifies one-cycle feasibility and the maximum common
  response on the achievable-response ellipsoid. Adaptive discovery,
  independent multiplicity-controlled screening, and independent final
  confirmation are separated so recommendation uncertainty is not estimated
  from the same simulation noise that selected the candidates.
- Desired gains can now be elicited as breeder-defined intervals with
  `define_desired_gain_intervals()`. `optimize_desired_gains()` restricts its
  search to the admissible interval cone and searches in genetic-standard-
  deviation units by default, while retaining representative requested gains
  and the trait-unit directions used by the simulator. A dedicated vignette
  explains signs, units, feasibility and what “optimum” means. Its new
  `mode = "interval"` recommends only actually simulated vectors, defines high
  gain as jointly reaching every breeder-specified lower bound, reports
  Jeffreys-binomial posterior probabilities and credible intervals per vector,
  and ranks conservatively by the lower joint-probability bound with balanced
  interval attainment as a tie-breaker.
- DGSI holdouts are split before optimisation; external validation data can
  select the replicate and never enter coefficient fitting.
- Mallard, Tallis and Kemp--Harville now use exact proportional-response
  contrasts. The undocumented Harville soft penalty and Mallard absolute-gain
  claim were removed.
- Covariance propagation passes sampled residual covariance directly to
  AlphaSimR, evaluates the point covariance separately, and computes rank churn
  only under a declared scalarisation.
- Optimizer checkpoints use complete SHA-256 fingerprints including objective
  and search inputs. Replicate outcomes are retained, scalar uncertainty is
  computed per replicate, paired shared-seed bootstrapping preserves objective
  covariance, and exhausted candidate pools stop explicitly.
- Diversity defaults off and requires verified marker--QTL disjointness unless
  an experimental override is recorded.
- Desired-gain families no longer invent a different aggregate merit for every
  direction. Merit-dependent criteria require common `aggregate_weights`.
- Multi-environment stability contrasts now feed `restricted_index()` through
  `constraint_matrix`; duplicate genotype--environment records are rejected.
- The default covariance repair preserves trait variances by repairing the
  correlation matrix. The aggressive ASRgenomics route remains explicit.

## Technical review of 31 July 2026 — remaining items

### New: covariance uncertainty reaches the recommendation

`optimize_desired_gains()` conditions on one covariance estimate, so a
direction that is best only because of how `G` happened to be estimated looked
identical to one that is robustly best.

`propagate_covariance_uncertainty()` re-evaluates a fixed set of directions
across draws from the sampling distribution of the covariance.
`draw_covariance_pairs()` draws the genetic and residual matrices separately
and reassembles `P* = G* + E*`, so every pair is admissible by construction
rather than by rejection.

Each draw **rebuilds the founder population**, because `G` is the architecture
the traits are built from rather than a matrix the simulation multiplies by.
That is what makes this expensive and why directions are supplied rather than
searched: the intended workflow is to take the Pareto set from
`optimize_desired_gains()` and pass it here.

The result separates the two sources of uncertainty. Monte Carlo error falls
with more replicates; covariance-estimation error does not and can only be
reduced by estimating `G` from more independent genetic units. When the
covariance component dominates, `variance_components$dominant_source` says so,
because running more replicates in that situation is wasted effort.

### Fixed: clonal calibration now reproduces the whole covariance

The previous release calibrated the genotypic *variances* and left the
correlations emergent, so the setup accepted a full covariance matrix as its
target and reproduced only its diagonal. Dominance perturbs the correlation
structure as well, and that perturbation is a smooth, near-identity function of
the additive correlations supplied, so the fixed-point iteration now corrects
both and projects the correction back onto the set of valid correlation
matrices. `calibration_error` records the largest remaining deviation in each
of the variances and correlations, and non-convergence warns rather than
passing silently.

### New: multi-environment and genotype-by-environment structure

`expand_environments()` builds the separable covariance
`Cov(g_je, g_kf) = G_jk * C_ef` over trait-environment combinations, so the
whole index layer applies to a multi-environment problem unchanged. It supports
environment-specific economic weights, and `include_stability = TRUE` returns
per-trait response contrasts that can be passed to `restricted_index()`.

The residual is expanded separately and `P` reassembled from the two, so the
expanded pair remains admissible; expanding `P` directly with the genetic
environmental correlation would not guarantee that. The reported
`interaction_share` says whether modelling the environments separately was
worth the complexity. `widen_environments()` reshapes long-format trial data.

### New: validation against the QGSI reference implementation

`run_qgsi()` implements Cerón-Rojas et al. (2026) and was previously checked
only against the equations as this package restates them, which establishes
that the code matches our reading and nothing more.

That paper's replication deposit
([doi:10.71682/10549385](https://doi.org/10.71682/10549385)) contains the
authors' own R scripts. The intermediate data files those scripts read were not
deposited, so their published figures cannot be recomputed — but the algorithm
can, and a reference implementation is a stronger comparison than a table,
because it runs on any input rather than only on theirs.

`tests/testthat/test-reference-implementation.R` transcribes their code
verbatim and runs both implementations on identical inputs. Compared: the
Smith–Hazel coefficients `b = P^-1 G w`, the accuracy `corrHI`, the aggregate
response `Rs`, the quadratic index variance
`Vpq = b'Pb + 2 tr(BB P BB P)` and the quadratic merit variance
`VHq = w'Gw + 2 tr(A G A G)`.

The accuracy comparison is the informative one. The reference computes a ratio
of standard deviations and this package computes a correlation; they coincide
only because `b'Gw = b'Pb` for the optimum index, an identity asserted as its
own test. Agreement establishes that both parties mean the same thing by
accuracy, which comparing coefficients alone would not.

### New: reproduction against published results

`vignette("DesiredGainR-reproduction")` checks the package against tables
computed in SAS PROC IML and published by Rahimi and Debnath (2023),
*Scientific Reports* 13:18977. The reference values are frozen in
`rahimi_debnath_2023()`.

The article's covariance matrices are not distributed, so coefficients cannot
be recomputed. What is checkable is stronger: Pesek–Baker gives expected
response *exactly* proportional to the desired gains, and dividing the
published gains by the published `d` vector gives the same constant, 0.3358, to
four significant figures across all seven traits. Seven ratios cannot coincide
by accident.

The vignette also explains the article's reported `R_HI = 0.0018` for
Pesek–Baker against 0.9887 for the optimum index. It is not an arithmetic
error: it follows from using the desired-gain vector as the aggregate weights,
which is what this package did until 0.5.0, and collapses because the trait
scales span a factor of 200.

### Release infrastructure

`cran-comments.md`, `.zenodo.json` for DOI minting on tag, `paper/paper.md`
with a bibliography for the method paper, and
`data-raw/release_checklist.R`, which reports every outstanding release step
rather than stopping at the first.

## Technical review of 31 July 2026 — Gates B and C

### New: restricted and proportional-gain indices

`restricted_index()` adds the closed-form constrained families:
Kempthorne–Nordskog, Tallis, Mallard, Harville, and the projection form of
restricted Smith–Hazel. Where a constraint is linear these give exactly what
`run_dgsi()` approaches by stochastic search, without the replicate dependence
or the selection optimism, so they also serve as a convergence check on the
search.

### New: upstream adapters, `summary()` and `coef()`

`import_covariance()` converts variance components and covariance matrices from
ASReml-R, sommer, breedR, BGLR and StageWise. It reorders traits to the
analysis order rather than assuming it, refuses an incomplete component set,
detects a correlation matrix passed as a covariance, and checks the imported
matrix against `P`. A silently transposed trait order produces an index that is
wrong in a way nothing downstream can detect, which is what this exists to
prevent.

`summary()` and `coef()` methods are added for fitted indices and for
`run_dgsi()` results, returning tidy tables rather than printed text.

### Fixed: best-of-replicate optimism is now removed, not just reported

`run_dgsi()` gains `replicate_selection`, defaulting to `"holdout"`: each
replicate's coefficients are scored on candidates held out of that comparison,
so the choice no longer depends on the candidates finally selected. The
previous behaviour remains as `replicate_selection = "training"` and is
labelled as biased.

### Fixed: the estimated-`P` exception is part of the API

`run_dgsi()` enforced admissibility for a supplied `P` but silently permitted an
inadmissible estimated one. It now stops unless
`allow_incompatible_estimated_P = TRUE`, and records the smallest residual
eigenvalue, the sampling threshold, the status and any override in
`covariance_provenance$compatibility`.

### Fixed: optimiser checkpoints identified only trait names

A checkpoint could resume evaluations produced under different founders, a
different genetic covariance, a different mating system or a different cycle
count and present them as this run's frontier. Checkpoints now carry a
fingerprint of every result-defining input, including package and AlphaSimR
versions, and a mismatch stops with the differing fields named rather than
resuming.

### Fixed: simulation error was computed and discarded

`.dgr_evaluate_direction()` returned Monte Carlo standard errors that the
optimiser threw away, so the surrogate fitted one homogeneous nugget and
treated a precisely evaluated direction and a noisy one alike. Replicate-level
outcomes, standard errors and their covariance are now retained and supplied to
the Gaussian process as a per-point nugget.

The frontier reports `pareto_probability`, the proportion of resampled
frontiers on which each direction remained non-dominated given its own Monte
Carlo error. This represents *simulation* error only; `G` and `P` are still
treated as fixed, and that limitation is stated in the result.

### Fixed: the optimiser could spend its budget re-evaluating one direction

Expected improvement at an already-observed point can stay positive under a
fitted nugget. Evaluated directions are now excluded from the acquisition pool,
and a repeat that does occur is pooled as additional replication with a
correspondingly reduced standard error rather than recorded as a new point.

### Changed: clonal replication

`simulate_selection_cycles()` gains `n_clonal_replicates`. A clonal trial
phenotypes several ramets per genotype and averages them, so selection is more
accurate than a single plot; the previous single-phenotype assumption
understated clonal response.

### Other

`MATURITY.md` states which components are stable, beta and experimental, and
what the package deliberately does not do. CI now covers R-devel, oldrel and
the declared minimum 4.1, treats warnings as failures, checks the source
tarball, verifies that roxygen output is committed, renders the Breeder Guide
as a build step, and runs one job with every optional dependency installed that
**fails if any test skips**. The legacy `inst/extdata` example CSVs are
excluded from the build.

## Technical review of 31 July 2026 — Gate A

### Fixed: DGSI desired gains were in the wrong space when traits were unscaled

`run_dgsi()` documents `dg` and `realised_response` in candidate
standard-deviation units, and its objective follows that definition. The Yamada
coefficient solve, however, operates in the units of `G` and `P`. With
`scale_traits = TRUE` the two coincide; with `scale_traits = FALSE` the
analysis space is raw trait units, so the standard-deviation vector was asking
for \eqn{d_j} *raw* units per trait.

For two traits with candidate standard deviations 10 and 1 and equal requested
gains, this targeted standardised response in the ratio 1:10 rather than 1:1.
The stochastic search could partly recover, since the map is invertible, but
the proposal coordinates were mislabelled, `optimised_d` was not in the
documented units, and the non-iterated Yamada comparator was simply wrong.

Desired gains are now converted to the analysis space before the solve, and the
index is invariant to the units the traits are recorded in. New tests assert
that invariance directly. **Re-run any `run_dgsi()` analysis that used
`scale_traits = FALSE` on traits with dissimilar scales.**

### New: theoretical transmitted response

`run_dgsi()` now reports `theoretical_response`, the classical
\eqn{i\,\mathbf{G}\mathbf{b}/\sqrt{\mathbf{b}'\mathbf{P}\mathbf{b}}}, in the
analysis space, in original trait units and standardised. `realised_response`
is a differential among the candidates that were selected; this is what the
next generation inherits. Keeping them apart is what allows DGSI to be compared
with the classical families on the same criterion.

### Fixed: release blockers

`index_uncertainty()`'s roxygen contained a raw `$$…$$` display that was copied
literally into the Rd file, where `\mathbf`, `\sim`, `\mathcal` and `\nu` parsed
as unknown macros. This was the sole source of both `R CMD check` warnings. It
is now `\deqn{}` with an ASCII fallback.

`inst/CITATION` reported version 0.3.1 while `DESCRIPTION` and `CITATION.cff`
said 0.5.0. It now derives the version from package metadata so it cannot drift
again.

`open_desiredgain_guide()` defaults to `"html"` rather than `"pdf"`. The PDF
edition needs LaTeX as well as pandoc, so defaulting to it meant the export
failed even where the guide had been built. A missing edition now names the one
that is present.

### Fixed: provenance that claimed more than it verified

The simulation recorded the *requested* thread count even when the assignment
failed on an older AlphaSimR, so a multithreaded run could carry a
single-thread reproducibility claim. It now reads back the effective count,
records both, and warns on a mismatch.

`.dgr_diversity_geno()` intersected a stored marker panel with whatever was
available and fell back to the full genotype matrix when nothing matched —
silently changing the basis of the diversity measure, and reintroducing the
QTL the panel exists to exclude. The panel is now required exactly, ordered by
identifier, and retained loci are stored by name rather than as a positional
logical vector.

### Changed: diversity optimisation requires a verified panel

`optimize_desired_gains(include_diversity = TRUE)` now errors when the setup
has no marker panel, rather than silently falling back to all segregating
sites. Without a panel the diversity axis partly measures genetic gain itself,
so a direction is penalised for working. Where the panel cannot be verified as
disjoint from the QTL on the installed AlphaSimR, it warns.

### Other

Seeds route through `.dgr_seed()` in `run_dgsi()` and
`optimize_desired_gains()` rather than `as.integer()`, so fractional seeds no
longer truncate silently, and a base seed too close to the integer limit for
its derived replicate seeds is rejected. `ridge_P`, `ridge_M` and
`plateau_tolerance` are validated for length, finiteness and sign.

## External review, August 2026 — release blockers

### Breaking: what a simulation cycle means has changed

`simulate_selection_cycles()` measured the *candidates* of each cycle rather
than the population its selection produced. Because random mating does not
shift a population mean, the cycle 1 gain was zero in expectation **for every
desired-gain direction**, so a one-cycle run could not distinguish directions
at all and `optimize_desired_gains(n_cycles = 1)` was comparing noise.

Cycle 0 is now the founder population, and cycle *t* records the population
produced by selecting parents from the cycle *t-1* candidates and crossing
them. Every row from cycle 1 onward is a transmitted selection response.
`cycle_table` gains a cycle 0 row, `parent_inbreeding`, `n_parents_selected`,
`qtl_segregating` and `genic_variance_ratio`. **Re-run any stored simulation
result.**

Genetic mean and diversity are now measured on the same population; previously
merit came from the whole evaluation population and relatedness from the
selected parents.

### Breaking: independent culling is a hard gate again

`selection_index(method = "independent_culling")` ranked candidates on a 0/1
pass indicator and then took the top `n_select`, so requesting 20 when 7 passed
selected 13 failures. `n_select` is now a maximum: 7 pass, 7 are returned, with
a warning. Selection intensity is computed from the number actually selected,
and `culling_report` names the gates each candidate failed. New `n_selected`
records the executed count; `n_select` still holds the request.

### New: covariance compatibility is validated

`P - G` must be positive semidefinite — a diagonal check is not enough, since a
linear combination of traits can exceed its phenotypic variance while no single
trait does. Enforced in `selection_index()`, `evaluate_index()`, `run_dgsi()`
and `index_uncertainty()`, and `true_G - Gamma` in `run_qgsi()` when accuracy
or MSPE is requested. Correlations are clamped only for floating-point
excursions, after validation. No valid run can now report a heritability above
one, an accuracy above one, a squared correlation above one, or a negative
MSPE.

`index_uncertainty()` with `P` fixed now discards draws where the drawn `G`
exceeds it, counts them in `n_inadmissible`, and describes itself as local
sensitivity analysis rather than joint resampling.

### Fixed: simulations mutated the caller's SimParam

`SimParam` is an R6 object, so every crossing call advanced the caller's
`lastId` and two runs from one setup were not independent. It is now deep-cloned
per call. `n_threads` defaults to 1 and warns above it, since thread count
affects the order random numbers are consumed. The result carries a
`provenance` record with the AlphaSimR, R and package versions, seed, thread
count and simulation settings.

### Fixed: diversity was measured on the loci under selection

Relatedness used `pullQtlGeno()`. Changing QTL frequencies is what genetic gain
*is*, so the metric rose whenever selection succeeded and a direction was
penalised for working. `founder_population()` gains
`n_markers_per_chromosome`, creating a neutral panel held disjoint from the
QTL, and diversity is measured there. QTL frequency change is reported
separately as `qtl_segregating` and `genic_variance_ratio`.

### Fixed: clonal mode was miscalibrated

AlphaSimR's `var` argument sets *additive* variance, so with dominance present
the dominance variance was added on top and the realised genotypic variance
exceeded the supplied `G` — the matrix a clonal programme selects on. Since
scaling QTL effects by *c* scales both components by *c²*, the correction is
exact in one pass and is now applied. `founder_population()` gains
`heritability = c("narrow", "broad")`, because the two diverge exactly when
dominance is simulated. `G_realised` records the achieved covariance; realised
correlations are emergent and warn above a 0.15 deviation. The sexual and
clonal phases are now explicit in `simulate_selection_cycles()`.

### Fixed: counts and defaults

`run_dgsi(n_select)` bypassed the strict validator, so `1.9` silently became
`1`. All counts now route through it; `NA`, `Inf`, zero, negative,
non-integral and out-of-range values fail with distinct messages.
`run_dgsi_qgsi_pipeline()` now defaults to `fallback_to_top_n = FALSE` and
`debug = FALSE`, `run_qgsi_desired_gain()` to `debug = FALSE`. Enabling the
fallback warns at request and again if it fires, and `$selection_rule` records
requested versus executed.

### Other

`evaluate_index()` reports `h2_index` and `accuracy_index`. `run_dgsi()` gains
a `print()` method and the additional class `desiredgainr_dgsi` alongside
`desired_gain_index`, which is retained for dependent packages. The Breeder's
Guide source ships in `inst/guide/` with a build script in `data-raw/`, and
tests verify the artifacts resolve through `system.file()`.

## Audit of 31 July 2026

`AUDIT-2026-07-31.md` records the full findings. The changes below follow from
it.

### Fixed: the diversity metric could not detect diversity loss

`.dgr_mean_relationship()` recomputed allele frequencies from the population
being measured. Centring on a sample's own frequencies forces every marker
column to sum to zero, so the whole genomic relationship matrix sums to zero and
its mean off-diagonal element equals $-1/(n-1)$ regardless of what has happened
to the germplasm. The reported value was a function of the number of parents and
almost nothing else.

This mattered beyond the reporting, because
`optimize_desired_gains(include_diversity = TRUE)` used it as the diversity axis
of a Pareto frontier. Any direction chosen for its diversity properties was
chosen on noise.

Allele frequencies are now fixed at the founder population and reused for every
subsequent measurement. `simulate_selection_cycles()` gains `mean_inbreeding`
per cycle and `final_inbreeding`, and `effective_size` is now derived from the
rate of inbreeding rather than from differencing the old metric.

**If you have run `optimize_desired_gains(include_diversity = TRUE)` on an
earlier version, re-run it.** The gain axis was unaffected; the diversity axis
was not meaningful.

### Fixed: `R_HI` for the desired-gain families

`evaluate_index()` was passed the desired-gain vector as the aggregate weights
for `pesek_baker` and `yamada`. There is no theory under which desired gains
define net merit, and the resulting quantity behaved erratically when trait
scales differed. The implied economic weights
$\mathbf{w} = \mathbf{G}^{-1}\mathbf{P}\mathbf{G}^{-1}\mathbf{d}$ are now used
instead, so $R_{HI}$ recovers its standard interpretation and the index families
become comparable on one criterion.

### Fixed: the reference set for estimating `P` was unguarded

`run_dgsi()` estimated `P` from `ref_data` with no check beyond two records. A
covariance matrix from fewer records than traits is singular, passes the
positive-semidefinite check unchanged, and is then made invertible by the ridge
term, so the failure was silent. It now stops below `length(trait_cols) + 1`
records, warns below five times the trait count, and records
`P_numerical_rank` in `covariance_provenance`.

### New: `bend_covariance()`

Multi-trait restricted maximum likelihood routinely returns a genetic covariance
matrix that is not positive definite. The package previously rejected it
outright, which sent users to an ad hoc repair outside the audit trail.
`.dgr_check_psd()` now names the repair in its rejection message.

`bend_covariance()` is a wrapper around `ASRgenomics::G.tuneup(bend = TRUE)`,
which bends by calling `Matrix::nearPD()` — Higham's (2002)
alternating-projections algorithm. That is the only bending route the package
offers, deliberately. A package that repairs covariance matrices two different
ways invites two different answers to the same question, and there is no
defensible basis on which a user would choose between them mid-analysis. A test
asserts that the wrapper returns exactly what `G.tuneup()` returns, so a result
obtained here is reproducible by anyone running the underlying function
directly.

Two of `G.tuneup()`'s fixed choices have real consequences for a trait
covariance, and the result surfaces both rather than letting them pass
unnoticed.

`keepDiag = FALSE`, so **the trait variances are not preserved**. The
heritabilities implied by the bent matrix are not the ones supplied.
`adjustment$max_variance_change` and `max_relative_variance_change` report how
far they moved, and the print method shows them.

`posd.tol = 1e-2`, which floors eigenvalues at a hundredth of the largest and
so caps the condition number near 100. That is the right target for the
\eqn{n \times n} genomic relationship matrix `G.tuneup()` was written for, which
is near-singular by construction. For a \eqn{p \times p} trait covariance it is
aggressive: condition numbers of \eqn{10^3} to \eqn{10^4} are ordinary and
genuine there, so bending will smooth away correlation structure the data
support. `eig_tol` is passed through but governs `nearPD()`'s `eig.tol`, not
`posd.tol`, so it does not lift that cap. Compare [matrix_diagnostics()] before
and after.

`ASRgenomics` is in `Suggests` rather than `Imports`, because it brings in
roughly ten further packages for plotting and pedigree handling and a single
function should not be able to make the whole package uninstallable.
`bend_covariance()` stops with an install instruction when it is absent.

### New: `index_uncertainty()`

Every quantity a selection index reports has treated `G` and `P` as known
constants. `index_uncertainty()` resamples them from a Wishart sampling
distribution, refits the index on each draw, and reports coefficient intervals,
sign stability, rank correlation, and retention of the selected set.

The genetic and residual matrices are resampled, and `P` reassembled as their
sum, rather than resampling `G` and `P` independently, because they are
estimated from the same records.

For the desired-gain families the reported cost of estimation error has an exact
reference point: since $\mathbf{G}\mathbf{b} = \mathbf{d}$ identically, an index
fitted on the true `G` delivers the desired gains in the requested proportions,
and the cosine between the achieved response and `d` isolates estimation error
with no confounding.

`genetic_df` is the user's judgement and governs every interval width. It counts
independent genetic units — families, parents, distinct clones — and not plots.

### New: `predict()` for a fitted index

Applies a fitted index to a new candidate set, reusing the centring and scaling
constants from fitting rather than recomputing them, so that scores remain
comparable across cycles.

### New: index heritability

`evaluate_index()` now reports `h2_index` ($\mathbf{b}'\mathbf{G}\mathbf{b} /
\mathbf{b}'\mathbf{P}\mathbf{b}$) and `accuracy_index`, its square root. This is
the only accuracy measure available for the desired-gain families before
economic weights are supplied.

### New: `scale_by` in `selection_index()`

`scale_traits` was documented as dividing by population standard deviations but
divided by the candidates' standard deviations. Those coincide only when the
candidates are an unselected random sample, which a late-stage trial is not.
`scale_by = "phenotypic"` divides by $\sqrt{\operatorname{diag}(\mathbf{P})}$
instead, giving a scaled `P` with a unit diagonal exactly and a scaled `G`
carrying the heritabilities. The default remains `"sample"`, so existing calls
are unchanged.

### Changed: the surrogate's kernel

The Gaussian process in `optimize_desired_gains()` fitted one lengthscale per
trait, estimated by marginal likelihood from however many objective evaluations
the budget allowed. Below roughly ten evaluations per trait those estimates are
dominated by noise, which then steers the acquisition function. A single shared
lengthscale is now used in that regime.

### Reported optimism in `run_dgsi()`

At that stage the winning replicate was chosen on the same candidates it then
selected, so its objective was biased downward. The result gained an
`optimism` element; the 1 August remediation above subsequently made a pre-fit
holdout the default while retaining the training rule explicitly.

## Example data

Nine datasets replace the previous examples, which were uncorrelated standard
normal variates. That earlier material made the genetic covariance matrix
effectively an identity matrix, and an identity matrix cannot illustrate a
multi-trait selection index, because the purpose of an index is to resolve the
tension between correlated traits measured on different scales.

The new data describe one simulated tropical maize programme with six traits:
grain yield, plant height, anthesis date, the anthesis-silking interval, ears
per plant, and grey leaf spot severity. They carry three properties the earlier
examples lacked.

1. **Genuine antagonism.** Grain yield correlates positively with anthesis
   date, yet the objective requires yield to rise while the cycle shortens, so
   the two cannot both be pushed freely. Yield is also strongly and negatively
   correlated with the anthesis-silking interval. Consequently
   `gain_feasibility()` returns informative answers rather than declaring every
   target attainable.
2. **Realistic scale differences.** Genetic standard deviations span a factor
   of roughly one hundred and twenty, from ears per plant to plant height,
   which is what gives the standardisation and conditioning diagnostics
   something to detect.
3. **Internal consistency.** The markers explain the trait values. Genetic
   values are generated from ninety quantitative trait loci drawn from the
   marker panel, then whitened and recoloured so that their realised covariance
   equals `dgr_G` to numerical precision. Heritabilities are recoverable
   exactly as `diag(dgr_G) / diag(dgr_P)`.

The datasets are `dgr_traits`, `dgr_G`, `dgr_P`, `dgr_candidates`, `dgr_gebv`,
`dgr_history`, `dgr_hap1`, `dgr_hap2` and `dgr_map`. The haplotype matrices and
map satisfy the `founder_haplotypes()` contract directly, and `dgr_history`
records a previous cycle's decision produced by an undisclosed weight vector,
so that `retrospective_weights()` has something real to recover.

`tests/testthat/test-example-data.R` asserts each of these properties, since
data that quietly stopped demonstrating them would make every example
misleading.

The demonstration script `inst/examples/demo_pipeline.R` has been rewritten
around the new data and now runs the whole package in the order the decisions
actually arise: inspect the covariance matrices, translate between weights and
desired gains, test feasibility, recover the historical objective, compare
index families, check the effective weights, measure sensitivity to the stated
weights, fit the optimised desired-gain index, and finally score genomic
estimated breeding values with the quadratic index.

The former `inst/extdata/example_pheno.csv` and `inst/extdata/example_gebv.csv`
are superseded and may be deleted.

## Corrections found by running the demonstration

- **`selection_index()` returned the wrong candidates in `$selected`.** The
  ranking table is sorted by score after the selection flag is computed, and
  the returned subset was taken by filtering that sorted table with the
  unsorted flag. The rows were therefore chosen positionally rather than by
  merit. The `selected` column of `$ranking` was always correct, as were the
  observed selection differentials, so the fault was confined to the
  `$selected` table itself — which is nonetheless the element most users read
  first.

  Nothing in a suite of 356 tests detected it, because the object remained
  internally plausible: the right number of rows, the right columns, and
  consistent values. It surfaced only from noticing that two numbers printed by
  the demonstration script could not both be true, a Spearman correlation of
  0.938 between two indices alongside an overlap of one candidate in twenty
  between their selected sets. Tests now assert the relationship between those
  two quantities, and assert that `$selected` matches an independently
  recomputed top-n.

Four defects survived a passing test suite and were exposed only by reading
the output of `demo_pipeline.R` on the new data. Each is a case where the code
did what it was told and the result was nonetheless misleading.

- **The index coefficient of variation diverged on a centred index.**
  Standardising the traits places the index mean at zero, so `sd / mean`
  reported values of the order of 10^17. The quantity is undefined there, and
  it is now withheld with `CV_I_note` explaining why.

- **The objective-setting functions had no `lower_is_better`.** Every other
  entry point in the package accepts it, so `implied_economic_weights()`,
  `implied_desired_gains()` and `gain_feasibility()` silently interpreted their
  inputs in the raw trait direction. A positive desired gain for plant height
  was therefore a request to grow taller. The author of the demonstration
  script made exactly this mistake. All three now accept `lower_is_better` and
  orient the covariance matrices internally.

- **`implied_desired_gains()` could not return the units it was given.** A
  round trip through `implied_economic_weights(gain_units = "genetic_sd")` came
  back in trait units, which looks like a failed inversion. It now takes
  `gain_units`, and a test asserts the round trip is exact in all three.

- **`run_qgsi()` gave no warning when trait scales were mismatched.** Economic
  weights multiply genomic estimated breeding values directly, so on the
  original scale a trait with a large standard deviation dominates whatever
  weight it was given. In the demonstration, plant height overwhelmed the index
  and the expected gain in grain yield came out *negative* while selecting for
  higher yield. `run_qgsi()` now issues the same conditional scale warning that
  `run_dgsi()` already carried, naming the dominating trait and its share.

- **`selection_index()` gains `center_traits`, separate from `scale_traits`.**
  Centring was previously unconditional, which placed the index mean at zero
  in every call and therefore made the coefficient of variation permanently
  undefined. The two transformations are now controlled independently, as they
  already were in `run_qgsi()`, so a published result reporting `CV_I` on an
  index built from raw trait values can be reproduced. Centring never changes
  a ranking, since it shifts every score by the same constant.

- **Units are now reported in the units the caller used.** `gain_feasibility()`
  returned the attainable response in original trait units regardless of the
  `gain_units` the target was stated in, so a request of one genetic standard
  deviation of yield came back as 0.41 tonnes per hectare. Comparing the two
  invites misreading. `attainable_response_input_units` and
  `shortfall_input_units` are now returned alongside the trait-unit versions.

- **Coefficients are reported on a comparable footing.**
  `retrospective_weights()` returns `coefficients_per_sd`, and the print method
  shows that rather than the raw vector. A coefficient carries inverse trait
  units, so ears per plant, whose standard deviation is 0.10, received the
  largest raw coefficient in the example programme while representing a modest
  emphasis.

- `run_qgsi()` no longer repeats the citation for the expected-gain equation on
  every row of `expected_gain_per_trait`, where six identical sentences swamped
  the printed table. It is carried as a `basis` attribute on that table and as
  `expected_gain_basis` in the result.

- Minor: `gain_feasibility()` no longer reports a required selection of one
  candidate when the required proportion is smaller than a single candidate,
  since rounding up overstated what the population could deliver. A required
  intensity beyond the reach of normal truncation is reported as such rather
  than converted into a proportion.

## Corrections found by the test suite

- **The two desired-gain formulations are one index, not two.** For a square
  invertible `G`, the Yamada form reduces algebraically to the Pesek-Baker
  form, since `P⁻¹G(GP⁻¹G)⁻¹d = P⁻¹G(G⁻¹PG⁻¹)d = G⁻¹d`. Both `method` values
  are retained as documented routes to one estimand. They remain worth
  distinguishing, however, because only the Yamada route is available when `G`
  is rank deficient, and the two are not numerically equivalent under an
  ill-conditioned matrix. The documentation previously claimed they gave
  different coefficients and rankings; that claim was wrong.

- **Expected per-trait gains rest on two approximations, not one.** The
  identity `Cov(gamma, I) = Gamma w` is exact, because the third central
  moments of a zero-mean multivariate normal vanish. Two further steps are not.
  First, the achieved selection differential on the index is assumed to equal
  the normal-theory value `i * sd(I)`; the index is a quadratic form and is
  right skewed, so the differential actually achieved is larger. Second,
  `E[gamma | I]` is assumed linear in `I`, which also fails for a quadratic
  index.

  Monte Carlo testing shows these two errors partly cancel. Correcting only the
  first, by substituting the observed differential, makes the prediction worse
  rather than better. The total index standard deviation is therefore used
  because it is the correct linear-regression denominator, and not because it
  demonstrably fits simulation best: at modest curvature the residual error of
  five to nine per cent exceeds the two per cent difference between the two
  candidate denominators, so simulation cannot arbitrate between them. Both
  assumptions are now reported alongside the gains.

- `simulate_selection_cycles()` tracked genetic gain with `bv()`, whose values
  are expressed relative to the supplied population and therefore have a mean
  of zero in every cycle. Cumulative gain was consequently always zero. Gain is
  now tracked with `gv()`, which is anchored by the trait mean set when the
  trait was added. Covariance estimation for rebuilding the index continues to
  use breeding values for additive programmes, since centring does not affect a
  covariance.

- `effective_weights()` failed on a rank-deficient covariance matrix, where the
  shares are undefined, and that failure propagated out of `run_dgsi()`.
  Shares are now returned as `NA` in that case, and `run_dgsi()` treats the
  diagnostic as non-fatal.

- `effective_weights()` no longer warns by default, and its threshold adapts to
  the number of traits. A concentrated effective weight can arise from a
  deliberately asymmetric objective as readily as from mismatched trait scales,
  so warning on every asymmetry was noise. The scale mismatch that Crosbie et
  al. (1980) describe is still reported, by the separate `run_dgsi()` warning
  that fires only when trait standard deviations differ by more than fivefold.

This release begins the multi-cycle simulation layer. Its purpose is narrow and
deliberately kept so: it exists to compare desired-gain directions over several
breeding cycles, because a direction that maximises response in the first cycle
may exhaust the genetic variation that response depends on within a few more.
It is not a crossing-plan tool, and parent, cross and mating allocation remain
outside the scope of this package.

## Founders come from real phased marker data

DesiredGainR does not simulate founder genomes. Therefore the linkage
disequilibrium, allele frequency spectrum and population structure of every
simulation are those of the breeder's own germplasm rather than those of an
assumed demography, which removes the principal scientific weakness of
simulation-based recommendations.

- `founder_haplotypes()` validates and packages phased founder haplotypes. It
  accepts variant-by-individual matrices `hap1` and `hap2` of 0 and 1, together
  with a variant map giving chromosome and genetic position. This is the shape
  returned by phased variant-call-format readers, including
  `HapBlockR::read_phased_vcf()`.

  Genotype dosage in 0, 1, 2 coding is rejected, and the error explains why. A
  haplotype is one chromosome copy, so at a biallelic locus it carries a single
  allele rather than an allele count; dosage is the sum of the homologues and
  therefore records how many alternative alleles an individual carries but not
  which copy carries them. Assigning heterozygous calls at random would destroy
  the linkage disequilibrium that determines how favourable alleles segregate
  together. Hence phase must be supplied rather than inferred.

- Missing calls are permitted on input and resolved under an explicit
  `missing_policy` of `"error"`, `"drop_variant"` or `"drop_individual"`.
  Phased variant call format output legitimately carries missing values,
  whereas AlphaSimR requires complete haplotypes, so they must be resolved
  rather than imputed silently. The rate observed before resolution, and
  everything removed, are recorded in the returned object.

- The contract is deliberately **diploid and biallelic**, and a polyploid
  dataset raises an explanatory error. Extending the two homologue matrices to
  further copies would not make the simulation polyploid, because
  autopolyploid and allopolyploid meiosis require an explicit crop-specific
  pairing and recombination model, including multivalent formation and double
  reduction. Applying a diploid meiotic model to polyploid homologues would
  produce confidently wrong results, so polyploid support must wait for a
  dedicated homologue-level representation.

- `dosage_diagnostics()` reports overall, per-individual and per-marker
  heterozygosity and missingness. Run it before any conversion, because a low
  overall rate can conceal a small number of badly affected individuals or
  markers, and that distribution determines which conversion policy is
  cheaper.

- `haplotypes_from_inbred_dosage()` derives phase from dosage for highly inbred
  diploid material. A dosage of 0 or 2 is unambiguous, so doubled haploids,
  recombinant inbred lines and advanced selfed generations can be converted
  without external phasing.

  No threshold is imposed on residual heterozygosity, because no universally
  appropriate level exists: it depends on the generation, mating history, crop,
  population type, genotyping error rate and quality-control procedure, and
  must therefore be measured from each dataset. Heterozygous calls are never
  resolved silently and are never assigned to a homologue. The available
  policies are (i) `"error"`, (ii) `"drop_variant"`, (iii)
  `"drop_individual"`, and (iv) `"mask"`, which sets the call to missing for
  resolution by `founder_haplotypes()`. Every rate observed and every call
  removed or masked is recorded in the returned `conversion` attribute.

- `founder_population()` builds the AlphaSimR founder population and calibrates
  the trait architecture to the genetic covariance matrix the breeder has
  already estimated: variances from the diagonal of `G` and correlations from
  the corresponding correlation matrix.

## Recurrent selection across cycles

- `simulate_selection_cycles()` runs a recurrent-selection programme forward
  under a fixed desired-gain direction and records, for every cycle, the
  genetic mean and variance per trait, the mean relationship among selected
  parents, and the implied effective population size.

  Three mating systems are supported: (i) `"self"` for self-pollinated line
  development, advancing by selfing or doubled haploidy, (ii) `"outcross"` for
  recurrent selection in a random-mating population, and (iii) `"clonal"` for
  clonally propagated crops, where selection acts on total genetic value
  because the clone inherits dominance intact. The clonal system therefore
  requires a setup built with `dominance_degree`, and `G` should then be the
  genotypic rather than the additive covariance.

  `reestimate_index` controls whether the index is rebuilt from each cycle's
  own simulated data, which is what a breeding programme actually does and
  which propagates estimation error across cycles, or held fixed at the cycle
  one solution, which isolates the effect of the desired-gain direction itself.
  A large divergence between the two settings means the recommendation is
  sensitive to estimation error, and the breeder should be told so.

## Searching for the best desired-gain direction

- `optimize_desired_gains()` searches for the desired-gain direction that
  performs best once several cycles have been simulated. Only the direction is
  optimised, because the attainable magnitude is fixed by selection intensity,
  so the domain is the unit sphere with one fewer free parameter than there are
  traits. That is small enough to cover densely for the trait numbers breeders
  use.

  **The simulation is the objective function and is never replaced by a cheaper
  approximation.** A Gaussian process is fitted to the accumulated results, but
  it decides only where the next simulation should be spent; it never filters,
  screens, or substitutes for an evaluation, and because the acquisition
  function retains an exploration term across the whole domain, no region can
  be permanently excluded by the surrogate. An earlier design used a
  deterministic infinitesimal-model recursion as the surrogate's prior mean;
  that was rejected, because it would embed a second genetic model whose
  assumptions could fail in precisely the situations where simulation is most
  needed.

  Four ranking modes are available: (i) `"pareto"`, the default, returning the
  non-dominated set of multi-cycle outcomes without requiring economic weights,
  (ii) `"economic"`, maximising a weighted sum, (iii) `"target"`, minimising
  the distance to stated absolute targets, and (iv) `"constrained"`, maximising
  a focal trait subject to floors on the remainder.

  Two properties matter for interpretation. Every direction is evaluated with
  the same sequence of seeds, so comparisons between directions share their
  stochasticity, which removes more comparison variance than raising the
  replicate count. The reported frontier is computed from posterior means
  rather than raw replicate averages, because a frontier built from raw draws
  is populated by fortunate runs.

  No single best direction is returned for the scalar modes either. The optimum
  is conditional on the supplied genetic covariance, the founder germplasm and
  the programme parameters, so the result reports every direction whose
  posterior outcome lies within `stability_tolerance` of the best.

- Diversity, measured as the negated mean relationship among selected parents,
  is included as an additional Pareto objective by default. Family balancing
  and coancestry control can cost a large share of nominal gain, so the
  trade-off is reported rather than assumed away.

- Searches are resumable. Supplying `checkpoint` writes the accumulated
  evaluations after every simulation and reloads them automatically on a
  later call, because a realistic budget takes hours and an interruption
  should not discard the work.

## Dependencies

- `AlphaSimR` is added to `Suggests`. Every simulation function fails with an
  explicit installation message when it is absent, and the tests skip rather
  than fail. DesiredGainR does not depend on, import, or call HapBlockR, which
  keeps the dependency arrow one-directional.

# DesiredGainR 0.4.0

This release adds the objective-setting layer. Selection indices are
straightforward to compute; deciding what to select for is not, and the
literature identifies the specification of economic weights as the principal
obstacle to routine index use. Hence the new functions do not compute better
indices. They help a breeder state, test, and defend the objective an index is
built from.

Every addition is additive. No existing function, argument, result element, or
column name has been renamed or removed.

## Translating between the two ways of stating an objective

The Smith-Hazel economic index and the Pesek-Baker desired-gain index are two
parameterisations of one linear index. Therefore each desired-gain vector
corresponds to exactly one economic-weight vector.

- `implied_economic_weights()` returns the weights implied by desired gains,
  through `w = G^-1 P G^-1 d`.
- `implied_desired_gains()` performs the reverse translation, `d = G P^-1 G w`.
  This is the more useful direction in practice: a breeder who proposes weights
  can immediately see what response those weights actually request.

## Testing whether an objective is attainable

- `gain_feasibility()` reports the selection intensity a target response
  requires, `i = sqrt(d' G^-1 P G^-1 d)`, the selected proportion that delivers
  it, whether that proportion is reachable in the available population, and the
  per-trait shortfall when it is not. Covarrubias-Pazaran (2021) directs
  breeders to external software for this check; it is now available in R.

## Recovering an objective from past decisions

- `retrospective_weights()` recovers the linear rule implicit in a programme's
  historical selections, through `b = P^-1 s`. This is the method used to
  introduce a selection index at the International Maize and Wheat Improvement
  Center (CIMMYT) and the International Rice Research Institute (IRRI). The
  recovered weights are a starting point and reproduce past bias, so they
  should be inspected and adjusted before adoption.

## Quantifying how much the objective matters

- `weight_sensitivity()` perturbs the stated weights and reports how often the
  selected set survives unchanged. A high stability proportion means further
  effort refining the weights will not change the decision; a low one means the
  objective, rather than the index, is the binding uncertainty.
- `effective_weights()` reports each trait's contribution as `b_j * sigma_gj`
  rather than `b_j`, and warns when one trait dominates. Crosbie et al. (1980)
  showed that equal weights on unstandardised traits concentrate selection on
  whichever trait carries the largest genetic variance.
- `matrix_diagnostics()` reports the condition number of a covariance matrix.
  An index built by inverting an ill-conditioned matrix is numerically
  meaningless while remaining superficially plausible.

## Classical index families

- `selection_index()` builds the Smith-Hazel, base, Pesek-Baker, Yamada, and
  Mulamba-Mock indices, together with independent culling and tandem selection,
  from one interface, so that families and objectives can be compared on
  identical data. Guimaraes et al. (2021) found the Mulamba-Mock rank-sum
  index, which requires neither economic weights nor covariance matrices, to
  give the most balanced multi-trait response of the methods they compared.

  **Two indices share the name Pesek-Baker.** Pesek and Baker (1969) give
  `b = G^-1 d`, which is `method = "pesek_baker"` and what Rahimi and Debnath
  (2023) implement. Yamada et al. (1975) give `b = P^-1 G (G P^-1 G)^-1 d`,
  which is `method = "yamada"` and what Joukhadar et al. (2024) and `run_dgsi()`
  use. Both satisfy `G'b = d`, but they yield different coefficients and
  different rankings. Hence any comparison with published results must state
  which formulation was used.

## Evaluating an index

- `evaluate_index()` reports the four criteria used by Rahimi and Debnath
  (2023) and by the established selection-index software: (i) `R_HI`, the
  correlation between the index and net merit, (ii) `delta_H`, the expected
  gain in aggregate merit, (iii) `RE`, efficiency relative to direct selection
  on a main trait, and (iv) `CV_I`, the coefficient of variation of the index,
  together with the expected response for each trait.

  Relative efficiency below 1 is the intended trade-off of multi-trait
  selection, not a defect; every value reported by Rahimi and Debnath was below
  1.

## Other changes

- `estimate_genetic_covariance()` gains `method = "adjusted_means_surrogate"`,
  which returns the covariance of across-environment adjusted means labelled as
  a working surrogate for genetic covariance. Covarrubias-Pazaran (2021)
  recommends this for operational use when no fitted genetic covariance matrix
  is available. Earlier releases refused this workflow; it is now supported
  explicitly, with its assumptions recorded, rather than left to be performed
  silently outside the package.
- `run_dgsi()` now returns `non_iterated`, the classical Yamada solution using
  the supplied desired gains without iteration. Joukhadar et al. (2024)
  benchmark against exactly this comparator, and it costs nothing because it is
  the search starting point.
- `run_dgsi()` now returns `effective_weights` and warns when trait standard
  deviations differ by more than fivefold while `scale_traits = FALSE`. The
  warning is conditional on the scales actually differing, so analyses on
  comparable scales remain quiet.

# DesiredGainR 0.3.1

No function, argument, result-element, or column name was removed or renamed.
Existing code continues to run unchanged. Two corrections do change returned
*values*; both superseded quantities remain available under new names.

## Corrected results (numerical change)

- `run_qgsi()`: `expected_gain_per_trait$Expected_Genetic_Gain` now divides
  `Gamma %*% w` by the **total** index standard deviation,
  `sqrt(w'Gamma w + 2 tr(W Gamma W Gamma))`. Releases up to 0.3.0 divided by
  `sqrt(w'Gamma w)`, the purely linear (LGSI) index standard deviation, which
  overstated every per-trait gain whenever `W` was non-zero. The previous
  values are still returned in the new column
  `Expected_Genetic_Gain_LinearSD`. A Monte Carlo test against simulated
  truncation selection now pins the corrected formula
  (`tests/testthat/test-theory-simulation.R`).
- `run_qgsi()`: `theoretical_parameters$squared_index_merit_correlation` now
  uses the general definition `Cov(H, I)^2 / (Var(I) Var(H))`. Releases up to
  0.3.0 returned `Var(I) / Var(H)`, which equals the squared correlation only
  when the index is the MSPE-optimal predictor of net merit — not guaranteed,
  because `run_qgsi()` accepts arbitrary user weights. The former ratio is
  retained as `theoretical_parameters$variance_ratio_index_to_merit`.
  Unchanged when `true_G` is not supplied (both remain `NA`), and unchanged
  when `true_G == Gamma` (both equal 1).
- `run_dgsi()`: ties in the index score are now broken by ascending candidate
  identifier instead of by input row order, so the selected set no longer
  depends on how the rows were sorted. Only tied candidates are affected.

## Fixed

- `run_dgsi()` no longer leaves the caller's random number generator reseeded.
  The RNG state is saved on entry and restored on exit, so `seed` gives
  reproducible results without altering the user's stream.
- `run_qgsi()` now warns when `center_traits = FALSE` and `Gamma` is estimated
  internally. In that combination the scores use uncentred GEBVs while the
  estimated `Gamma` is a centred covariance, so the reported index mean,
  variances, and expected gains describe a different index from the one in
  `ranked_geno`. Supply `Gamma` explicitly, or use `center_traits = TRUE`.

## Documentation

- `run_dgsi()`: `dg` is now documented as being in **candidate
  standard-deviation units**, which is what the objective has always compared
  against, regardless of `scale_traits`. A new Details section gives the
  realised-response and objective formulas, explains the `0.25` and `0.05`
  floors that were previously undocumented constants, and states that
  best-of-`n_rep` training selection is optimistically biased. The current
  default avoids that rule with a pre-fit holdout.
- Added a worked `@examples` block to `run_dgsi()`.

## Internal

- `globalVariables()` is now declared once in `R/globals.R` instead of in two
  files with overlapping and partly obsolete entries.
- `.desiredgainr_z()` and `.validate_square_matrix()` are unused by any
  exported function. They are retained for now in case external code calls
  them via `:::`, and are marked for removal in a future release.
- `inst/CITATION` uses `bibentry()` instead of the deprecated `citEntry()`.
- `DESCRIPTION`: added `Language: en-GB` and `Depends: R (>= 4.1)`; dropped the
  hand-written `Author:`/`Maintainer:` fields, which duplicated `Authors@R` and
  could drift out of sync.
- `run_dgsi()` hoists the loop-invariant column means and standard deviations
  out of the search loop. Results are identical; the search is substantially
  faster on large candidate sets.
- Added `tests/testthat/test-theory-simulation.R`: Monte Carlo validation of
  the QGSI index mean, variance, per-trait gain and squared correlation, plus
  a Pesek-Baker consistency check, an RNG-restoration check, and a row-order
  invariance check.
- Added `tests/testthat/test-downstream-contract.R`, which pins the call
  surface and the result elements used by `hapblockr::build_selection_index()`
  (`coefficients`, `realised_response`, `ranked_geno$SelectionIndex`,
  `ranked_geno$QGSI`, `replicate_diagnostics$Chosen`, `component_summary`, the
  nine argument names that `dgsi_control` protects, and the
  `center_traits = TRUE` default). None of the 0.3.1 corrections touch these.

# DesiredGainR 0.3.0

- Removed the author-named DGSI wrapper; `run_dgsi()` is now the sole public
  desired-gain index entry point.
- Added `estimate_genetic_covariance()` with clearly distinguished prediction
  covariance, PEV-corrected, SE-diagonal-corrected, and
  relationship-adjusted estimators, together with assumptions and provenance.
- Implemented QGSI genomic covariance estimation from reference GEBVs using
  Supplementary Equations 19.1 and 19.2 of Cerón-Rojas et al. (2026).
- Added model-based QGSI mean, linear and quadratic variance, normal-selection
  intensity, expected net-merit response, and expected per-trait gains.
- Added optional simulation-only squared accuracy and mean squared prediction
  error calculations when the true genetic covariance matrix is supplied.
- Added exact top-ranked selection and clearly separated observed GEBV
  differentials from expected or realised genetic gain.
- Separated DGSI desired gains from QGSI linear and quadratic economic weights
  throughout the combined pipeline.
- Removed the unsupported constructor that fabricated quadratic economic
  weights from desired gains and empirical correlations.
- Added equation-level tests for scores, contributions, covariance estimators,
  response parameters, relationship-matrix alignment, and selection.

# DesiredGainR 0.2.0

- Renamed the package from DGQGSI to DesiredGainR.
- Added canonical `run_dgsi()` and `run_qgsi()` interfaces.
- Required an explicit genetic covariance matrix for DGSI and recorded the
  provenance of `P` and `G`.
- Changed threshold selection to eligibility followed by index ranking.
- Added automatic best-replicate selection, convergence, coefficient, ranking
  and selected-set stability diagnostics.
- Made missing-value handling explicit.
- Required an explicit symmetric quadratic-weight matrix for QGSI.
- Added candidate-specific linear, squared and cross-product contributions.

# DGQGSI 0.1.0

- Initial package structure.
