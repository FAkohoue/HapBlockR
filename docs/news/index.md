# Changelog

## HapBlockR 0.3.12.9000 (development)

### Scientific correctness

- Fixed
  [`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md)’s
  `optimal_solution_found` quality gate, caught by `R CMD check` failing
  to rebuild the *full pipeline* vignette: `identical(sol$status, 0)`
  compared [`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html)’s
  returned status against a double `0` literal using
  [`identical()`](https://rdrr.io/r/base/identical.html), which treats
  integer and double as different types even when numerically equal – so
  the gate failed on every genuinely optimal solve whenever lpSolve
  returned an integer status. Replaced with `isTRUE(sol$status == 0)`,
  which compares by value regardless of that type distinction.

- Fixed a regression in
  [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)’s
  family-size eligibility rule caught by the test suite: under
  `family_select_mode = "count"` with a named-vector `n_per_family`, an
  earlier revision in this same development version used each family’s
  own named quota as its exclusion threshold when available. That broke
  the documented and tested behaviour – `n_per_family` is only
  guaranteed to be defined for whichever families end up *chosen*, which
  isn’t decided until after ranking, so a name already present at the
  pre-ranking exclusion stage cannot yet be trusted as that family’s
  real quota. `rank_k` is now the exclusion threshold for every family
  under a named-vector quota, named or not, matching the original design
  and `test-family-selection.R`.

- Parent shortlisting now separates coverage-only
  [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  from the explicit joint
  [`select_parents_ga_ts()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga_ts.md)
  objective. The hybrid tool requires a positive whole-genome merit
  contribution, records its own result method and cannot silently fall
  back to coverage-only selection. Merit floors, `merit_priority` and
  `merit_weight` now belong only to the hybrid interface.

- [`select_parents_ga()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
  now evaluates the complete documented objective and enforces exact
  founder count and relatedness limits as hard postconditions.
  Replicated searches report stability and automatically return the
  feasible replicate with the greatest complete objective.

- `select_parents_ocs(engine = "optisel")` no longer silently
  substitutes a heuristic plan after a failed solver call. Parent-count
  limits trigger a genuine constrained re-solve.

- Cross-validation now returns row-level out-of-fold predictions and
  pooled predictive ability, RMSE, MAE, bias, and calibration. Random,
  grouped, and forward designs are available.

- Epistasis fine mapping removes degenerate interaction terms, uses
  deterministic folds, and reports the actual number of tests.

- Association testing no longer labels an unadjusted result as adjusted
  when the requested mixed model cannot be fitted. Fallback behaviour is
  explicit.

- Core-collection selection validates metric-distance assumptions and
  uses square-root Euclidean distances derived from genomic
  relationships.

- Fixed `get_V_inv_sqrt(method = "chol")` (the default): the Cholesky
  whitening factor was returned untransposed, so it did not satisfy
  `A V t(A) = I` for a non-diagonal relationship matrix. This silently
  defeated `rV2`’s kinship whitening for structured/related populations
  – exactly the case the metric exists for – while remaining invisible
  on the diagonal matrices the previous test suite used.
  `method = "eigen"` was unaffected. A non-diagonal regression test now
  covers this.

- Merit-floor interfaces now provide `min_sel_mode = "sd_above_mean"`
  for intuitive selection of directionally superior candidates. A value
  of zero retains candidates at or above the mean and positive values
  require the stated SD superiority. Deliberately broad candidate pools
  use the clearer `"relaxed_pool"` name. The former `"sd_below_mean"`
  name remains only as a deprecated compatibility alias.

- Fixed
  [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md)’s
  `too_small` family-exclusion check under
  `family_select_mode = "count"`: when `n_per_family` is a per-family
  named vector, the check now compares each family’s eligible membership
  against ITS OWN quota. It previously borrowed a single scalar
  threshold (the largest quota across all families) for every family,
  which could wrongly exclude a family whose own, smaller quota it would
  have comfortably satisfied.

- Epistasis fine mapping’s `p_bonf` now uses the number of interaction
  tests that actually produced a valid fit, not the pre-loop
  candidate-pair count (which counted pairs later skipped for
  degenerate/rank-deficient OLS fits). This makes the multiple-testing
  correction match the number of tests actually reported, as already
  documented above.

- [`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md)’s
  single-lead-SNP-per-block comparison now computes a genuine
  2-population, 1-df Cochran’s Q heterogeneity test – the same test its
  sibling function
  [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md)
  already used for a single shared allele – instead of declaring Q
  “undefined” and falling back to pooled-effect significance. A QTL with
  a large effect in one population and a much smaller (even same-signed)
  effect in the other is now correctly flagged as not replicating.

- [`cv_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/cv_haplotype_prediction.md)’s
  quality gate no longer requires every individual/trait/repetition
  combination to be populated. It now only flags an individual actually
  tested MORE than once for the same trait within the same repetition
  (real test-set duplication). The previous all-cells-populated check
  produced spurious failures under multi-trait `"forward"` validation
  whenever traits had different phenotyping coverage.

- [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)’s
  strict phasing-reliability gate now applies to
  `variance_model = "linked"` as well as `"phased"` – both build progeny
  variance from phased haplotype blocks and depend equally on phasing
  accuracy.

- [`assess_phasing_accuracy()`](https://FAkohoue.github.io/HapBlockR/reference/assess_phasing_accuracy.md)’s
  per-chromosome switch-error computation now explicitly sorts variants
  by position before counting orientation transitions. It previously
  assumed the truth set’s stored row order was already position-sorted
  within each chromosome.

- Fixed a block-overlap-resolution gap in the C++
  `resolve_overlap_cpp()` kernel used by
  [`Big_LD()`](https://FAkohoue.github.io/HapBlockR/reference/Big_LD.md):
  a single batch pass could leave a chain of 3+ mutually overlapping
  blocks incorrectly resolved – two independently computed boundary
  splits could squeeze the shared middle block into an invalid interval
  (silently dropped) while the two outer blocks, whose direct overlap
  was never itself checked, survived as still-overlapping in the output.
  The resolver now repeats the detect/resolve/clean pass on the current
  block state until a pass finds no remaining overlaps.

- [`select_parents_by_family()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_by_family.md),
  [`validate_crosses_exact()`](https://FAkohoue.github.io/HapBlockR/reference/validate_crosses_exact.md),
  [`estimate_diplotype_effects()`](https://FAkohoue.github.io/HapBlockR/reference/estimate_diplotype_effects.md),
  [`compare_block_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_block_effects.md),
  and
  [`compare_gwas_effects()`](https://FAkohoue.github.io/HapBlockR/reference/compare_gwas_effects.md)
  now return versioned `hapblockr_result` objects (parameters,
  identifiers, transformations, quality gates, decision table, and
  uncertainty), matching every other breeder-facing decision function.
  [`validate()`](https://FAkohoue.github.io/HapBlockR/reference/validate.md)
  and
  [`build_breeding_exchange()`](https://FAkohoue.github.io/HapBlockR/reference/build_breeding_exchange.md)
  previously errored on these five results because they were plain
  lists.

- [`certify_mating_plan()`](https://FAkohoue.github.io/HapBlockR/reference/certify_mating_plan.md)’s
  `cross_rule` and `no_repeated_cross` violation rows now identify the
  actual parent pair rather than a row index, and the repeated-cross
  violation reports the true observed repeat count instead of a
  hardcoded value of 2.

### Identity, provenance, and reproducibility

- Genotype import preserves immutable physical sample and variant IDs,
  separates display labels, verifies strict chunk identity, and uses
  SHA-256 cache manifests.
- VCF parsing has explicit multiallelic policies and stricter genotype
  semantics.
- Breeder-facing outputs use a versioned `hapblockr_result` contract
  with methods for validation, printing, summaries, plots, and
  data-frame conversion.
- The contract records parameters, seeds, identifiers, hashes,
  transformations, software, quality-control gates, warnings, fallbacks,
  exclusions, decisions, uncertainty, and validation status.

### Uncertainty and decisions

- Prediction output reports GEBV prediction error variance and
  reliability where defensible, together with a minimum-reliability
  recommendation gate.
- [`compute_local_gebv()`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md)
  can propagate supplied SNP-effect standard errors under an explicit
  independent-marker-error approximation.
- [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
  reports predicted progeny standard deviation, downside quantiles,
  uncertainty intervals, conservative cross reliability, phasing
  reliability where required, and recommendation reasons.
- [`assess_decision_stability()`](https://FAkohoue.github.io/HapBlockR/reference/assess_decision_stability.md)
  compares threshold, environment, population, and model scenarios
  through selection frequency and Jaccard overlap.

### Multi-trait, environment, and operations

- [`fit_multitrait_gblup()`](https://FAkohoue.github.io/HapBlockR/reference/fit_multitrait_gblup.md)
  fits a genuine multivariate genomic BLUP with genetic and residual
  covariance diagnostics, PEV, and reliability.
- [`build_selection_index()`](https://FAkohoue.github.io/HapBlockR/reference/build_selection_index.md)
  constructs Smith-Hazel and Pesek-Baker indices internally and
  delegates DGSI and QGSI to DesiredGainR. Trait directions are
  mandatory; declared units are optional.
- [`fit_gxe_gblup()`](https://FAkohoue.github.io/HapBlockR/reference/fit_gxe_gblup.md)
  fits genomic main and reaction-norm effects using an optional
  environmental kernel;
  [`build_environment_kernel()`](https://FAkohoue.github.io/HapBlockR/reference/build_environment_kernel.md)
  derives a kernel from validated environmental covariates.
- [`prepare_breeding_targets()`](https://FAkohoue.github.io/HapBlockR/reference/prepare_breeding_targets.md)
  replaces field-trial analysis. HapBlockR now accepts externally
  analysed genotype-level adjusted means, BLUEs, identity BLUPs,
  pedigree BLUPs, breeding values, general combining ability or total
  genetic values with declared uncertainty, and propagates their
  precision through internal marker, haplotype, multi-trait and GxE
  models.
- [`screen_candidate_crosses()`](https://FAkohoue.github.io/HapBlockR/reference/screen_candidate_crosses.md)
  and
  [`certify_mating_plan()`](https://FAkohoue.github.io/HapBlockR/reference/certify_mating_plan.md)
  cover sex and role, flowering, fertility, reciprocal effects,
  availability, parent and period capacity, required and forbidden
  pairs, heterotic groups, subpopulation quotas, quarantine, family
  size, violations, and binding constraints.
- Canonical metadata validation and checksum-protected
  BrAPI/MIAPPE-oriented exchange tables support traceable integration
  without embedding remote service calls in modelling functions.

### External tools and compiled code

- Beagle JAR archives are no longer distributed. Beagle 5.x must be
  explicitly configured by argument, option, environment variable, or a
  JAR beside the output prefix.
- [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md)
  records Java, JAR, command, reference, map, input, and output
  provenance; verifies the reported Beagle version; and checks sample,
  variant, allele, genotype, and imputation quality.
- [`assess_phasing_accuracy()`](https://FAkohoue.github.io/HapBlockR/reference/assess_phasing_accuracy.md)
  adds truth-set dosage accuracy, phased-allele concordance, call rate,
  and per-sample and per-chromosome switch-error estimates.
  [`phase_with_beagle()`](https://FAkohoue.github.io/HapBlockR/reference/phase_with_beagle.md)
  can enforce these metrics through explicit thresholds.
- Compiled parallel methods use at most two threads, check for user
  interruption, and build with optional OpenMP through portable
  configure logic.

### Documentation and repository

- Package-level claims now distinguish file-backed from in-memory
  pathways and state the limits of diploid, biallelic haplotype support.
- The breeder’s guide is generated from version-controlled source as a
  tagged PDF and HTML edition. It now includes detailed tool variants, a
  worked decision, uncertainty and feasibility gates,
  result-interpretation guidance, references, change history, and a
  sign-off template.
- `tools/build_breeder_guide_accessible.cjs` now actually builds both
  shipped editions from the real `.Rmd` source: it previously pointed at
  a non-existent `.md` file and, even if redirected, would have parsed
  the document’s YAML frontmatter and inline R date expression as
  literal text (no R/knitr/pandoc is invoked). The script now strips and
  evaluates the frontmatter, builds a title block, reproduces
  `number_sections`-style chapter numbering, derives the in-page
  Contents list from the actual chapters instead of a hand-maintained
  anchor list, and writes the HTML edition in addition to the PDF.
- British English is declared for package prose. Citation, contribution,
  conduct, security, and support files have been added.
- The README is now a short, task-oriented front door: install, one
  example, the capability map, and links to the vignettes and Breeder’s
  Guide, rather than a complete method catalogue. Its Breeder’s Guide
  reference previously named the wrong source path and claimed a Word
  format that
  [`open_breeder_guide()`](https://FAkohoue.github.io/HapBlockR/reference/open_breeder_guide.md)
  has never offered (`c("pdf", "html")` only); both are now corrected
  package-wide (`R/breeder_guide.R`, `man/open_breeder_guide.Rd`,
  `README.md`).
- The pkgdown navbar now has explicit Tutorial, Vignettes, and Breeder
  guide tabs (`_pkgdown.yml`), and the GitHub Actions pkgdown workflow
  publishes the built Breeder’s Guide alongside the site so the new
  navbar link resolves to a real page.

## HapBlockR 0.3.12

- Added family- and genetic-cluster-quota selection, local-GEBV tools,
  genomic mating, optimum-contribution selection, Pareto selection,
  exact cross-plan validation, and core-collection support.
- Expanded haplotype association, population comparison, epistasis,
  prediction, phasing, and large-data pathways.

## HapBlockR 0.3.0

- Introduced the integrated LD-block, haplotype, prediction, and
  breeding decision workflow.

## HapBlockR 0.2.0

- Added compiled LD calculations and multi-chromosome block detection.

## HapBlockR 0.1.0

- Initial public development release.

Detailed historical changes remain available in the repository history
and release notes.
