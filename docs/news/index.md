# Changelog

## HapBlockR 0.3.12.9000 (development)

### Scientific correctness

- Selection-index objectives now use one sign convention: `directions`
  controls increase or decrease and weights or desired gains must be
  non-negative favourable-direction magnitudes. Silent
  [`abs()`](https://rdrr.io/r/base/MathFun.html) coercion has been
  removed.
- QGSI now receives `n_select`, reports linear and quadratic weights
  without fabricating one global coefficient vector, and reports
  model-expected gains from the total QGSI variance. DGSI now reports
  both its model-expected genetic response and its separately labelled
  standardised realised selected-set differential. Coefficient tables
  have one schema across all four methods.
- The DesiredGainR integration now requires version 0.5.0 and consumes
  its authoritative DGSI theoretical response and QGSI expected-gain
  outputs. Responses are converted from the engine’s analysis scale to
  original trait units, including when scaling is requested, and
  selection intensity reflects the number actually selected.
  `qgsi_control` exposes non-structural `run_qgsi()` controls without
  allowing HapBlockR’s data, objective, direction or selection-count
  arguments to be overridden.
- Integration tests now compare delegated outputs directly with the
  updated DesiredGainR engines, including reference scaling, holdout
  selection, eligibility thresholds, explicit genomic covariance and
  relationship-aware covariance estimation. Control names must match the
  installed API exactly. DGSI returns `coefficients_original_units` and
  `score_intercept` for correct original-unit effect propagation and
  score reconstruction, while preserving the engine’s coefficients and
  complete result. Result provenance records the dependency version and
  hashes the objective and engine controls.
- Multivariate and GxE REML likelihoods now use residual degrees of
  freedom. Prediction error variances include fixed-effect estimation
  uncertainty.
- Full sampling covariance uses exact target record keys and checked
  diagonal precision. `sampling_covariance_mode` explicitly selects
  sampling-only or sampling-plus-residual models, preventing silent
  covariance loss or unintended residual duplication. The sampling-only
  optimiser no longer allocates a dense zero matrix during every
  objective evaluation.
- GxE across-environment predictions again select one row explicitly per
  genotype, avoiding BLAS-dependent duplicate IDs caused by
  final-decimal differences. Prepared-target validation now requires
  `record_key` and uses the same trait-by-environment precision grouping
  during preparation and covariance consistency checks.
- Kernel diagonal validation is independent of the spectral PSD
  tolerance, and selection intensity requires one positive finite real
  scalar.
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
- Merit-floor interfaces now provide `min_sel_mode = "sd_above_mean"`
  for intuitive selection of directionally superior candidates. A value
  of zero retains candidates at or above the mean and positive values
  require the stated SD superiority. Deliberately broad candidate pools
  use the clearer `"relaxed_pool"` name. The former `"sd_below_mean"`
  name remains only as a deprecated compatibility alias.

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
  tagged PDF and editable Word document. It now includes detailed tool
  variants, a worked decision, uncertainty and feasibility gates,
  result-interpretation guidance, references, change history, and a
  sign-off template.
- British English is declared for package prose. Citation, contribution,
  conduct, security, and support files have been added.
- The README is now a task-oriented introduction rather than a complete
  method catalogue.

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
