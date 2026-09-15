# Evidence and readiness

DesiredGainR combines an operational quantitative-genetic core with advanced
recommendation, uncertainty and simulation tools. This document shows the
evidence supporting each layer and how it adds value to breeding decisions.

## Operationally ready core

- Objective definition and feasibility: `implied_economic_weights()`,
  `implied_desired_gains()`, `gain_feasibility()`, `retrospective_weights()`,
  `weight_sensitivity()`, and `effective_weights()`.
- Covariance diagnostics and validation: `matrix_diagnostics()`,
  `bend_covariance()`, and `import_covariance()`.

`gain_feasibility()` answers the single-cycle question exactly and without
simulation. Prefer it whenever it is sufficient.

## Validated index workflows

- Classical index families: `selection_index()`, `evaluate_index()`,
  `predict()`, and `summary()`.
- Restricted and proportional-gain indices: `restricted_index()`.
- Desired-gain index: `run_dgsi()`.
- Quadratic genomic index: `run_qgsi()` and `compare_dg_and_qgsi()`.

Version 0.5.0 strengthened DGSI with explicit desired-gain unit handling when
traits differ in scale. Earlier analyses should be rerun to benefit from the
verified implementation.

## Advanced uncertainty and environment analysis

- Covariance uncertainty: `index_uncertainty()`, `draw_covariance_pairs()`, and
  `propagate_covariance_uncertainty()`.
- Multi-environment expansion: `expand_environments()` and
  `widen_environments()`.

These functions use a Wishart approximation whose degrees of freedom must count
independent genetic units, not plots. Multi-environment functions expand an
already fitted covariance structure; they do not fit the trial model.

## Advanced recommendation and simulation tools

These tools support structured scenario comparison and breeder decisions when
their stated population, covariance and simulation assumptions match the
programme.

- Recurrent simulation: `founder_population()` and
  `simulate_selection_cycles()`.
- Direction optimisation: `optimize_desired_gains()`.
- Population-driven recommendation: `suggest_desired_gains()`.
- Diversity optimisation: `optimize_desired_gains(include_diversity = TRUE)`.

The simulation layer is credible for comparing directions under its stated
assumptions, offers phenotypic or cross-fitted RR-BLUP selection, and reports
the realised covariance, calibration error and held-out prediction accuracy.
`suggest_desired_gains()` has exact conditional one-cycle mathematics checked
against independent Clarabel and OSQP oracles. It separates multi-start
discovery, independent screening, multiple-finalist confirmation, model
calibration, architecture draws and optional covariance draws, and refuses to
call unresolved cases recommendations.

Real-programme evidence ships under `inst/validation/` and is explained in
`vignette("DesiredGainR-empirical-validation")`. It includes common-check and
year-adjusted recurrent-cycle trends, frozen-index temporal validation,
Genomes-to-Fields transport, the CIMMYT wheat anchor and the Zhang et al. maize
RCGS programme. The suite reproduces published gain, recovers raw genomic
prediction, verifies cross-file linkage and tests recommendations in held-out
environments. A future prospective vector-comparison experiment can extend
this evidence by quantifying incremental field response among strategies.

Diversity defaults off and requires a neutral marker panel known to be
disjoint from QTL, unless an experimental override for unknown overlap is
recorded. Known overlap is rejected.

DesiredGainR connects fitted genetic evaluations to transparent candidate
ranking, response prediction and objective choice. It complements mating-plan
and optimal-contribution tools, and provides independent solver and reference
checks suitable for high-value decisions.
