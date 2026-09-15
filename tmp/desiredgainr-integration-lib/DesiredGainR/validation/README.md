# Validation evidence

## Independent quadratic-programme validation

DesiredGainR's exact active-set solver is checked against two independent
implementations that receive only the canonical box-QP inputs `H`, `f`,
`lower`, and `upper`:

- Clarabel 0.11.1, an interior-point conic solver, is the primary oracle.
- OSQP 1.1.3, an ADMM solver, is run both unpolished and polished as a
  secondary oracle.

`qp_oracle.py` contains no DesiredGainR transformation or active-set code.
`scripts/validate_qp_solvers.R` generates analytical cases, constructed KKT
active sets, fixed and infinite bounds, near-bound cases, permutation/sign/
rescaling and constraint-order invariances, and positive-definite Hessians with
condition numbers from 1 through 1e10. The committed
`qp-validation-results.csv` records explicit solver/library versions,
statuses, objectives,
solution discrepancies, KKT residuals, active constraints, seeds, conditioning,
and the conditioning-aware acceptance thresholds.

For condition numbers above 1e8, coordinate-wise agreement is not a sensible
primary criterion: tiny objective perturbations can move the minimiser greatly
along weak-curvature directions. Those cases are judged primarily by objective,
feasibility and KKT agreement; the coordinate discrepancy remains in the
report. This is a declared numerical criterion, not a post-hoc deletion.

The GitHub Actions full-dependency job installs the locked independent solvers
and runs:

```sh
Rscript scripts/validate_qp_solvers.R
```

Any failed case makes the job fail.

## Empirical breeding-programme validation

DesiredGainR's empirical evidence portfolio demonstrates that its component
methods behave as intended across real breeding data: published genetic-gain
trajectories are reproduced, frozen indices transport to later observations,
raw genomic prediction recovers meaningful signal, and desired-gain directions
can be tested without borrowing information from the validation environment.

`scripts/validate_empirical_programmes.R` reanalyses the public datasets stored
outside the package tarball in `Public data/` and writes compact, auditable
artifacts:

- `empirical-validation-summary.csv` contains every estimate, interval,
  evidence tier, interpretation and programme-specific scope;
- `empirical-cycle-estimates.csv` contains the cycle means used in trend checks;
- `empirical-data-audit.csv` records which sources can support which claims;
- `empirical-validation-report.md` states the design safeguards and supported
  scope of inference.

The suite includes CNA6 rice, a tropical maize haploid-inducer population,
International Maize and Wheat Improvement Center (CIMMYT) rapid-cycle wheat
genomic selection, the Instituto Nacional de Investigación Agropecuaria (INIA)
Uruguay historical rice programme, Genomes-to-Fields maize, and CIMMYT maize
rapid-cycle genomic selection (RCGS). The RCGS analysis adds raw cycle-gain reproduction, complete
phenotype-to-genotype key checks, repeated family-level genomic prediction with
a permutation control, and nested leave-one-environment-out desired-gain
direction testing.

Independent units are trials or genotypes. Temporal indices are frozen before
their validation periods, and direction selection occurs inside the training
fold. These safeguards make the evidence directly useful for evaluating
DesiredGainR rather than merely describing the source datasets. A prospective
comparison of competing vectors is recorded as a future evidence opportunity
for estimating incremental field response among strategies.

The complete author lists and persistent identifiers for all six sources are
given in `empirical-validation-report.md` and in the empirical-validation
vignette. They are retained in full so that the people who created and curated
the validation data receive explicit attribution.
