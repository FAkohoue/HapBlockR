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
