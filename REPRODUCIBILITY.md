# Reproducibility scope and source provenance

This package is a new, lightweight implementation of the algorithmic workflow described in the manuscript. It was assembled by examining the research archive `codes(3).zip` and the revision of *Propagation-guided lasso for sparse intervention in complex networks* (12 November 2026).

## Supported operations

- Network normalization using absolute row- and column-degree and zero-degree inverses set to zero.
- PRINCE sensitivity matrix `(1-alpha) (I-alpha G')^{-1}`.
- Sparse interventions using `glmnet` with nonnegative-state constraints through lower coefficient bounds.
- Two elbow selectors (geometric vs derivative heuristic); note their distinct scopes.
- Separate NORSP-style forward-selection baseline with constrained active-set refitting.
- Seeded signed directed network simulation, optional observation noise and edge rewiring, support metrics.

## Deviations / limitations

- The original BioModels script `2 - LASSO_estimate_D_revised.R` uses `glmnet` with default intercept and feature-standardization settings. The new package explicitly uses **no intercept** and **no feature standardization** to reflect the model equation `y = S Delta`. Therefore, it is *not guaranteed to match historical coefficients* from the archived pipeline. Confirm or change these options if exact replication is needed.
- The original BioModels derivative elbow and the geometric elbow described for the new synthetic benchmark are not identical. The derivative selector included here is a practical reconstruction, not an exact line-by-line copy of the original local-extrema code.
- The original MATLAB SBML dynamic experiments, raw BioModels data, and plots are not distributed with this minimal package.
- The NORSP-style algorithm is a reconstructed forward-selection comparator, **not** the original NORSP code or an independently verified reproduction.
- Random-network simulation follows the manuscript's broad design but does not claim to regenerate the paper's reported numerical summary statistics.
- The R package must pass `R CMD check` in an R-enabled environment; creation in a runtime without R did not permit execution of that check locally. GitHub Actions is included for automated verification after upload.

## Suggested archival organization

- `PGLasso/`: portable R package and example.
- `reproduction/biomodels/`: original published-source pipeline and scripts, after license review and removal of external dependencies/binaries.
- `reproduction/simulated_networks/`: exact scripts, output tables, and seeds used for paper's comparative benchmark.
- GitHub release linked to Zenodo for an immutable DOI after the release is approved.

Never present illustrative simulation results as reproducing manuscript accuracy percentages without independent verification.
