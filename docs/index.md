# PGLasso

## Propagation-guided lasso for sparse intervention in complex networks

PGLasso uses network propagation and penalized regression to identify a small set of candidate interventions that moves a biological network from an initial state toward a target state. It requires an adjacency matrix, not a fully specified dynamic system.

**Core functions:** `pglasso_sensitivity()`, `pglasso_fit()`, and `norsp_forward()`.

### Installation

```r
install.packages("glmnet")
install.packages("remotes")
remotes::install_github("YOUR-GITHUB-USERNAME/PGLasso")
```

### Example

```r
library(PGLasso)
sim <- pglasso_simulate(n = 50, n_true = 3, seed = 2026)
result <- pglasso_compare(sim)
result$summary
```

The comparison function labeled `norsp_forward()` is a reconstructed **NORSP-style** forward-selection baseline, not the original NORSP implementation. The simulation is an illustrative example, not a regeneration of the paper's full benchmark.

[Read the package README](../README.md) | [Reproducibility notes](../REPRODUCIBILITY.md)

**Manuscript:** Shi, T., Yu, H., and Hageman Blair, R. *Propagation-guided lasso for sparse intervention in complex networks*, under revision (2026).

**Contact:** hageman@buffalo.edu
