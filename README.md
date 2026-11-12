# PGLasso — propagation-guided sparse network intervention

**Research R package (v0.1.0)**

PGLasso estimates a sparse vector of interventions that moves a network from an initial state toward a specified target. It uses a PRINCE-derived network sensitivity matrix and constrained lasso. A separate **NORSP-style forward-selection** function provides an approximate comparison method.

> **Publication:** Shi, T., Yu, H., and Hageman Blair, R. *Propagation-guided lasso for sparse intervention in complex networks*. Manuscript under revision (2026). Add a DOI/journal citation when available.

## Installation

Install R (version 4.1 or later), and then:

```r
install.packages("glmnet")
install.packages("remotes")
remotes::install_github("hageman/PGLasso")
```

Or download the repository and install from the extracted package directory:

```r
remotes::install_local("path/to/PGLasso")
```

## Five-minute example

```r
library(PGLasso)
sim <- pglasso_simulate(n = 50, n_true = 3, seed = 2026)
results <- pglasso_compare(sim, max_nodes = 10)
print(results$summary)
results$pglasso$selected_names
results$norsp_forward$selected_names
```

See [`inst/examples/simulation_example.R`](inst/examples/simulation_example.R) for the longer simulation, including a comparison plot.

## Function reference

| Function | Purpose |
|---|---|
| `pglasso_sensitivity(G, alpha)` | PRINCE normalized-adjacency sensitivity matrix |
| `pglasso_fit(S, z0, zT, ...)` | Constrained lasso and elbow-based selection |
| `norsp_forward(S, z0, zT, ...)` | Reconstructed NORSP-style forward selection |
| `pglasso_simulate(...)` | Reproducible directed signed network + target |
| `pglasso_compare(sim, ...)` | Support recovery and error-ratio comparison |
| `pglasso_example()` | Small example runnable from the R console |

### Input convention

`G[i,j]` represents an edge **from** node `i` **to** node `j`. Signed edges can be -1 or +1. For BioModels Jacobians, the original analysis formed `G = t(sign(J))` and removed self-loops. The sensitivity matrix is:

`S = (1-alpha) * solve(I - alpha * D_out^(-1/2) * G * D_in^(-1/2))`.

The intervention satisfies `z0 + Delta >= 0` when `nonnegative=TRUE`. `glmnet` uses a scaled quadratic objective, so its lambda values should **not** be compared numerically with lambda values from an unscaled optimization statement. The outputs include the graph-estimated target and graph error ratio; testing an intervention in a nonlinear biochemical model requires a separate dynamic simulator.

### Elbow selection and the paper

The package provides `selector="geometric"` for the simulated-network benchmark and `selector="derivative"` for an approximation of the original BioModels derivative-based heuristic. These are **different rules**; be explicit about which is used in reports. The simulation example is illustrative and **not a reproduction of all published experiments or summary figures**.

### What is—and is not—NORSP here

`norsp_forward()` is a **reconstructed NORSP-style greedy forward-selection baseline**, using the same network sensitivity matrix and a geometric elbow to choose the intervention set size. It is **not an exact implementation of the original NORSP publication/software**. Neither function estimates edges or solves an optimal-control problem over time.

### Reproducing the BioModels experiments

The original research archive includes MATLAB/SBML pipelines and R analysis scripts for 82 mechanistic BioModels. This lightweight R package does **not** bundle those large source data, MATLAB dependencies, or model outputs. To reproduce the paper's BioModels numerical results, the original pipeline/data must be deposited separately and the precise preprocessing and tuning configuration documented. See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## Development and testing

The GitHub Actions workflow runs `R CMD check` and `testthat` tests on commits and pull requests. Run locally with:

```r
install.packages("testthat")
testthat::test_local("PGLasso")
```

## License

Released under the [MIT License](LICENSE).

## Contact

For research correspondence: **hageman@buffalo.edu**.
