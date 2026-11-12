# PGLasso

### Propagation-guided lasso for sparse intervention in complex networks

**Open-source R software for selecting sparse interventions in complex networks.**

[GitHub repository](https://github.com/hageman/PGLasso) · [Read the README](https://github.com/hageman/PGLasso#readme) · [Simulation example](https://github.com/hageman/PGLasso/blob/main/inst/examples/simulation_example.R)

---

## What is PGLasso?

PGLasso identifies candidate network nodes to perturb—and estimates their perturbation magnitudes—to move a network from an observed initial state toward a specified target. The method combines a **PRINCE network-propagation sensitivity matrix** with **constrained lasso regression**, without requiring a fully specified nonlinear dynamic model.

### What the software includes

- **PGLasso estimation:** `pglasso_sensitivity()` and `pglasso_fit()`
- **NORSP-style forward selection:** `norsp_forward()`; a reconstructed comparator, not the original NORSP implementation
- **Synthetic-network demonstrations:** `pglasso_simulate()` and `pglasso_compare()`
- **Installation instructions, function documentation, and automated tests**

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github("hageman/PGLasso")
library(PGLasso)
```

## Try an example

```r
sim <- pglasso_simulate(n = 50, n_true = 3, seed = 2026)
results <- pglasso_compare(sim, max_nodes = 10)
print(results$summary)
```

This example is illustrative and does not reproduce all numerical results in the associated manuscript. The full MATLAB/SBML BioModels pipeline is not included in the R package; see [reproducibility notes](https://github.com/hageman/PGLasso/blob/main/REPRODUCIBILITY.md).

## Cite PGLasso

Shi, T., Yu, H., and Hageman Blair, R. (2026). *Propagation-guided lasso for sparse intervention in complex networks.* **Manuscript submitted.**

```bibtex
@unpublished{shi2026pglasso,
  author = {Shi, Tiange and Yu, Han and Hageman Blair, Rachael},
  title  = {Propagation-guided lasso for sparse intervention in complex networks},
  year   = {2026},
  note   = {Manuscript submitted},
  url    = {https://github.com/hageman/PGLasso}
}
```

**License:** [MIT](https://github.com/hageman/PGLasso/blob/main/LICENSE.md) · **Contact:** [hageman@buffalo.edu](mailto:hageman@buffalo.edu)
