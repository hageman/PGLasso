<div align="center">

# PGLasso

### Propagation-guided lasso for sparse intervention in complex networks

**An open-source R package for identifying small, interpretable sets of network interventions.**

[Project website](https://hageman.github.io/PGLasso/) · [Installation](#installation) · [Quick start](#quick-start) · [Methods](#how-it-works) · [Citation](#citation)

</div>

---

## Overview

How can we identify a small number of perturbations that move a complex network from its initial state toward a desired target, when its full dynamics are unknown?

**PGLasso** combines network propagation with sparse regression. Starting from a network adjacency matrix and initial and target node states, it constructs a propagation-derived sensitivity matrix and uses constrained lasso regression to identify candidate intervention nodes and estimate their perturbation magnitudes.

The package is designed for research on biological and other complex networks in which a static network representation is available but a complete mechanistic model may not be.

### Highlights

- **Sparse intervention selection:** identify a manageable set of candidate nodes and estimate how strongly to perturb them.
- **Network-based sensitivity:** use PRINCE propagation to capture direct and indirect effects encoded by network topology.
- **Multiple network types:** work with undirected, directed, and directed signed adjacency matrices.
- **Feasibility-aware fitting:** support nonnegative post-intervention states through constrained regression.
- **Transparent comparison:** compare PGLasso with a separate reconstructed **NORSP-style forward-selection** baseline.
- **Reproducible demonstrations:** generate synthetic networks and evaluate support recovery and graph-based target error.

## Installation

PGLasso requires **R 4.1 or later**. Install the current development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("hageman/PGLasso")
library(PGLasso)
```

The package uses `glmnet`, which will be installed as a dependency if needed. PGLasso is currently distributed through GitHub; it is **not yet available on CRAN**.

## Quick start

This example generates a **50-node directed signed network** with three true intervention nodes and compares PGLasso with forward selection.

```r
library(PGLasso)

# Simulate a network, target state, and three true perturbations
sim <- pglasso_simulate(n = 50, n_true = 3, seed = 2026)

# Estimate sparse interventions with both methods
results <- pglasso_compare(sim, max_nodes = 10)

# Compare exact support recovery, precision, recall, F1, and target error
print(results$summary)

# See the selected intervention nodes
results$pglasso$selected
results$norsp_forward$selected
```

For this particular seed, both methods recover all three true intervention nodes. The example illustrates how to run the software; it is **not** a reproduction of the manuscript's aggregate simulation results.

A longer example, including a visualization, is available in [`inst/examples/simulation_example.R`](inst/examples/simulation_example.R).

## How it works

PGLasso begins with an adjacency matrix **G**, an initial state **z₀**, a target state **zₜ**, and a propagation parameter **α**. PRINCE propagation produces a sensitivity matrix **S**:

```text
S = (1 − α) (I − α G′)⁻¹
```

where **G′** is the degree-normalized adjacency matrix. The change from initial to target state is modeled as **zₜ − z₀ ≈ SΔ**. Lasso estimates a sparse perturbation vector **Δ**, optionally constrained so that post-intervention node states remain nonnegative.

The resulting coefficients describe **which nodes to perturb and by how much**. A solution is selected along the regularization path using an elbow criterion.

PGLasso does **not** infer network edges or estimate a time-dependent control signal. Its predicted response is based on the network-propagation approximation; biological validation in a nonlinear dynamical system requires a separate simulator.

## Main functions

| Function | Purpose |
| --- | --- |
| `pglasso_sensitivity()` | Compute the PRINCE-derived sensitivity matrix from a network. |
| `pglasso_fit()` | Estimate a sparse intervention vector using constrained lasso. |
| `norsp_forward()` | Run the reconstructed NORSP-style forward-selection comparator. |
| `pglasso_simulate()` | Generate a reproducible synthetic directed signed network experiment. |
| `pglasso_compare()` | Compare intervention recovery and graph-based target approximation. |
| `pglasso_example()` | Run a smaller built-in demonstration. |

**Network convention:** `G[i, j]` represents an edge directed **from node i to node j**. For signed networks, `+1` represents activation and `-1` inhibition.

### About the NORSP comparison

`norsp_forward()` implements a **reconstructed NORSP-style greedy forward-selection baseline** for comparison. It is **not** the original published NORSP software or a verified exact reproduction of that implementation. Please retain this distinction when citing or reporting results.

## Reproducibility and scope

The repository includes the standalone R package, automated tests, and illustrative synthetic-network examples. The original study additionally analyzed **82 mechanistic models from BioModels** using MATLAB/SBML and R scripts. Those original MATLAB pipelines and full experimental outputs are **not bundled in this repository**. BioModels accession numbers are reported in the manuscript's supplementary materials.

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for methodological details and distinctions between the illustrative examples and the reported experiments.

The R package has been tested on macOS with **11 passing automated tests** and successfully installed directly from GitHub. A single timestamp-verification note remained during the local R package check.

## Citation

If you use PGLasso in research, please cite the associated manuscript:

> Shi, T., Yu, H., and Hageman Blair, R. (2026). **Propagation-guided lasso for sparse intervention in complex networks.** *Manuscript submitted.*

**BibTeX** (no journal specified):

```bibtex
@unpublished{shi2026pglasso,
  author = {Shi, Tiange and Yu, Han and Hageman Blair, Rachael},
  title  = {Propagation-guided lasso for sparse intervention in complex networks},
  year   = {2026},
  note   = {Manuscript submitted},
  url    = {https://github.com/hageman/PGLasso}
}
```

From R, you can also access the package citation with `citation("PGLasso")` after installation.

## License and contact

PGLasso is open-source software released under the [MIT License](LICENSE.md).

**Research correspondence:** Rachael Hageman Blair · [hageman@buffalo.edu](mailto:hageman@buffalo.edu)
