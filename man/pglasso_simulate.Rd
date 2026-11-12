\name{pglasso_simulate}
\alias{pglasso_simulate}
\title{Simulate a directed signed network intervention}
\usage{pglasso_simulate(n = 50L, expected_degree = 4, n_true = 3L,
  alpha = 0.9, noise_sd = 0.05, rewire_frac = 0, seed = 123L)}
\arguments{
\item{n}{Number of network nodes.}
\item{expected_degree}{Expected outgoing edge count.}
\item{n_true}{Number of perturbed nodes.}
\item{alpha}{Propagation parameter.}
\item{noise_sd}{Observation noise relative to standard deviation of the target change.}
\item{rewire_frac}{Fraction of edges rewired in the observed network.}
\item{seed}{Random seed.}
}
\value{A list containing original and observed networks, the observed sensitivity matrix, initial and target states, true perturbation and selected nodes.}
\description{Generates a reproducible illustrative signed directed network and a sparse target state.}
