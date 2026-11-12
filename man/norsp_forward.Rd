\name{norsp_forward}
\alias{norsp_forward}
\title{Reconstructed NORSP-style forward selection}
\usage{norsp_forward(S, z0, zT, max_nodes = 10L,
  nonnegative = TRUE, selector = "geometric")}
\arguments{
\item{S}{Square sensitivity matrix.}
\item{z0}{Initial state vector.}
\item{zT}{Target state vector.}
\item{max_nodes}{Maximum number of selected candidate nodes.}
\item{nonnegative}{Enforce a nonnegative perturbed initial state.}
\item{selector}{Currently \code{"geometric"} only.}
}
\value{A list with the chosen intervention \code{delta}, \code{selected} nodes, error ratio and candidate path.}
\description{Greedily adds candidate nodes using active-set least-squares fitting and an elbow criterion. This comparator is not the original NORSP software.}
