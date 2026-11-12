\name{pglasso_fit}
\alias{pglasso_fit}
\title{Fit a sparse network intervention using constrained lasso}
\usage{pglasso_fit(S, z0, zT, nlambda = 300L, max_nodes = 10L,
  selector = c("geometric", "derivative"), nonnegative = TRUE)}
\arguments{
\item{S}{Square numeric sensitivity matrix returned by \code{pglasso_sensitivity}.}
\item{z0}{Initial node-state vector.}
\item{zT}{Target node-state vector.}
\item{nlambda}{Number of regularization values passed to \code{glmnet}.}
\item{max_nodes}{Maximum number of nonzero perturbations permitted when selecting the regularization level. Set to \code{NULL} for the complete path.}
\item{selector}{Geometric elbow or derivative-based elbow heuristic.}
\item{nonnegative}{Whether to require the perturbed initial state to remain nonnegative.}
}
\value{A list containing \code{delta}, \code{selected}, \code{selected_names}, \code{predicted_state}, \code{error_ratio}, \code{lambda}, \code{path}, and \code{fit}.}
\description{Uses \code{glmnet} to fit lasso solutions to the network-derived model matrix, then selects a sparse intervention from the regularization path.}
\examples{
G <- matrix(c(0, 1, -1, 0), 2, 2)
S <- pglasso_sensitivity(G, 0.5)
pglasso_fit(S, z0 = c(1, 1), zT = c(0.9, 1.1), max_nodes = 2)
}
