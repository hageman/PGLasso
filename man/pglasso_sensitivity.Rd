\name{pglasso_sensitivity}
\alias{pglasso_sensitivity}
\title{Compute a propagation-derived network sensitivity matrix}
\usage{pglasso_sensitivity(G, alpha = 0.9)}
\arguments{
\item{G}{A square numeric adjacency matrix; rows denote source nodes and columns targets.}
\item{alpha}{Propagation fraction, strictly between 0 and 1.}
}
\value{A numeric sensitivity matrix, of the same dimensions as \code{G}.}
\description{Uses row and column absolute-degree normalization and PRINCE propagation. Nodes with zero degree receive zero inverse-degree weight.}
\examples{
G <- matrix(c(0, 1, -1, 0), 2, 2)
pglasso_sensitivity(G, 0.5)
}
