\name{pglasso_compare}
\alias{pglasso_compare}
\title{Compare sparse intervention solutions}
\usage{pglasso_compare(sim, max_nodes = 10L)}
\arguments{
\item{sim}{A simulated network returned by \code{pglasso_simulate}.}
\item{max_nodes}{Maximum selected interventions.}
}
\value{A list with \code{pglasso}, \code{norsp_forward} and \code{summary} results.}
\description{Compares PGLasso and a reconstructed NORSP-style forward-selection baseline with exact support recovery, precision, recall, F1 score and error ratio.}
