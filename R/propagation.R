#' Compute a PRINCE propagation sensitivity matrix
#'
#' @param G Square adjacency matrix. Rows are sources and columns are targets.
#'   Signed edges may be -1, 0, or 1. Undirected matrices should be symmetric.
#' @param alpha Propagation parameter, strictly between 0 and 1.
#' @return A numeric matrix S = (1-alpha) (I-alpha G')^{-1}.
#' @export
pglasso_sensitivity <- function(G, alpha = 0.9) {
  G <- .pgl_matrix(G)
  if (nrow(G) != ncol(G) || !nrow(G)) stop("G must be a nonempty square matrix")
  if (any(!is.finite(G))) stop("G must contain finite values")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  degree_out <- rowSums(abs(G))
  degree_in <- colSums(abs(G))
  inv_out <- ifelse(degree_out > 0, 1 / sqrt(pmax(degree_out, 1)), 0)
  inv_in <- ifelse(degree_in > 0, 1 / sqrt(pmax(degree_in, 1)), 0)
  normalized <- sweep(sweep(G, 1L, inv_out, `*`), 2L, inv_in, `*`)
  n <- nrow(G)
  S <- tryCatch((1 - alpha) * solve(diag(n) - alpha * normalized),
                error = function(e) stop("Cannot compute sensitivity: ", conditionMessage(e)))
  dimnames(S) <- dimnames(G)
  S
}

.pgl_matrix <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("Input must be a numeric matrix")
  storage.mode(x) <- "double"
  x
}

.pgl_inputs <- function(S, z0, zT, max_nodes) {
  S <- .pgl_matrix(S)
  n <- nrow(S)
  if (n == 0 || ncol(S) != n) stop("S must be a nonempty square matrix")
  if (any(!is.finite(S))) stop("S must contain finite values")
  if (!is.numeric(z0) || !is.numeric(zT) || length(z0) != n ||
      length(zT) != n || any(!is.finite(z0)) || any(!is.finite(zT))) {
    stop("z0 and zT must be finite numeric vectors with nrow(S) entries")
  }
  if (!is.null(max_nodes) && (!is.numeric(max_nodes) || length(max_nodes) != 1L ||
       !is.finite(max_nodes) || max_nodes < 1 || max_nodes != floor(max_nodes))) {
    stop("max_nodes must be a positive integer")
  }
  list(S = S, z0 = as.numeric(z0), zT = as.numeric(zT),
       y = as.numeric(zT - z0), n = n)
}

.pgl_error <- function(y, fitted) {
  denom <- sum(y^2)
  if (denom < .Machine$double.eps) return(if (sum((y - fitted)^2) < 1e-12) 0 else Inf)
  sum((y - fitted)^2) / denom
}

#' Select a sparse intervention by a geometric elbow
#'
#' Selects the point with largest (normalized residual improvement - normalized
#' intervention count) from candidate solutions. This is the reconstructed
#' synthetic benchmark rule; it is not the derivative-based original BioModels
#' elbow rule.
.pgl_elbow <- function(errors, sizes) {
  stopifnot(length(errors) == length(sizes), length(errors) > 0L)
  improvement <- max(errors) - errors
  if (max(improvement) > min(improvement)) {
    improvement <- (improvement - min(improvement)) / diff(range(improvement))
  } else improvement[] <- 0
  if (max(sizes) > min(sizes)) {
    complexity <- (sizes - min(sizes)) / diff(range(sizes))
  } else complexity <- rep(0, length(sizes))
  which.max(improvement - complexity)
}

#' Choose an elbow based on the original BioModels derivative heuristic
#'
#' Attempts the earliest local maximum or sign change in the finite-difference
#' error derivative with respect to log(lambda), falling back to lowest error.
.pgl_derivative_elbow <- function(errors, lambdas) {
  n <- length(errors)
  if (n < 4L) return(which.min(errors))
  d <- diff(sqrt(pmax(errors, 0))) / diff(log(lambdas))
  idx <- integer()
  for (i in seq.int(2L, length(d) - 1L)) {
    if (is.finite(d[i]) && is.finite(d[i - 1L]) && is.finite(d[i + 1L])) {
      if (d[i - 1L] < d[i] && d[i] > d[i + 1L] ||
          d[i - 1L] * d[i + 1L] <= 0) idx <- c(idx, i)
    }
  }
  if (length(idx)) idx[1L] else which.min(errors)
}
