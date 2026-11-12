#' NORSP-style forward-selection comparison
#'
#' Reconstructed forward subset selection baseline, NOT the published NORSP
#' software or an exact replication. Candidate nodes are added greedily to
#' minimize least squares residual. For nonnegative states, each active-set fit
#' is solved using constrained quadratic optimization (L-BFGS-B).
#'
#' @param S Propagation sensitivity matrix.
#' @param z0,zT Initial and target network states.
#' @param max_nodes Maximum number of selected nodes, default 10.
#' @param nonnegative Constrain z0 + Delta >= 0.
#' @param selector Currently only the geometric elbow rule.
#' @return List with selected nodes, perturbations, error and forward path.
#' @export
norsp_forward <- function(S, z0, zT, max_nodes = 10L,
                          nonnegative = TRUE, selector = "geometric") {
  if (!identical(selector, "geometric")) stop("Only the geometric selector is supported")
  inp <- .pgl_inputs(S, z0, zT, max_nodes)
  nsteps <- min(max_nodes, inp$n)
  chosen <- integer()
  deltas <- matrix(0, nrow = inp$n, ncol = nsteps + 1L)
  errors <- numeric(nsteps + 1L)
  errors[1L] <- .pgl_error(inp$y, rep(0, inp$n))
  sizes <- 0:nsteps
  for (step in seq_len(nsteps)) {
    candidates <- setdiff(seq_len(inp$n), chosen)
    attempts <- lapply(candidates, function(candidate) {
      active <- c(chosen, candidate)
      X <- inp$S[, active, drop = FALSE]
      if (nonnegative) {
        lower <- -inp$z0[active]
        objective <- function(b) sum((inp$y - as.numeric(X %*% b))^2)
        gradient <- function(b) -2 * as.numeric(crossprod(X, inp$y - as.numeric(X %*% b)))
        start <- pmax(lower, rep(0, length(active)))
        opt <- optim(start, objective, gradient, method = "L-BFGS-B",
                     lower = lower, control = list(maxit = 1000L, factr = 1e7))
        coef <- opt$par
      } else {
        coef <- as.numeric(qr.coef(qr(X), inp$y))
        coef[is.na(coef)] <- 0
      }
      d <- numeric(inp$n)
      d[active] <- coef
      list(delta = d, error = .pgl_error(inp$y, as.numeric(inp$S %*% d)))
    })
    winner <- which.min(vapply(attempts, `[[`, numeric(1), "error"))
    chosen <- c(chosen, candidates[winner])
    deltas[, step + 1L] <- attempts[[winner]]$delta
    errors[step + 1L] <- attempts[[winner]]$error
  }
  k <- .pgl_elbow(errors, sizes)
  delta <- deltas[, k]
  names(delta) <- if (!is.null(colnames(inp$S))) colnames(inp$S) else paste0("node", seq_len(inp$n))
  selected <- which(abs(delta) > 1e-7)
  structure(list(delta = delta, selected = selected,
                 selected_names = names(delta)[selected],
                 predicted_state = inp$z0 + as.numeric(inp$S %*% delta),
                 error_ratio = errors[k], index = k, method = "NORSP-style forward selection",
                 path = data.frame(size = sizes, error_ratio = errors)),
            class = "pglasso_result")
}
