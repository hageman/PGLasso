#' Fit PGLasso using a sensitivity matrix
#'
#' @param S Precomputed propagation sensitivity matrix.
#' @param z0,zT Initial and target network state vectors.
#' @param alpha Unused if S is provided; compute S separately with
#'   pglasso_sensitivity().
#' @param nlambda Number of lambda values for glmnet.
#' @param max_nodes Maximum selected nodes when choosing an elbow. For the
#'   benchmark, default 10; NULL permits the complete path.
#' @param selector "geometric" (benchmark) or "derivative" (BioModels).
#' @param nonnegative Enforce z0 + Delta >= 0 via coefficient lower bounds.
#' @return List including Delta, selected nodes, predicted state, error ratio,
#'   lambda, full path diagnostics and the glmnet fit.
#' @export
pglasso_fit <- function(S, z0, zT, nlambda = 300L, max_nodes = 10L,
                        selector = c("geometric", "derivative"),
                        nonnegative = TRUE) {
  selector <- match.arg(selector)
  inp <- .pgl_inputs(S, z0, zT, max_nodes)
  if (nlambda < 4L) stop("nlambda must be at least 4")
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install glmnet to fit PGLasso")
  # glmnet scales the squared-loss objective by 1/(2n). Lambda is therefore
  # reported on glmnet's scale, not the unscaled manuscript convention.
  lower <- if (nonnegative) -inp$z0 else rep(-Inf, inp$n)
  fit <- glmnet::glmnet(x = inp$S, y = inp$y, family = "gaussian",
                        alpha = 1, nlambda = as.integer(nlambda),
                        intercept = FALSE, standardize = FALSE,
                        lower.limits = lower)
  betas <- as.matrix(fit$beta)
  sizes <- colSums(abs(betas) > 1e-7)
  predictions <- inp$S %*% betas
  errors <- vapply(seq_len(ncol(betas)), function(k) .pgl_error(inp$y, predictions[, k]), numeric(1))
  eligible <- which(sizes <= if (is.null(max_nodes)) inp$n else max_nodes)
  if (!length(eligible)) stop("No lambda values meet max_nodes")
  local <- if (selector == "geometric") .pgl_elbow(errors[eligible], sizes[eligible])
           else .pgl_derivative_elbow(errors[eligible], fit$lambda[eligible])
  k <- eligible[local]
  delta <- as.numeric(betas[, k])
  names(delta) <- if (!is.null(colnames(inp$S))) colnames(inp$S) else paste0("node", seq_len(inp$n))
  selected <- which(abs(delta) > 1e-7)
  structure(list(delta = delta, selected = selected,
                 selected_names = names(delta)[selected],
                 predicted_state = inp$z0 + as.numeric(inp$S %*% delta),
                 error_ratio = errors[k], lambda = fit$lambda[k], index = k,
                 selector = selector, method = "PGLasso",
                 path = data.frame(lambda = fit$lambda, size = sizes, error_ratio = errors),
                 fit = fit), class = "pglasso_result")
}
