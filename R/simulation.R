#' Generate a reproducible directed signed network intervention
#'
#' @param n Number of nodes (at least 5).
#' @param expected_degree Expected number of outgoing edges per node.
#' @param n_true Number of truly perturbed nodes.
#' @param alpha Propagation parameter.
#' @param noise_sd Gaussian noise SD relative to the SD of the true state change.
#' @param rewire_frac Fraction of edges to rewire in the observed network.
#' @param seed Random seed, used locally.
#' @return List with true and observed adjacency, S, z0, zT, and true Delta.
#' @export
pglasso_simulate <- function(n = 50L, expected_degree = 4, n_true = 3L,
                             alpha = 0.9, noise_sd = 0.05,
                             rewire_frac = 0, seed = 123L) {
  if (n < 5 || n_true < 1 || n_true > n || expected_degree <= 0 ||
      noise_sd < 0 || rewire_frac < 0 || rewire_frac > 1) stop("Invalid simulation parameters")
  if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  prob <- min(1, expected_degree / (n - 1))
  G <- matrix(sample(c(-1, 0, 1), n * n, replace = TRUE,
                     prob = c(prob / 2, 1 - prob, prob / 2)), n, n)
  diag(G) <- 0
  dimnames(G) <- list(paste0("node", seq_len(n)), paste0("node", seq_len(n)))
  z0 <- runif(n, 1, 2)
  true_nodes <- sort(sample.int(n, n_true))
  true_delta <- numeric(n)
  true_delta[true_nodes] <- runif(n_true, -0.75, -0.25)
  S_true <- pglasso_sensitivity(G, alpha = alpha)
  signal <- as.numeric(S_true %*% true_delta)
  zT <- z0 + signal + rnorm(n, sd = noise_sd * sd(signal))
  G_obs <- G
  edges <- which(G_obs != 0, arr.ind = TRUE)
  nr <- round(nrow(edges) * rewire_frac)
  if (nr) {
    rows <- sample.int(nrow(edges), nr)
    for (r in rows) {
      source <- edges[r, 1]
      target <- edges[r, 2]
      sign <- G_obs[source, target]
      G_obs[source, target] <- 0
      choices <- setdiff(which(G_obs[source, ] == 0), source)
      if (length(choices)) G_obs[source, choices[sample.int(length(choices), 1L)]] <- sign
    }
  }
  list(G_true = G, G_observed = G_obs,
       S = pglasso_sensitivity(G_obs, alpha = alpha),
       z0 = z0, zT = zT, delta_true = true_delta,
       true_nodes = true_nodes, alpha = alpha, seed = seed)
}

#' Run and compare PGLasso with NORSP-style forward selection
#' @param sim Output of pglasso_simulate().
#' @param max_nodes Intervention budget.
#' @return Named list with results and a support-recovery summary.
#' @export
pglasso_compare <- function(sim, max_nodes = 10L) {
  pg <- pglasso_fit(sim$S, sim$z0, sim$zT, max_nodes = max_nodes)
  fs <- norsp_forward(sim$S, sim$z0, sim$zT, max_nodes = max_nodes)
  true <- sim$true_nodes
  summarise <- function(x) {
    tp <- length(intersect(x$selected, true))
    precision <- if (length(x$selected)) tp / length(x$selected) else 0
    recall <- tp / length(true)
    data.frame(method = x$method, n_selected = length(x$selected),
               exact_recovery = setequal(x$selected, true),
               precision = precision, recall = recall,
               F1 = if (precision + recall) 2 * precision * recall / (precision + recall) else 0,
               error_ratio = x$error_ratio)
  }
  list(pglasso = pg, norsp_forward = fs,
       summary = rbind(summarise(pg), summarise(fs)))
}

#' Run a small reproducible comparison
#' @export
pglasso_example <- function() {
  sim <- pglasso_simulate(n = 30, n_true = 3, seed = 42)
  comparison <- pglasso_compare(sim, max_nodes = 6)
  print(comparison$summary)
  invisible(list(simulation = sim, comparison = comparison))
}
