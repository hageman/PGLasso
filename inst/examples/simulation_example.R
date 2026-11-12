# Run after installing package with install.packages(".", repos = NULL, type = "source")
library(PGLasso)

set.seed(2026)
sim <- pglasso_simulate(n = 50, expected_degree = 4,
                        n_true = 3, alpha = 0.9,
                        noise_sd = 0.05, rewire_frac = 0.05,
                        seed = 2026)
cat("True intervention nodes:", paste(sim$true_nodes, collapse = ", "), "\n")
results <- pglasso_compare(sim, max_nodes = 10)
print(results$summary)
cat("PGLasso selected:", paste(results$pglasso$selected, collapse = ", "), "\n")
cat("Forward selected:", paste(results$norsp_forward$selected, collapse = ", "), "\n")

# Visualize error against model size (method-specific candidate paths).
plot(results$norsp_forward$path$size,
     results$norsp_forward$path$error_ratio,
     type = "b", pch = 1, xlab = "Candidate number of selected nodes",
     ylab = "Graph-based error ratio", ylim = c(0, 1),
     main = "PGLasso vs NORSP-style forward selection")
with(results$pglasso$path,
     points(size, error_ratio, pch = 16))
legend("topright", legend = c("Forward selection", "PGLasso"),
       pch = c(1, 16), bty = "n")

# To reproduce manuscript benchmark statistics, additional decisions regarding
# network generation, perturbation distributions, and replicates must be
# matched to the original analysis scripts. This illustrative example does not
# claim to reproduce the paper's published 87.3% versus 67.0% result.
