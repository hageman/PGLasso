test_that("PRINCE matches the identity network case", {
  G <- matrix(0, 4, 4)
  expect_equal(pglasso_sensitivity(G, 0.5), 0.5 * diag(4))
})

test_that("undirected and signed sensitivity matrices are finite", {
  G <- matrix(c(0, 1, 0, -1, 0, 1, 0, -1, 0), 3, 3)
  S <- pglasso_sensitivity(G, 0.9)
  expect_equal(dim(S), c(3L, 3L))
  expect_true(all(is.finite(S)))
})

test_that("PGLasso and baseline return feasible estimates", {
  skip_if_not_installed("glmnet")
  sim <- pglasso_simulate(n = 15, n_true = 2, seed = 4)
  pg <- pglasso_fit(sim$S, sim$z0, sim$zT, max_nodes = 5)
  fs <- norsp_forward(sim$S, sim$z0, sim$zT, max_nodes = 5)
  expect_true(all(sim$z0 + pg$delta >= -1e-6))
  expect_true(all(sim$z0 + fs$delta >= -1e-6))
  expect_lte(length(pg$selected), 5)
  expect_lte(length(fs$selected), 5)
  expect_true(is.finite(pg$error_ratio))
  expect_true(is.finite(fs$error_ratio))
})

test_that("simulation deterministic by seed", {
  a <- pglasso_simulate(seed = 13)
  b <- pglasso_simulate(seed = 13)
  expect_identical(a$G_true, b$G_true)
  expect_identical(a$delta_true, b$delta_true)
})
