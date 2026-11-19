test_that("simulations are reproducible and leave R's random-number stream unchanged", {
  set.seed(99)
  before <- .Random.seed
  a <- simulate_toy_data(list(n_id = 40, seed = 3))
  expect_identical(before, .Random.seed)
  b <- simulate_toy_data(list(n_id = 40, seed = 3))
  tr <- toy_mapping(a)$trait
  expect_identical(a$age, b$age)
  expect_identical(a[[tr]], b[[tr]])
})

test_that("the paper simulation covers the Methods' 27 correlation scenarios, with coefficients uncorrelated with each other", {
  vals <- c(-0.5, 0, 0.5)
  for (r0 in vals) for (r1 in vals) for (r2 in vals) {
    d <- simulate_toy_data(list(trait = "paper", n_id = 60, seed = 2, rho_int = r0, rho_slope = r1, rho_shape = r2))
    expect_gt(nrow(d), 0)
    expect_equal(paper_rho(attr(d, "truth")$cfg), c(r0, r1, r2))
  }
  # without the grid, the correlations follow the selective-disappearance choice
  expect_equal(paper_rho(list(sd_type = "both", sd_dir = -1)), c(-0.5, -0.5, 0))
  expect_equal(paper_rho(list(sd_type = "dependent", sd_dir = 1)), c(0, 0.5, 0))
  # joint normal as in Eq. 2: given correlations with lifespan, the coefficients do not correlate with each other
  rho <- c(0.5, 0.5, -0.5); set.seed(1); n <- 50000
  z <- stats::rnorm(n); U <- matrix(stats::rnorm(n * 3), n, 3) %*% chol(diag(3) - outer(rho, rho))
  C <- stats::cor(cbind(z, sapply(1:3, function(j) rho[j] * z + U[, j])))
  expect_equal(unname(C[1, 2:4]), rho, tolerance = 0.02)
  expect_lt(max(abs(C[2:4, 2:4][upper.tri(diag(3))])), 0.02)
})
