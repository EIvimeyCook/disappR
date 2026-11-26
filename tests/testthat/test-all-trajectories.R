# Step 4's "All individuals' trajectories" figure: every individual drawn, simpler functions for few records

test_that("individuals with too few records get a simpler function, and a single age a cross", {
  set.seed(1)
  d <- data.frame(id = rep(c("a", "b", "c", "d"), c(6, 2, 1, 3)), age = c(1:6, 1:2, 4, 2:4))
  d$trait <- 10 + 2 * d$age - 0.3 * d$age^2 + stats::rnorm(nrow(d), 0, 0.2)
  q <- all_individual_trajectories(d, "Quadratic")
  used <- tapply(q$lines$fun_used, q$lines$id, unique)
  expect_identical(used[["a"]], "Quadratic"); expect_identical(used[["d"]], "Quadratic")   # 3 records: an exact quadratic
  expect_identical(used[["b"]], "Linear")                                                    # 2 records: a line
  expect_identical(unique(q$points$id), "c")                                                 # 1 record: a cross
  expect_equal(as.numeric(q$counts), c(2, 1, 1))
  k <- all_individual_trajectories(d, "Cubic")
  ku <- tapply(k$lines$fun_used, k$lines$id, unique)
  expect_identical(c(ku[["a"]], ku[["d"]], ku[["b"]]), c("Cubic", "Quadratic", "Linear"))
  # each line spans the individual's own ages and passes the 2-record individual's points exactly
  b <- q$lines[q$lines$id == "b", ]
  expect_equal(range(b$age), c(1, 2)); expect_equal(b$fitted[c(1, nrow(b))], d$trait[d$id == "b"], tolerance = 1e-8)
})

test_that("count traits are drawn on the log scale and the figure leaves the step-4 fits unchanged", {
  set.seed(2)
  d <- data.frame(id = rep(1:30, each = 5), age = rep(1:5, 30)); d$trait <- stats::rpois(nrow(d), exp(1 + 0.2 * d$age - 0.03 * d$age^2))
  before <- fit_individual_function(d, "Quadratic", link = "log")
  z <- all_individual_trajectories(d, "Logarithmic", link = "log")
  expect_true(all(z$lines$fitted > 0)); expect_true(all(z$lines$fun_used %in% c("Logarithmic", "Linear")))
  expect_identical(fit_individual_function(d, "Quadratic", link = "log"), before)
})
