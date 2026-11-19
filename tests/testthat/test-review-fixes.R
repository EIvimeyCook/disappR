# Fixes from an independent review of 0.25.2 (November 2026)

test_that("the random-effects description names every age term of the ageing function", {
  m <- list(map = list(id = "ID"), has_group = FALSE)
  expect_identical(random_display(m, "correlated_all", basis_label("Quadratic")), "(1 + age + age² | ID)")
  expect_match(random_display(m, "uncorrelated_all", basis_label("Cubic")), "(0 + age³ | ID)", fixed = TRUE)
  expect_identical(random_display(m, "correlated", basis_label("Quadratic")), "(1 + age | ID)")
  info <- model_info_content("M2", list(age_function = "Quadratic", random_slope = "correlated_all"), m)
  expect_match(info$random, "age²", fixed = TRUE)
})

test_that("continuous proportions get beta advice; successes out of trials keep binomial", {
  set.seed(1)
  expect_identical(suggest_family(stats::rbeta(200, 40, 50))$family, "beta")
  expect_identical(suggest_family(stats::rbinom(200, 5, 0.4) / 5)$family, "binomial")
  n <- sample(2:12, 200, TRUE)
  expect_identical(suggest_family(stats::rbinom(200, n, 0.3) / n)$family, "binomial")
  z <- suggest_family(c(0, 1, stats::rbeta(200, 40, 50)))   # exact 0 or 1: beta cannot be fitted as is
  expect_identical(z$family, "gaussian"); expect_match(z$text, "Smithson")
  expect_identical(suggest_family(stats::rbinom(200, 1, 0.4))$family, "binomial")
})

test_that("the exported GAMM script computes concurvity as the app does, without stopping on rank-deficient fits", {
  skip_if_not_installed("mgcv")
  d <- do.call(rbind, lapply(1:40, function(i) { a <- 1:6; data.frame(id = i, age = a, trait = sin(a / 2) + stats::rnorm(1) + stats::rnorm(6, 0, .3)) }))
  b <- standardise_data(d, list(id = "id", age = "age", trait = "trait", alr = "__AUTO_LAST__", life = "", entry = "__AUTO_FIRST__",
                                covars = character(0), cov_factor = character(0), cov_int = character(0), group = "", nested = TRUE,
                                random = character(0), censor = "", censor_value = "", condition = "", trials = "", start_mode = "afr",
                                start_age = NA_real_, age_round = NA_real_, cov_age = character(0)))
  r <- gamm_fit(gamm_data(b$data, b$meta), utils::modifyList(gamm_defaults(), list(k = 5)))
  code <- paste(gamm_code(r), collapse = "\n")
  expect_match(code, "try(print(concurvity(fit, full = FALSE)))", fixed = TRUE)
  expect_false(grepl("concurvity(fit, full = TRUE)", code, fixed = TRUE))
})

test_that("non-linear (nlme) fits skip R2 quietly", {
  skip_if_not_installed("nlme")
  f <- nlme::nlme(height ~ SSasymp(age, Asym, R0, lrc), data = as.data.frame(datasets::Loblolly), fixed = Asym + R0 + lrc ~ 1, groups = ~ Seed,
                  random = Asym ~ 1, start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
  out <- utils::capture.output(v <- r2_pair(f), type = "message")
  expect_true(all(is.na(v))); expect_length(out, 0)
})
