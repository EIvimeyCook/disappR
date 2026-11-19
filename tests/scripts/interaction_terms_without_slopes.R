# Companion to interaction_false_positives.R: when there is NO selective disappearance but individuals differ in their
# ageing trajectories, are the age x ALR interaction terms of Model 4 themselves 'significant', in which direction,
# and how far are Model 4's and Model 2's predicted trajectories from the true one? Paper design (250 individuals,
# complete sampling), random intercept only vs random slope and shape. Run from the package root:
#   Rscript tests/scripts/interaction_terms_without_slopes.R  [> inst/validation/random_slopes/interaction_terms_without_slopes.txt]
root <- getwd(); setwd(file.path(root, "inst", "app")); suppressMessages(source("global.R"))
reps <- as.integer(Sys.getenv("DISAPPR_IFP_REPS", "25")); cores <- as.integer(Sys.getenv("DISAPPR_IFP_CORES", "3"))
structs <- c("random intercept" = "none", "+ random slope and shape" = "correlated_all")
grid <- expand.grid(rep = seq_len(reps), st = names(structs), stringsAsFactors = FALSE)
one <- function(i) {
  g <- grid[i, ]
  d <- simulate_toy_data(list(trait = "paper", n_id = 250, seed = 1000 + g$rep, rho_int = 0, rho_slope = 0, rho_shape = 0))
  tr <- attr(d, "truth"); b <- standardise_data(d, toy_mapping(d))
  f <- fit_model_suite(b$data, b$meta, c("M2", "M4"), "Quadratic", "gaussian", random_slope = structs[[g$st]], include_invalid = TRUE)
  a <- stats::setNames(f$aic$AIC, sub("^Model ", "M", f$aic$Model))
  m4 <- f$fits$M4; cf <- lme4::fixef(m4); V <- as.matrix(stats::vcov(m4))
  it <- grep(":ALR$|^ALR:", names(cf), value = TRUE)
  z <- cf[it] / sqrt(diag(V)[it])
  W <- as.numeric(t(cf[it]) %*% solve(V[it, it]) %*% cf[it])           # joint Wald test of the interaction terms
  m2 <- f$fits$M2
  lrt <- 2 * (as.numeric(stats::logLik(m4)) - as.numeric(stats::logLik(m2)))
  ages <- 1:30
  dev <- function(m) { p <- predict_population_curve(m, f, ages); tc <- toy_true_curve(tr, p$age); mean(abs(100 * (p$fitted - tc) / tc)) }
  data.frame(g, dAIC = a[["M2"]] - a[["M4"]], p_joint_wald = stats::pchisq(W, length(it), lower.tail = FALSE),
             p_lrt = stats::pchisq(lrt, length(it), lower.tail = FALSE), any_term_p05 = any(2 * stats::pnorm(-abs(z)) < 0.05),
             z_lin = unname(z[1]), z_quad = unname(z[2]), dev_M2 = dev(m2), dev_M4 = dev(m4))
}
R <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), one, mc.cores = cores, mc.preschedule = FALSE))
options(width = 220)
S <- do.call(rbind, lapply(split(R, R$st), function(x) data.frame(`Random effects` = x$st[1], n = nrow(x),
  `M4 >2 AIC better (%)` = round(100 * mean(x$dAIC > 2)), `joint Wald p<0.05 (%)` = round(100 * mean(x$p_joint_wald < 0.05)),
  `LRT p<0.05 (%)` = round(100 * mean(x$p_lrt < 0.05)), `either term p<0.05 (%)` = round(100 * mean(x$any_term_p05)),
  `age x ALR z>0 (%)` = round(100 * mean(x$z_lin > 0)), `age2 x ALR z>0 (%)` = round(100 * mean(x$z_quad > 0)),
  `median |z| (age, age2 x ALR)` = paste(round(stats::median(abs(x$z_lin)), 1), round(stats::median(abs(x$z_quad)), 1), sep = ", "),
  `mean |D|% M2` = round(mean(x$dev_M2), 2), `mean |D|% M4` = round(mean(x$dev_M4), 2), check.names = FALSE)))
cat("No selective disappearance (all lifespan correlations 0); individuals differ in intercept, slope and shape (CV 0.2).\n\n")
print(S[order(match(S$`Random effects`, names(structs))), ], row.names = FALSE)
cat("\n|D|%: mean absolute deviation of the model's predicted population trajectory from the true one, ages 1-30.\n")
