# How often do the interaction models (Model 4: age terms x ALR; Model 5: mean age x delta-age terms) beat their
# additive counterparts (Models 2 and 3) by more than 2 AIC, with and without random slope (and shape) terms?
# Design of Sanghvi et al. 2026 (complete-sampling 'CS' simulation, step 1's paper simulation): 250 individuals,
# lifespan 25 (SD 5), quadratic trajectories whose intercept, slope and shape vary among individuals (CV 0.2) and are
# correlated with lifespan as set. Run from the package root:
#   Rscript tests/scripts/interaction_false_positives.R  [> inst/validation/random_slopes/interaction_false_positives.txt]
# DISAPPR_IFP_REPS (default 25) replicates per scenario; DISAPPR_IFP_CORES (default 3) parallel workers.
root <- getwd(); setwd(file.path(root, "inst", "app")); suppressMessages(source("global.R"))
reps <- as.integer(Sys.getenv("DISAPPR_IFP_REPS", "25")); cores <- as.integer(Sys.getenv("DISAPPR_IFP_CORES", "3"))
scen <- list("none (no selective disappearance)" = c(0, 0, 0), "intercept only (age-independent)" = c(0.5, 0, 0),
             "slope (age-dependent)" = c(0, 0.5, 0), "shape (age-dependent)" = c(0, 0, 0.5))
structs <- c("random intercept" = "none", "+ random slope" = "correlated", "+ random slope and shape" = "correlated_all")
grid <- expand.grid(rep = seq_len(reps), sc = names(scen), st = names(structs), stringsAsFactors = FALSE)
one <- function(i) {
  g <- grid[i, ]; r <- scen[[g$sc]]
  d <- simulate_toy_data(list(trait = "paper", n_id = 250, seed = 1000 + g$rep, rho_int = r[1], rho_slope = r[2], rho_shape = r[3]))
  b <- standardise_data(d, toy_mapping(d))
  f <- tryCatch(fit_model_suite(b$data, b$meta, c("M2", "M3", "M4", "M5"), "Quadratic", "gaussian", random_slope = structs[[g$st]], include_invalid = TRUE),
                error = function(e) NULL)
  if (is.null(f) || !isTRUE(f$ok)) return(data.frame(g, M4_vs_M2 = NA, M5_vs_M3 = NA, invalid = NA))
  a <- stats::setNames(f$aic$AIC, sub("^Model ", "M", f$aic$Model))
  inv <- if (!is.null(f$aic$Eligible)) sum(!f$aic$Eligible) else NA
  data.frame(g, M4_vs_M2 = a[["M2"]] - a[["M4"]], M5_vs_M3 = a[["M3"]] - a[["M5"]], invalid = inv)
}
t0 <- Sys.time()
R <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), one, mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("%d datasets x 3 random structures in %.0f min (%d replicates per scenario)\n\n", reps * length(scen), as.numeric(difftime(Sys.time(), t0, units = "mins")), reps))
S <- do.call(rbind, lapply(split(R, list(R$sc, R$st), drop = TRUE), function(x) data.frame(
  Scenario = x$sc[1], `Random effects` = x$st[1], n = sum(is.finite(x$M4_vs_M2)),
  `Model 4 beats Model 2 by >2 AIC (%)` = round(100 * mean(x$M4_vs_M2 > 2, na.rm = TRUE)),
  `median dAIC (M2 - M4)` = round(stats::median(x$M4_vs_M2, na.rm = TRUE), 1),
  `Model 5 beats Model 3 by >2 AIC (%)` = round(100 * mean(x$M5_vs_M3 > 2, na.rm = TRUE)),
  `fits flagged invalid` = sum(x$invalid, na.rm = TRUE), check.names = FALSE)))
options(width = 250)
S <- S[order(match(S$Scenario, names(scen)), match(S$`Random effects`, names(structs))), ]
print(S, row.names = FALSE)
cat("\nIn the first two scenarios a declared age-dependent selective disappearance is a false positive; in the last two it is a true one (power).\n")
utils::write.csv(R, file.path(root, "inst", "validation", "random_slopes", "interaction_false_positives_replicates.csv"), row.names = FALSE)
