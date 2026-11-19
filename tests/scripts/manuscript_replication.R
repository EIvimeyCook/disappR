# The empirical results of Sanghvi et al. 2026 (Methods in Ecology and Evolution; Results, 'Empirical implementation',
# and Supplementary section S5) reproduced through the app's engine with the bundled examples' settings: the fruit fly
# (example 'fly') and great tit laying-date (example 'bouwhuis') datasets, Models 1-5. Run from the package root:
#   Rscript tests/scripts/manuscript_replication.R  [> inst/validation/manuscript_replication.txt]
# Reference values: the authors' analysis code ('sel dis code.Rmd', OSF kevnm, section 'Empirical data') run on the
# authors' data files (identical to the bundled ones): glmmTMB, ML. The manuscript text rounds from an earlier run
# (fly: 38.4, 24.1, 28.0, 2.7; great tit: 11.8, 9.5, 53.7, 48.7); these are the values the code gives today.
# Exits with status 1 if any dAIC differs from the authors' code by more than 0.05 or Model 4 is not best.
app_dir <- file.path("inst", "app")
if (!dir.exists(app_dir)) stop("Run from the package root.")
owd <- setwd(app_dir); on.exit(setwd(owd), add = TRUE)
suppressMessages(source("global.R"))
code <- list(fly = c(M1 = 38.34, M2 = 24.05, M3 = 27.27, M5 = 2.67), bouwhuis = c(M1 = 11.77, M2 = 9.43, M3 = 53.29, M5 = 48.62))
text <- list(fly = c(M1 = 38.4, M2 = 24.1, M3 = 28.0, M5 = 2.7), bouwhuis = c(M1 = 11.8, M2 = 9.5, M3 = 53.7, M5 = 48.7))
bad <- character(0)
for (key in names(code)) {
  ex <- EXAMPLES[[key]]
  b <- standardise_data(load_example_file(ex$file), ex$mapping)
  r <- fit_model_suite(b$data, b$meta, c("M1", "M2", "M3", "M4", "M5"), ex$age_function, ex$family)
  a <- stats::setNames(r$aic$AIC, sub("^Model ", "M", r$aic$Model))
  d <- a - a[["M4"]]
  cat(sprintf("%s (%s, %s; %d records, %d individuals): best %s\n", key, ex$family, ex$age_function, sum(is.finite(b$data$trait)), length(unique(b$data$id)), names(which.min(a))))
  for (m in names(code[[key]])) {
    cat(sprintf("  Model %s: dAIC app %6.2f | authors' code %6.2f | manuscript text %5.1f\n", sub("M", "", m), d[[m]], code[[key]][[m]], text[[key]][[m]]))
    if (abs(d[[m]] - code[[key]][[m]]) > 0.05) bad <- c(bad, sprintf("%s %s: %.2f vs %.2f", key, m, d[[m]], code[[key]][[m]]))
  }
  if (!identical(names(which.min(a)), "M4")) bad <- c(bad, paste(key, "best model is not Model 4"))
}
if (length(bad)) { cat("\nMISMATCH:\n ", paste(bad, collapse = "\n  "), "\n"); quit(status = 1) }
cat("\nThe app reproduces the manuscript's empirical model comparisons (to 0.05 AIC of the authors' code).\n")
