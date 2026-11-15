# Does the GAMM's age-dependent proxy form (s(age, by = ALR)) pick up individual differences in ageing as selective
# disappearance? Paper simulation (step 1) with no selective disappearance, 250 individuals; the proxy comparison of
# step 7 with a random intercept, a random intercept and slope, and a random smooth per individual. Run from the
# package root:  Rscript tests/scripts/gamm_proxy_false_positives.R   (the random smooths take about 20 minutes)
root <- getwd(); setwd(file.path(root, "inst", "app")); suppressMessages(source("global.R"))
reps <- c(ri = 12L, rs = 12L, fs = 8L)
one <- function(i, ind) {
  d <- simulate_toy_data(list(trait = "paper", n_id = 250, seed = 1000 + i, rho_int = 0, rho_slope = 0, rho_shape = 0))
  b <- standardise_data(d, toy_mapping(d)); gd <- gamm_data(b$data, b$meta)
  s <- utils::modifyList(gamm_defaults(), list(individual = ind, covars = FALSE, random = FALSE))
  pc <- tryCatch(gamm_proxy_comparison(gd, s), error = function(e) NULL)
  data.frame(rep = i, individual = ind, verdict = if (is.null(pc)) NA_character_ else pc$verdicts[["ALR"]] %||% NA_character_)
}
g <- do.call(rbind, lapply(names(reps), function(k) data.frame(i = seq_len(reps[[k]]), ind = k)))
R <- do.call(rbind, parallel::mclapply(seq_len(nrow(g)), function(k) one(g$i[k], g$ind[k]), mc.cores = 4, mc.preschedule = FALSE))
print(table(Individual = R$individual, `ALR verdict` = R$verdict, useNA = "ifany"))
cat("\nWith no selective disappearance, 'age_dependent' is a false positive.\n")
