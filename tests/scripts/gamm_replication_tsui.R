# Replication of the 12 Tsui et al. 2026 (Evol. Lett.) GAMMs with data in the deposit (the telomere data are not), on
# the bundled files: (a) the authors' model G formula refitted with mgcv::gam (REML), (b) the same analysis through the
# disappR engine as a user would run it with the example's settings. Compares the age smooth (edf, statistic, p) and
# the lifespan, ALR or AFR term. Run from the package root:
#   Rscript tests/scripts/gamm_replication_tsui.R  [> inst/validation/gamm/replication_results_tsui.txt]
suppressPackageStartupMessages(library(mgcv))
if (!file.exists("DESCRIPTION")) stop("Run from the package root.")
for (f in c("utils.R", "provenance.R", "simulation.R", "import.R", "validation.R", "data-preparation.R", "formulas.R", "model-fit.R",
            "model-comparison.R", "diagnostics.R", "trajectories.R", "disappearance.R", "missingness.R", "interpretation.R",
            "effects.R", "evidence.R", "ui-helpers.R", "export.R", "gamm.R", "api.R")) suppressMessages(source(file.path("R", f), local = TRUE))
errors <- character(0)
line <- function(...) cat(..., "\n", sep = "")
rd <- function(f) utils::read.csv(file.path("inst", "app", "data", f), stringsAsFactors = FALSE)
fac <- function(d, v) { for (x in v) d[[x]] <- factor(d[[x]]); d }
age_row <- function(fit) { st <- summary(fit)$s.table; st[grep("^s\\(age", rownames(st))[1], , drop = FALSE] }
term_row <- function(fit, lab) {
  sm <- summary(fit)
  if (lab %in% rownames(sm$p.table)) return(sprintf("%s %s (p %s)", lab, format_num(sm$p.table[lab, 1], 3), format_p(sm$p.table[lab, 4])))
  i <- grep(paste0("^s\\(", lab, "\\)"), rownames(sm$s.table))
  if (length(i)) return(sprintf("s(%s) edf %s (p %s)", lab, format_num(sm$s.table[i[1], 1], 3), format_p(sm$s.table[i[1], 4])))
  ""
}
show <- function(label, fit, n, extra) {
  a <- age_row(fit)
  line(sprintf("  %-34s %5d rows  age smooth edf %6.3f  stat %8.2f  p %-8s  %s", label, n, a[1, 1], a[1, 3], format_p(a[1, 4]), extra))
}
authors <- function(label, formula, data, family, proxy) {
  t0 <- Sys.time()
  fit <- tryCatch(gam(stats::as.formula(formula), data = data, method = "REML", family = family), error = function(e) { errors <<- c(errors, paste(label, conditionMessage(e))); NULL })
  if (!is.null(fit)) show("(a) authors' formula, mgcv::gam", fit, nrow(fit$model), term_row(fit, proxy))
  invisible(fit)
}
app <- function(label, id, trait = NULL, gamm = list(), proxy_label = "") {
  ex <- disappr_example(id)
  mp <- ex$mapping; if (!is.null(trait)) mp$trait <- trait
  b <- standardise_data(ex$data, mp)
  gd <- gamm_data(b$data, b$meta)
  s <- utils::modifyList(utils::modifyList(gamm_defaults(), EXAMPLES[[id]]$gamm), c(gamm, list(family = gamm_suggest_family(EXAMPLES[[id]]$family))))
  r <- tryCatch(gamm_fit(gd, s), error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  if (!isTRUE(r$ok)) { errors <<- c(errors, paste(label, r$message)); return(invisible(NULL)) }
  show(sprintf("(b) disappR engine (%s)", r$engine), r$fit, r$n, term_row(r$fit, proxy_label))
  invisible(r)
}
run <- function(title, ...) { line(""); line(title) }

physio <- fac(rd("tsui_2026_warbler_body_mass.csv"), c("Observer", "BirdID", "SexEstimate", "newstat", "birthyear"))
physio$newlifespan <- physio$newlifespan
re3 <- "+ s(birthyear, bs = 're') + s(Observer, bs = 're') + s(BirdID, bs = 're')"
run("Model 2: wing length (Gaussian; lifespan smooth)")
authors("wing", paste("WingLength_z ~ s(age_year, k = 15) + s(newlifespan) + RightTarsus_z + newstat + SexEstimate + avg_invert", re3), physio, gaussian(), "newlifespan")
app("wing", "tsui_mass", "WingLength_z", list(k = 15, proxy_form = "smooth"), "LS")
run("Model 3: tarsus length (Gaussian; lifespan smooth)")
authors("tarsus", paste("RightTarsus_z ~ s(age_year, k = 10) + s(newlifespan) + SexEstimate + avg_invert + newstat", re3), physio, gaussian(), "newlifespan")
app("tarsus", "tsui_mass", "RightTarsus_z", list(k = 10, proxy_form = "smooth"), "LS")
run("Model 4: fat score (Poisson; lifespan linear)")
authors("fat", paste("FatScore ~ s(age_year, k = 15) + newlifespan + SexEstimate + newstat + avg_invert + CatchTime_mins", re3), physio, poisson(), "newlifespan")
app("fat", "tsui_mass", "FatScore", list(k = 15, proxy_form = "linear", family = "poisson"), "LS")

oxy <- fac(rd("tsui_2026_warbler_oxidative.csv"), c("BirdID", "SexEstimate", "newstat", "birthyear"))
re2 <- "+ s(birthyear, bs = 're') + s(BirdID, bs = 're')"
run("Model 5: oxidative damage, ROMs (Gaussian; lifespan smooth)")
authors("roms", paste("sum_roms_z ~ s(age_year, k = 10) + s(newlifespan) + SexEstimate + newstat + avg_invert + CatchTime_mins", re2), oxy, gaussian(), "newlifespan")
app("roms", "tsui_oxidative", NULL, list(), "LS")
run("Model 6: antioxidant capacity, OXY (Gaussian; lifespan linear)")
authors("oxy", paste("sum_oxy_z ~ s(age_year, k = 10) + newlifespan + SexEstimate + newstat + avg_invert + CatchTime_mins", re2), oxy, gaussian(), "newlifespan")
app("oxy", "tsui_oxidative", "sum_oxy_z", list(proxy_form = "linear"), "LS")

blood <- fac(rd("tsui_2026_warbler_blood.csv"), c("Observer", "BirdID", "SexEstimate", "newstat", "birthyear"))
run("Model 7: haematocrit (Gaussian; lifespan smooth)")
authors("haematocrit", paste("Cells_z ~ s(age, k = 15) + s(newlifespan) + SexEstimate + newstat + avg_invert + CatchTime_mins", re3), blood, gaussian(), "newlifespan")
app("haematocrit", "tsui_blood", NULL, list(), "LS")
run("Model 8: white blood cells (Gaussian; lifespan linear)")
authors("wbc", paste("wbc_z ~ s(age, k = 15) + newlifespan + SexEstimate + newstat + avg_invert + CatchTime_mins", re3), blood, gaussian(), "newlifespan")
app("wbc", "tsui_blood", "wbc_z", list(proxy_form = "linear"), "LS")

mal <- fac(rd("tsui_2026_warbler_malaria.csv"), c("BirdID", "SexEstimate", "newstat", "birthyear"))
run("Model 9: malaria infection (binomial; lifespan linear)")
authors("malaria", paste("TestResult ~ s(age, k = 8) + SexEstimate + newlifespan + newstat + avg_invert", re2), mal, binomial(), "newlifespan")
app("malaria", "tsui_malaria", NULL, list(), "LS")

prov <- fac(rd("tsui_2026_warbler_provisioning.csv"), c("Observer", "BirdID", "SexEstimate", "newstat", "WatchType", "birthyear"))
run("Model 11: provisioning rate (Gaussian; lifespan smooth)")
authors("provisioning", "prate_z ~ s(age, k = 10) + s(newlifespan) + SexEstimate + BroodSize + nr_helpers + newstat + avg_invert + s(birthyear, bs = 're') + s(WatchType, bs = 're') + s(BirdID, bs = 're') + s(Observer, bs = 're')",
        prov, gaussian(), "newlifespan")
app("provisioning", "tsui_provisioning", NULL, list(), "LS")

rep <- fac(rd("tsui_2026_warbler_reproduction.csv"), c("BirdID", "SexEstimate", "newstat", "birthyear"))
run("Model 12: annual reproductive success (Poisson; first-breeding year linear, last-breeding year smooth)")
authors("ars", paste("ars ~ s(parent_age, k = 15) + minrepro + s(maxrepro) + SexEstimate + flno + newstat", re2), rep, poisson(), "maxrepro")
authors("ars, ages instead of years", paste("ars ~ s(parent_age, k = 15) + minreproage + s(maxreproage) + SexEstimate + flno + newstat", re2), rep, poisson(), "maxreproage")
app("ars", "tsui_reproduction", NULL, list(), "ALR")

surv <- fac(rd("tsui_2026_warbler_survival.csv"), c("SexEstimate", "statbefore", "occasionyear", "birthyear"))
run("Model 13: annual survival (binomial; no bird random effect, no proxy)")
authors("survival", "survival ~ s(age, k = 19) + SexEstimate + avg_bug + yearlytq + statbefore + s(occasionyear, bs = 're') + s(birthyear, bs = 're')", surv, binomial(), "")
app("survival", "tsui_survival", NULL, list(), "")

if (length(errors)) { cat("\nErrors:\n", paste("-", errors, collapse = "\n"), "\n"); quit(status = 1) }
cat("\nReplication script finished without errors.\n")
