# Replication of the age-dependent analyses of Sanghvi et al. 2021 (Evol. Ecol. 35: 367-385) on the bundled seed-beetle
# data: does a model with group-specific age curves beat one shared curve (with a group intercept)? (a) mgcv::gam by ML
# with the published structure (beetle and container random effects, adult lifespan as a fixed term); (b) disappR's
# group comparison (G, GS, GI, S, I) with the example's settings. Run from the package root:
#   Rscript tests/scripts/gamm_replication_beetles.R  [> inst/validation/gamm/replication_results_beetles.txt]
# Published: fecundity, density-specific curves best (dAIC 27.5 to the next model); flying ability, sex-specific
# trajectories best (dAIC 21.8; ordinal 'ocat' family in the paper, Gaussian here); propensity, sex-specific (dAIC 118).
suppressPackageStartupMessages(library(mgcv))
if (!file.exists("DESCRIPTION")) stop("Run from the package root.")
for (f in c("utils.R", "provenance.R", "simulation.R", "import.R", "validation.R", "data-preparation.R", "formulas.R", "model-fit.R",
            "model-comparison.R", "diagnostics.R", "trajectories.R", "disappearance.R", "missingness.R", "interpretation.R",
            "effects.R", "evidence.R", "ui-helpers.R", "export.R", "gamm.R", "api.R")) suppressMessages(source(file.path("R", f), local = TRUE))
errors <- character(0)
line <- function(...) cat(..., "\n", sep = "")
rd <- function(f) utils::read.csv(file.path("inst", "app", "data", f), stringsAsFactors = TRUE)
authors <- function(label, trait, age, grp, d, fam, k = 8, lifespan = TRUE) {
  d$id <- factor(d$id); d$container <- factor(d$container)
  re <- paste(if (lifespan) "+ lifespan", "+ s(id, bs = 're') + s(container, bs = 're')")
  f_shared <- sprintf("%s ~ s(%s, k = %d) + %s %s", trait, age, k, grp, re)
  f_group <- sprintf("%s ~ s(%s, by = %s, k = %d) + %s %s", trait, age, grp, k, grp, re)
  fits <- lapply(list(shared = f_shared, group = f_group), function(f)
    tryCatch(gam(stats::as.formula(f), data = d, family = fam, method = "ML"), error = function(e) { errors <<- c(errors, paste(label, conditionMessage(e))); NULL }))
  if (any(vapply(fits, is.null, logical(1)))) return(invisible(NULL))
  a <- vapply(fits, stats::AIC, numeric(1))
  line(sprintf("  (a) mgcv::gam (ML, %s%s), %d records: one shared curve AIC %.1f; %s-specific curves AIC %.1f; dAIC %.1f in favour of %s",
               fam$family, if (lifespan) ", lifespan term" else ", no lifespan term", nrow(fits$shared$model), a[["shared"]], grp, a[["group"]], abs(diff(a)), if (a[["group"]] < a[["shared"]]) paste0(grp, "-specific curves") else "one shared curve"))
}
app <- function(label, id) {
  ex <- disappr_example(id)
  b <- standardise_data(if (!is.null(ex$subset)) ex$data else ex$data, ex$mapping)
  gd <- gamm_data(b$data, b$meta)
  s <- utils::modifyList(utils::modifyList(gamm_defaults(), EXAMPLES[[id]]$gamm), list(family = gamm_suggest_family(EXAMPLES[[id]]$family)))
  gc <- tryCatch(gamm_group_comparison(gd, s), error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  if (!isTRUE(gc$ok)) { errors <<- c(errors, paste(label, gc$message)); return(invisible(NULL)) }
  tab <- gc$table[is.finite(gc$table$AIC), , drop = FALSE]
  line(sprintf("  (b) disappR group comparison (%s, %s family): %s", gc$engine, s$family,
               paste(sprintf("%s %+.1f", tab$Model, tab$dAIC)[order(tab$AIC)], collapse = ", ")))
  line("      ", gc$verdict)
}
beetles <- rd("sanghvi_2021_beetle_flight.csv")
virgins <- beetles[beetles$mating == "Virgin", , drop = FALSE]
fec <- rd("sanghvi_2021_beetle_fecundity.csv")

line("==== Fecundity of mated females by developmental density (published: density-specific curves best, dAIC 27.5) ====")
authors("fecundity", "eggs", "day", "density", fec, gaussian(), lifespan = FALSE)
authors("fecundity", "eggs", "day", "density", fec, gaussian())
authors("fecundity", "eggs", "day", "density", fec, nb())
app("fecundity", "sanghvi2021_fecundity")
line("")
line("==== Flying ability of virgins by sex (published: sex-specific trajectories best, dAIC 21.8, ordinal family) ====")
authors("ability", "ability", "age", "sex", virgins, gaussian())
app("ability", "sanghvi2021_flight")
line("")
line("==== Flying propensity of virgins by sex (published: sex-specific trajectories best, dAIC 118) ====")
authors("propensity", "propensity", "age", "sex", virgins, nb())
app("propensity", "sanghvi2021_propensity")

if (length(errors)) { cat("\nErrors:\n", paste("-", errors, collapse = "\n"), "\n"); quit(status = 1) }
cat("\nReplication script finished without errors.\n")
