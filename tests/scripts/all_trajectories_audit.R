# Audit of step 4's "All individuals' trajectories" figure (all_individual_trajectories): every bundled example (each
# trait of multi-trait files) and a spread of simulations, under every ageing function on both scales. Checks that it
# runs, that every individual with a trait value is drawn once (line or crosses), that curves are finite, that the
# function used never exceeds what the records allow, and the time taken. Run from the package root:
#   Rscript tests/scripts/all_trajectories_audit.R      (exits with status 1 on any problem)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
funs <- c(AGE_FUNCTIONS, A3_NONLINEAR)
# distinct ages a function needs (its number of parameters); the raw-scale fallback is a line
need <- c(Linear = 2, Quadratic = 3, Cubic = 4, Logarithmic = 2, "Asymptotic exponential" = 2, "Linear (raw scale)" = 2)
need[[A3_NONLINEAR]] <- 2
cap <- function(d, n = 2000L) { ids <- unique(d$id); if (length(ids) <= n) d else d[d$id %in% ids[unique(round(seq(1, length(ids), length.out = n)))], , drop = FALSE] }
check_one <- function(label, d) {
  d <- cap(d)
  dd <- d[is.finite(d$trait) & is.finite(d$age), , drop = FALSE]
  ids <- unique(as.character(dd$id)); n_age <- tapply(dd$age, as.character(dd$id), function(a) length(unique(a)))
  out <- list()
  for (lk in c("identity", "log")) for (fn in funs) {
    if (identical(lk, "log") && any(dd$trait < 0)) next
    t0 <- Sys.time()
    z <- tryCatch(all_individual_trajectories(d, fn, link = lk), error = function(e) e)
    sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    prob <- character(0)
    if (inherits(z, "error")) prob <- paste("error:", conditionMessage(z)) else {
      drawn <- c(unique(z$lines$id), unique(z$points$id))
      if (anyDuplicated(drawn)) prob <- c(prob, "an individual drawn twice")
      miss <- setdiff(ids, drawn); if (length(miss)) prob <- c(prob, sprintf("%d individuals not drawn", length(miss)))
      if (nrow(z$lines) && !all(is.finite(z$lines$fitted))) prob <- c(prob, "non-finite curve values")
      if (nrow(z$lines)) {
        u <- tapply(z$lines$fun_used, z$lines$id, `[`, 1)
        over <- names(u)[n_age[names(u)] < need[u]]
        if (length(over)) prob <- c(prob, sprintf("%d individuals given a function their records cannot support", length(over)))
      }
      pt_multi <- if (nrow(z$points)) sum(n_age[unique(z$points$id)] >= 2) else 0
      if (sum(z$counts) != length(ids)) prob <- c(prob, "counts do not add up")
    }
    out[[length(out) + 1]] <- data.frame(Data = label, Scale = lk, Function = fn, Individuals = length(ids),
      Chosen = if (inherits(z, "error")) NA else z$counts[[fn]], Simpler = if (inherits(z, "error")) NA else sum(z$counts) - z$counts[[fn]] - z$counts[["points"]],
      Crosses = if (inherits(z, "error")) NA else z$counts[["points"]], `Crosses with 2+ ages` = if (inherits(z, "error")) NA else pt_multi,
      Seconds = round(sec, 1), Problems = paste(prob, collapse = "; "), check.names = FALSE, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
res <- list()
# every bundled example, and every extra trait of the multi-trait files
for (k in names(EXAMPLES)) {
  ex <- EXAMPLES[[k]]; x <- disappr_example(k)
  traits <- unique(c(ex$mapping$trait, names(ex$traits)))
  for (tr in traits) {
    m <- ex$mapping; tp <- example_trait_preset(ex, tr); m$trait <- tr
    if (!is.null(tp$covars)) { m$covars <- tp$covars; m$cov_factor <- tp$cov_factor %||% m$cov_factor }
    b <- tryCatch(standardise_data(x$data, m), error = function(e) NULL)
    if (is.null(b)) { res[[length(res) + 1]] <- data.frame(Data = paste(k, tr), Scale = NA, Function = NA, Individuals = NA, Chosen = NA, Simpler = NA, Crosses = NA, `Crosses with 2+ ages` = NA, Seconds = NA, Problems = "data did not prepare", check.names = FALSE); next }
    res[[length(res) + 1]] <- check_one(paste0(k, if (!identical(tr, ex$mapping$trait)) paste0(":", tr) else ""), b$data)
    cat(".", sep = "")
  }
}
cat("\n")
# simulations: every trait type, sampling design and ageing form, the paper simulation, and small or sparse data
sims <- list(
  "mass, quadratic, complete" = list(), "Poisson counts, MCAR" = list(trait = "count_pois", missingness = "mcar"),
  "NB counts, missing when old" = list(trait = "count_nb", missingness = "mwo"), "ZINB counts, missing when young" = list(trait = "count_zinb", missingness = "mwy"),
  "cubic, trait-dependent missing" = list(form = "Cubic", missingness = "trait"), "logarithmic, condition missing" = list(form = "Logarithmic", missingness = "condition"),
  "asymptotic, individual AFR" = list(form = "Asymptotic exponential", afr_mode = "individual", sa_type = "both"),
  "linear, short-lived (3 occasions)" = list(form = "Linear", mean_ls = 3), "paper simulation" = list(trait = "paper"),
  "families and diet" = list(groups = TRUE, diet = TRUE), "30 individuals" = list(n_id = 30), "large rate variation" = list(rate_var = "high"))
for (nm in names(sims)) {
  d <- simulate_toy_data(utils::modifyList(list(seed = 7), sims[[nm]])); b <- standardise_data(d, toy_mapping(d))
  res[[length(res) + 1]] <- check_one(paste("sim:", nm), b$data); cat(".", sep = "")
}
# edge cases: every individual seen once; two ages each; zeros and negative values; one individual
edge <- list("all single records" = data.frame(id = 1:50, age = sample(1:10, 50, TRUE), trait = stats::rnorm(50, 10)),
             "two ages each" = data.frame(id = rep(1:40, each = 2), age = rep(1:2, 40), trait = stats::rpois(80, 4)),
             "negative trait values" = data.frame(id = rep(1:30, each = 5), age = rep(1:5, 30), trait = stats::rnorm(150)),
             "one individual" = data.frame(id = 1, age = 1:6, trait = c(3, 5, 6, 6, 5, 3)),
             "all-zero counts" = data.frame(id = rep(1:20, each = 4), age = rep(1:4, 20), trait = 0))
for (nm in names(edge)) {
  e <- edge[[nm]]
  b <- standardise_data(e, list(id = "id", age = "age", trait = "trait", alr = "__AUTO_LAST__", life = "", entry = "__AUTO_FIRST__", covars = character(0),
                                cov_factor = character(0), cov_int = character(0), group = "", nested = TRUE, random = character(0), censor = "",
                                censor_value = "", condition = "", trials = "", start_mode = "afr", start_age = NA_real_, age_round = NA_real_, cov_age = character(0)))
  res[[length(res) + 1]] <- check_one(paste("edge:", nm), b$data); cat(".", sep = "")
}
cat("\n")
R <- do.call(rbind, res)
bad <- R[nzchar(R$Problems), , drop = FALSE]
options(width = 220)
cat(sprintf("%d dataset x function x scale runs on %d datasets; %d with problems; slowest %.1f s (%s)\n", nrow(R), length(unique(R$Data)), nrow(bad),
            max(R$Seconds, na.rm = TRUE), paste(R[which.max(R$Seconds), c("Data", "Function", "Scale")], collapse = ", ")))
if (nrow(bad)) print(bad[, c("Data", "Scale", "Function", "Problems")], row.names = FALSE)
cat(sprintf("Crosses for individuals with 2+ ages (every function failed): %d runs, %d individuals\n", sum(R$`Crosses with 2+ ages` > 0, na.rm = TRUE), sum(R$`Crosses with 2+ ages`, na.rm = TRUE)))
print(utils::head(R[order(-R$`Crosses with 2+ ages`), c("Data", "Scale", "Function", "Individuals", "Crosses with 2+ ages")], 8), row.names = FALSE)
out_csv <- Sys.getenv("DISAPPR_AUDIT_CSV", file.path(tempdir(), "all_trajectories_audit.csv")); utils::write.csv(R, out_csv, row.names = FALSE)
if (nrow(bad)) quit(status = 1)
