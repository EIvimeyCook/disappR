# Out-of-sample tests of the GAMM step (step 7): data the GAMM layer was not built or tuned on.
# Run from the package root:   Rscript tests/scripts/gamm_out_of_sample.R  [> inst/validation/gamm/out_of_sample_results.txt]
#   A. The app's own simulator (step 1's simulated data, with a known true trajectory): every trait type and age
#      function, selective disappearance (none, age-independent, age-dependent, both; both directions), selective
#      appearance, missingness and grouped data. Each dataset is fitted as a user would (ALR linear; AFR linear when
#      entry ages differ), and compared with the truth: how far the GAMM curve is from the true within-individual
#      trajectory (mean |D| %, as the app's 'deviation from truth' table), with and without the proxy term, against the
#      best step-5 GLMM (and, under age-dependent selection, the GAMM with the age-dependent form ti(age, ALR)); the sign of
#      the ALR term; and, for a subset, whether the proxy comparison names the simulated kind of selective disappearance.
#   B. The 13 bundled examples that were analysed with parametric models (steps 1-6), through every GAMM option:
#      proxy forms, AFR, known lifespan, individual random slopes or none, observation-level random effect, terminal
#      record, numeric covariates as smooths, group structures G-I on a categorical covariate, the gamm4 and bam
#      engines, and the group, proxy and shape comparisons; with the direction of the ALR effect compared with the
#      step-5 GLMM (Model 2) and the correlation between the GAMM and GLMM population curves.
# Every fit must succeed or return a plain message; every result goes through every post-fit function (tables, model
# output, guide, checks, curves, derivative, peak, residuals, cautions, R code). Fails (status 1) only on R errors or
# broken invariants (non-finite curves, intervals that do not contain the fit, R code that does not parse).
# DISAPPR_OOS_MAX = n runs only the first n of each part. DISAPPR_OOS_CORES sets the number of parallel workers (default 2); DISAPPR_OOS_PART = "A" or "B" runs one part.
root <- getwd()
app_dir <- file.path(root, "inst", "app")
if (!dir.exists(app_dir)) stop("Run from the package root.")
setwd(app_dir); on.exit(setwd(root), add = TRUE)
suppressMessages(source("global.R"))
if (!requireNamespace("mgcv", quietly = TRUE)) { cat("mgcv not installed: nothing to test.\n"); quit(status = 0) }
cores <- max(1L, as.integer(Sys.getenv("DISAPPR_OOS_CORES", "2")))
part <- toupper(Sys.getenv("DISAPPR_OOS_PART", "AB"))
n_max <- as.integer(Sys.getenv("DISAPPR_OOS_MAX", "0"))   # > 0: only the first n scenarios and examples (a quick check)
first_n <- function(x) if (n_max > 0) utils::head(x, n_max) else x
pmap <- function(x, f) if (cores > 1 && .Platform$OS.type == "unix") parallel::mclapply(x, f, mc.cores = cores, mc.preschedule = FALSE) else lapply(x, f)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- every post-fit function, with invariants ----------------------------------------------------------------------
post_fit <- function(label, r, gd) {
  err <- character(0); bad <- character(0)
  run <- function(nm, expr) { x <- tryCatch(expr, error = function(e) e); if (inherits(x, "error")) err <<- c(err, sprintf("%s | %s: %s", label, nm, conditionMessage(x))); x }
  kc <- run("gamm_kcheck", gamm_kcheck(r)); cc <- run("gamm_concurvity", gamm_concurvity(r))
  qr <- run("gamm_quantile_residuals", gamm_quantile_residuals(r))
  kc <- if (is.data.frame(kc)) kc else data.frame(); cc <- if (is.data.frame(cc)) cc else data.frame()
  o <- run("gamm_output", gamm_output(r))
  if (is.list(o) && !inherits(o, "error")) run("gamm_output_reading", gamm_output_reading(o))
  run("gamm_re_test", gamm_re_test(r))
  if (is.list(o) && !inherits(o, "error")) run("gamm_guide", gamm_guide(r, o))
  run("gamm_check_summary", gamm_check_summary(r, kc, cc, if (inherits(qr, "error")) NULL else qr))
  for (fn in c("gamm_tables", "gamm_difference", "gamm_vcomp", "gamm_residuals", "gamm_support")) run(fn, get(fn)(r))
  run("gamm_cautions", gamm_cautions(r, kc, cc, gd = gd))
  for (sc in c("typical", "average")) {
    cv <- run(paste("gamm_curve", sc), gamm_curve(r, scale = sc))
    if (is.data.frame(cv)) {
      if (!nrow(cv) || any(!is.finite(cv$fit))) bad <- c(bad, sprintf("%s | %s curve has non-finite values", label, sc))
      ok <- is.finite(cv$lower_sim) & is.finite(cv$upper_sim)
      if (any(cv$lower_sim[ok] > cv$fit[ok] + 1e-8 | cv$upper_sim[ok] < cv$fit[ok] - 1e-8)) bad <- c(bad, sprintf("%s | %s curve outside its interval", label, sc))
    }
  }
  dv <- run("gamm_derivative", gamm_derivative(r))
  if (is.data.frame(dv)) {
    run("gamm_periods", gamm_periods(dv, simultaneous = TRUE))
    po <- run("gamm_peak_onset", gamm_peak_onset(r, n_draws = 200))
    run("gamm_reading", gamm_reading(r, dv, if (inherits(po, "error")) NULL else po))
  }
  code <- tryCatch(paste(gamm_code(r), collapse = "\n"), error = function(e) "")
  if (!nzchar(code) || inherits(tryCatch(parse(text = code), error = function(e) e), "error")) bad <- c(bad, paste(label, "| R code does not parse"))
  list(errors = err, broken = bad)
}
fit_one <- function(label, gd, s) {
  t0 <- Sys.time()
  r <- tryCatch(gamm_fit(gd, s), error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) return(list(label = label, status = "ERROR", errors = paste(label, "| gamm_fit:", conditionMessage(r)), secs = secs))
  if (!isTRUE(r$ok)) return(list(label = label, status = "message", message = r$message, secs = secs,
                                 errors = if (!is.character(r$message) || !nzchar(r$message)) paste(label, "| failed without a message") else character(0)))
  pf <- post_fit(label, r, gd)
  list(label = label, status = "ok", r = r, secs = secs, errors = pf$errors, broken = pf$broken)
}
alr_term <- function(r) {
  pt <- tryCatch(gamm_tables(r)$parametric, error = function(e) NULL)
  if (is.null(pt)) return(c(est = NA, p = NA))
  i <- which(pt$Term == "ALR")
  if (!length(i)) return(c(est = NA, p = NA))
  c(est = pt$Estimate[i[1]], p = pt$P[i[1]])
}
all_errors <- character(0); all_broken <- character(0)
line <- function(...) cat(sprintf(...), "\n", sep = "")

# =====================================================================================================================
# A. Simulated data with a known truth
# =====================================================================================================================
if (grepl("A", part)) {
line("==== A. The app's simulator: GAMM against the true trajectory ====")
sc <- list()
add <- function(...) sc[[length(sc) + 1]] <<- list(...)
for (f in AGE_FUNCTIONS) for (sd in list(c("none", 1), c("independent", 1), c("independent", -1), c("dependent", 1), c("both", 1)))
  add(trait = "mass", form = f, sd_type = sd[1], sd_dir = as.numeric(sd[2]), cmp = f %in% c("Quadratic", "Logarithmic"))
for (f in AGE_FUNCTIONS) for (sd in list(c("none", 1), c("independent", 1), c("both", -1)))
  add(trait = "count_pois", form = f, sd_type = sd[1], sd_dir = as.numeric(sd[2]), cmp = identical(f, "Quadratic"))
for (sd in c("none", "independent", "both")) add(trait = "count_nb", form = "Quadratic", sd_type = sd, sd_dir = 1, cmp = TRUE)
for (sd in c("none", "both")) add(trait = "count_zinb", form = "Quadratic", sd_type = sd, sd_dir = 1, cmp = FALSE)
for (sd in c("none", "both")) add(trait = "paper", form = "Quadratic", sd_type = sd, sd_dir = 1, cmp = FALSE)
for (sa in c("independent", "dependent")) add(trait = "mass", form = "Quadratic", sd_type = "independent", sd_dir = 1, afr_mode = "individual", sa_type = sa, cmp = TRUE)
add(trait = "count_pois", form = "Quadratic", sd_type = "both", sd_dir = 1, afr_mode = "individual", sa_type = "both", cmp = FALSE)
for (m in c("mcar", "mwo", "mwy", "trait")) add(trait = "mass", form = "Quadratic", sd_type = "independent", sd_dir = 1, missingness = m, cmp = FALSE)
add(trait = "mass", form = "Quadratic", sd_type = "both", sd_dir = 1, groups = TRUE, diet = TRUE, cmp = FALSE)
add(trait = "count_pois", form = "Cubic", sd_type = "independent", sd_dir = -1, groups = TRUE, diet = TRUE, strength = "subtle", cmp = FALSE)

sim_one <- function(i) {
  cfg <- sc[[i]]; cmp <- isTRUE(cfg$cmp); cfg$cmp <- NULL
  cfg <- utils::modifyList(list(n_id = 200, seed = 100 + i), cfg)
  label <- sprintf("%-10s %-22s SD %-11s %+d%s%s%s", cfg$trait, cfg$form, cfg$sd_type, cfg$sd_dir, if (identical(cfg$afr_mode, "individual")) paste(" SA", cfg$sa_type) else "",
                   if (!is.null(cfg$missingness)) paste(" miss", cfg$missingness) else "", if (isTRUE(cfg$groups)) " groups+diet" else "")
  out <- list(label = label, errors = character(0), broken = character(0))
  d <- tryCatch(simulate_toy_data(cfg), error = function(e) e)
  if (inherits(d, "error")) { out$errors <- paste(label, "| simulate:", conditionMessage(d)); return(out) }
  tr <- attr(d, "truth")
  b <- tryCatch(standardise_data(d, toy_mapping(d)), error = function(e) e)
  if (inherits(b, "error")) { out$errors <- paste(label, "| standardise:", conditionMessage(b)); return(out) }
  gd <- gamm_data(b$data, b$meta)
  fam <- toy_family(cfg$trait)
  afr <- identical(cfg$afr_mode, "individual")
  s <- utils::modifyList(gamm_defaults(), list(family = gamm_suggest_family(fam), proxy = "ALR", proxy_form = "linear", afr_form = if (afr) "linear" else "none"))
  a <- fit_one(paste(label, "ALR"), gd, s)
  n <- fit_one(paste(label, "no proxy"), gd, utils::modifyList(s, list(proxy = "none", afr_form = "none")))
  out$errors <- c(a$errors, n$errors); out$broken <- c(a$broken, n$broken)
  if (!identical(a$status, "ok") || !identical(n$status, "ok")) { out$status <- paste(a$status, n$status, a$message %||% "", n$message %||% ""); return(out) }
  # age-dependent selective disappearance: also the age-dependent form the proxy comparison points to, ti(age, ALR)
  t <- if (cfg$sd_type %in% c("dependent", "both")) fit_one(paste(label, "ALR ti"), gd, utils::modifyList(s, list(proxy_form = "ti"))) else NULL
  if (!is.null(t)) { out$errors <- c(out$errors, t$errors); out$broken <- c(out$broken, t$broken) }
  obs <- stats::aggregate(id ~ age, data = b$data[is.finite(b$data$trait), ], FUN = function(z) length(unique(z)))
  names(obs) <- c("age", "n")
  curves <- rbind(data.frame(age = (ca <- gamm_curve(a$r, by_proxy = FALSE))$age, fitted = ca$fit, Method = "GAMM with proxy"),
                  data.frame(age = (cn <- gamm_curve(n$r, by_proxy = FALSE))$age, fitted = cn$fit, Method = "GAMM without proxy"))
  if (identical(t$status, "ok")) curves <- rbind(curves, data.frame(age = (ct <- gamm_curve(t$r, by_proxy = FALSE))$age, fitted = ct$fit, Method = "GAMM with ti(age, ALR)"))
  # the step-5 benchmark: the best of Models 1, 2 and 4 (with Model 7 when entry ages differ), with the true age function
  gl <- tryCatch(fit_model_suite(b$data, b$meta, c("M1", "M2", "M4", if (afr) "M7"), cfg$form, fam), error = function(e) e)
  if (inherits(gl, "error") || !isTRUE(gl$ok)) out$glmm <- "failed" else {
    best <- sub("^Model ", "M", gl$aic$Model[[1]])   # the AIC table holds labels ("Model 2"), the fits ids ("M2")
    pc <- tryCatch(predict_population_curve(gl$fits[[best]], gl, sort(unique(c(ca$age)))), error = function(e) NULL)
    if (is.data.frame(pc) && nrow(pc)) curves <- rbind(curves, data.frame(age = pc$age, fitted = pc$fitted, Method = paste0("GLMM (", best, ")")))
    out$glmm <- best
  }
  dev <- deviation_from_truth(curves, tr, obs)
  out$dev <- stats::setNames(round(dev$Mean_abs_D, 1), dev$Method)
  at <- alr_term(a$r)
  out$alr <- at
  out$sd_level <- cfg$sd_type %in% c("independent", "both"); out$sd_dir <- cfg$sd_dir; out$sd_type <- cfg$sd_type; out$trait <- cfg$trait
  out$engine <- a$r$engine; out$n <- a$r$n; out$secs <- round(a$secs + n$secs)
  if (cmp) {
    pc <- tryCatch(gamm_proxy_comparison(gd, s), error = function(e) e)
    if (inherits(pc, "error")) out$errors <- c(out$errors, paste(label, "| gamm_proxy_comparison:", conditionMessage(pc)))
    else if (isTRUE(pc$ok)) { out$verdict <- paste(names(pc$verdicts), pc$verdicts, sep = ": ", collapse = "; "); out$alr_verdict <- pc$verdicts[["ALR"]] }
    nv <- tryCatch(gamm_naive_vs_corrected(gd, a$r), error = function(e) e)
    if (inherits(nv, "error")) out$errors <- c(out$errors, paste(label, "| gamm_naive_vs_corrected:", conditionMessage(nv)))
  }
  out$status <- "ok"
  out
}
t0 <- Sys.time()
RA <- pmap(first_n(seq_along(sc)), sim_one)
for (x in RA) {
  if (inherits(x, "try-error")) { all_errors <- c(all_errors, paste("worker:", as.character(x))); next }
  all_errors <- c(all_errors, x$errors); all_broken <- c(all_broken, x$broken)
  if (!identical(x$status, "ok")) { line("  %s  NOT FITTED: %s", x$label, x$status %||% ""); next }
  dv <- x$dev
  line("  %s | %s %4d rec %3ds | mean |D|%%: with proxy %5s, without %5s%s, %s %5s | ALR %+.3g (p %s)%s",
       x$label, x$engine, x$n, x$secs, dv[["GAMM with proxy"]] %||% NA, dv[["GAMM without proxy"]] %||% NA,
       if (!is.null(dv[["GAMM with ti(age, ALR)"]])) sprintf(", ti %5s", dv[["GAMM with ti(age, ALR)"]]) else "", sprintf("GLMM %-3s", x$glmm),
       if (length(gi <- grep("^GLMM", names(dv)))) dv[[gi[1]]] else NA, x$alr[["est"]], format_p(x$alr[["p"]]),
       if (!is.null(x$verdict)) paste0(" | proxy comparison: ", x$verdict) else "")
}
ok <- Filter(function(x) identical(x$status, "ok") && length(x$dev), RA)
sel <- Filter(function(x) isTRUE(x$sd_level), ok)
sgn <- Filter(function(x) identical(x$sd_type, "independent") && x$trait %in% c("mass", "count_pois", "count_nb"), ok)
better <- vapply(sel, function(x) isTRUE(x$dev[["GAMM with proxy"]] <= x$dev[["GAMM without proxy"]]), logical(1))
right <- vapply(sgn, function(x) isTRUE(sign(x$alr[["est"]]) == x$sd_dir && x$alr[["p"]] < 0.05), logical(1))
gam_v_glmm <- vapply(ok, function(x) { g <- grep("^GLMM", names(x$dev)); if (!length(g)) NA_real_ else x$dev[["GAMM with proxy"]] - x$dev[[g[1]]] }, numeric(1))
line("")
line("  %d of %d simulated datasets fitted (%.0f min).", length(ok), length(RA), as.numeric(difftime(Sys.time(), t0, units = "mins")))
line("  With selective disappearance through lifespan (%d datasets): the ALR term brings the GAMM closer to the true curve in %d.", length(sel), sum(better))
line("  Age-independent selective disappearance (%d datasets): ALR term in the simulated direction and p < 0.05 in %d.", length(sgn), sum(right))
line("  GAMM with proxy minus best GLMM, mean |D| percentage points: median %.2f, range %.2f to %.2f (negative: the GAMM is closer).",
     stats::median(gam_v_glmm, na.rm = TRUE), min(gam_v_glmm, na.rm = TRUE), max(gam_v_glmm, na.rm = TRUE))
dep <- Filter(function(x) x$sd_type %in% c("dependent", "both") && !is.null(x$dev[["GAMM with ti(age, ALR)"]]), ok)
if (length(dep)) {
  dl <- vapply(dep, function(x) x$dev[["GAMM with proxy"]], numeric(1)); dt <- vapply(dep, function(x) x$dev[["GAMM with ti(age, ALR)"]], numeric(1))
  dg <- vapply(dep, function(x) { g <- grep("^GLMM", names(x$dev)); if (length(g)) x$dev[[g[1]]] else NA_real_ }, numeric(1))
  line("  Age-dependent selective disappearance (%d datasets), median mean |D|%%: GAMM with linear ALR %.1f, with ti(age, ALR) %.1f, best GLMM %.1f.",
       length(dep), stats::median(dl), stats::median(dt), stats::median(dg, na.rm = TRUE))
}
vv <- Filter(function(x) !is.null(x$alr_verdict), ok)
if (length(vv)) {
  want <- function(x) switch(x$sd_type, none = "none", independent = "age_independent", "age_dependent")
  hit <- vapply(vv, function(x) identical(x$alr_verdict, want(x)), logical(1))
  line("  Proxy comparison (marginal-likelihood AIC) gave the simulated kind of selective disappearance in %d of %d datasets%s.", sum(hit), length(vv),
       if (any(!hit)) paste0(" (other: ", paste(vapply(vv[!hit], function(x) sprintf("%s %s -> %s", x$sd_type, x$trait, x$alr_verdict), character(1)), collapse = "; "), ")") else "")
}
line("")
}

# =====================================================================================================================
# B. Bundled examples analysed with parametric models
# =====================================================================================================================
if (grepl("B", part)) {
line("==== B. Parametric-study examples through every GAMM option ====")
keys <- c("fly", "bichet", "bichet_marmot", "moullec_swift", "pasztor_apollo", "wynn", "sanghvi_female", "allain", "bouwhuis",
          "warner", "mckennaell_breeding", "mckennaell_weight", "szejnersigal_activity")
ex_one <- function(key) {
  ex <- EXAMPLES[[key]]
  out <- list(key = key, errors = character(0), broken = character(0), rows = character(0))
  raw <- load_example_file(ex$file)
  if (!is.null(ex$subset)) raw <- raw[!is.na(raw[[ex$subset$var]]) & as.character(raw[[ex$subset$var]]) %in% ex$subset$levels, , drop = FALSE]
  b <- tryCatch(standardise_data(raw, ex$mapping), error = function(e) e)
  if (inherits(b, "error")) { out$errors <- paste(key, "| standardise:", conditionMessage(b)); return(out) }
  gd <- gamm_data(b$data, b$meta)
  fam <- ex$family %||% "gaussian"
  base <- utils::modifyList(gamm_defaults(), list(family = gamm_suggest_family(fam)))
  grp <- gamm_group_choices(b$meta, gd)
  has_ls <- any(is.finite(gd$LS)); has_afr <- any(is.finite(gd$AFR)) && length(unique(gd$AFR[is.finite(gd$AFR)])) >= 3
  num_cov <- any(vapply(b$meta$covars %||% character(0), function(v) identical(b$meta$cov_types[[v]], "numeric"), logical(1)))
  variants <- list("default (ALR linear)" = list(),
                   "ALR smooth" = list(proxy_form = "smooth"), "ALR vc s(age, by = ALR)" = list(proxy_form = "vc"),
                   "ALR ti(age, ALR)" = list(proxy_form = "ti"), "no proxy" = list(proxy = "none"),
                   "random slopes" = list(individual = "rs"), "no individual effect" = list(individual = "none"),
                   "terminal record" = list(terminal = TRUE), "observation-level RE" = list(olre = TRUE),
                   "k = 5" = list(k = 5), "engine bam" = list(engine = "bam"), "engine gamm4" = list(engine = "gamm4"),
                   "no covariates or random intercepts" = list(covars = FALSE, random = FALSE))
  if (has_ls) variants[["known lifespan (LS) smooth"]] <- list(proxy = "LS", proxy_form = "smooth")
  if (has_afr) variants[["AFR linear + ALR linear"]] <- list(afr_form = "linear")
  if (num_cov) variants[["numeric covariates as smooths"]] <- list(cov_smooth = TRUE)
  if (length(grp)) for (st in c("GS", "GI", "S", "I")) variants[[sprintf("group %s, %s", names(grp)[1], st)]] <- list(group = grp[[1]], structure = st)
  res <- list()
  for (v in names(variants)) {
    f <- fit_one(paste(key, "|", v), gd, utils::modifyList(base, variants[[v]]))
    out$errors <- c(out$errors, f$errors); out$broken <- c(out$broken, f$broken)
    out$rows <- c(out$rows, sprintf("    %-40s %-7s %6.1fs %s", v, f$status, f$secs,
                                    if (identical(f$status, "ok")) sprintf("%s, %d rec, AIC %.1f", f$r$engine, f$r$n, gamm_aic(f$r)) else substr(f$message %||% "", 1, 110)))
    if (identical(v, "default (ALR linear)") && identical(f$status, "ok")) res$default <- f$r
  }
  # the comparisons the app offers
  cmp <- list(shape = function() gamm_parametric_reference(gd, base), proxy = function() gamm_proxy_comparison(gd, base))
  if (length(grp)) cmp$groups <- function() gamm_group_comparison(gd, utils::modifyList(base, list(group = grp[[1]])))
  for (nm in names(cmp)) {
    t1 <- Sys.time(); x <- tryCatch(cmp[[nm]](), error = function(e) e); secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(x, "error")) { out$errors <- c(out$errors, sprintf("%s | %s comparison: %s", key, nm, conditionMessage(x))); next }
    txt <- if (is.data.frame(x)) paste("best:", x$Model[[1]]) else if (isTRUE(x$ok)) (if (!is.null(x[["verdict"]])) x[["verdict"]] else paste(x$texts, collapse = " ")) else paste("message:", x$message)
    out$rows <- c(out$rows, sprintf("    %-40s %6.1fs %s", paste(nm, "comparison"), secs, substr(txt, 1, 160)))
  }
  if (!is.null(res$default)) {
    nv <- tryCatch(gamm_naive_vs_corrected(gd, res$default), error = function(e) e)
    if (inherits(nv, "error")) out$errors <- c(out$errors, paste(key, "| naive vs corrected:", conditionMessage(nv)))
  }
  # agreement with step 5: Model 2 (age function + ALR) against the default GAMM (smooth of age + ALR)
  gl <- tryCatch(fit_model_suite(b$data, b$meta, "M2", ex$age_function %||% "Quadratic", fam), error = function(e) e)
  if (!inherits(gl, "error") && isTRUE(gl$ok) && !is.null(res$default)) {
    cf <- tryCatch(glmmTMB::fixef(gl$fits$M2)$cond, error = function(e) tryCatch(lme4::fixef(gl$fits$M2), error = function(e) NULL))
    g_alr <- if (!is.null(cf) && "ALR" %in% names(cf)) cf[["ALR"]] else NA_real_
    a_alr <- alr_term(res$default)
    cg <- gamm_curve(res$default, by_proxy = FALSE)
    pc <- tryCatch(predict_population_curve(gl$fits$M2, gl, cg$age), error = function(e) NULL)
    rho <- if (is.data.frame(pc) && nrow(pc) > 3) stats::cor(cg$fit[match(pc$age, cg$age)], pc$fitted, use = "complete.obs") else NA_real_
    out$agree <- c(glmm = g_alr, gamm = a_alr[["est"]], p = a_alr[["p"]], rho = rho)
  } else out$agree <- c(glmm = NA, gamm = NA, p = NA, rho = NA)
  out
}
t0 <- Sys.time()
RB <- pmap(first_n(keys), ex_one)
agree <- list()
for (x in RB) {
  if (inherits(x, "try-error")) { all_errors <- c(all_errors, paste("worker:", as.character(x))); next }
  all_errors <- c(all_errors, x$errors); all_broken <- c(all_broken, x$broken)
  line("  %s", x$key); for (rw in x$rows) line("%s", rw)
  ag <- x$agree
  line("    agreement with step 5 (Model 2): ALR coefficient GLMM %+.3g, GAMM %+.3g (p %s), %s; population curves correlate r = %.3f",
       ag[["glmm"]], ag[["gamm"]], format_p(ag[["p"]]), if (isTRUE(sign(ag[["glmm"]]) == sign(ag[["gamm"]]))) "same direction" else "DIFFERENT DIRECTION", ag[["rho"]])
  agree[[x$key]] <- ag
}
A <- do.call(rbind, agree)
same <- sign(A[, "glmm"]) == sign(A[, "gamm"])
sig_diff <- !same & A[, "p"] < 0.05
line("")
line("  %d examples (%.0f min). ALR effect in the same direction as the GLMM in %d of %d; opposite and significant in the GAMM in %d.",
     length(RB), as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(same, na.rm = TRUE), sum(is.finite(A[, "glmm"]) & is.finite(A[, "gamm"])), sum(sig_diff, na.rm = TRUE))
line("  GAMM vs GLMM population curves: median correlation %.3f (lowest %.3f, %s).", stats::median(A[, "rho"], na.rm = TRUE),
     min(A[, "rho"], na.rm = TRUE), rownames(A)[which.min(A[, "rho"])])
}

if (length(all_broken)) { cat("\nBROKEN INVARIANTS (", length(all_broken), "):\n  ", paste(all_broken, collapse = "\n  "), "\n", sep = "") }
if (length(all_errors)) { cat("\nERRORS (", length(all_errors), "):\n  ", paste(all_errors, collapse = "\n  "), "\n", sep = ""); }
if (length(all_errors) || length(all_broken)) quit(status = 1)
cat("\nOut-of-sample GAMM tests finished: no errors, no broken invariants.\n")
