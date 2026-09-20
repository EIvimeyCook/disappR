# ============================================================================
# disappR 0.9.4 - global definitions
# Workflow follows Figure 8 of the associated paper:
#   A1 binned age-trait trajectories, A2 proxy-trait association at each age,
#   A3 individual parametric fits; B1 random effects, B2 ageing function,
#   B3 proxy choice (ALR vs mean age), B4-B5 comparative Models 1-8.
# ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

HAS_GLMMTMB  <- requireNamespace("glmmTMB", quietly = TRUE)
HAS_LMERTEST <- requireNamespace("lmerTest", quietly = TRUE)
HAS_DHARMA   <- requireNamespace("DHARMa", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

AGE_FUNCTIONS <- c("Linear", "Quadratic", "Cubic", "Logarithmic", "Asymptotic exponential")
# A3 only: a function that is non-linear in its parameters, a * exp(b * age), fitted by nls per individual
A3_NONLINEAR <- "Exponential (a \u00b7 exp(b \u00b7 age))"
A3_FUNCTIONS <- c(AGE_FUNCTIONS, A3_NONLINEAR)
# Colour-blind-safe, clearly distinct colours for ageing functions (Okabe-Ito)
FUNCTION_COLOURS <- stats::setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"), A3_FUNCTIONS)
MODEL_IDS <- paste0("M", 1:10)
options(shiny.maxRequestSize = 200 * 1024^2)   # uploads up to 200 MB
model_label <- function(m) sub("^M", "Model ", m)

# Colours follow the manuscript figures (Models 1-6) where possible.
METHOD_COLOURS <- c(
  "Model 1" = "darkgrey", "Model 2" = "hotpink1", "Model 3" = "darkgreen",
  "Model 4" = "orange", "Model 5" = "steelblue2", "Model 6" = "red",
  "Model 7" = "#8C510A", "Model 8" = "#01665E", "Model 9" = "#762A83", "Model 10" = "#1B7837",
  "Observed" = "grey35", "Decomposition" = "purple", "True (simulated)" = "black",
  "Individual fits: mean of coefficients" = "#E7298A", "Individual fits: mean of functions" = "#66A61E"
)

# Clearly distinguishable colours: Okabe-Ito for up to 7 levels, then an HCL qualitative palette.
distinct_colours <- function(n) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}
# Ordered but clearly distinguishable colours (for ordered bins): the viridis HCL palette without its palest end.
ordered_colours <- function(n) {
  n <- max(1L, as.integer(n))
  grDevices::hcl.colors(n + 1, "Viridis")[seq_len(n)]
}
warm_palette <- c("#A85F4A", "#C8944C", "#7C8060", "#B77C75", "#6D5A4A", "#D4AA7D")
bin_palette <- c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2CBB7", "#EE8468", "#B40426")

theme_disappR <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#FBF7F1", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#FBF7F1", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E7DDD0", linewidth = 0.35),
      axis.title = ggplot2::element_text(colour = "#4B4037", face = "bold"),
      axis.text = ggplot2::element_text(colour = "#62564C"),
      plot.title = ggplot2::element_text(colour = "#453A32", face = "bold", size = ggplot2::rel(1.05)),
      plot.subtitle = ggplot2::element_text(colour = "#75665A", size = ggplot2::rel(0.85)),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold", colour = "#4B4037"),
      legend.text = ggplot2::element_text(colour = "#62564C")
    )
}

empty_plot <- function(msg) {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = msg, colour = "#75665A", size = 5) +
    ggplot2::theme_void()
}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

mode_or_na <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

finite_mean <- function(v) {
  v <- v[is.finite(v)]
  if (length(v)) mean(v) else NA_real_
}

cor_safe <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 4) return(NA_real_)
  if (!isTRUE(stats::sd(a[ok]) > 0) || !isTRUE(stats::sd(b[ok]) > 0)) return(NA_real_)
  stats::cor(a[ok], b[ok])
}

format_p <- function(x) {
  if (length(x) != 1 || !is.finite(x)) return("NA")
  if (x < 0.001) "<0.001" else sprintf("%.3f", x)
}

format_num <- function(x, digits = 3) format(signif(x, digits), trim = TRUE, scientific = FALSE)

# Least-squares line of y on x with its 95% confidence band on a grid (trend lines drawn in figures).
# A line through only two distinct x values has no residual degrees of freedom and therefore no band (NA).
lm_band <- function(x, y, n = 50) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(unique(x)) < 2) return(data.frame())
  fit <- stats::lm(y ~ x)
  grid <- seq(min(x), max(x), length.out = n)
  pr <- suppressWarnings(stats::predict(fit, newdata = data.frame(x = grid), interval = "confidence", level = 0.95))
  out <- data.frame(x = grid, fit = as.numeric(pr[, "fit"]), lo = as.numeric(pr[, "lwr"]), hi = as.numeric(pr[, "upr"]))
  out$lo[!is.finite(out$lo)] <- NA_real_
  out$hi[!is.finite(out$hi)] <- NA_real_
  out
}

# Sampling step: the smallest interval that is common (>= 5% of intervals) among within-individual
# intervals between consecutive records. Staggered cohorts, floating-point noise in ages and occasional
# extra records therefore do not shrink it. Without repeated records it falls back to gaps between
# distinct ages; when no interval is common (irregular sampling) it is the median interval. A smallest
# common interval that is close to (>= 0.75 times) but not a divisor of a clearly dominant interval (>= 60%)
# is an unequal first or last interval, e.g. a day-1 record followed by weekly sampling on days 7, 14, 21,
# not a finer schedule: the dominant interval is then the step.
infer_age_step <- function(age, id = NULL) {
  ok <- is.finite(age)
  if (!is.null(id)) ok <- ok & !is.na(id)
  if (sum(ok) < 2) return(NA_real_)
  a <- age[ok]
  tol <- 1e-6 * max(1, diff(range(a)))
  d <- numeric(0)
  if (!is.null(id)) {
    g <- as.character(id[ok])
    o <- order(g, a)
    a2 <- a[o]
    g2 <- g[o]
    same <- g2[-1] == g2[-length(g2)]
    d <- diff(a2)[same]
  }
  d <- d[is.finite(d) & d > tol]
  if (!length(d)) {
    u <- sort(unique(a))
    d <- diff(u)
    d <- d[d > tol]
    if (!length(d)) return(NA_real_)
  }
  d <- signif(d, 6)
  tab <- table(d)
  vals <- as.numeric(names(tab))
  share <- as.numeric(tab) / length(d)
  common <- vals[share >= 0.05]
  if (!length(common)) return(stats::median(d))
  step <- min(common)
  dom <- vals[which.max(share)]
  if (max(share) >= 0.6 && step / dom >= 0.75 && abs(dom / step - round(dom / step)) > 0.01) step <- dom
  step
}

# ---------------------------------------------------------------------------
# Teaching simulator
# ---------------------------------------------------------------------------
# Individual ageing trajectories F_i(age) = sum_k b_ki g_k(u), with u = age / (mean lifespan / 20) so
# that the trajectory shape is the same for 3 or 30 sampling occasions, and g_k the basis of the
# chosen form (linear, quadratic, cubic, logarithmic, exponential; the same bases as the models).
# Each individual's coefficients = population means
#   + an individual component (heterogeneity in level and ageing rate),
#   + a lifespan-linked component (selective disappearance): on the intercept for age-independent,
#     on the first ageing coefficient ("rate") for age-dependent selective disappearance,
#   + an AFR-linked component (selective appearance), defined in the same way.
# The lifespan and AFR latent variables are independent, so AFR and ALR are uncorrelated by
# construction; with individual-specific AFR, individuals are included only if they outlive the
# latest possible AFR, which depends on lifespan alone. Count traits use a log link.
# "paper" reproduces the manuscript's multivariate-normal quadratic parameterisation exactly.
TOY_TRAITS <- c(
  "Body mass \u2014 continuous (Gaussian)" = "mass",
  "Fecundity \u2014 counts (Poisson)" = "count_pois",
  "Fecundity \u2014 overdispersed counts (negative binomial)" = "count_nb",
  "Fecundity \u2014 zero-inflated overdispersed counts (ZINB)" = "count_zinb",
  "Paper simulation \u2014 continuous fecundity (Methods Eq. 1\u20134)" = "paper"
)

SELECTION_TYPES <- c("None" = "none", "Age-independent" = "independent", "Age-dependent" = "dependent",
                     "Age-independent and age-dependent" = "both")

TOY_DEFAULTS <- list(trait = "mass", form = "Quadratic", shape = "default", strength = "dramatic", sd_type = "both", sd_dir = 1,
                     rate_var = "low", mean_ls = 20, missingness = "complete",
                     afr_mode = "same", sa_type = "none", sa_dir = 1, diet = FALSE, groups = FALSE,
                     n_id = 300, seed = 1)

# Parameters of each simulated form on the scaled time u. mu: population coefficients (intercept
# first); s_rate: lifespan- or AFR-linked SD of the first ageing coefficient at dramatic strength;
# u_small / u_large: individual SDs of the ageing coefficients (the intercept SD is fixed).
# The small values were tuned so that, without selective disappearance, random-intercept models
# usually do not favour the interaction models (verified by simulation for every form).
toy_form_spec <- function(form, count = FALSE) {
  if (!isTRUE(count)) {
    switch(form,
      Linear      = list(mu = c(50, -0.8), s_rate = 0.425, u_small = 0.04, u_large = 0.5),
      Cubic       = list(mu = c(40, 1.2, 0.06, -0.0035), s_rate = 0.425,
                         u_small = c(0.04, 0.0012, 0.00003), u_large = c(0.5, 0.015, 0.0005)),
      Logarithmic = list(mu = c(30, 8), s_rate = 2.8, u_small = 0.27, u_large = 3.3),
      "Asymptotic exponential" = list(mu = c(28, 6), s_rate = 2.0, u_small = 0.17, u_large = 2.3),
      list(mu = c(40, 2, -0.06), s_rate = 0.425, u_small = c(0.04, 0.0012), u_large = c(0.5, 0.015)))
  } else {
    switch(form,
      Linear      = list(mu = c(log(40), -0.06), s_rate = 0.0255, u_small = 0.003, u_large = 0.03),
      Cubic       = list(mu = c(log(20), 0.08, 0.004, -0.0003), s_rate = 0.0255,
                         u_small = c(0.003, 0.00015, 0.000005), u_large = c(0.03, 0.0015, 0.00005)),
      Logarithmic = list(mu = c(log(10), 0.5), s_rate = 0.17, u_small = 0.02, u_large = 0.2),
      "Asymptotic exponential" = list(mu = c(log(12), 0.35), s_rate = 0.12, u_small = 0.012, u_large = 0.14),
      list(mu = c(log(30), 0.12, -0.007), s_rate = 0.0255, u_small = c(0.003, 0.00015), u_large = c(0.03, 0.0015)))
  }
}

# Basis of a simulated form (without intercept). The exponential basis is exp(-(u - m) / s) (a decline that is
# steep early in life and flattens later, with a positive coefficient), with m and s
# the mean and SD of all expected (pre-missingness) scaled ages, mirroring the app's asymptotic exponential function.
toy_form_basis <- function(form, u, ref = c(0, 1)) {
  switch(form,
    Linear = cbind(u),
    Cubic = cbind(u, u^2, u^3),
    Logarithmic = cbind(log(pmax(u, 1e-9))),
    "Asymptotic exponential" = cbind(exp(-(u - ref[1]) / ref[2])),
    cbind(u, u^2))
}

# Biologically motivated shapes within each simulated ageing form: the population mean coefficients (intercept
# first) on the scaled time u (u = 20 at the mean lifespan), for continuous traits (mu) and on the log scale for
# counts (mu_count). "default" is the form's original shape (toy_form_spec()), so default simulations are
# unchanged. Shapes change the mean coefficients only: the lifespan- and AFR-linked components (selection) and the
# individual variation stay the same, so selective disappearance and appearance keep their strength. Continuous
# shapes stay positive to about 1.75 times the mean lifespan.
TOY_SHAPES <- list(
  Linear = list(
    default = list(label = "Moderate senescence (steady decline)"),
    slow = list(label = "Slow senescence", mu = c(50, -0.4), mu_count = c(log(40), -0.03)),
    fast = list(label = "Fast senescence", mu = c(50, -1.25), mu_count = c(log(40), -0.1)),
    none = list(label = "No senescence (flat)", mu = c(50, 0), mu_count = c(log(40), 0)),
    improve = list(label = "Improvement with age (e.g. growth or experience)", mu = c(50, 0.5), mu_count = c(log(40), 0.03))),
  Quadratic = list(
    default = list(label = "Improvement, then senescence"),
    early = list(label = "Early peak, then fast senescence", mu = c(40, 1.3, -0.065), mu_count = c(log(30), 0.1, -0.008)),
    late = list(label = "Late peak, then slow senescence", mu = c(40, 1.6, -0.04), mu_count = c(log(30), 0.112, -0.004)),
    ushape = list(label = "Decline, then improvement (U-shape)", mu = c(50, -1.4, 0.04), mu_count = c(log(40), -0.1, 0.0028))),
  Cubic = list(
    default = list(label = "Improvement, plateau, then accelerating late-life decline"),
    twophase = list(label = "Early decline, mid-life plateau, then late-life decline",
                    mu = c(45, -0.864, 0.072, -0.002), mu_count = c(log(40), -0.0864, 0.0072, -0.0002)),
    terminal = list(label = "Stable, then late-life increase (e.g. terminal investment)",
                    mu = c(40, 0.45, -0.045, 0.0015), mu_count = c(log(20), 0.03, -0.003, 0.0001))),
  Logarithmic = list(
    default = list(label = "Rapid early improvement, then levelling off"),
    gradual_up = list(label = "Gradual improvement, then levelling off", mu = c(40, 4), mu_count = c(log(20), 0.25)),
    rapid_down = list(label = "Rapid early decline, then slower decline", mu = c(55, -8), mu_count = c(log(40), -0.5)),
    gradual_down = list(label = "Gradual decline, then levelling off", mu = c(50, -4), mu_count = c(log(30), -0.25))),
  "Asymptotic exponential" = list(
    default = list(label = "Rapid early decline, then plateau"),
    gradual_down = list(label = "Gradual decline, then plateau", mu = c(34, 3), mu_count = c(log(14), 0.18)),
    rapid_up = list(label = "Rapid early improvement, then plateau (growth to an asymptote)", mu = c(55, -6), mu_count = c(log(30), -0.35)),
    gradual_up = list(label = "Gradual improvement, then plateau", mu = c(47, -3), mu_count = c(log(24), -0.18)))
)
toy_shape_choices <- function(form) {
  sh <- TOY_SHAPES[[form %||% "Quadratic"]] %||% TOY_SHAPES[["Quadratic"]]
  stats::setNames(names(sh), vapply(names(sh), function(k) paste0(sh[[k]]$label, if (identical(k, "default")) " (default)" else ""), character(1)))
}
# Mean coefficients of a non-default shape; NULL for the default (or an unknown) shape.
toy_shape_mu <- function(form, shape, count = FALSE) {
  shape <- shape %||% "default"
  sh <- TOY_SHAPES[[form %||% ""]][[shape]]
  if (is.null(sh) || identical(shape, "default")) return(NULL)
  if (isTRUE(count)) sh$mu_count else sh$mu
}
toy_shape_label <- function(form, shape) {
  sh <- TOY_SHAPES[[form %||% ""]]
  sh[[shape %||% "default"]]$label %||% sh[["default"]]$label %||% ""
}

simulate_toy_data <- function(cfg = list()) {
  cfg <- utils::modifyList(TOY_DEFAULTS, cfg)
  for (nm in c("sd_dir", "sa_dir", "mean_ls", "n_id", "seed")) cfg[[nm]] <- as.numeric(cfg[[nm]])
  set.seed(as.integer(cfg$seed))
  type <- cfg$trait
  paper <- identical(type, "paper")
  is_count <- grepl("^count", type)
  form <- if (paper || !isTRUE(cfg$form %in% AGE_FUNCTIONS)) "Quadratic" else cfg$form
  dramatic <- identical(cfg$strength, "dramatic") && !paper
  n <- max(30L, as.integer(cfg$n_id))
  mean_ls <- if (paper) 25 else min(30, max(3, round(cfg$mean_ls)))
  cc <- mean_ls / 20
  individual_afr <- identical(cfg$afr_mode, "individual")
  sd_level <- cfg$sd_type %in% c("independent", "both")
  sd_rate <- cfg$sd_type %in% c("dependent", "both")
  sa_level <- individual_afr && cfg$sa_type %in% c("independent", "both")
  sa_rate <- individual_afr && cfg$sa_type %in% c("dependent", "both")
  sd_dir <- if (isTRUE(cfg$sd_dir < 0)) -1 else 1
  sa_dir <- if (isTRUE(cfg$sa_dir < 0)) -1 else 1

  z_ls <- stats::rnorm(n)
  z_afr <- stats::rnorm(n)
  if (paper) {
    mu_b <- c(600, 9, -0.3)
    sd_b <- abs(mu_b) * 0.2
    rho <- c(if (sd_level) 0.5 * sd_dir else 0, if (sd_rate) 0.5 * sd_dir else 0, 0)
    E <- matrix(stats::rnorm(n * 3), n, 3)
    B <- sapply(1:3, function(j) mu_b[j] + sd_b[j] * (rho[j] * z_ls + sqrt(1 - rho[j]^2) * E[, j]))
    B[, 1] <- B[, 1] + (if (sa_level) sa_dir * 0.5 * sd_b[1] else 0) * z_afr
    B[, 2] <- B[, 2] + (if (sa_rate) sa_dir * 0.5 * sd_b[2] else 0) * z_afr
    LS_raw <- 25 + 5 * z_ls
    sigma_e <- 5
    shift_diet <- -60
    u_sd <- 40
  } else {
    sp <- toy_form_spec(form, is_count)
    # biologically motivated shape within the form (population mean coefficients only; default = toy_form_spec)
    shape_mu <- toy_shape_mu(form, cfg$shape, is_count)
    if (length(shape_mu) == length(sp$mu)) sp$mu <- shape_mu
    k_str <- if (dramatic) 1 else 0.35
    s_level <- if (is_count) 0.2125 else 5.1
    np <- length(sp$mu)
    s_U <- c(if (is_count) 0.25 else 4, if (identical(cfg$rate_var, "high")) sp$u_large else sp$u_small)
    E <- matrix(stats::rnorm(n * np), n, np)
    B <- matrix(sp$mu, n, np, byrow = TRUE) + sweep(E, 2, s_U, "*")
    B[, 1] <- B[, 1] + k_str * s_level * ((if (sd_level) sd_dir * z_ls else 0) + (if (sa_level) sa_dir * z_afr else 0))
    B[, 2] <- B[, 2] + k_str * sp$s_rate * ((if (sd_rate) sd_dir * z_ls else 0) + (if (sa_rate) sa_dir * z_afr else 0))
    mu_b <- sp$mu
    LS_raw <- mean_ls + 0.3 * mean_ls * z_ls
    sigma_e <- 1.5
    shift_diet <- if (is_count) -0.3 else -5
    u_sd <- if (is_count) 0.25 else 4
  }

  diet <- if (isTRUE(cfg$diet)) stats::rbinom(n, 1, 0.5) else rep(0, n)
  n_groups <- 25L
  grp <- if (isTRUE(cfg$groups)) sample(seq_len(n_groups), n, replace = TRUE) else rep(1L, n)
  u_tr <- if (isTRUE(cfg$groups)) stats::rnorm(n_groups, 0, u_sd)[grp] else rep(0, n)
  B[, 1] <- B[, 1] + shift_diet * diet + u_tr
  LS <- pmax(if (paper) 1 else 2, round(LS_raw))

  if (individual_afr) {
    span <- max(2, round(0.4 * mean_ls))
    AFR <- 1 + floor(stats::pnorm(0.9 * z_afr + sqrt(1 - 0.81) * stats::rnorm(n)) * span)
    alive <- which(LS >= span)            # inclusion depends on lifespan only, not on AFR
  } else {
    AFR <- rep(1, n)
    alive <- seq_len(n)
  }
  cond <- as.numeric(scale(0.7 * z_ls + sqrt(1 - 0.49) * stats::rnorm(n)))

  n_age <- as.integer(LS[alive] - AFR[alive] + 1)
  idx <- rep(alive, n_age)
  age <- AFR[idx] + sequence(n_age) - 1
  u <- age / cc
  ref <- c(mean(u), stats::sd(u))
  if (!is.finite(ref[2]) || ref[2] <= 0) ref[2] <- 1
  X <- if (paper) cbind(1, age, age^2) else cbind(1, toy_form_basis(form, u, ref))
  eta <- rowSums(B[idx, , drop = FALSE] * X)
  m <- length(age)
  trait <- switch(type,
    paper = pmax(0, eta) + stats::rnorm(m, 0, sigma_e),
    mass = eta + stats::rnorm(m, 0, sigma_e),
    count_pois = stats::rpois(m, exp(pmin(eta, log(2000)))),
    count_nb = stats::rnbinom(m, size = 3, mu = exp(pmin(eta, log(2000)))),
    count_zinb = ifelse(stats::runif(m) < 0.25, 0, stats::rnbinom(m, size = 3, mu = exp(pmin(eta, log(2000)))))
  )

  s <- if (dramatic || paper) 1 else 0.6
  zt <- as.numeric(scale(if (is_count) log1p(trait) else trait))
  keep_p <- switch(cfg$missingness,
    complete = rep(1, m),
    mcar = rep(if (paper) 0.5 else 0.75, m),
    mwo = if (paper) 1 / (1 + exp(0.25 * (age - 40))) else 1 / (1 + exp(0.5 * s / cc * (age - 0.7 * mean_ls))),
    mwy = if (paper) 1 - 1 / (1 + exp(0.25 * (age - 10))) else 1 - 1 / (1 + exp(0.5 * s / cc * (age - 0.35 * mean_ls))),
    trait = stats::plogis(1.6 + 2 * s * zt),
    condition = stats::plogis(1.6 + 2 * s * cond[idx]),
    rep(1, m)
  )
  keep <- stats::runif(m) < keep_p

  trait_col <- if (identical(type, "mass")) "body_mass" else "fecundity"
  out <- data.frame(ID = sprintf("ID%04d", idx), age = age, trait = trait, lifespan = LS[idx],
                    condition = round(cond[idx], 3), stringsAsFactors = FALSE)
  names(out)[3] <- trait_col
  if (individual_afr) out$AFR <- AFR[idx]
  if (isTRUE(cfg$diet)) out$diet <- ifelse(diet[idx] == 1, "Poor", "Standard")
  if (isTRUE(cfg$groups)) out$family <- sprintf("F%02d", grp[idx])
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "truth") <- list(type = type, form = form, paper = paper, mu_b = mu_b, cc = cc, ref = ref,
                             shift_T = shift_diet, mean_T = mean(diet[alive]),
                             mean_expT = mean(exp(shift_diet * diet[alive])),
                             zi = if (identical(type, "count_zinb")) 0.25 else 0,
                             shape = if (paper) "default" else (cfg$shape %||% "default"), cfg = cfg)
  out
}

# Typical-individual trajectory implied by the simulation (population mean coefficients),
# marginalised over the diet covariate in the same way as the app's model predictions.
toy_true_curve <- function(truth, age) {
  X <- if (isTRUE(truth$paper) || is.null(truth$form)) cbind(1, age, age^2) else cbind(1, toy_form_basis(truth$form, age / truth$cc, truth$ref))
  base <- as.numeric(X %*% truth$mu_b)
  if (truth$type %in% c("paper", "mass")) {
    base + truth$shift_T * truth$mean_T
  } else {
    (1 - truth$zi) * exp(base) * truth$mean_expT
  }
}

toy_mapping <- function(df) {
  cols <- names(df)
  afr_ind <- "AFR" %in% cols
  list(id = "ID", age = "age", trait = if ("body_mass" %in% cols) "body_mass" else "fecundity",
       alr = "__AUTO_LAST__", life = "lifespan", entry = if (afr_ind) "AFR" else "__AUTO_FIRST__",
       condition = "condition", covars = intersect("diet", cols), cov_factor = intersect("diet", cols),
       group = if ("family" %in% cols) "family" else "", nested = TRUE, random = character(0),
       censor = "", censor_value = "", start_mode = if (afr_ind) "afr" else "same", start_age = 1, age_round = NA_real_,
       cov_age = character(0), cov_int = character(0))
}

toy_family <- function(type) {
  switch(type, mass = "gaussian", paper = "gaussian", count_pois = "poisson", count_nb = "nbinom2",
         count_zinb = "zinb", "gaussian")
}

# Expectations verified by simulation (random-intercept models, several seeds per scenario).
toy_truth_text <- function(cfg) {
  cfg <- utils::modifyList(TOY_DEFAULTS, cfg)
  cfg$mean_ls <- as.numeric(cfg$mean_ls)
  paper <- identical(cfg$trait, "paper")
  form <- if (paper) "Quadratic" else cfg$form
  trait <- names(TOY_TRAITS)[TOY_TRAITS == cfg$trait]
  mag <- if (paper) "|\u03c1| = 0.5 as in the Methods" else if (identical(cfg$strength, "dramatic")) "dramatic strength" else "moderate strength"
  L <- if (paper) 25 else cfg$mean_ls
  dir_txt <- function(d) if (isTRUE(as.numeric(d) < 0)) "negative" else "positive"
  sel_lab <- function(x) tolower(names(SELECTION_TYPES)[SELECTION_TYPES == x])
  individual_afr <- identical(cfg$afr_mode, "individual")
  sim <- paste0("Simulated: ", trait, "; ", tolower(form), " ageing",
                if (!paper) paste0(" (", tolower(toy_shape_label(form, cfg$shape)), ")") else "",
                "; selective disappearance: ", sel_lab(cfg$sd_type),
                if (cfg$sd_type != "none") paste0(" (", dir_txt(cfg$sd_dir), ")") else "",
                "; AFR ", if (individual_afr) paste0("individual-specific, selective appearance: ", sel_lab(cfg$sa_type),
                                                     if (cfg$sa_type != "none") paste0(" (", dir_txt(cfg$sa_dir), ")") else "")
                else "the same for everyone",
                "; ", mag, "; mean lifespan ", L, " occasions",
                if (!paper) paste0("; ", if (identical(cfg$rate_var, "high")) "large" else "small", " individual differences in ageing rates") else "",
                ".")
  expect <- character(0)
  if (cfg$sd_type %in% c("dependent", "both")) {
    expect <- c(expect,
      "Trajectory plot: lifespan (ALR/LS) bins diverge with age, so bin differences and their trend lines change with age.",
      "Lifespan plot: the LS\u2013trait regression coefficient changes across age bins (see the coefficient table).",
      "Models 4, 5 and 6 should beat Models 1\u20133 by a wide margin; the decomposition drifts from the truth at late ages.")
  } else if (identical(cfg$sd_type, "independent")) {
    expect <- c(expect,
      "Trajectory plot: bins are offset but parallel, so bin differences are similar at every age (flat trend lines away from zero).",
      "Lifespan plot: similar non-zero LS\u2013trait coefficients in every age bin.",
      "Model 2 should clearly beat Model 1; the interaction models add little.")
  } else {
    expect <- c(expect,
      "Trajectory plot: bins overlap and their differences scatter around zero; lifespan plot: coefficients near zero in every age bin.",
      "Model 1 should usually be best or within about 2 AIC of the best; occasional runs favour interactions by chance.")
  }
  flat <- !paper && identical(form, "Linear") && identical(cfg$shape, "none")
  expect <- c(expect, if (flat) {
    "No senescence: the simulated typical trajectory is flat, so any change of the observed means with age comes from selective disappearance (or appearance). A flat trajectory is a special case of every ageing function, so the ageing-function comparison need not single out one of them."
  } else paste0("The population-level ageing-function comparison should usually favour the ",
    form, " function", if (identical(form, "Asymptotic exponential")) " (approximately: the exponential basis is scaled by the SD of the sampled ages, which changes under age-dependent missingness)" else "",
    "; under selective disappearance Model 1 can prefer a more flexible function because the sampled records are distorted."))
  if (!paper && !identical(cfg$shape %||% "default", "default")) {
    expect <- c(expect, "These expectations were verified by simulation for the default shape of each ageing form. Other shapes change only the population mean coefficients, not the strength of selection or the individual variation, so the patterns should be similar, but AIC margins can differ.")
  }
  if (identical(cfg$rate_var, "high") || paper) {
    expect <- c(expect, "Individuals differ strongly in ageing rates: with random intercepts only, Models 4\u20135 can be favoured even without selective disappearance, because the interaction terms absorb among-individual differences in ageing (as in the manuscript). Choose a random slope under 'Random effects for individuals' on the Modelling tab and compare.")
  }
  miss <- switch(cfg$missingness,
    complete = "Complete sampling: ALR equals lifespan and mean age tracks it closely, so Models 4 and 5 should give similar AICs.",
    mcar = paste0("MCAR (", if (paper) "retention 0.5, as in the Methods" else "about 75% of occasions sampled", "): with age-dependent selection, Model 4 should beat Model 5 because ALR is the less error-prone lifespan proxy."),
    mwo = "Missing when old: late ages are thinly sampled and ALR underestimates lifespan; Model 6 (true LS) should be clearly best, and Model 4 should beat Model 5.",
    mwy = "Missing when young: first records are delayed, so AFR (first record) varies; with age-dependent selection Model 4 should beat Model 5.",
    trait = "Trait-dependent missingness (low values more often missing): the Sampling tab should flag the prior-trait association, and observed means are inflated.",
    condition = "Condition-dependent missingness (condition covaries with lifespan): the Sampling tab should flag condition and LS.")
  expect <- c(expect, miss)
  if (individual_afr) {
    expect <- c(expect, switch(cfg$sa_type,
      none = "AFR varies at random: AFR bins overlap in the trajectory plot (selective appearance question) and Models 7\u20138 should not improve on Models 2 and 4 (within about 2 AIC).",
      independent = "Age-independent selective appearance: AFR bins are offset but parallel in the trajectory plot, and Models 7 and 8 should clearly beat Models 2 and 4.",
      dependent = "Age-dependent selective appearance: AFR bins diverge with age in the trajectory plot, and Model 8 should clearly beat Model 4.",
      "Age-independent and age-dependent selective appearance: AFR bins are offset and diverge with age, and Models 7 and 8 should clearly beat Models 2 and 4."),
      "AFR is independent of lifespan by construction, so AFR and ALR should be nearly uncorrelated.")
    if (L < 8) expect <- c(expect, "With short lifespans and individual-specific AFR, individuals that die before the latest possible AFR are not sampled at all.")
  }
  if (isTRUE(cfg$diet)) expect <- c(expect, "Diet: the 'Poor' diet lowers the trait at all ages but does not affect lifespan. It is mapped as a covariate (and can be used to split the lifespan-group figures into panels); removing it adds unexplained variation but should not create selective-disappearance signals.")
  if (isTRUE(cfg$groups)) expect <- c(expect, "Families: 25 families share trait effects, mapped as a nested random effect (1 | family) + (1 | family:ID).")
  if (grepl("^count", cfg$trait)) expect <- c(expect, "Counts: use the matching family (Poisson / negative binomial / ZINB); compare with Gaussian to see how rankings can change.")
  if (L < 6) expect <- c(expect, "Few occasions per individual: individual fits are very imprecise and interaction models have little within-individual information.")
  list(sim = sim, expect = expect)
}

load_fly_example <- function() {
  f <- file.path("data", "fly_fecundity.csv")
  if (!file.exists(f)) return(NULL)
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

FLY_MAPPING <- list(
  id = "F1_ID", age = "F1_Age", trait = "F2_Count", alr = "ALR", life = "LS",
  entry = "__AUTO_FIRST__", condition = "", covars = c("Rep", "Paternal_age", "Paternal_sperm_age"),
  cov_factor = c("Rep", "Paternal_age", "Paternal_sperm_age"), cov_int = character(0),
  group = "F0_ID", nested = TRUE, random = character(0), censor = "Censored", censor_value = "0",
  start_mode = "same", start_age = 4, age_round = NA_real_, cov_age = character(0)
)

# ---------------------------------------------------------------------------
# Bundled empirical examples. Each entry carries the published data, the mapping
# and the model settings that reproduce (or come close to) the published analysis.
# ---------------------------------------------------------------------------
load_example_file <- function(file) {
  f <- file.path("data", file)
  if (!file.exists(f)) return(NULL)
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

example_map <- function(...) {
  base <- list(id = "", age = "", trait = "", alr = "__AUTO_LAST__", life = "", entry = "__AUTO_FIRST__",
               condition = "", covars = character(0), cov_factor = character(0), cov_int = character(0),
               group = "", nested = TRUE, random = character(0), censor = "", censor_value = "",
               start_mode = "afr", start_age = NA_real_, age_round = NA_real_, cov_age = character(0))
  utils::modifyList(base, list(...))
}

EXAMPLES <- list(
  fly = list(
    label = "Sanghvi et al. 2025, American Naturalist \u2014 Drosophila melanogaster (daily fecundity)",
    file = "fly_fecundity.csv", mapping = FLY_MAPPING,
    family = "zinb", age_function = "Quadratic", models = c("M1", "M2", "M3", "M4", "M5"),
    note = paste("Lifetime fecundity of laboratory female Drosophila melanogaster, assayed every 14 days until death, with lifespan,",
                 "replicate and paternal age treatments. The paper fitted zero-inflated negative binomial mixed models with",
                 "offspring nested in fathers, and found that accounting for lifespan changed the estimated ageing trajectory,",
                 "with the interaction models supported over the additive ones. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  bichet = list(
    label = "Bichet et al. 2022, Journal of Animal Ecology \u2014 common tern (immune parameters)",
    file = "bichet_2022_tern_immunity.csv",
    mapping = example_map(id = "ID", age = "age", trait = "HA",
                          covars = c("sex", "storage_time_HAHL", "initial_lysis", "assay_batch_HAHL"),
                          cov_factor = c("sex", "assay_batch_HAHL"), random = "year_sampling"),
    family = "gaussian", age_function = "Linear", models = c("M1", "M2", "M3", "M4"),
    note = paste("Innate immune measures (haemagglutination titre and haptoglobin) in common terns of known age, sampled over many",
                 "years. Age here is derived as year of sampling minus year of birth. The paper partitioned age into",
                 "among- and within-individual components in linear mixed models and found a within-individual increase in",
                 "haemagglutination with age, no change in haptoglobin, and no evidence of selective (dis)appearance in either.",
                 "Switch the trait to hapto for the second parameter. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  wynn = list(
    label = "Wynn et al. 2025, Journal of Animal Ecology \u2014 common tern (navigational efficiency)",
    file = "wynn_2025_tern_navigation.csv",
    mapping = example_map(id = "individual", age = "age", trait = "deflection",
                          covars = "season", cov_factor = "season", random = "track_id"),
    family = "gaussian", age_function = "Linear", models = c("M1", "M2", "M3", "M4", "M5"),
    note = paste("Navigational efficiency of migrating common terns, measured as the deflection of each track from the direction of",
                 "the goal. This response is not in the archived file and is derived here from the archived geolocator positions",
                 "(the angle between the bearing to the next fix and the bearing to the goal, with short steps treated as",
                 "stopovers and removed; see data/PROVENANCE.md). The paper fitted linear mixed models partitioning age into",
                 "among- and within-individual components, with track nested in individual, and found that older birds navigate",
                 "more efficiently among individuals but not within them, which they interpret as selective disappearance.",
                 "Track identity is mapped as an extra random intercept rather than as the individual. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  sanghvi_female = list(
    label = "Sanghvi et al. 2022, Evolution \u2014 seed beetle (female fecundity)",
    file = "sanghvi_2022_beetle_female_fecundity.csv",
    mapping = example_map(id = "ID", age = "Adult_age", trait = "Daily_Eggs", life = "Adult_lifespan",
                          covars = c("DevT", "AdultT", "Block"), cov_factor = c("DevT", "AdultT", "Block"),
                          cov_int = "DevT|||AdultT", cov_age = c("DevT", "AdultT"),
                          group = "Family", nested = TRUE),
    family = "nbinom2", age_function = "Quadratic", models = c("M1", "M2", "M3", "M4", "M5", "M6"),
    note = paste("Daily egg counts of female seed beetles from the same 2 \u00d7 2 temperature experiment. The paper fitted negative",
                 "binomial mixed models with the temperature treatments interacting with age and age\u00b2, adult lifespan as a",
                 "fixed effect for selective disappearance, and random slopes of age; it found that hot developmental and hot",
                 "adult temperatures each accelerated the decline in fecundity with age. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  allain = list(
    label = "Allain et al. 2023, Oikos \u2014 eastern chipmunk (reproduction)",
    file = "allain_2023_chipmunk_reproduction.csv",
    mapping = example_map(id = "ID", age = "age", trait = "nb_juv", entry = "AFR", life = "lifespan",
                          covars = c("season", "sex", "site"), cov_factor = c("season", "sex", "site"),
                          cov_age = "season"),
    family = "poisson", age_function = "Quadratic", models = c("M1", "M2", "M4", "M6", "M7", "M8", "M9", "M10"),
    note = paste("Reproductive events of wild eastern chipmunks, with age and age at first reproduction recorded in months and",
                 "lifespan known for individuals that died. The paper fitted generalised linear mixed models (binomial for the",
                 "probability of weaning, Poisson for the number of juveniles) comparing candidate models by AICc, and found",
                 "that age at first reproduction and its interaction with age improved the models substantially: reproductive",
                 "ageing depended on when individuals started breeding. weaned_juv is a binary alternative trait for the",
                 "binomial family. Lifespan is unknown for some individuals, so keeping Model 6 ticked restricts every model to",
                 "the records with a known lifespan, because AIC comparisons need the same rows throughout. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  bouwhuis = list(
    label = "Bouwhuis et al. 2009, Proc. R. Soc. B \u2014 great tit (recruit production)",
    file = "bouwhuis_2009_great_tit_recruitment.csv",
    mapping = example_map(id = "female", age = "f_min_age", trait = "LD", alr = "f_ALR",
                          covars = c("YR_FL", "loc_density", "f_status", "pred"),
                          cov_factor = c("f_status", "pred"), random = c("year", "area")),
    family = "poisson", age_function = "Quadratic", models = c("M1", "M2", "M3", "M4", "M5"),
    note = paste("Nearly five decades of breeding by female great tits in Wytham Woods: annual recruit production together with",
                 "clutch size, brood size and fledgling number, with year quality, breeding density, female status and nest-box",
                 "type as covariates and year and wood sector as random effects. The paper fitted cross-classified mixed models",
                 "and showed that selective disappearance of poorer breeders masks part of the within-individual decline, so the",
                 "onset of senescence is earlier than a population-level analysis suggests, while the component traits show age",
                 "effects without a lifespan effect. Age at last reproduction is supplied in the data and deliberately differs",
                 "from the last recorded age for some females, which the integrity table flags. The trait opens on laying date",
                 "(LD); recruits, clutch size (CS), brood size (BS) and fledgling number (FL), the traits analysed in the",
                 "paper, can be selected from the trait dropdown on the Data tab. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  warner = list(
    label = "Warner et al. 2016, PNAS \u2014 painted turtle (reproduction)",
    file = "warner_2016_turtle_reproduction.csv",
    mapping = example_map(id = "Female ID", age = "Reproductive Age", trait = "Avg Egg Mass (g)",
                          covars = "Plastron Length (mm)", random = "Year"),
    family = "gaussian", age_function = "Quadratic", models = c("M1", "M2", "M3", "M4", "M5"),
    note = paste("More than two decades of nesting by wild painted turtles: egg mass, clutch size and hatching success by",
                 "reproductive age, with plastron length as a covariate and maternal identity as a random effect. The paper",
                 "fitted linear mixed models and found that egg mass increases with reproductive age while clutch size does not,",
                 "alongside separate mark-recapture analyses of mortality that are outside the scope of this app. Some records",
                 "are recorded as UNK and are read as missing, and females nest several times per season, so duplicate",
                 "ID \u00d7 age records are expected. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  mckennaell_breeding = list(
    label = "McKenna-Ell et al. 2023, Biology Letters \u2014 Soay sheep (breeding probability, offspring survival)",
    file = "mckennaell_2023_soay_breeding_survival.csv",
    mapping = example_map(id = "FemaleID", age = "Age", trait = "Fecundity", alr = "AgeLastObs",
                          covars = c("BredYearling", "EarlyLifeRec"), cov_factor = "BredYearling",
                          cov_age = c("BredYearling", "EarlyLifeRec"), random = c("ObsYear", "FemaleCohort")),
    family = "binomial", age_function = "Linear", models = c("M1", "M2", "M3", "M4", "M5"),
    note = paste("Annual reproduction of known-age female Soay sheep on St Kilda from age 5 onwards (later life), with two binary",
                 "traits: whether a female gave birth to a live lamb (Fecundity) and, for females that did, whether at least one",
                 "lamb survived its first winter (OffspringRecruitment). The paper fitted binomial generalised linear mixed models",
                 "with a linear effect of age, age at last observation for selective disappearance, whether the female bred as a",
                 "yearling (a factor) and her early-life recruitment (lambs raised in her first years), each early-life measure",
                 "also interacting with age, and random intercepts for female, year and birth cohort. It found senescent declines",
                 "and selective disappearance in both traits, a faster decline in breeding probability in females that bred as",
                 "yearlings, and no effect of early-life reproduction on the rate of ageing in offspring survival. Model 2 with",
                 "these covariates and interactions is the published model. Switch the trait to OffspringRecruitment for offspring",
                 "survival: records without a lamb are then blank and are not modelled. With 'Standardise' ticked, the early-life",
                 "main effects refer to the mean age rather than to age 0; untick it to compare every coefficient with the",
                 "published table. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  mckennaell_weight = list(
    label = "McKenna-Ell et al. 2023, Biology Letters \u2014 Soay sheep (offspring birth weight)",
    file = "mckennaell_2023_soay_offspring_weight.csv",
    mapping = example_map(id = "FemaleID", age = "Age", trait = "OffspringBirthWt", alr = "AgeLastObs",
                          covars = c("OffspringCaptureAge", "OffspringSex", "OffspringTwinStatus", "BredYearling", "EarlyLifeRec"),
                          cov_factor = c("OffspringSex", "OffspringTwinStatus", "BredYearling"),
                          cov_age = c("BredYearling", "EarlyLifeRec"), random = c("ObsYear", "FemaleCohort")),
    family = "gaussian", age_function = "Linear", models = c("M1", "M2", "M3", "M4", "M5"),
    note = paste("Birth weight of lambs born to known-age female Soay sheep aged 5 and older, one row per lamb. The paper fitted a",
                 "Gaussian linear mixed model with the lamb's age at capture, its sex and twin status, the mother's age (linear),",
                 "her age at last observation for selective disappearance, whether she bred as a yearling and her early-life",
                 "recruitment, both also interacting with age, and random intercepts for mother, year and the mother's birth",
                 "cohort. It found that birth weight declined with maternal age, that mothers observed to older ages had heavier",
                 "lambs (selective disappearance), and that early-life reproduction did not change the rate of this decline.",
                 "Model 2 with these covariates and interactions is the published model. Twins give their mother two records in",
                 "the same year: keep 'Keep all' so that each lamb is an observation, as in the paper. Age at last observation is",
                 "later than the last weighed lamb for mothers that stopped breeding, which the integrity table flags. With",
                 "'Standardise' ticked, the early-life main effects refer to the mean age; untick it to compare every coefficient",
                 "with the published table. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses.")),
  szejnersigal_activity = list(
    label = "Szejner-Sigal et al. 2025, Proc. R. Soc. B \u2014 alfalfa leafcutting bee (locomotor activity)",
    file = "szejnersigal_2025_bee_activity.csv",
    mapping = example_map(id = "id", age = "age", trait = "total.act", life = "age.death"),
    family = "gaussian", age_function = "Quadratic", models = c("M1", "M2", "M4", "M6"),
    extra = list(M1 = "LS"), subset = list(var = "sex", levels = "f"),
    note = paste("Locomotor activity (beam breaks in four hours) of individually marked alfalfa leafcutting bees (Megachile",
                 "rotundata) measured weekly in the laboratory from emergence until death (days 1, 7, 14, 21 and so on), with",
                 "each bee's age at death. The paper fitted linear mixed models separately for females and males, with age as a",
                 "quadratic (a negative parabola), lifespan as a covariate and a random intercept for each bee, and found that",
                 "activity rose to a peak in mid-life and then declined, earlier and at lower levels in males, with little link",
                 "between early-life activity and lifespan. The example opens on females (Subset the data: sex = f); choose m for",
                 "males. Age at death is mapped as the known lifespan (LS) and LS is added to Model 1, so Model 1 is the published",
                 "model (age, age\u00b2 and lifespan); Model 2 uses the age at the last weekly record instead, and Model 6 lets",
                 "lifespan interact with age. The column age.group20 holds the paper's short-, average- and long-lived groups and",
                 "can be used for panels. The first interval is six days (day 1 to day 7), so the integrity table reports ages off",
                 "the weekly schedule; the sampling grid assigns each record to the nearest weekly occasion. ",
                 "The settings in the app are defaulted to be closely aligned to the original fit reported in the study for",
                 "that trait, and users can change any of them. The defaults generally reproduce the published pattern, but",
                 "results may differ for several reasons: how the data were subset, model terms or covariates that were not",
                 "fully specified in the paper, different random-effect structures, and differences in software or estimation.",
                 "These examples are exploratory rather than exact re-runs of the published analyses."))
)
EXAMPLE_CHOICES <- stats::setNames(names(EXAMPLES), vapply(EXAMPLES, function(e) e$label, character(1)))

# ---------------------------------------------------------------------------
# Column guessing for uploaded data
# ---------------------------------------------------------------------------
guess_mapping <- function(df) {
  cols <- names(df)
  low <- gsub("[^a-z0-9]+", "_", tolower(cols))
  n <- nrow(df)
  is_num <- vapply(df, function(x) {
    ch <- as.character(x)
    present <- !is.na(x) & nzchar(ch)
    if (!any(present)) return(FALSE)
    v <- suppressWarnings(as.numeric(ch[present]))
    isTRUE(mean(is.finite(v)) > 0.95)
  }, logical(1))
  nu <- vapply(df, function(x) length(unique(x[!is.na(x)])), integer(1))
  find_col <- function(patterns, pool) {
    for (p in patterns) {
      h <- which(pool & low == p)
      if (length(h)) return(cols[[h[[1]]]])
    }
    for (p in patterns) {
      h <- which(pool & grepl(p, low, fixed = TRUE))
      if (length(h)) return(cols[[h[[1]]]])
    }
    ""
  }

  id_hits <- which(nu >= 2 & nu < n & grepl("(^|_)id($|_)|individual|animal|subject|ring|(^|_)ind($|_)", low))
  id <- if (length(id_hits)) cols[[id_hits[[which.max(nu[id_hits])]]]] else cols[[1]]

  excl_age <- grepl("mean|delta|centr|age2|age_2|_sq|sq_|squared|paternal|maternal|sperm|alr|afr|first|last|entry|death|lifespan", low)
  age <- find_col(c("age", "time", "occasion", "year"), is_num & !excl_age & cols != id)
  if (!nzchar(age)) age <- cols[[min(2, length(cols))]]

  excl_trait <- grepl("(^|_)ls($|_)|lifespan|alr|afr|mean|delta|age|censor|(^|_)rep($|_)|obs|(^|_)id($|_)|year", low)
  trait <- find_col(c("trait", "count", "fecund", "fertil", "offspring", "egg", "clutch", "value", "phenotype",
                      "mass", "weight", "size", "date", "score"), is_num & !excl_trait & cols != id & cols != age)
  if (!nzchar(trait)) {
    h <- which(is_num & !excl_trait & cols != id & cols != age & nu > 2)
    trait <- if (length(h)) cols[[h[[1]]]] else cols[[min(3, length(cols))]]
  }

  life <- find_col(c("ls", "lifespan", "life_span", "longevity", "age_at_death", "death_age"), is_num)
  alr <- find_col(c("alr", "age_last_record", "last_record", "age_at_last"), is_num)
  afr <- find_col(c("afr", "age_first_record", "first_record", "age_at_first", "entry_age"), is_num)
  cond <- find_col(c("condition", "state", "body_condition", "quality"), rep(TRUE, length(cols)))
  list(
    id = id, age = age, trait = trait,
    alr = if (nzchar(alr)) alr else "__AUTO_LAST__",
    life = life,
    entry = if (nzchar(afr)) afr else "__AUTO_FIRST__",
    condition = cond, covars = character(0), cov_factor = character(0), cov_int = character(0), group = "", nested = TRUE, random = character(0),
    censor = "", censor_value = "", start_mode = "afr", start_age = NA_real_, age_round = NA_real_, cov_age = character(0)
  )
}

cov_name <- function(nm) paste0("cv_", make.names(nm))
re_name <- function(nm) paste0("re_", make.names(nm))
display_term <- function(x) gsub("(^|[^A-Za-z0-9_])(cv|re)_", "\\1", x)

# UI helpers shared by ui.R and server.R (objects defined in ui.R are not visible to the server)
sign_choices <- c("\u2212 (negative)" = "-1", "0 (none)" = "0", "+ (positive)" = "1")
box_note <- function(...) div(class = "small-note", ...)

# Human-readable random-effect structure using the user's column names.
# Individual-level random-effect structures. TRUE/FALSE (older settings) map to correlated/none.
RANDOM_STRUCTURES <- c("Random intercept: (1 | ID)" = "none",
                       "Automatic: chosen from the data support" = "auto",
                       "Uncorrelated intercept and slope: (1 | ID) + (0 + age | ID)" = "uncorrelated",
                       "Correlated intercept and slope: (1 + age | ID)" = "correlated")
normalise_slope <- function(x) {
  if (isTRUE(x)) return("correlated")
  if (is.null(x) || !length(x) || isFALSE(x) || is.na(x[[1]])) return("none")
  x <- as.character(x[[1]])
  if (x %in% c("none", "auto", "uncorrelated", "correlated")) x else "none"
}
slope_text <- function(x) {
  switch(normalise_slope(x), correlated = "correlated random intercept and slope", uncorrelated = "uncorrelated random intercept and slope",
         auto = "random slope chosen automatically", "random intercept only")
}

random_display <- function(meta, random_slope = FALSE, slope_label = "age") {
  map <- meta$map
  g2 <- isTRUE(meta$has_group2)
  grp_lab <- if (g2 && isTRUE(meta$nested)) paste0(map$group2, ":", map$group) else map$group
  id_lab <- if (isTRUE(meta$has_group) && isTRUE(meta$nested)) paste0(grp_lab, ":", map$id) else map$id
  rs <- normalise_slope(random_slope)
  terms <- c(
    if (g2) paste0("(1 | ", map$group2, ")"),
    if (isTRUE(meta$has_group)) paste0("(1 | ", grp_lab, ")"),
    switch(rs,
      correlated = paste0("(1 + ", slope_label, " | ", id_lab, ")"),
      uncorrelated = paste0("(1 | ", id_lab, ") + (0 + ", slope_label, " | ", id_lab, ")"),
      auto = paste0("(1 | ", id_lab, ") with a random slope if the data support it"),
      paste0("(1 | ", id_lab, ")")),
    if (length(meta$random_terms)) paste0("(1 | ", unname(meta$random_labels[meta$random_terms]), ")")
  )
  paste(terms, collapse = " + ")
}

# ---------------------------------------------------------------------------
# Robust reading of uploaded delimited files
# ---------------------------------------------------------------------------
# Detects the separator (comma, semicolon, tab, pipe) and decimal commas, strips a UTF-8 byte-order
# mark, treats common missing-value codes as NA, trims white space, keeps every column as text (so IDs
# such as "007" are not turned into 7; numbers are converted where they are used), repairs blank or
# duplicated column names and drops empty rows and columns. Falls back to latin1 for non-UTF-8 files.
read_user_csv <- function(path, file_name = "") {
  ext <- tolower(tools::file_ext(if (nzchar(file_name %||% "")) file_name else path))
  if (ext %in% c("xls", "xlsx", "xlsm", "ods")) {
    return(list(data = NULL, note = "This is a spreadsheet workbook: save the sheet as CSV (comma- or semicolon-separated) and upload that file."))
  }
  magic <- tryCatch(readBin(path, "raw", n = 8), error = function(e) raw(0))
  if (length(magic) >= 4 && (identical(magic[1:4], as.raw(c(0x50, 0x4b, 0x03, 0x04))) ||
                             identical(magic[1:4], as.raw(c(0xd0, 0xcf, 0x11, 0xe0))))) {
    return(list(data = NULL, note = "This file is a spreadsheet workbook (Excel), not a CSV, even if it is named .csv: open it and save the sheet as CSV."))
  }
  first <- tryCatch(suppressWarnings(readLines(path, n = 60, warn = FALSE, encoding = "UTF-8")), error = function(e) character(0))
  # invalid UTF-8 (e.g. latin1 files) would make the regular expressions below fail: replace bad bytes for sniffing only
  first <- iconv(first, from = "UTF-8", to = "UTF-8", sub = "?")
  first[is.na(first)] <- ""
  first <- first[nzchar(trimws(first))]
  if (!length(first)) return(list(data = NULL, note = "The file is empty or could not be read."))
  hdr <- sub("^\ufeff", "", first[[1]])
  seps <- c(",", ";", "\t", "|")
  cnt <- vapply(seps, function(s) {
    m <- gregexpr(s, hdr, fixed = TRUE)[[1]]
    as.numeric(sum(m > 0))
  }, numeric(1))
  sep <- if (max(cnt) > 0) seps[[which.max(cnt)]] else ","
  dec <- "."
  if (!identical(sep, ",") && length(first) > 1) {
    fields <- trimws(gsub("\"", "", unlist(strsplit(first[-1], sep, fixed = TRUE))))
    if (sum(grepl("^-?[0-9]+,[0-9]+$", fields)) > sum(grepl("^-?[0-9]+\\.[0-9]+$", fields))) dec <- ","
  }
  na_codes <- c("NA", "", ".", "-", "NaN", "N/A", "n/a", "na", "#N/A", "NULL", "null")
  read_with <- function(enc) {
    bad <- FALSE
    val <- tryCatch(withCallingHandlers(
      utils::read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, na.strings = na_codes,
                      strip.white = TRUE, colClasses = "character", fileEncoding = enc, comment.char = "",
                      blank.lines.skip = TRUE),
      warning = function(w) {
        if (grepl("invalid|incomplete|embedded nul", conditionMessage(w))) bad <<- TRUE
        invokeRestart("muffleWarning")
      }), error = function(e) NULL)
    if (bad) NULL else val
  }
  x <- read_with("UTF-8-BOM")
  if (is.null(x)) x <- read_with("latin1")
  if (is.null(x)) return(list(data = NULL, note = "The file could not be parsed as a delimited text (CSV) file."))
  nm <- trimws(sub("^\ufeff", "", names(x)))
  blank <- is.na(nm) | !nzchar(nm)
  nm[blank] <- paste0("column_", which(blank))
  names(x) <- make.unique(nm)
  x[] <- lapply(x, function(col) {
    col <- trimws(col)
    col[!is.na(col) & col %in% na_codes] <- NA
    if (identical(dec, ",")) {
      ok <- !is.na(col)
      if (any(ok) && mean(grepl("^-?[0-9]+(,[0-9]+)?$", col[ok])) > 0.9) col <- gsub(",", ".", col, fixed = TRUE)
    }
    col
  })
  if (ncol(x)) x <- x[, vapply(x, function(col) any(!is.na(col)), logical(1)), drop = FALSE]
  if (ncol(x)) x <- x[rowSums(!is.na(x)) > 0, , drop = FALSE]
  rownames(x) <- NULL
  sep_name <- switch(sep, `;` = "semicolon", `\t` = "tab", `|` = "pipe", "comma")
  note <- sprintf("Read %d rows and %d columns (%s-separated%s).", nrow(x), ncol(x), sep_name,
                  if (identical(dec, ",")) ", decimal commas converted" else "")
  list(data = x, note = note)
}

# Numeric if at least 95% of the non-missing values are numbers, otherwise text.
maybe_numeric <- function(x) {
  if (is.numeric(x) || is.logical(x)) return(x)
  ch <- trimws(as.character(x))
  present <- !is.na(ch) & nzchar(ch)
  if (!any(present)) return(x)
  v <- suppressWarnings(as.numeric(ch))
  if (mean(is.finite(v[present])) >= 0.95) {
    v[!is.finite(v)] <- NA_real_
    v
  } else {
    ch[!present] <- NA_character_
    ch
  }
}

# ---------------------------------------------------------------------------
# Standardise the user's data into the internal format
# ---------------------------------------------------------------------------
# The column of the second grouping level (the group containing the group, e.g. family above father), or "" when it
# is not mapped, missing from the data, or duplicates the group, ID, age or trait column.
group2_column <- function(map, cols) {
  g2 <- map$group2 %||% ""
  g1 <- map$group %||% ""
  ok <- length(g2) == 1 && !is.na(g2) && nzchar(g2) && g2 %in% cols &&
    length(g1) == 1 && !is.na(g1) && nzchar(g1) && g1 %in% cols && !g2 %in% c(g1, map$id, map$age, map$trait)
  if (isTRUE(ok)) g2 else ""
}

# Individual identifiers as the models see them, built from the raw columns exactly as in standardise_data():
# with nesting, group + ID, and top-level group + group + ID when a second grouping level is mapped.
individual_ids <- function(df, map) {
  ids <- trimws(as.character(df[[map$id]]))
  clean <- function(nm) {
    v <- trimws(as.character(df[[nm]]))
    v[!is.na(v) & !nzchar(v)] <- NA_character_
    v
  }
  g1 <- map$group %||% ""
  if (!isTRUE(map$nested) || !(length(g1) == 1 && !is.na(g1) && nzchar(g1) && g1 %in% names(df))) return(ids)
  g <- clean(g1)
  g2 <- group2_column(map, names(df))
  if (nzchar(g2)) {
    h <- clean(g2)
    g <- ifelse(is.na(g) | is.na(h), g, paste(h, g, sep = "/"))
  }
  ifelse(is.na(g), ids, paste(g, ids, sep = "/"))
}

# Checks shown under the covariate boxes of the Data tab when a column seems to be in the wrong box: text in a
# CONTINUOUS covariate (values that are not numbers become missing), or a CATEGORICAL covariate that looks
# continuous (numeric with decimals or many distinct values, so each value would become its own level).
covariate_type_warnings <- function(df, num = character(0), fac = character(0)) {
  num <- intersect(num %||% character(0), names(df))
  fac <- intersect(fac %||% character(0), names(df))
  msgs <- character(0)
  for (nm in intersect(num, fac)) {
    msgs <- c(msgs, sprintf("Warning: '%s' is in both boxes and is used as CATEGORICAL: remove it from one of them.", nm))
  }
  for (nm in setdiff(num, fac)) {
    ch <- trimws(as.character(df[[nm]]))
    present <- !is.na(ch) & nzchar(ch)
    if (!any(present)) next
    v <- suppressWarnings(as.numeric(ch[present]))
    bad <- !is.finite(v)
    ex <- paste0("'", utils::head(unique(ch[present][bad]), 3), "'", collapse = ", ")
    if (mean(bad) > 0.2) {
      msgs <- c(msgs, sprintf("Warning: '%s' is mapped as CONTINUOUS but %.0f%% of its values are not numbers (e.g. %s), so it looks CATEGORICAL. As a continuous covariate those values are treated as missing: move it to the CATEGORICAL box.",
                              nm, 100 * mean(bad), ex))
    } else if (any(bad)) {
      msgs <- c(msgs, sprintf("Note: %d values of the CONTINUOUS covariate '%s' are not numbers (e.g. %s) and are treated as missing.", sum(bad), nm, ex))
    } else {
      u <- sort(unique(v))
      # (whole numbers counting from zero, such as numbers of offspring, are left alone: they are counts, not codes)
      if (length(u) >= 3 && length(u) <= 5 && all(abs(u - round(u)) < 1e-8) && !(u[[1]] == 0 && all(diff(u) == 1))) {
        msgs <- c(msgs, sprintf("Note: the CONTINUOUS covariate '%s' has only %d distinct whole-number values (%s). If these are codes for groups (e.g. treatments or blocks), move it to the CATEGORICAL box.",
                                nm, length(u), paste(format(u, trim = TRUE), collapse = ", ")))
      }
    }
  }
  for (nm in setdiff(fac, num)) {
    ch <- trimws(as.character(df[[nm]]))
    present <- !is.na(ch) & nzchar(ch)
    if (!any(present)) next
    v <- suppressWarnings(as.numeric(ch[present]))
    if (mean(is.finite(v)) < 0.95) next
    u <- unique(v[is.finite(v)])
    decimals <- any(abs(u - round(u)) > 1e-8)
    if ((decimals && length(u) > 5) || length(u) > 20) {
      msgs <- c(msgs, sprintf("Warning: '%s' is mapped as CATEGORICAL but looks CONTINUOUS (%d distinct numeric values%s): each value would become its own level. Move it to the CONTINUOUS box unless the numbers are group codes.",
                              nm, length(u), if (decimals) ", with decimals" else ""))
    }
  }
  msgs
}

# Message for the column mapped as age: age must be numeric; values that are not numbers are dropped.
age_type_message <- function(x, col = "age") {
  ch <- trimws(as.character(x))
  present <- !is.na(ch) & nzchar(ch)
  if (!any(present)) return(sprintf("Warning: '%s' has no values.", col))
  v <- suppressWarnings(as.numeric(ch[present]))
  bad <- !is.finite(v)
  if (!any(bad)) return("")
  ex <- paste0("'", utils::head(unique(ch[present][bad]), 3), "'", collapse = ", ")
  sprintf("Warning: %d of %d values in '%s' are not numbers (e.g. %s). Age must be numeric, so these rows are dropped.%s",
          sum(bad), sum(present), col, ex,
          if (mean(bad) > 0.5) " The column looks categorical: map a numeric age column instead (convert dates or age classes to numbers first)." else "")
}

standardise_data <- function(df, map, dup_action = "keep") {
  n_raw <- nrow(df)
  has_col <- function(nm) length(nm) == 1 && !is.na(nm) && nzchar(nm) && nm %in% names(df)
  out <- data.frame(
    id = trimws(as.character(df[[map$id]])),
    age = safe_numeric(df[[map$age]]),
    trait = safe_numeric(df[[map$trait]]),
    .row = seq_len(n_raw),
    stringsAsFactors = FALSE
  )
  # optional rounding of ages to a sampling resolution (e.g. irregular field ages)
  age_res <- suppressWarnings(as.numeric(map$age_round %||% NA_real_))
  rounding <- length(age_res) == 1 && is.finite(age_res) && age_res > 0
  round_age <- function(a) if (rounding) round(a / age_res) * age_res else a
  out$age <- round_age(out$age)
  keep <- !is.na(out$id) & nzchar(out$id) & is.finite(out$age)
  out <- out[keep, , drop = FALSE]
  rows <- out$.row
  n_dropped <- n_raw - nrow(out)

  has_group <- has_col(map$group)
  # Optional second grouping level above the group (e.g. individuals within fathers within families).
  group2_col <- group2_column(map, names(df))
  has_group2 <- has_group && nzchar(group2_col)
  out$group <- if (has_group) trimws(as.character(df[[map$group]][rows])) else rep(NA_character_, nrow(out))
  out$group[!is.na(out$group) & !nzchar(out$group)] <- NA_character_
  out$group2 <- if (has_group2) trimws(as.character(df[[group2_col]][rows])) else rep(NA_character_, nrow(out))
  out$group2[!is.na(out$group2) & !nzchar(out$group2)] <- NA_character_
  n_multi_group <- 0L
  n_multi_group2 <- 0L
  if (has_group && nrow(out)) {
    if (has_group2) {
      ng2 <- tapply(out$group2, out$group, function(g) length(unique(g[!is.na(g)])))
      n_multi_group2 <- sum(ng2 > 1)
      # Three levels: a group is identified by top-level group + group (equivalent to (1 | group2/group)).
      if (isTRUE(map$nested)) out$group <- ifelse(is.na(out$group) | is.na(out$group2), out$group, paste(out$group2, out$group, sep = "/"))
    }
    ng <- tapply(out$group, out$id, function(g) length(unique(g[!is.na(g)])))
    n_multi_group <- sum(ng > 1)
    # Nested design: an individual is identified by group + ID (equivalent to (1 | group/ID)).
    if (isTRUE(map$nested)) out$id <- ifelse(is.na(out$group), out$id, paste(out$group, out$id, sep = "/"))
  }

  if (nrow(out)) {
    last_rec <- stats::ave(out$age, out$id, FUN = max)
    first_rec <- stats::ave(out$age, out$id, FUN = min)
  } else {
    last_rec <- numeric(0)
    first_rec <- numeric(0)
  }

  alr_mapped <- has_col(map$alr) && !identical(map$alr, "__AUTO_LAST__")
  out$alr <- if (alr_mapped) round_age(safe_numeric(df[[map$alr]])[rows]) else last_rec

  life_auto <- identical(map$life, "__AUTO_LAST__")
  out$life <- if (life_auto) last_rec else if (has_col(map$life)) safe_numeric(df[[map$life]])[rows] else rep(NA_real_, nrow(out))
  # Censored individuals (e.g. escaped or still alive at the end of the study): their recorded
  # lifespan is not a death age, so LS is set to missing for every row of those individuals.
  n_censored <- 0L
  censored_ids <- character(0)
  has_censor <- has_col(map$censor) && !life_auto && nrow(out) > 0
  if (has_censor) {
    cens_row <- as.character(df[[map$censor]][rows]) == as.character(map$censor_value %||% "")
    cens_row[is.na(cens_row)] <- FALSE
    cens_ids <- unique(out$id[cens_row])
    n_censored <- length(cens_ids)
    censored_ids <- cens_ids
    out$life[out$id %in% cens_ids] <- NA_real_
  }

  # AFR only. The age at first trait expression is never used here: it affects the missingness window alone.
  entry_mapped <- has_col(map$entry) && !identical(map$entry, "__AUTO_FIRST__")
  out$entry <- if (entry_mapped) round_age(safe_numeric(df[[map$entry]])[rows]) else first_rec

  # Number of binomial trials behind each proportion (optional; used as prior weights).
  has_trials <- has_col(map$trials)
  out$.trials <- if (has_trials) safe_numeric(df[[map$trials]])[rows] else rep(NA_real_, nrow(out))
  n_bad_trials <- if (has_trials) sum(!is.finite(out$.trials) | out$.trials <= 0) else 0L
  out$condition <- rep(NA_real_, nrow(out))
  out$condition_label <- rep(NA_character_, nrow(out))
  if (has_col(map$condition)) {
    c0 <- df[[map$condition]][rows]
    cn <- safe_numeric(c0)
    if (sum(is.finite(cn)) >= 0.8 * sum(!is.na(c0) & nzchar(as.character(c0)))) {
      out$condition <- cn
    } else {
      out$condition_label <- as.character(c0)
    }
  }

  covars <- character(0)
  cov_labels <- character(0)
  cov_types <- character(0)
  cov_age <- character(0)
  cov_pairs <- character(0)
  skipped <- character(0)
  # map$cov_factor lists categorical covariates; the others are continuous. Without it (older presets),
  # a covariate is continuous when at least 95% of its values are numbers.
  declared <- !is.null(map$cov_factor)
  for (nm in unique(map$covars %||% character(0))) {
    if (!has_col(nm) || nm %in% c(map$id, map$age, map$trait, map$group, if (has_group2) group2_col)) next
    raw_vals <- df[[nm]][rows]
    vals <- if (!declared) {
      maybe_numeric(raw_vals)
    } else if (nm %in% map$cov_factor) {
      v <- trimws(as.character(raw_vals))
      v[!is.na(v) & !nzchar(v)] <- NA_character_
      v
    } else {
      safe_numeric(raw_vals)
    }
    if (!any(!is.na(vals) & nzchar(as.character(vals)))) {
      skipped <- c(skipped, nm)
      next
    }
    cn <- cov_name(nm)
    while (cn %in% names(out)) cn <- paste0(cn, "_")
    out[[cn]] <- vals
    if (nm %in% (map$cov_age %||% character(0))) cov_age <- c(cov_age, cn)
    covars <- c(covars, cn)
    cov_labels[cn] <- nm
    cov_types[cn] <- if (is.numeric(vals)) "continuous" else "categorical"
  }
  # interactions chosen as "a|||b" (two covariates) or "a|||age" (covariate x ageing terms)
  # (with no covariates, names(cov_labels) is NULL and setNames() would fail on it)
  internal_of <- if (length(cov_labels)) stats::setNames(names(cov_labels), unname(cov_labels)) else character(0)
  for (pr in unique(map$cov_int %||% character(0))) {
    parts <- strsplit(pr, "|||", fixed = TRUE)[[1]]
    if (length(parts) != 2 || !parts[[1]] %in% names(internal_of)) next
    if (identical(parts[[2]], "age")) {
      cov_age <- union(cov_age, internal_of[[parts[[1]]]])
    } else if (parts[[2]] %in% names(internal_of) && !identical(parts[[1]], parts[[2]])) {
      cov_pairs <- union(cov_pairs, paste(internal_of[[parts[[1]]]], internal_of[[parts[[2]]]], sep = ":"))
    }
  }

  random_terms <- character(0)
  random_labels <- character(0)
  for (nm in map$random %||% character(0)) {
    if (!has_col(nm) || nm %in% c(map$id, map$age, map$trait, map$group, if (has_group2) group2_col)) next
    vals <- trimws(as.character(df[[nm]][rows]))
    vals[!is.na(vals) & !nzchar(vals)] <- NA_character_
    if (all(is.na(vals))) {
      skipped <- c(skipped, nm)
      next
    }
    rn <- re_name(nm)
    while (rn %in% names(out)) rn <- paste0(rn, "_")
    out[[rn]] <- vals
    random_terms <- c(random_terms, rn)
    random_labels[rn] <- nm
  }

  n_dup_raw <- if (nrow(out)) sum(duplicated(paste(out$id, out$age, sep = "\r"))) else 0L
  if (identical(dup_action, "mean") && n_dup_raw > 0) {
    key <- paste(out$id, out$age, sep = "\r")
    tr_mean <- tapply(out$trait, key, finite_mean)
    first <- !duplicated(key)
    out <- out[first, , drop = FALSE]
    out$trait <- as.numeric(tr_mean[key[first]])
  }

  out$.row <- NULL
  out <- out[order(out$id, out$age), , drop = FALSE]
  rownames(out) <- NULL
  meta <- list(
    n_raw = n_raw, n_dropped = n_dropped, n_dup = n_dup_raw, dup_action = dup_action,
    alr_mapped = alr_mapped, life_auto = life_auto, has_life = life_auto || has_col(map$life),
    entry_mapped = entry_mapped, has_group = has_group, nested = isTRUE(map$nested),
    n_multi_group = n_multi_group, has_group2 = has_group2, n_multi_group2 = n_multi_group2, covars = covars, cov_labels = cov_labels,
    random_terms = random_terms, random_labels = random_labels,
    has_trials = has_trials, n_bad_trials = n_bad_trials, trials_col = if (has_trials) map$trials else "",
    has_censor = has_censor, n_censored = n_censored, age_round = if (rounding) age_res else NA_real_,
    skipped_columns = unique(skipped), cov_age = cov_age, cov_types = cov_types, cov_pairs = cov_pairs,
    censored_ids = censored_ids, map = map
  )
  list(data = out, meta = meta)
}

# ---------------------------------------------------------------------------
# Individual-level summaries (vectorised)
# ---------------------------------------------------------------------------
individual_metrics <- function(dat) {
  if (!nrow(dat)) return(data.frame())
  ids <- unique(dat$id)
  f <- factor(dat$id, levels = ids)
  okt <- is.finite(dat$trait)
  fo <- factor(dat$id[okt], levels = ids)
  num_mean <- function(v) as.numeric(tapply(v, f, finite_mean))
  out <- data.frame(
    id = ids,
    n_rows = tabulate(as.integer(f), nbins = length(ids)),
    n_trait = tabulate(as.integer(fo), nbins = length(ids)),
    first_recorded = as.numeric(tapply(dat$age, f, min)),
    last_recorded = as.numeric(tapply(dat$age, f, max)),
    stringsAsFactors = FALSE
  )
  if (any(okt)) {
    out$first_observed <- as.numeric(tapply(dat$age[okt], fo, min))
    out$last_observed <- as.numeric(tapply(dat$age[okt], fo, max))
    out$mean_age <- as.numeric(tapply(dat$age[okt], fo, mean))
    out$trait_mean <- as.numeric(tapply(dat$trait[okt], fo, mean))
  } else {
    out$first_observed <- NA_real_
    out$last_observed <- NA_real_
    out$mean_age <- NA_real_
    out$trait_mean <- NA_real_
  }
  out$alr <- num_mean(dat$alr)
  out$lifespan <- num_mean(dat$life)
  out$entry <- num_mean(dat$entry)
  out$condition <- num_mean(dat$condition)
  out$condition_label <- as.character(tapply(dat$condition_label, f, mode_or_na))
  out$group <- as.character(tapply(dat$group, f, function(g) {
    g <- g[!is.na(g)]
    if (length(g)) g[[1]] else NA_character_
  }))
  out
}

id_age_means <- function(dat) {
  d <- dat[is.finite(dat$trait) & is.finite(dat$age), c("id", "age", "trait"), drop = FALSE]
  if (!nrow(d)) return(data.frame(id = character(0), age = numeric(0), trait = numeric(0)))
  z <- stats::aggregate(trait ~ id + age, data = d, FUN = mean)
  z[order(z$id, z$age), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Data integrity diagnostics (new in 0.7.0)
# ---------------------------------------------------------------------------
# Starting suggestion only: a Poisson GLM on a quadratic in age gives a Pearson dispersion
# and the ratio of observed to expected zeros. Individual heterogeneity also inflates
# both, so confirm with the error-family check (GLMM AICs) on the Model comparison tab.
suggest_family <- function(trait, age = NULL) {
  ok <- is.finite(trait) & (if (is.null(age)) TRUE else is.finite(age))
  v <- trait[ok]
  if (length(v) < 10) return(list(family = "gaussian", kind = "Unknown", text = "Too few trait values to assess the distribution."))
  is_count <- all(v >= 0) && all(abs(v - round(v)) < 1e-8)
  # Proportions (values within 0-1 that are not all 0 or 1): binomial with a trials column.
  if (!is_count && all(v >= 0) && all(v <= 1)) {
    return(list(family = "binomial", kind = "Proportion (0-1)", zeros = mean(v == 0), dispersion = NA_real_,
                text = "Proportion trait (values between 0 and 1): binomial mixed models. Choose the number of trials (for example clutch size or the number of eggs laid) under 'Weights' below the error family on the Modelling tab, so that each proportion is weighted by its sample size; without it the proportions cannot be fitted as binomial and the Gaussian option is the fallback. If the proportions are more variable than binomial sampling allows, use the beta-binomial family."))
  }
  if (!is_count) return(list(family = "gaussian", kind = "Continuous", text = "Continuous trait: Gaussian mixed models (lme4)."))
  # A binary (0/1) trait is not a count: the count families would be wrong, and there is no binomial
  # family here, so Gaussian is offered as a linear-probability approximation.
  if (all(v %in% c(0, 1))) {
    return(list(family = "binomial", kind = "Binary (0/1)", zeros = mean(v == 0), dispersion = NA_real_,
                text = sprintf("Binary trait (only 0 and 1, %.0f%% zeros): binomial mixed models with a logit link (glmmTMB). The count families do not apply. The Gaussian option fits a linear-probability model instead: it does not keep predictions within 0 and 1 and, because the variance of a binary trait changes with the predicted probability, it can favour the interaction models (4, 6) even when there is no selective disappearance at all, so prefer the binomial family here.",
                               100 * mean(v == 0))))
  }
  zeros <- mean(v == 0)
  disp <- NA_real_
  zero_ratio <- NA_real_
  if (!is.null(age) && length(unique(age[ok])) >= 3 && mean(v) > 0) {
    a <- age[ok]
    g <- tryCatch(suppressWarnings(stats::glm(v ~ a + I(a^2), family = stats::poisson())), error = function(e) NULL)
    if (!is.null(g)) {
      mu <- stats::fitted(g)
      disp <- sum((v - mu)^2 / mu) / max(1, length(v) - 3)
      expected_zeros <- mean(stats::dpois(0, mu))
      zero_ratio <- if (expected_zeros > 0) zeros / expected_zeros else Inf
    }
  }
  if (!is.finite(disp) && mean(v) > 0) disp <- stats::var(v) / mean(v)
  # Underdispersed counts (clutch or litter size, for example) are not Poisson: Poisson would
  # overstate the residual variance, so Gaussian is the better starting point.
  under <- isTRUE(disp < 0.7)
  fam <- if (under) "gaussian" else if (isTRUE(zero_ratio > 1.5) && zeros > 0.1) "zinb" else if (isTRUE(disp > 2)) "nbinom2" else "poisson"
  zr_txt <- if (is.finite(zero_ratio)) {
    if (zero_ratio > 999) " (more than 999\u00d7 the Poisson expectation)" else sprintf(" (%.1f\u00d7 the Poisson expectation)", zero_ratio)
  } else ""
  txt <- sprintf("Non-negative integer trait: %.0f%% zeros%s, Pearson dispersion %.1f.%s Starting suggestion: %s; confirm with the error-family check.",
                 100 * zeros, zr_txt, disp,
                 if (under) " The counts are underdispersed (dispersion well below 1), so a Poisson model would overstate the error." else "",
                 family_label(fam))
  list(family = fam, kind = if (under) "Counts (underdispersed)" else "Counts", text = txt, zeros = zeros, dispersion = disp)
}

family_label <- function(f) {
  switch(f, gaussian = "Gaussian (lme4)", poisson = "Poisson", nbinom2 = "negative binomial (nbinom2)",
         nbinom1 = "negative binomial (nbinom1)", zip = "zero-inflated Poisson",
         zinb = "zero-inflated negative binomial (nbinom2)", zinb1 = "zero-inflated negative binomial (nbinom1)", f)
}

data_integrity <- function(dat, meta) {
  rows <- list()
  add <- function(check, result, status, advice) {
    rows[[length(rows) + 1]] <<- data.frame(Check = check, Result = result, Status = status, Advice = advice,
                                            stringsAsFactors = FALSE)
  }
  if (!nrow(dat)) return(list(table = data.frame(), family = suggest_family(numeric(0))))
  im <- individual_metrics(dat)
  step <- infer_age_step(dat$age, dat$id)
  schedule <- list(irregular = FALSE, n_schedules = 1, step = step, off_share = 0)
  add("Rows used", sprintf("%d of %d rows", meta$n_raw - meta$n_dropped, meta$n_raw),
      if (meta$n_dropped > 0) "Note" else "OK",
      if (meta$n_dropped > 0) "Rows without an ID or a numeric age were removed." else "")
  n_na <- sum(!is.finite(dat$trait))
  add("Missing trait values", sprintf("%d rows", n_na), if (n_na > 0) "Note" else "OK",
      if (n_na > 0) "Rows without a trait value count as recorded ages (for ALR/AFR) but are not modelled." else "")
  add("Duplicate ID \u00d7 age records", sprintf("%d rows", meta$n_dup), if (meta$n_dup > 0) "Warning" else "OK",
      if (meta$n_dup > 0) paste0("Kept as separate observations (", meta$dup_action,
                                 "). Check for ID or age typos; choose 'average duplicates' if they are technical replicates. If the same ID is reused in different groups (families, vials, sites), do not average: map the grouping column and keep nesting, so that group + ID identifies the individual.") else "")
  n_single <- sum(im$n_trait == 1)
  add("Individuals with one trait record", sprintf("%d of %d (%.0f%%)", n_single, nrow(im), 100 * n_single / nrow(im)),
      if (n_single / nrow(im) > 0.25) "Note" else "OK",
      "They inform among-individual terms but contribute no within-individual change or decomposition pairs.")
  if (isTRUE(meta$alr_mapped)) {
    mism <- sum(abs(im$alr - im$last_recorded) > 1e-8, na.rm = TRUE)
    add("Mapped ALR vs last recorded age", sprintf("%d individuals differ", mism), if (mism > 0) "Warning" else "OK",
        if (mism > 0) "ALR should equal the last recorded age; check the column or use automatic calculation." else "")
  }
  if (isTRUE(meta$has_censor)) {
    add("Censored individuals", sprintf("%d individuals (%s = %s)", meta$n_censored, meta$map$censor, meta$map$censor_value),
        if (meta$n_censored > 0) "Note" else "OK",
        "LS is set to missing for censored individuals: their recorded lifespan is not an age at death.")
  }
  if (isTRUE(meta$has_life) && !isTRUE(meta$life_auto)) {
    n_ls_na <- sum(!is.finite(im$lifespan))
    n_ls_bad <- sum(is.finite(im$lifespan) & im$lifespan < im$last_recorded - 1e-8)
    add("Lifespan (LS) missing", sprintf("%d individuals", n_ls_na), if (n_ls_na > 0) "Note" else "OK",
        if (n_ls_na > 0) "Model 6 drops these individuals; by default Model 6 is then left unselected so Models 1-5 use all data." else "")
    add("LS earlier than last recorded age", sprintf("%d individuals", n_ls_bad), if (n_ls_bad > 0) "Warning" else "OK",
        if (n_ls_bad > 0) "Impossible values: check lifespan or age records for these IDs (listed below)." else "")
    if (is.finite(step)) {
      gap <- (im$lifespan - im$last_recorded) / step
      gap <- gap[is.finite(gap)]
      if (length(gap)) {
        add("LS \u2212 ALR (time steps)", sprintf("median %.2f; %.0f%% \u2265 1 step", stats::median(gap), 100 * mean(gap >= 1 - 1e-8)),
            if (mean(gap >= 1 - 1e-8) > 0.1) "Note" else "OK",
            "Gaps of \u2265 1 step mean an individual was alive at a sampling occasion without a record (terminal missingness), or that ages and LS use different day conventions. See the lifespan margin on the Sampling tab.")
      }
    }
  }
  if (isTRUE(meta$has_group)) {
    nbad <- meta$n_multi_group %||% 0L
    add("IDs linked to >1 higher-level group", sprintf("%d IDs", nbad), if (nbad > 0) "Warning" else "OK",
        if (nbad > 0) paste("Possible typo in the grouping column.", if (isTRUE(meta$nested)) "With nesting, each group/ID combination is treated as a separate individual." else "Without nesting these rows share one individual.") else "Nested: individuals are identified by group + ID.")
  }
  if (isTRUE(meta$has_group2)) {
    nb2 <- meta$n_multi_group2 %||% 0L
    add("Groups linked to >1 top-level group", sprintf("%d groups", nb2), if (nb2 > 0) (if (isTRUE(meta$nested)) "Note" else "Warning") else "OK",
        if (nb2 > 0) paste("The same group label occurs in more than one top-level group.", if (isTRUE(meta$nested)) "With nesting, each top-level group/group combination is a separate group (e.g. father 1 of family A and father 1 of family B are different fathers)." else "Without nesting these rows share one group-level random intercept.") else "Three levels: individuals within groups within top-level groups.")
  }
  if (is.finite(step)) {
    fa <- stats::ave(dat$age, dat$id, FUN = min)
    rel <- (dat$age - fa) / step
    off <- abs(rel - round(rel)) > 0.01
    first_ages <- fa[!duplicated(dat$id)]
    phase <- round(((first_ages - min(dat$age)) / step) %% 1, 2)
    phase[phase >= 1] <- 0
    n_sched <- length(unique(phase))
    irregular <- mean(off) > 0.2
    schedule <- list(irregular = irregular, n_schedules = n_sched, step = step, off_share = mean(off))
    add("Sampling schedule",
        sprintf("%d distinct ages; step = %s; %d off-schedule rows%s", length(unique(dat$age)), format_num(step), sum(off),
                if (!irregular && n_sched > 1) sprintf("; %d offset schedules", n_sched) else ""),
        if (irregular) "Warning" else if (sum(off) > 0 || n_sched > 1) "Note" else "OK",
        if (irregular) {
          sprintf("Ages are irregular (not on a common sampling schedule), so the missingness grid, per-age means and the decomposition are approximate. Consider 'Round ages to multiples of' on the Data tab (e.g. %s).", format_num(signif(step, 2)))
        } else if (n_sched > 1) {
          "Individuals follow schedules offset from each other (e.g. staggered cohorts): expected occasions are anchored on each individual's first record."
        } else if (sum(off) > 0) {
          "Some ages are not on the schedule: they are assigned to the nearest expected occasion."
        } else "")
  }
  if (length(meta$skipped_columns)) {
    add("Covariates or random effects without values", paste(meta$skipped_columns, collapse = ", "), "Warning",
        "These columns contain no values and were not used.")
  }
  tv <- dat$trait[is.finite(dat$trait)]
  if (length(tv) >= 10) {
    spread <- stats::IQR(tv)
    if (!isTRUE(spread > 0)) spread <- stats::sd(tv)
    if (isTRUE(spread > 0)) {
      n_ext <- sum(abs(tv - stats::median(tv)) > 20 * spread)
      if (n_ext > 0) {
        add("Extreme trait values", sprintf("%d values more than 20 interquartile ranges from the median (maximum %s)", n_ext, format_num(max(abs(tv)))),
            "Note", "Check whether these are genuine (heavy-tailed counts) or data-entry errors or unit changes: a few extreme values can dominate models and plots.")
      }
    }
  }
  n_code <- sum(dat$trait %in% c(-99, -999, -9999)) + sum(dat$age %in% c(-99, -999, -9999)) + sum(dat$life %in% c(-99, -999, -9999))
  if (n_code > 0) {
    add("Possible missing-value codes", sprintf("%d values equal to -99, -999 or -9999", n_code), "Warning",
        "Recode missing-value codes as blank or NA before uploading; otherwise they are analysed as real values.")
  }
  if (any(dat$age <= 0)) {
    add("Ages at or below zero", sprintf("%d rows", sum(dat$age <= 0)), "Note",
        "Allowed. The logarithmic function then uses log(age \u2212 minimum age + 1).")
  }
  for (cv in meta$covars %||% character(0)) {
    x <- dat[[cv]]
    lab <- meta$cov_labels[[cv]]
    if (is.numeric(x)) {
      if (!isTRUE(stats::sd(x, na.rm = TRUE) > 0)) add(paste("Covariate", lab), "no variation", "Warning", "Omitted from the models.")
    } else {
      nl <- length(unique(x[!is.na(x)]))
      if (nl < 2) {
        add(paste("Covariate", lab), sprintf("%d level", nl), "Warning", "A factor needs at least two levels: omitted from the models.")
      } else if (nl > 50) {
        add(paste("Covariate", lab), sprintf("%d levels", nl), "Note", "Many levels: consider a random intercept instead of a fixed effect (slow fits, many parameters).")
      }
    }
    if (any(is.na(x))) add(paste("Covariate", lab, "missing"), sprintf("%d rows", sum(is.na(x))), "Note", "Rows with a missing covariate are dropped from every model.")
  }
  if (isTRUE(meta$has_trials)) {
    nb <- meta$n_bad_trials %||% 0L
    add("Number of binomial trials", sprintf("%s: %d rows without a positive number of trials", meta$trials_col, nb),
        if (nb > 0) "Warning" else "OK",
        if (nb > 0) "Rows without a positive number of trials are dropped from binomial models." else "Used as prior weights in the binomial families.")
  }
  add("Records per individual", sprintf("median %s (range %d\u2013%d)", format_num(stats::median(im$n_trait)), min(im$n_trait), max(im$n_trait)),
      if (stats::median(im$n_trait) < 4) "Note" else "OK",
      if (stats::median(im$n_trait) < 4) "Few records per individual: individual parametric fits will be restricted to long-lived individuals." else "")
  n_afr <- length(unique(im$entry[is.finite(im$entry)]))
  add("AFR variation", sprintf("%d distinct values", n_afr), if (n_afr < 2) "Warning" else "OK",
      if (n_afr < 2) "Fewer than two AFR values: selective appearance cannot be assessed and Models 7-10 are unavailable. Fitted anyway they would repeat Models 1-4, because their AFR terms are dropped as rank deficient." else "")
  if (n_afr >= 2) {
    r_aa <- cor_safe(im$entry, im$alr)
    add("AFR\u2013ALR correlation", if (is.finite(r_aa)) sprintf("r = %.2f", r_aa) else "NA",
        if (isTRUE(abs(r_aa) > 0.3)) "Note" else "OK",
        "Strongly correlated AFR and ALR make appearance (AFR) and disappearance (ALR) terms in Models 7-8 hard to separate.")
  }
  fam <- suggest_family(dat$trait, dat$age)
  add("Trait distribution", fam$text, if (fam$family != "gaussian") "Warning" else "OK",
      if (fam$family != "gaussian") "Model rankings can change with the error family: choose it on the Model comparison tab." else "")
  bad_ids <- if (isTRUE(meta$has_life) && !isTRUE(meta$life_auto)) im$id[is.finite(im$lifespan) & im$lifespan < im$last_recorded - 1e-8] else character(0)
  list(table = do.call(rbind, rows), family = fam, bad_ls_ids = bad_ids, schedule = schedule)
}

# ---------------------------------------------------------------------------
# Expected sampling grid and missingness (B3 inputs)
# ---------------------------------------------------------------------------
# Each individual's expected window starts at its mapped/automatic AFR and ends at
# the last scheduled occasion at or before (LS - margin) when LS is known, otherwise at
# its last record. Occasions after death cannot be sampled, so the end uses floor().
# start_mode = "afr": expected occasions start at each individual's AFR (mapped or first record).
# start_mode = "same": trait expression starts at start_age for everyone, so unrecorded occasions
# between start_age and an individual's first record count as missing.
build_missing_grid <- function(dat, margin = 0, start_mode = "afr", start_age = NA_real_, max_cells = 300000) {
  d <- dat[is.finite(dat$age), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  step0 <- infer_age_step(d$age, d$id)
  if (!is.finite(step0) || step0 <= 0) step0 <- 1
  ids <- unique(d$id)
  f <- factor(d$id, levels = ids)
  fi <- as.integer(f)
  first_age <- as.numeric(tapply(d$age, f, min))
  entry <- as.numeric(tapply(d$entry, f, function(v) if (any(is.finite(v))) min(v[is.finite(v)]) else NA_real_))
  life <- as.numeric(tapply(d$life, f, finite_mean))
  amax <- max(d$age)
  common_start <- identical(start_mode, "same") && isTRUE(is.finite(start_age))
  # An age at first trait expression (AFE), where supplied, opens the window earlier than entry into the
  # dataset: occasions between expression and entry then count as missed rather than as "not yet expressing".
  start_ref <- if (common_start) rep(start_age, length(ids)) else ifelse(is.finite(entry), entry, first_age)
  build <- function(step) {
    # Each individual's schedule is anchored on its first record, so individuals sampled on schedules
    # offset from each other (staggered cohorts) are not given spurious missed occasions.
    n_before <- pmax(0, floor((first_age - start_ref) / step + 1e-9))
    anchor <- first_age - n_before * step
    k_rec <- round((d$age - anchor[fi]) / step)
    last_k <- as.numeric(tapply(k_rec, f, max))
    end_k <- ifelse(is.finite(life), floor((life - margin - anchor) / step + 1e-9), last_k)
    end_k <- pmin(end_k, floor((amax - anchor) / step + 1e-9))
    end_k <- pmax(end_k, last_k)
    end_k[!is.finite(end_k)] <- last_k[!is.finite(end_k)]
    start_k <- pmax(0, end_k - 4999)
    list(anchor = anchor, k_rec = k_rec, start_k = start_k, end_k = end_k,
         n_cells = as.integer(pmax(1, end_k - start_k + 1)))
  }
  step <- step0
  b <- build(step)
  coarsen <- 1
  while (sum(as.numeric(b$n_cells)) > max_cells && coarsen < 1e6) {
    coarsen <- coarsen * max(2, ceiling(sum(as.numeric(b$n_cells)) / max_cells))
    step <- step0 * coarsen
    b <- build(step)
  }

  gi <- rep(seq_along(ids), b$n_cells)
  grid <- data.frame(id = ids[gi], stringsAsFactors = FALSE)
  grid$k <- rep(b$start_k, b$n_cells) + sequence(b$n_cells) - 1
  grid$age <- b$anchor[gi] + grid$k * step

  okt <- is.finite(d$trait)
  dkey <- paste(d$id, b$k_rec, sep = "\r")
  gkey <- paste(grid$id, grid$k, sep = "\r")
  grid$trait <- NA_real_
  if (any(okt)) {
    cell_mean <- tapply(d$trait[okt], dkey[okt], mean)
    grid$trait <- as.numeric(cell_mean[gkey])
  }
  grid$trait_observed <- is.finite(grid$trait)
  grid$missing <- !grid$trait_observed

  fo <- rep(NA_real_, length(ids))
  lo <- rep(NA_real_, length(ids))
  if (any(okt)) {
    fok <- factor(d$id[okt], levels = ids)
    fo <- as.numeric(tapply(b$k_rec[okt], fok, min))
    lo <- as.numeric(tapply(b$k_rec[okt], fok, max))
  }
  grid$pattern <- ifelse(!grid$missing, "Observed",
                  ifelse(!is.finite(fo[gi]), "Missed: no trait records",
                  ifelse(grid$k < fo[gi], "Missed: before first record",
                  ifelse(grid$k > lo[gi], "Missed: after last record", "Missed: between records"))))

  im <- individual_metrics(d)
  mm <- match(grid$id, im$id)
  grid$life <- im$lifespan[mm]
  grid$entry <- im$entry[mm]
  grid$alr <- im$alr[mm]
  grid$condition <- im$condition[mm]
  grid$condition_label <- im$condition_label[mm]

  n <- nrow(grid)
  tr <- grid$trait
  loi <- cummax(ifelse(is.finite(tr), seq_len(n), 0))
  prev <- c(0, loi[-n])
  okp <- prev > 0
  if (any(okp)) okp[okp] <- grid$id[prev[okp]] == grid$id[okp]
  grid$prior_trait <- NA_real_
  grid$prior_trait[okp] <- tr[prev[okp]]
  grid$anchor <- b$anchor[gi]
  grid$step <- step
  attr(grid, "step_note") <- if (coarsen > 1) {
    sprintf("The sampling step was coarsened from %s to %s to keep the grid below %s cells; consider rounding ages on the Data tab.",
            format_num(step0), format_num(step), format(max_cells, big.mark = ","))
  } else ""
  grid
}

missingness_drivers <- function(grid, alpha = 0.05) {
  empty <- data.frame(Predictor = character(0), Effect = character(0), P_value = numeric(0), Signal = character(0))
  if (!nrow(grid) || length(unique(grid$missing)) < 2 || sum(grid$missing) < 10) return(empty)
  specs <- c(Age = "age", LS = "life", AFR = "entry", `Prior observed trait` = "prior_trait", Condition = "condition")
  ans <- lapply(names(specs), function(lbl) {
    x <- grid[[specs[[lbl]]]]
    ok <- is.finite(x)
    if (sum(ok) < 20 || length(unique(x[ok])) < 2) return(NULL)
    z <- data.frame(missing = as.integer(grid$missing[ok]), x = as.numeric(scale(x[ok])))
    if (length(unique(z$missing)) < 2) return(NULL)
    fit <- tryCatch(suppressWarnings(stats::glm(missing ~ x, data = z, family = stats::binomial())), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)$coefficients
    if (!"x" %in% rownames(sm)) return(NULL)
    p <- sm["x", 4]
    data.frame(Predictor = lbl, Effect = sprintf("OR %.2f per SD", exp(sm["x", 1])), P_value = p,
               Signal = if (is.finite(p) && p < alpha) "Associated" else "Not detected", stringsAsFactors = FALSE)
  })
  ans <- Filter(Negate(is.null), ans)
  if (!length(ans)) return(empty)
  do.call(rbind, ans)
}

missingness_summary <- function(dat, grid, alpha = 0.05, life_known = FALSE) {
  im <- individual_metrics(dat)
  out <- list(percent = NA_real_, p = NA_real_, step = NA_real_,
              r_mean_alr = cor_safe(im$mean_age, im$alr), r_alr_ls = NA_real_, r_mean_ls = NA_real_,
              gap_steps = numeric(0), pattern = "Unknown", drivers = missingness_drivers(grid, alpha),
              guidance = "Insufficient information.", type = "Could not construct the sampling grid.")
  if (!nrow(grid)) return(out)
  out$step <- grid$step[[1]]
  out$step_note <- attr(grid, "step_note") %||% ""
  out$percent <- 100 * mean(grid$missing)
  out$p <- 1 - mean(grid$missing)
  if (isTRUE(life_known)) {
    out$r_alr_ls <- cor_safe(im$alr, im$lifespan)
    out$r_mean_ls <- cor_safe(im$mean_age, im$lifespan)
    g <- (im$lifespan - im$alr) / out$step
    out$gap_steps <- g[is.finite(g)]
  }
  tab <- table(grid$pattern[grid$missing])
  out$pattern <- if (!length(tab)) "No missing cells" else
    paste(paste0(names(tab), ": ", as.integer(tab), " (", round(100 * as.numeric(tab) / sum(tab)), "%)"), collapse = "; ")
  drv <- out$drivers
  sig <- if (nrow(drv)) drv$Predictor[is.finite(drv$P_value) & drv$P_value < alpha] else character(0)
  out$type <- if (out$percent < 1) {
    "Sampling is effectively complete within the expected windows."
  } else if (!length(sig)) {
    "No association between missingness and the tested observed variables (compatible with, but not proof of, MCAR)."
  } else {
    paste0("Missingness is associated with: ", paste(sig, collapse = ", "), ".")
  }
  rtxt <- if (is.finite(out$r_mean_alr)) sprintf(" (r mean age\u2013ALR = %.3f)", out$r_mean_alr) else ""
  ls_txt <- if (is.finite(out$r_alr_ls) && is.finite(out$r_mean_ls)) {
    sprintf(" With known lifespan: r(ALR, LS) = %.2f, r(mean age, LS) = %.2f.", out$r_alr_ls, out$r_mean_ls)
  } else ""
  out$guidance <- if (out$p >= 0.95) {
    paste0("Detection is near complete (p = ", sprintf("%.2f", out$p), "): ALR and mean age are near-equivalent proxies", rtxt,
           ". Expect Models 4 and 5 to receive similar support; differences then mainly reflect their different design spaces, not proxy quality.", ls_txt)
  } else {
    paste0("Detection is incomplete (p = ", sprintf("%.2f", out$p), "): intermittent gaps shift mean age but not ALR, so ALR is usually the less error-prone lifespan proxy", rtxt,
           ". In the manuscript's MCAR, MWO and MWY simulations Model 4 outperformed Model 5; compare both rather than defaulting to mean-age centring.", ls_txt)
  }
  if (!isTRUE(life_known)) {
    out$guidance <- paste(out$guidance, "Without known LS the window ends at each individual's last record, so terminal missingness is invisible and p is an upper bound.")
  }
  if (nzchar(out$step_note %||% "")) out$type <- paste(out$type, out$step_note)
  out
}

missingness_by_variable <- function(dat, grid, variable = "Age") {
  if (!nrow(grid)) return(data.frame())
  if (identical(variable, "Age")) {
    z <- stats::aggregate(as.integer(grid$missing), list(value = grid$age), mean)
    names(z)[2] <- "missing_rate"
    z$n <- as.numeric(table(grid$age)[as.character(z$value)])
    z$type <- "continuous"
    return(z)
  }
  im <- individual_metrics(dat)
  val <- switch(variable, LS = im$lifespan, ALR = im$alr, Trait = im$trait_mean,
                Condition = if (any(is.finite(im$condition))) im$condition else im$condition_label, im$mean_age)
  if (is.null(val) || length(val) != nrow(im)) val <- rep(NA_real_, nrow(im))
  map <- stats::setNames(val, im$id)
  x <- grid
  x$x <- unname(map[x$id])
  if (is.numeric(x$x)) {
    ok <- is.finite(x$x)
    if (!any(ok)) return(data.frame())
    if (length(unique(x$x[ok])) > 8) {
      q <- unique(stats::quantile(x$x[ok], probs = seq(0, 1, length.out = 7), type = 7))
      if (length(q) >= 3) {
        x <- x[ok, , drop = FALSE]
        x$group <- cut(x$x, breaks = q, include.lowest = TRUE)
        z <- stats::aggregate(as.integer(x$missing), list(group = x$group), mean)
        names(z)[2] <- "missing_rate"
        mids <- tapply(x$x, x$group, mean)
        z$value <- as.numeric(mids[as.character(z$group)])
        z$type <- "binned"
        return(z[order(z$value), , drop = FALSE])
      }
    }
    x <- x[ok, , drop = FALSE]
    z <- stats::aggregate(as.integer(x$missing), list(value = x$x), mean)
    names(z)[2] <- "missing_rate"
    z$type <- "continuous"
    return(z)
  }
  x$group <- as.character(x$x)
  x$group[is.na(x$group) | !nzchar(x$group)] <- "Unknown"
  z <- stats::aggregate(as.integer(x$missing), list(group = x$group), mean)
  names(z)[2] <- "missing_rate"
  z$type <- "categorical"
  z
}

# Sampling grid for display: every shown individual x expected occasion across the displayed age range
# (on that individual's own schedule), classified as observed, missed (expected but unrecorded), final
# expected occasion (death if LS is known, else ALR), or not expected (before AFR or after death / ALR).
# The number of tiles is capped so that the plot stays drawable; attr "n_shown" gives the individuals kept.
grid_display <- function(grid, ids, max_tiles = 60000) {
  g <- grid[grid$id %in% ids, , drop = FALSE]
  if (!nrow(g)) return(data.frame())
  step <- g$step[[1]]
  ids <- ids[ids %in% g$id]
  anc <- tapply(g$anchor, g$id, function(v) v[[1]])
  a_i <- as.numeric(anc[ids])
  amin <- min(g$age)
  amax <- max(g$age)
  kmin <- ceiling((amin - a_i) / step - 1e-9)
  kmax <- floor((amax - a_i) / step + 1e-9)
  n <- as.integer(pmax(1, kmax - kmin + 1))
  if (sum(as.numeric(n)) > max_tiles) {
    keep <- unique(round(seq(1, length(ids), length.out = max(2, floor(max_tiles / mean(n))))))
    ids <- ids[keep]; a_i <- a_i[keep]; kmin <- kmin[keep]; n <- n[keep]
  }
  full <- data.frame(id = rep(ids, n), stringsAsFactors = FALSE)
  full$k <- rep(kmin, n) + sequence(n) - 1
  full$age <- rep(a_i, n) + full$k * step
  m <- match(paste(full$id, full$k, sep = "\r"), paste(g$id, g$k, sep = "\r"))
  status <- ifelse(is.na(m), "Not expected", ifelse(g$missing[m], "Missed", "Observed"))
  last_k <- tapply(g$k, g$id, max)
  status[!is.na(m) & full$k == as.numeric(last_k[full$id])] <- "Death or ALR"
  full$status <- status
  attr(full, "n_shown") <- length(ids)
  full
}

coverage_by_age <- function(grid) {
  if (!nrow(grid)) return(data.frame())
  ages <- sort(unique(grid$age))
  n_obs <- as.numeric(tapply(grid$trait_observed, grid$age, sum)[as.character(ages)])
  n_exp <- as.numeric(table(grid$age)[as.character(ages)])
  data.frame(Age = ages, Expected = n_exp, Observed = n_obs, Percent_observed = round(100 * n_obs / n_exp, 1))
}

# ---------------------------------------------------------------------------
# Proxies: individual-level values used in A1/A2
# ---------------------------------------------------------------------------
proxy_choices <- function(im, meta, target = "disappearance") {
  if (identical(target, "all")) {
    ch <- c("ALR", "Mean age")
    if (isTRUE(meta$has_life) && !isTRUE(meta$life_auto) && sum(is.finite(im$lifespan)) >= 4) ch <- c(ch, "LS")
    if (length(unique(im$entry[is.finite(im$entry)])) >= 2) ch <- c(ch, "AFR")
    return(ch)
  }
  if (identical(target, "appearance")) {
    if (length(unique(im$entry[is.finite(im$entry)])) < 2) return(character(0))
    return(c("AFR"))
  }
  ch <- c("ALR", "Mean age")
  if (isTRUE(meta$has_life) && !isTRUE(meta$life_auto) && sum(is.finite(im$lifespan)) >= 4) ch <- c(ch, "LS")
  ch
}

proxy_values <- function(im, proxy) {
  v <- switch(proxy, ALR = im$alr, `Mean age` = im$mean_age, LS = im$lifespan, AFR = im$entry, im$alr)
  if (is.null(v) || length(v) != nrow(im)) v <- rep(NA_real_, nrow(im))
  stats::setNames(v, im$id)
}

# ---------------------------------------------------------------------------
# A1: trait trajectories within proxy bins (individual-level binning)
# ---------------------------------------------------------------------------
# With no more distinct proxy values than requested bins, every distinct value is its own bin.
make_proxy_bins <- function(values, n_bins = 4, method = "equal") {
  x <- values[is.finite(values)]
  n_bins <- max(2L, as.integer(n_bins))
  ux <- sort(unique(x))
  if (length(ux) < 2) return(NULL)
  if (length(ux) <= n_bins) {
    labs <- make.unique(paste0(seq_along(ux), ": ", format_num(ux)))
    return(factor(labs[match(values, ux)], levels = labs))
  }
  br <- if (identical(method, "quantile")) {
    unique(as.numeric(stats::quantile(x, probs = seq(0, 1, length.out = n_bins + 1), type = 7)))
  } else {
    seq(min(x), max(x), length.out = n_bins + 1)
  }
  if (length(br) < 3) br <- seq(min(x), max(x), length.out = n_bins + 1)
  if (length(unique(br)) < 3) return(NULL)
  labs <- make.unique(paste0(seq_len(length(br) - 1), ": ", format_num(br[-length(br)]), "\u2013", format_num(br[-1])))
  cut(values, breaks = br, include.lowest = TRUE, labels = labs)
}

# Upper limit for the number of bins: distinct values (snapped to the sampling step) shared by at
# least n_min individuals. Individuals with a single record count at their only age.
max_bins_for <- function(values, n_min = 3, step = NA_real_, floor_bins = 3L) {
  v <- values[is.finite(values)]
  if (!length(v)) return(as.integer(floor_bins))
  if (is.finite(step) && step > 0) v <- round(v / step) * step
  tab <- table(v)
  as.integer(max(floor_bins, sum(tab >= max(1, n_min))))
}

# facet: optional named vector (individual id -> panel label); NULL puts everyone in one panel.
binned_trajectory <- function(dat, pvals, n_bins = 4, method = "equal", n_min = 3, facet = NULL) {
  b <- make_proxy_bins(pvals, n_bins, method)
  if (is.null(b)) return(data.frame())
  bmap <- stats::setNames(as.character(b), names(pvals))
  ia <- id_age_means(dat)
  if (!nrow(ia)) return(data.frame())
  ia$bin <- unname(bmap[ia$id])
  ia$facet <- if (is.null(facet)) "All" else unname(facet[ia$id])
  ia <- ia[!is.na(ia$bin) & !is.na(ia$facet), , drop = FALSE]
  if (!nrow(ia)) return(data.frame())
  key <- paste(ia$facet, ia$bin, ia$age, sep = "\r")
  first <- !duplicated(key)
  s <- ia[first, c("facet", "bin", "age"), drop = FALSE]
  kk <- key[first]
  s$mean <- as.numeric(tapply(ia$trait, key, mean)[kk])
  s$n <- as.numeric(tapply(ia$trait, key, length)[kk])
  s$se <- as.numeric(tapply(ia$trait, key, function(v) if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_)[kk])
  s$bin <- factor(s$bin, levels = levels(b))
  s <- s[s$n >= n_min, , drop = FALSE]
  s[order(s$facet, s$bin, s$age), , drop = FALSE]
}

# Differences in mean trait between proxy bins at each age (input: binned_trajectory output), per facet.
# mode = "successive": each bin minus the bin below it; mode = "all": every pair of bins.
# Differences are signed (higher bin minus lower bin) so that crossing trajectories stay visible.
bin_differences <- function(s, mode = "successive") {
  empty <- data.frame(age = numeric(0), difference = numeric(0), pair = character(0), facet = character(0))
  if (!nrow(s)) return(empty)
  if (!"facet" %in% names(s)) s$facet <- "All"
  lev <- levels(s$bin)
  if (length(lev) < 2) return(empty)
  prs <- if (identical(mode, "all")) {
    utils::combn(seq_along(lev), 2, simplify = FALSE)
  } else {
    lapply(seq_len(length(lev) - 1), function(j) c(j, j + 1))
  }
  num <- sub(":.*$", "", lev)
  pair_levels <- vapply(prs, function(pr) paste0("Bin ", num[[pr[2]]], " \u2212 bin ", num[[pr[1]]]), character(1))
  rows <- list()
  for (fv in unique(as.character(s$facet))) {
    sf <- s[as.character(s$facet) == fv, , drop = FALSE]
    mean_of <- stats::setNames(sf$mean, paste(as.character(sf$bin), sf$age, sep = "\r"))
    ages <- sort(unique(sf$age))
    for (j in seq_along(prs)) {
      pr <- prs[[j]]
      a_lo <- mean_of[paste(lev[[pr[1]]], ages, sep = "\r")]
      a_hi <- mean_of[paste(lev[[pr[2]]], ages, sep = "\r")]
      ok <- !is.na(a_lo) & !is.na(a_hi)
      if (!any(ok)) next
      rows[[length(rows) + 1]] <- data.frame(age = ages[ok], difference = as.numeric(a_hi[ok] - a_lo[ok]),
                                             pair = pair_levels[[j]], facet = fv, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  out$pair <- factor(out$pair, levels = intersect(pair_levels, unique(out$pair)))
  out
}

# Least-squares trend of each bin difference against age (lines in the A1 difference plot).
bin_difference_trends <- function(dz) {
  if (!nrow(dz)) return(data.frame())
  if (!"facet" %in% names(dz)) dz$facet <- "All"
  grp <- paste(dz$facet, as.character(dz$pair), sep = "\r")
  rows <- lapply(split(dz, grp), function(x) {
    if (length(unique(x$age)) < 2) return(NULL)
    fit <- stats::lm(difference ~ age, data = x)
    cf <- suppressWarnings(summary(fit)$coefficients)
    if (!"age" %in% rownames(cf)) return(NULL)
    p <- if (nrow(x) > 2 && ncol(cf) >= 4) cf["age", 4] else NA_real_
    data.frame(Facet = x$facet[[1]], Pair = as.character(x$pair[[1]]), N_ages = nrow(x),
               Age_from = min(x$age), Age_to = max(x$age), Intercept = cf[1, 1], Slope_per_age = cf["age", 1],
               P = p, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$Pair <- factor(out$Pair, levels = levels(dz$pair))
  out <- out[order(out$Facet, out$Pair), , drop = FALSE]
  out$Pair <- as.character(out$Pair)
  rownames(out) <- NULL
  out
}

# A2: each individual's mean trait within age bins, against an individual-level proxy such as LS.
# With no more distinct ages than bins, every age is its own bin. age_mid = mean sampled age in the bin.
trait_by_age_bins <- function(dat, pvals, n_age_bins = 4, facet = NULL) {
  ia <- id_age_means(dat)
  if (!nrow(ia)) return(data.frame())
  ia$proxy <- unname(pvals[ia$id])
  ia$facet <- if (is.null(facet)) "All" else unname(facet[ia$id])
  ia <- ia[is.finite(ia$trait) & is.finite(ia$proxy) & !is.na(ia$facet), , drop = FALSE]
  if (!nrow(ia)) return(data.frame())
  ua <- sort(unique(ia$age))
  n_age_bins <- max(2L, as.integer(n_age_bins))
  if (length(ua) <= n_age_bins) {
    ia$age_bin <- factor(as.character(ia$age), levels = as.character(ua))
  } else {
    br <- seq(min(ua), max(ua), length.out = n_age_bins + 1)
    labs <- make.unique(paste0(format_num(br[-length(br)]), "\u2013", format_num(br[-1])))
    ia$age_bin <- cut(ia$age, breaks = br, include.lowest = TRUE, labels = labs)
  }
  z <- stats::aggregate(trait ~ id + facet + age_bin, data = ia, FUN = mean)
  z$age_bin <- factor(as.character(z$age_bin), levels = levels(ia$age_bin))
  z$proxy <- unname(pvals[as.character(z$id)])
  mids <- tapply(ia$age, ia$age_bin, mean)
  z$age_mid <- as.numeric(mids[as.character(z$age_bin)])
  z
}

# A2 coefficient table: regression coefficient of the trait on the proxy within each age bin (per
# facet), its 95% CI, and the change from the previous age bin.
a2_bin_slopes <- function(z) {
  if (!nrow(z)) return(data.frame())
  if (!"facet" %in% names(z)) z$facet <- "All"
  lev <- levels(z$age_bin)
  rows <- list()
  for (fv in unique(as.character(z$facet))) {
    for (ab in lev) {
      x <- z[as.character(z$facet) == fv & as.character(z$age_bin) == ab, , drop = FALSE]
      if (!nrow(x)) next
      est <- se <- lo <- hi <- r <- p <- NA_real_
      if (nrow(x) >= 3 && length(unique(x$proxy)) >= 2) {
        fit <- stats::lm(trait ~ proxy, data = x)
        cf <- suppressWarnings(summary(fit)$coefficients)
        if ("proxy" %in% rownames(cf)) {
          est <- cf["proxy", 1]
          se <- cf["proxy", 2]
          p <- if (ncol(cf) >= 4) cf["proxy", 4] else NA_real_
          tq <- stats::qt(0.975, df = nrow(x) - 2)
          lo <- est - tq * se
          hi <- est + tq * se
        }
        r <- suppressWarnings(stats::cor(x$proxy, x$trait))
      }
      rows[[length(rows) + 1]] <- data.frame(Facet = fv, Age_bin = ab, Mean_age = x$age_mid[[1]], N = nrow(x),
                                             Coefficient = est, SE = se, Lower_95 = lo, Upper_95 = hi, r = r, P = p,
                                             stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$Change_from_previous <- stats::ave(out$Coefficient, out$Facet, FUN = function(v) c(NA_real_, diff(v)))
  out
}

# Weighted trend of the A2 coefficients across age bins (per facet): change per unit age.
a2_slope_trend <- function(tab) {
  if (!nrow(tab)) return(data.frame())
  rows <- lapply(split(tab, tab$Facet), function(x) {
    x <- x[is.finite(x$Coefficient) & is.finite(x$SE) & x$SE > 0, , drop = FALSE]
    if (nrow(x) < 3) return(NULL)
    fit <- stats::lm(Coefficient ~ Mean_age, data = x, weights = 1 / x$SE^2)
    cf <- suppressWarnings(summary(fit)$coefficients)
    if (!"Mean_age" %in% rownames(cf)) return(NULL)
    data.frame(Facet = x$Facet[[1]], Change_per_age = cf["Mean_age", 1], P = cf["Mean_age", 4], N_bins = nrow(x),
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# A2 (summary statistic): association between the proxy and the trait at each age
# ---------------------------------------------------------------------------
proxy_slopes_by_age <- function(dat, pvals, metric = "slope", max_ages = 30, min_n = 5) {
  empty <- list(table = data.frame(), trend_est = NA_real_, trend_p = NA_real_, mean_est = NA_real_,
                mean_p = NA_real_, metric = metric)
  ia <- id_age_means(dat)
  if (!nrow(ia)) return(empty)
  ia$proxy <- unname(pvals[ia$id])
  ia <- ia[is.finite(ia$trait) & is.finite(ia$proxy), , drop = FALSE]
  if (!nrow(ia)) return(empty)
  if (length(unique(ia$age)) > max_ages) {
    br <- unique(as.numeric(stats::quantile(ia$age, probs = seq(0, 1, length.out = 9))))
    grp <- cut(ia$age, breaks = br, include.lowest = TRUE)
    ia$age_grp <- as.numeric(stats::ave(ia$age, grp, FUN = mean))
  } else {
    ia$age_grp <- ia$age
  }
  res <- lapply(split(ia, ia$age_grp), function(x) {
    n <- nrow(x)
    if (n < min_n || length(unique(x$proxy)) < 2 || length(unique(x$trait)) < 2) return(NULL)
    if (identical(metric, "r")) {
      r <- stats::cor(x$proxy, x$trait)
      if (!is.finite(r) || n < 4) return(NULL)
      zr <- atanh(max(min(r, 0.999999), -0.999999))
      se <- 1 / sqrt(n - 3)
      data.frame(age = x$age_grp[[1]], estimate = r, lo = tanh(zr - 1.96 * se), hi = tanh(zr + 1.96 * se),
                 tz = zr, tse = se, n = n)
    } else {
      fit <- stats::lm(trait ~ proxy, data = x)
      cf <- summary(fit)$coefficients
      if (!"proxy" %in% rownames(cf)) return(NULL)
      est <- cf["proxy", 1]
      se <- cf["proxy", 2]
      if (!is.finite(se) || se <= 0) return(NULL)
      data.frame(age = x$age_grp[[1]], estimate = est, lo = est - 1.96 * se, hi = est + 1.96 * se,
                 tz = est, tse = se, n = n)
    }
  })
  res <- Filter(Negate(is.null), res)
  if (!length(res)) return(empty)
  tab <- do.call(rbind, res)
  tab <- tab[order(tab$age), , drop = FALSE]
  rownames(tab) <- NULL
  out <- empty
  out$table <- tab
  w <- 1 / tab$tse^2
  if (nrow(tab) >= 2) {
    m0 <- stats::lm(tz ~ 1, data = tab, weights = w)
    s0 <- summary(m0)$coefficients
    out$mean_est <- s0[1, 1]
    out$mean_p <- s0[1, 4]
  }
  if (nrow(tab) >= 3) {
    m1 <- stats::lm(tz ~ age, data = tab, weights = w)
    s1 <- summary(m1)$coefficients
    if ("age" %in% rownames(s1)) {
      out$trend_est <- s1["age", 1]
      out$trend_p <- s1["age", 4]
    }
  }
  out
}

a2_interpretation <- function(a2, proxy, alpha = 0.05) {
  if (!nrow(a2$table)) return("Not estimable: need at least 5 individuals with trait data and proxy variation at two or more ages.")
  lab <- if (identical(a2$metric, "r")) "correlation" else "regression slope"
  main <- if (is.finite(a2$mean_p) && a2$mean_p < alpha) {
    sprintf("The %s between %s and the trait is on average %s (weighted mean p = %s), consistent with an age-independent component of selection.",
            lab, proxy, if (a2$mean_est > 0) "positive" else "negative", format_p(a2$mean_p))
  } else {
    sprintf("No consistent average %s between %s and the trait (p = %s).", lab, proxy, format_p(a2$mean_p))
  }
  trend <- if (is.finite(a2$trend_p) && a2$trend_p < alpha) {
    sprintf("The %s changes systematically with age (weighted trend p = %s): consistent with age-dependent selective disappearance/appearance.",
            lab, format_p(a2$trend_p))
  } else if (is.finite(a2$trend_p)) {
    sprintf("No systematic change of the %s with age (trend p = %s).", lab, format_p(a2$trend_p))
  } else {
    "Too few ages to test for a trend with age."
  }
  paste(main, trend)
}

# ---------------------------------------------------------------------------
# Ageing-function bases
# ---------------------------------------------------------------------------
# With standardise = TRUE, polynomial bases use z = (age - mean)/sd and powers of z;
# log and exponential bases are centred and scaled. These are linear
# reparameterisations of the raw bases (given an intercept), so likelihoods and AICs
# are unchanged while numerical conditioning improves (cf. the NaN/Hessian problems
# of raw age x ALR interactions in glmmTMB).
make_age_params <- function(age, fun, standardise = TRUE) {
  age <- age[is.finite(age)]
  p <- list(fun = fun, standardise = isTRUE(standardise), mean_age = mean(age), sd_age = stats::sd(age),
            min_age = min(age), all_positive = min(age) > 0)
  if (!is.finite(p$sd_age) || p$sd_age <= 0) p$sd_age <- 1
  lx <- if (p$all_positive) log(age) else log(age - p$min_age + 1)
  p$mean_log <- mean(lx)
  p$sd_log <- stats::sd(lx)
  if (!is.finite(p$sd_log) || p$sd_log <= 0) p$sd_log <- 1
  ex <- exp(-(age - p$mean_age) / p$sd_age)
  p$mean_exp <- mean(ex)
  p$sd_exp <- stats::sd(ex)
  if (!is.finite(p$sd_exp) || p$sd_exp <= 0) p$sd_exp <- 1
  p
}

age_basis <- function(age, p) {
  age <- as.numeric(age)
  z <- (age - p$mean_age) / p$sd_age
  if (p$fun %in% c("Linear", "Quadratic", "Cubic")) {
    x <- if (isTRUE(p$standardise)) z else age
    out <- data.frame(f1 = x)
    if (p$fun %in% c("Quadratic", "Cubic")) out$f2 <- x^2
    if (identical(p$fun, "Cubic")) out$f3 <- x^3
  } else if (identical(p$fun, "Logarithmic")) {
    lx <- if (isTRUE(p$all_positive)) log(pmax(age, 1e-12)) else log(pmax(age - p$min_age + 1, 1e-12))
    out <- data.frame(f1 = if (isTRUE(p$standardise)) (lx - p$mean_log) / p$sd_log else lx)
  } else {
    ex <- exp(-z)
    out <- data.frame(f1 = if (isTRUE(p$standardise)) (ex - p$mean_exp) / p$sd_exp else ex)
  }
  out
}

MODEL_FUNCTIONS <- c(AGE_FUNCTIONS, A3_NONLINEAR)

basis_label <- function(fun) {
  if (identical(fun, A3_NONLINEAR)) return("a\u00b7exp(b\u00b7z_age)")
  switch(fun, Linear = "age", Quadratic = "age + age\u00b2", Cubic = "age + age\u00b2 + age\u00b3",
         Logarithmic = "log(age)", "Asymptotic exponential" = "exp(\u2212z_age)", "age")
}

# ---------------------------------------------------------------------------
# A3: individual parametric fits
# ---------------------------------------------------------------------------
# Minimum records (at distinct ages) to fit a function to one individual: Linear 2, Quadratic 3, Cubic 4 (fits
# through every point are allowed), Logarithmic and both exponential functions 3.
A3_MIN_RECORDS <- c(Linear = 2L, Quadratic = 3L, Cubic = 4L, Logarithmic = 3L, "Asymptotic exponential" = 3L)
min_records_for <- function(fun) {
  if (identical(fun, A3_NONLINEAR)) return(3L)
  v <- unname(A3_MIN_RECORDS[fun])
  if (length(v) != 1 || is.na(v)) 3L else v
}

# min_resid_df is kept for compatibility; the minimum number of records now depends on the function.
individual_fit_one <- function(age, trait, p, min_resid_df = 1) {
  b <- age_basis(age, p)
  k <- ncol(b) + 1
  n <- length(trait)
  if (n < max(k, min_records_for(p$fun)) || length(unique(age)) < k) return(NULL)
  X <- cbind(`(Intercept)` = 1, as.matrix(b))
  fit <- tryCatch(stats::lm.fit(X, trait), error = function(e) NULL)
  if (is.null(fit) || fit$rank < k) return(NULL)
  cf <- fit$coefficients
  if (any(!is.finite(cf))) return(NULL)
  rss <- sum(fit$residuals^2)
  tss <- sum((trait - mean(trait))^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  df_res <- n - k
  adj <- if (df_res > 0 && is.finite(r2)) 1 - (1 - r2) * (n - 1) / df_res else NA_real_
  kk <- k + 1
  aic <- if (rss > 1e-12) n * log(rss / n) + 2 * kk else NA_real_
  aicc <- if (is.finite(aic) && (n - kk - 1) > 0) aic + 2 * kk * (kk + 1) / (n - kk - 1) else NA_real_
  list(fit = fit, coefficients = cf, stat = data.frame(n = n, R2 = r2, Adjusted_R2 = adj, AICc = aicc))
}

# Non-linear individual function a * exp(b * age), fitted by nls with log-linear starting values.
individual_fit_nonlinear <- function(age, trait, min_resid_df = 1) {
  n <- length(trait)
  k <- 2
  if (n < min_records_for(A3_NONLINEAR) || length(unique(age)) < k) return(NULL)
  pos <- trait > 0
  st <- if (sum(pos) >= 2 && length(unique(age[pos])) >= 2) {
    c0 <- stats::coef(stats::lm(log(trait[pos]) ~ age[pos]))
    list(a = exp(c0[[1]]), b = c0[[2]])
  } else {
    list(a = max(mean(trait), 1e-6), b = 0)
  }
  dd <- data.frame(age = age, trait = trait)
  fit <- tryCatch(suppressWarnings(stats::nls(trait ~ a * exp(b * age), data = dd, start = st,
                                              control = stats::nls.control(maxiter = 200, warnOnly = TRUE))),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$convInfo$isConv)) return(NULL)
  cf <- stats::coef(fit)
  if (any(!is.finite(cf))) return(NULL)
  rss <- sum(stats::residuals(fit)^2)
  tss <- sum((trait - mean(trait))^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  df_res <- n - k
  adj <- if (df_res > 0 && is.finite(r2)) 1 - (1 - r2) * (n - 1) / df_res else NA_real_
  kk <- k + 1
  aic <- if (rss > 1e-12) n * log(rss / n) + 2 * kk else NA_real_
  aicc <- if (is.finite(aic) && (n - kk - 1) > 0) aic + 2 * kk * (kk + 1) / (n - kk - 1) else NA_real_
  list(fit = fit, coefficients = c(a = cf[["a"]], b = cf[["b"]]), stat = data.frame(n = n, R2 = r2, Adjusted_R2 = adj, AICc = aicc))
}

# Two population-level reconstructions from individual fits f(age; theta_i), i = 1..N:
#   mean of coefficients:  f(age; theta_bar), theta_bar = (1/N) sum_i theta_i
#   mean of functions:     f_bar(age) = (1/N) sum_i f(age; theta_i)
# They coincide when f is linear in its parameters (all AGE_FUNCTIONS) and differ otherwise.
# 95% bands: delta method with S/N for the coefficient mean; SD_i[f(age; theta_i)] / sqrt(N) for the function mean.
fit_individual_function <- function(dat, fun, min_resid_df = 1, draw_ids = character(0)) {
  d <- dat[is.finite(dat$trait) & is.finite(dat$age), , drop = FALSE]
  empty <- list(curves = data.frame(), coefs = data.frame(), fit_stats = data.frame(), mean_curve = data.frame(),
                n_total = length(unique(d$id)), fitted_ids = character(0), nonlinear = identical(fun, A3_NONLINEAR))
  if (!nrow(d)) return(empty)
  nonlin <- identical(fun, A3_NONLINEAR)
  p <- if (nonlin) NULL else make_age_params(d$age, fun, standardise = FALSE)
  design <- function(ages) cbind(`(Intercept)` = 1, as.matrix(age_basis(ages, p)))
  curve_of <- function(cf, ages) {
    if (nonlin) cf[["a"]] * exp(cf[["b"]] * ages) else as.numeric(design(ages)[, names(cf), drop = FALSE] %*% cf)
  }
  sp <- split(seq_len(nrow(d)), d$id)
  coef_rows <- list(); stat_rows <- list(); curve_rows <- list()
  for (z in names(sp)) {
    ii <- sp[[z]]
    one <- if (nonlin) individual_fit_nonlinear(d$age[ii], d$trait[ii], min_resid_df) else individual_fit_one(d$age[ii], d$trait[ii], p, min_resid_df)
    if (is.null(one)) next
    cr <- as.data.frame(as.list(one$coefficients), check.names = FALSE)
    cr$id <- z
    coef_rows[[z]] <- cr
    st <- one$stat
    st$id <- z
    stat_rows[[z]] <- st
    if (z %in% draw_ids) {
      ag <- seq(min(d$age[ii]), max(d$age[ii]), length.out = 40)
      curve_rows[[z]] <- data.frame(id = z, age = ag, fitted = curve_of(one$coefficients, ag))
    }
  }
  if (!length(coef_rows)) return(empty)
  coefs <- do.call(rbind, coef_rows)
  fit_stats <- do.call(rbind, stat_rows)
  curves <- if (length(curve_rows)) do.call(rbind, curve_rows) else data.frame()
  cn <- setdiff(names(coefs), "id")
  cm <- as.matrix(coefs[cn])
  N <- nrow(cm)
  mcf <- stats::setNames(colMeans(cm), cn)
  # The average curves are drawn only across ages observed in >= 10 individuals, because beyond
  # that they extrapolate few individuals' fits.
  n_at_age <- table(id_age_means(d)$age)
  ok_age <- as.numeric(names(n_at_age))[as.numeric(n_at_age) >= 10]
  rng <- if (length(ok_age)) range(ok_age) else range(d$age)
  ag <- seq(rng[1], rng[2], length.out = 120)
  pred <- curve_of(mcf, ag)
  if (nonlin) {
    Fm <- outer(cm[, "a"], rep(1, length(ag))) * exp(outer(cm[, "b"], ag))
  } else {
    Xg <- design(ag)[, cn, drop = FALSE]
    Fm <- cm %*% t(Xg)
  }
  fun_mean <- colMeans(Fm)
  se <- se_fun <- rep(NA_real_, length(ag))
  if (N >= 3) {
    S <- stats::cov(cm)
    G <- if (nonlin) cbind(exp(mcf[["b"]] * ag), mcf[["a"]] * ag * exp(mcf[["b"]] * ag)) else Xg
    se <- sqrt(pmax(0, rowSums((G %*% (S / N)) * G)))
    se_fun <- apply(Fm, 2, stats::sd) / sqrt(N)
  }
  list(curves = curves, coefs = coefs, fit_stats = fit_stats,
       mean_curve = data.frame(age = ag, fitted = pred, lo = pred - 1.96 * se, hi = pred + 1.96 * se,
                               fitted_fun = fun_mean, lo_fun = fun_mean - 1.96 * se_fun, hi_fun = fun_mean + 1.96 * se_fun),
       n_total = length(sp), fitted_ids = fit_stats$id, nonlinear = nonlin)
}

compare_individual_functions <- function(dat, min_resid_df = 1, max_individuals = 2000) {
  ids <- unique(dat$id)
  if (length(ids) > max_individuals) {
    keep <- ids[unique(round(seq(1, length(ids), length.out = max_individuals)))]
    dat <- dat[dat$id %in% keep, , drop = FALSE]
  }
  n_considered <- length(unique(dat$id))
  per <- lapply(A3_FUNCTIONS, function(fn) {
    z <- fit_individual_function(dat, fn, min_resid_df)
    if (!nrow(z$fit_stats)) return(NULL)
    s <- z$fit_stats
    s$Function <- fn
    s
  })
  per <- Filter(Negate(is.null), per)
  if (!length(per)) return(data.frame())
  per <- do.call(rbind, per)
  funs <- unique(per$Function)
  fin <- per[is.finite(per$AICc), , drop = FALSE]
  common <- character(0)
  if (nrow(fin)) {
    nf <- tapply(fin$Function, fin$id, function(v) length(unique(v)))
    common <- names(nf)[nf == length(funs)]
  }
  per$Delta_AICc <- NA_real_
  cc <- per$id %in% common
  if (any(cc)) per$Delta_AICc[cc] <- stats::ave(per$AICc[cc], per$id[cc], FUN = function(v) v - min(v))
  out <- do.call(rbind, lapply(funs, function(fn) {
    z <- per[per$Function == fn, , drop = FALSE]
    zc <- z[z$id %in% common, , drop = FALSE]
    data.frame(Function = fn, N_considered = n_considered, N_fitted = nrow(z),
               Mean_R2 = mean(z$R2, na.rm = TRUE), Mean_adj_R2 = mean(z$Adjusted_R2, na.rm = TRUE),
               Median_adj_R2 = stats::median(z$Adjusted_R2, na.rm = TRUE),
               N_common = length(common),
               Mean_dAICc = if (nrow(zc)) mean(zc$Delta_AICc) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  out[order(out$Mean_dAICc, -out$Mean_adj_R2, na.last = TRUE), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Decomposition (Rebke et al. 2010)
# ---------------------------------------------------------------------------
# Pairs are records of the same individual on two successive sampling occasions, and only individuals sampled at both
# occasions contribute (survivor-restricted). On a common schedule two records are successive occasions when their
# interval rounds to one sampling step (0.5 to 1.5 steps), so unequal intervals between occasions (e.g. day 1, then
# weekly on days 7, 14, 21) are linked while a missed occasion (about two steps) is not, and staggered schedules follow
# their own occasions. On an irregular schedule the records are first placed on a common grid of occasions
# (decomposition_on_grid()) and the result carries a caution.
# Start: the first two consecutive sampling ages at which at least one individual was sampled at both. The starting
# value is the mean trait at the first of these ages of the individuals sampled at both; each later value adds the
# mean within-individual change between the next pair of consecutive ages (individuals sampled at both ages only).
# The curve stops at the first later interval with no individual sampled at both ages. On a regular schedule this
# is exactly the chain of one-step differences.
decomposition_trajectory <- function(dat) {
  ia <- id_age_means(dat)
  if (nrow(ia) < 2) return(data.frame())
  step <- infer_age_step(ia$age, ia$id)
  if (!is.finite(step) || step <= 0) return(data.frame())
  tol <- 0.01 * step
  ia <- ia[order(ia$id, ia$age), , drop = FALSE]
  # an irregular schedule (more than a fifth of the records are not a whole number of steps after the individual's first
  # record, the 'Sampling schedule' rule of the integrity table) is decomposed on a common grid of occasions instead
  rel <- (ia$age - stats::ave(ia$age, ia$id, FUN = min)) / step
  if (mean(abs(rel - round(rel)) > 0.01) > 0.2) return(decomposition_on_grid(ia))
  n <- nrow(ia)
  gap <- diff(ia$age)
  j <- which(ia$id[-1] == ia$id[-n] & gap >= 0.5 * step & gap < 1.5 * step)
  if (!length(j)) return(data.frame())
  # group floating-point variants of the same ages, but report the actual ages
  pairs <- data.frame(key = round(ia$age[j] / tol), end_key = round(ia$age[j + 1] / tol), a = ia$age[j], b = ia$age[j + 1],
                      diff = ia$trait[j + 1] - ia$trait[j], start = ia$trait[j])
  pairs$link <- paste(pairs$key, pairs$end_key)
  # one link per starting occasion: its most common next occasion (the only one on a common schedule)
  lk <- pairs[!duplicated(pairs$link), c("link", "key", "end_key"), drop = FALSE]
  lk$n <- as.numeric(table(pairs$link)[lk$link])
  lk <- lk[order(lk$key, -lk$n, lk$end_key), , drop = FALSE]
  lk <- lk[!duplicated(lk$key), , drop = FALSE]
  lk$a <- as.numeric(tapply(pairs$a, pairs$link, mean)[lk$link])
  lk$b <- as.numeric(tapply(pairs$b, pairs$link, mean)[lk$link])
  lk$inc <- as.numeric(tapply(pairs$diff, pairs$link, mean)[lk$link])
  lk$start <- as.numeric(tapply(pairs$start, pairs$link, mean)[lk$link])
  # walk from the first starting occasion through consecutive links and stop at the first gap
  i <- 1L
  current <- lk$start[[1]]
  out <- data.frame(age = lk$a[[1]], fitted = current, n_pairs = lk$n[[1]])
  repeat {
    current <- current + lk$inc[i]
    out <- rbind(out, data.frame(age = lk$b[i], fitted = current, n_pairs = lk$n[i]))
    nxt <- which(abs(lk$a - lk$b[i]) <= tol)
    if (!length(nxt)) break
    i <- nxt[[1]]
  }
  out
}

# Decomposition on an irregular schedule: every record is placed on the nearest occasion of a common grid, spaced by the
# median interval between an individual's successive records and starting at the earliest age; repeated records of an
# individual within one occasion are averaged; the mean within-individual change is then taken between successive
# occasions over the individuals sampled at both, and chained exactly as on a regular schedule. Plotted ages are the
# mean ages of the records on each occasion. When records of different ages share an occasion, the result carries
# attr(, "caution").
decomposition_on_grid <- function(ia) {
  n0 <- nrow(ia)
  if (n0 < 2) return(data.frame())
  gaps <- diff(ia$age)[ia$id[-1] == ia$id[-n0]]
  gaps <- gaps[is.finite(gaps) & gaps > 1e-6 * max(1, diff(range(ia$age)))]
  if (!length(gaps)) return(data.frame())
  grid <- stats::median(gaps)
  occ <- round((ia$age - min(ia$age)) / grid)
  occ_age <- tapply(ia$age, occ, mean)
  merged <- any(tapply(ia$age, occ, function(a) diff(range(a)) > 0.01 * grid))
  key <- paste(ia$id, occ, sep = "\r")
  one <- !duplicated(key)
  ob <- data.frame(id = ia$id[one], occ = occ[one], trait = as.numeric(tapply(ia$trait, key, mean)[key[one]]),
                   stringsAsFactors = FALSE)
  ob <- ob[order(ob$id, ob$occ), , drop = FALSE]
  n <- nrow(ob)
  if (n < 2) return(data.frame())
  j <- which(ob$id[-1] == ob$id[-n] & diff(ob$occ) == 1)
  if (!length(j)) return(data.frame())
  k_of <- ob$occ[j]
  change <- ob$trait[j + 1] - ob$trait[j]
  inc <- tapply(change, k_of, mean)
  npr <- tapply(change, k_of, length)
  k <- min(k_of)
  current <- mean(ob$trait[j][k_of == k])
  out <- data.frame(age = as.numeric(occ_age[as.character(k)]), fitted = current, n_pairs = as.numeric(npr[as.character(k)]))
  repeat {
    current <- current + as.numeric(inc[as.character(k)])
    out <- rbind(out, data.frame(age = as.numeric(occ_age[as.character(k + 1)]), fitted = current,
                                 n_pairs = as.numeric(npr[as.character(k)])))
    k <- k + 1
    if (!as.character(k) %in% names(inc)) break
  }
  if (isTRUE(merged)) {
    attr(out, "caution") <- paste0("The ages are not on a common sampling schedule, so each record was placed on the nearest occasion of a grid ",
                                   format_num(grid), " age units apart (the median interval between an individual's successive records), ",
                                   "and records of different ages were merged: the decomposition is approximate.")
  }
  out
}

observed_trajectory <- function(dat) {
  ia <- id_age_means(dat)
  if (!nrow(ia)) return(data.frame())
  z <- stats::aggregate(trait ~ age, data = ia, FUN = mean)
  z$n <- as.numeric(table(ia$age)[as.character(z$age)])
  names(z)[names(z) == "trait"] <- "fitted"
  z
}

# ---------------------------------------------------------------------------
# Comparative Models 1-8
# ---------------------------------------------------------------------------
# among = "linear": ALR, LS, AFR and mean age enter linearly (default, as in the manuscript).
# among = "same": for polynomial ageing functions they enter with the same polynomial order
# (e.g. ALR + ALR2 with a quadratic; mean_f1 + mean_f2 in the centring models).
# cov_age: covariates that interact with the ageing terms (covariate x age), in addition to their main effect.
model_formula_strings <- function(b, covars = character(0), among = "linear", cov_age = character(0), cov_pairs = character(0)) {
  delta <- paste0("delta_", b)
  poly_among <- identical(among, "same") && length(b) > 1
  mean_terms <- if (poly_among) paste0("mean_", b) else "mean_f1"
  if (poly_among) {
    prox <- function(v) paste(c(v, paste0(v, seq_along(b)[-1])), collapse = " + ")
    bs <- paste(b, collapse = " + ")
    core <- c(
      M1 = bs,
      M2 = paste(bs, "+", prox("ALR")),
      M3 = paste(c(mean_terms, delta), collapse = " + "),
      M4 = paste0("(", bs, ") * (", prox("ALR"), ")"),
      M5 = paste0("(", paste(mean_terms, collapse = " + "), ") * (", paste(delta, collapse = " + "), ")"),
      M6 = paste0("(", bs, ") * (", prox("LS"), ")"),
      M7 = paste(bs, "+", prox("ALR"), "+", prox("AFR")),
      M8 = paste0("(", bs, ") * (", prox("ALR"), ") + (", bs, ") * (", prox("AFR"), ")"),
      M9 = paste0("(", bs, ") * (", prox("ALR"), ") + ", prox("AFR")),
      M10 = paste0("(", bs, ") * (", prox("AFR"), ") + ", prox("ALR"))
    )
    return(add_covariate_terms(core, b, covars, cov_age, mean_terms, cov_pairs))
  }
  core <- c(
    M1 = paste(b, collapse = " + "),
    M2 = paste(c(b, "ALR"), collapse = " + "),
    M3 = paste(c("mean_f1", delta), collapse = " + "),
    M4 = paste(paste0(b, " * ALR"), collapse = " + "),
    M5 = paste0("mean_f1 * (", paste(delta, collapse = " + "), ")"),
    M6 = paste(paste0(b, " * LS"), collapse = " + "),
    M7 = paste(c(b, "ALR", "AFR"), collapse = " + "),
    M8 = paste(c(paste0(b, " * ALR"), paste0(b, " * AFR")), collapse = " + "),
    M9 = paste(c(paste0(b, " * ALR"), "AFR"), collapse = " + "),
    M10 = paste(c(paste0(b, " * AFR"), "ALR"), collapse = " + ")
  )
  add_covariate_terms(core, b, covars, cov_age, mean_terms, cov_pairs)
}

# Extra terms added to chosen models: extra = list(models = c("M5"), terms = c("ALR", "AFR_x_age", "cv_diet")).
# Keys: ALR, AFR, LS (additive), mean_age (additive), ALR_x_age, AFR_x_age, LS_x_age (proxy x ageing terms,
# with main effects), a covariate's internal name (additive) or "<covariate>:age" (covariate x ageing terms).
# Built terms (term builder): "term:" + up to three tokens (age, ALR, AFR, LS, mean_age or a covariate's internal
# name) joined by "+" (additive) or "*" (interaction with main effects, as in R), e.g. "term:ALR*AFR*age".
EXTRA_TERM_KEYS <- c("ALR (additive)" = "ALR", "AFR (additive)" = "AFR", "LS (additive)" = "LS", "Mean age (additive)" = "mean_age",
                     "ALR \u00d7 age" = "ALR_x_age", "AFR \u00d7 age" = "AFR_x_age", "LS \u00d7 age" = "LS_x_age")
BUILT_TERM_TOKENS <- c("age", "ALR", "AFR", "LS", "mean_age")
# Additive components of a built term, each a vector of interacting tokens: "term:ALR*AFR+cv_diet" gives
# list(c("ALR", "AFR"), "cv_diet"). NULL when the key is not a built term.
parse_built_term <- function(k) {
  if (!is.character(k) || length(k) != 1 || is.na(k) || !startsWith(k, "term:")) return(NULL)
  comps <- strsplit(strsplit(substring(k, 6), "+", fixed = TRUE)[[1]], "*", fixed = TRUE)
  comps <- lapply(comps, function(x) unique(trimws(x[nzchar(trimws(x))])))
  comps[vapply(comps, length, integer(1)) > 0]
}
built_term_ok <- function(k, covars = character(0)) {
  comps <- parse_built_term(k)
  toks <- unlist(comps)
  length(comps) > 0 && length(toks) <= 3 && all(toks %in% c(BUILT_TERM_TOKENS, covars))
}
built_term_age_interaction <- function(k) {
  comps <- parse_built_term(k)
  length(comps) > 0 && any(vapply(comps, function(cp) length(cp) > 1 && "age" %in% cp, logical(1)))
}
# Fixed-effect terms of a built key: age stands for every ageing term of the function (for the centring Models 3 and
# 5, the individual mean and within-individual deviation terms), ALR/AFR/LS for their linear or polynomial
# among-individual terms; components are joined as written (A * B = A + B + A:B).
built_term_formula <- function(k, b, among = "linear", centred = FALSE) {
  comps <- parse_built_term(k)
  if (!length(comps)) return("")
  poly_among <- identical(among, "same") && length(b) > 1
  prox <- function(v) if (poly_among) paste(c(v, paste0(v, seq_along(b)[-1])), collapse = " + ") else v
  mean_terms <- if (poly_among) paste(paste0("mean_", b), collapse = " + ") else "mean_f1"
  age_terms <- if (isTRUE(centred)) paste(c(strsplit(mean_terms, " + ", fixed = TRUE)[[1]], paste0("delta_", b)), collapse = " + ") else paste(b, collapse = " + ")
  wrap <- function(x) if (grepl("+", x, fixed = TRUE)) paste0("(", x, ")") else x
  piece <- function(tok) wrap(switch(tok, age = age_terms, ALR = prox("ALR"), AFR = prox("AFR"), LS = prox("LS"), mean_age = mean_terms, tok))
  paste(vapply(comps, function(cp) paste(vapply(cp, piece, character(1)), collapse = " * "), character(1)), collapse = " + ")
}
# Plain-language label of an extra-term key (menus, summaries).
extra_term_display <- function(k, meta = NULL) {
  labs <- meta$cov_labels %||% character(0)
  tok <- function(t) if (t %in% names(labs)) unname(labs[[t]]) else switch(t, age = "age", mean_age = "mean age", display_term(t))
  comps <- parse_built_term(k)
  if (length(comps)) return(paste(vapply(comps, function(cp) paste(vapply(cp, tok, character(1)), collapse = " \u00d7 "), character(1)), collapse = " + "))
  if (k %in% EXTRA_TERM_KEYS) return(names(EXTRA_TERM_KEYS)[match(k, EXTRA_TERM_KEYS)])
  if (grepl(":age$", k)) return(paste0(tok(sub(":age$", "", k)), " \u00d7 age"))
  if (k %in% names(labs)) return(paste0(labs[[k]], " (additive)"))
  k
}
apply_extra_terms <- function(fs, b, extra = NULL, among = "linear") {
  ex <- per_model_extra(extra)
  if (!length(ex)) return(fs)
  poly_among <- identical(among, "same") && length(b) > 1
  prox <- function(v) if (poly_among) paste(c(v, paste0(v, seq_along(b)[-1])), collapse = " + ") else v
  mean_terms <- if (poly_among) paste(paste0("mean_", b), collapse = " + ") else "mean_f1"
  for (m in intersect(names(ex), names(fs))) {
    age_terms <- if (m %in% c("M3", "M5")) paste(c(strsplit(mean_terms, " + ", fixed = TRUE)[[1]], paste0("delta_", b)), collapse = " + ") else paste(b, collapse = " + ")
    add <- vapply(ex[[m]], function(k) {
      if (k %in% c("ALR", "AFR", "LS")) return(prox(k))
      if (identical(k, "mean_age")) return(mean_terms)
      if (k %in% c("ALR_x_age", "AFR_x_age", "LS_x_age")) return(paste0("(", age_terms, ") * (", prox(sub("_x_age$", "", k)), ")"))
      if (startsWith(k, "term:")) return(built_term_formula(k, b, among, m %in% c("M3", "M5")))
      if (grepl(":age$", k)) return(paste0(sub(":age$", "", k), " * (", age_terms, ")"))
      k
    }, character(1))
    fs[[m]] <- paste(fs[[m]], "+", paste(add, collapse = " + "))
  }
  fs
}

# Extra terms are stored per model: list(M5 = "ALR", M6 = c("ALR", "cv_diet:age")). The older form
# list(models = ..., terms = ...) (the same terms for several models) is converted.
per_model_extra <- function(extra) {
  if (is.null(extra) || !length(extra)) return(list())
  if (!is.null(extra$models) && !is.null(extra$terms)) {
    return(stats::setNames(rep(list(unique(extra$terms)), length(extra$models)), extra$models))
  }
  extra[intersect(names(extra), MODEL_IDS)]
}
# Keep only extra-term keys that are known or refer to covariates in the analysed data.
clean_extra_terms <- function(extra, covars) {
  ex <- per_model_extra(extra)
  out <- list()
  for (m in names(ex)) {
    k <- unique(as.character(ex[[m]]))
    keep <- k[k %in% EXTRA_TERM_KEYS | k %in% covars | (grepl(":age$", k) & sub(":age$", "", k) %in% covars) |
              vapply(k, built_term_ok, logical(1), covars = covars, USE.NAMES = FALSE)]
    if (length(keep)) out[[m]] <- keep
  }
  if (length(out)) out else NULL
}

# TRUE if every term of the smaller model is also in the larger one (so a likelihood-ratio test is valid).
formula_nested <- function(small, large) {
  tl <- function(s) {
    lab <- tryCatch(attr(stats::terms(stats::as.formula(paste("y ~", s))), "term.labels"), error = function(e) NULL)
    if (is.null(lab)) return(NULL)
    vapply(strsplit(lab, ":", fixed = TRUE), function(p) paste(sort(p), collapse = ":"), character(1))
  }
  a <- tl(small)
  b <- tl(large)
  !is.null(a) && !is.null(b) && all(a %in% b) && length(b) > length(a)
}

add_covariate_terms <- function(core, b, covars, cov_age, mean_terms = "mean_f1", cov_pairs = character(0)) {
  cov_age <- intersect(cov_age, covars)
  cov_pairs <- cov_pairs[vapply(strsplit(cov_pairs, ":", fixed = TRUE), function(p) length(p) == 2 && all(p %in% covars), logical(1))]
  if (length(cov_pairs)) core[] <- paste(core, "+", paste(cov_pairs, collapse = " + "))
  if (length(cov_age)) {
    age_default <- paste(b, collapse = " + ")
    age_centred <- paste(c(mean_terms, paste0("delta_", b)), collapse = " + ")
    for (m in names(core)) {
      at <- if (m %in% c("M3", "M5")) age_centred else age_default
      core[[m]] <- paste0(core[[m]], " + ", paste0(cov_age, ":(", at, ")", collapse = " + "))
    }
  }
  if (length(covars)) core[] <- paste0(paste(covars, collapse = " + "), " + ", core)
  core
}

# Data support for individual ageing slopes (A3 diagnostics and random-slope advice).
individual_data_support <- function(dat) {
  d <- dat[is.finite(dat$trait) & is.finite(dat$age), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  sp <- split(d$age, d$id)
  n_obs <- vapply(sp, length, integer(1))
  n_ages <- vapply(sp, function(a) length(unique(a)), integer(1))
  span <- vapply(sp, function(a) diff(range(a)), numeric(1))
  q <- stats::quantile(n_obs, c(0.25, 0.75), names = FALSE)
  list(individuals = length(sp), observations = nrow(d), median_obs = stats::median(n_obs), iqr_low = q[[1]], iqr_high = q[[2]],
       pct_2 = 100 * mean(n_ages >= 2), pct_3 = 100 * mean(n_ages >= 3), pct_4 = 100 * mean(n_ages >= 4),
       n_3 = sum(n_ages >= 3), distinct_ages = length(unique(d$age)), median_span = stats::median(span))
}

# Heuristic advice: slope variance is separable from residual variance only through individuals with
# >= 3 distinct ages, and lme4 requires more observations than individual-level random effects.
# Warn when the data cannot support individual slopes at all.
data_support_warning <- function(sup) {
  if (is.null(sup)) return("")
  msgs <- character(0)
  if (isTRUE(sup$individuals < 30)) msgs <- c(msgs, sprintf("only %s individuals", format(sup$individuals, big.mark = ",")))
  if (isTRUE(sup$n_3 < 20)) msgs <- c(msgs, sprintf("only %s individuals with 3 or more distinct ages", format(sup$n_3, big.mark = ",")))
  if (isTRUE(sup$median_obs < 3)) msgs <- c(msgs, sprintf("a median of %s records per individual", format_num(sup$median_obs)))
  if (isTRUE(sup$distinct_ages < 4)) msgs <- c(msgs, sprintf("only %s distinct ages in the data", format(sup$distinct_ages, big.mark = ",")))
  if (!length(msgs)) return("")
  paste0("Warning: these data give weak support for individual-level slopes (", paste(msgs, collapse = "; "),
         "). Individual fits, the function comparison and random slopes will all be imprecise here, and model rankings that depend on among-individual differences in ageing rate should be treated as tentative.")
}

random_slope_advice <- function(sup) {
  if (is.null(sup)) return(list(level = "none", recommended = "none", text = "No trait data."))
  n_re <- 2 * sup$individuals
  facts <- sprintf("%d of %d individuals (%.0f%%) have \u2265 3 distinct ages; %s observations for %s individual-level random effects; %d distinct ages; median age span %s.",
                   sup$n_3, sup$individuals, sup$pct_3, format(sup$observations, big.mark = ","), format(n_re, big.mark = ","),
                   sup$distinct_ages, format_num(sup$median_span))
  if (sup$distinct_ages < 3 || sup$n_3 < 10 || sup$observations <= n_re) {
    list(level = "none", recommended = "none", facts = facts,
         text = paste("Random slopes are not supported by these data:", facts, "Use a random intercept."))
  } else if (sup$pct_3 >= 50 && sup$n_3 >= 30) {
    list(level = "good", recommended = "correlated", facts = facts,
         text = paste("Random slopes are supported:", facts, "A correlated intercept and slope can be estimated; if it has boundary or convergence notes, compare with the uncorrelated version."))
  } else {
    list(level = "limited", recommended = "uncorrelated", facts = facts,
         text = paste("Support for random slopes is limited:", facts, "The uncorrelated structure (no intercept\u2013slope correlation) is the more stable choice."))
  }
}

# Individual IDs are already made unique within groups when nesting is selected
# (standardise_data), so (1 | group) + (1 | id) is the nested (1 | group/id) structure;
# without nesting the same terms give crossed random effects. With a second grouping level the groups are also made
# unique within top-level groups, so (1 | group2) + (1 | group) + (1 | id) is (1 | group2/group/id).
random_term_string <- function(b1 = "f1", random_slope = FALSE, has_group = FALSE, extra = character(0), has_group2 = FALSE) {
  id_term <- switch(normalise_slope(random_slope),
    correlated = paste0("(1 + ", b1, " | id)"),
    uncorrelated = paste0("(1 | id) + (0 + ", b1, " | id)"),
    "(1 | id)")
  terms <- c(if (isTRUE(has_group2)) "(1 | group2)", if (isTRUE(has_group)) "(1 | group)", id_term, if (length(extra)) paste0("(1 | ", extra, ")"))
  paste(terms, collapse = " + ")
}

MODEL_FAMILIES <- c(
  "Gaussian (lme4::lmer)" = "gaussian",
  "Poisson (glmmTMB)" = "poisson",
  "Negative binomial, quadratic variance (nbinom2)" = "nbinom2",
  "Negative binomial, linear variance (nbinom1)" = "nbinom1",
  "Zero-inflated Poisson" = "zip",
  "Zero-inflated negative binomial (nbinom2)" = "zinb",
  "Zero-inflated negative binomial (nbinom1)" = "zinb1",
  "Binomial \u2014 binary 0/1, or a proportion with a trials column" = "binomial",
  "Beta-binomial \u2014 overdispersed proportions (glmmTMB)" = "betabinomial"
)
COUNT_FAMILIES <- c("poisson", "nbinom2", "nbinom1", "zip", "zinb", "zinb1")
# Families with a zero-inflation component (glmmTMB ziformula).
ZI_FAMILIES <- c("zip", "zinb", "zinb1")
# Binomial families: the response is a binary outcome (0/1) or a proportion of successes, in which
# case the number of trials must be mapped on the Data tab and is passed as prior weights.
BINOMIAL_FAMILIES <- c("binomial", "betabinomial")

model_definition_table <- function(age_function, covars = character(0), random_str = "(1 | id)", cov_labels = character(0),
                                   among = "linear", cov_age = character(0), cov_pairs = character(0)) {
  raw <- basis_label(age_function)
  poly <- identical(among, "same") && age_function %in% c("Quadratic", "Cubic")
  pw <- if (identical(age_function, "Cubic")) c("", "\u00b2", "\u00b3") else c("", "\u00b2")
  px <- function(v) if (poly) paste0("(", paste0(v, pw, collapse = " + "), ")") else v
  first <- switch(age_function, Logarithmic = "log(age)", "Asymptotic exponential" = "exp(\u2212z_age)", "age")
  deltas <- switch(age_function,
    Linear = "\u0394age", Quadratic = "\u0394age + \u0394age\u00b2", Cubic = "\u0394age + \u0394age\u00b2 + \u0394age\u00b3",
    Logarithmic = "\u0394log(age)", "Asymptotic exponential" = "\u0394exp(\u2212z_age)", "\u0394age")
  cv_txt <- if (length(covars)) paste0(paste(unname(cov_labels[covars] %||% display_term(covars)), collapse = " + "), " + ") else ""
  cv_age_txt <- if (length(intersect(cov_age, covars))) paste0(" + ", paste(unname(cov_labels[intersect(cov_age, covars)] %||% display_term(cov_age)), collapse = ", "), " \u00d7 age terms") else ""
  if (length(cov_pairs)) cv_age_txt <- paste0(cv_age_txt, " + ", paste(readable_terms(cov_pairs, age_function, cov_labels), collapse = " + "))
  mean_txt <- if (poly) paste0("mean(", paste0(first, pw, collapse = ") + mean("), ")") else paste0("mean(", first, ")")
  data.frame(
    Model = model_label(MODEL_IDS),
    Name = c("Naive", "Additive ALR", "Mean-centring", "ALR interaction", "Centring interaction",
             "LS interaction (positive control)", "Additive ALR + AFR", "Interactive ALR + AFR",
             "Interactive ALR + additive AFR", "Additive ALR + interactive AFR"),
    Fixed_effects = paste0(cv_txt, c(
      raw, paste(raw, "+", px("ALR")), paste0(mean_txt, " + ", deltas), paste0("(", raw, ") \u00d7 ", px("ALR")),
      paste0(if (poly) paste0("(", mean_txt, ")") else mean_txt, " \u00d7 (", deltas, ")"), paste0("(", raw, ") \u00d7 ", px("LS")),
      paste(raw, "+", px("ALR"), "+", px("AFR")), paste0("(", raw, ") \u00d7 ", px("ALR"), " + (", raw, ") \u00d7 ", px("AFR")),
      paste0("(", raw, ") \u00d7 ", px("ALR"), " + ", px("AFR")), paste0("(", raw, ") \u00d7 ", px("AFR"), " + ", px("ALR"))), cv_age_txt),
    Random = random_str,
    Question = c("Negative control: no selective disappearance term.",
                 "Age-independent selective disappearance (van de Pol & Verhulst 2006).",
                 "Within/among-individual separation (van de Pol & Wright 2009; Fay et al. 2022).",
                 "Age-dependent selective disappearance through ALR \u00d7 ageing terms.",
                 "Age-dependent selective disappearance through mean age \u00d7 \u0394age terms.",
                 "Known latent lifespan as the interacting covariate (simulation benchmark).",
                 "Age-independent disappearance (ALR) and appearance (AFR).",
                 "Age-dependent disappearance and appearance.",
                 "Age-dependent disappearance (ALR) with age-independent appearance (AFR).",
                 "Age-independent disappearance (ALR) with age-dependent appearance (AFR)."),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Per-model "i" help: meaning of each model and its exact specification under the current settings
# ---------------------------------------------------------------------------
MODEL_MEANING <- list(
  M1 = list(name = "Naive", proxy = "none", dis = "not accounted for", app = "not accounted for",
            text = "No selective disappearance or appearance term. The trajectory mixes within-individual ageing with changes in which individuals are still being sampled, so it is biased when individuals that disappear early differ from the others (negative control)."),
  M2 = list(name = "Additive ALR", proxy = "ALR (age at last record)", dis = "age-independent", app = "not accounted for",
            text = "ALR-based. ALR enters as a main effect, so longer- and shorter-lived individuals may differ by the same amount at every age (van de Pol & Verhulst 2006)."),
  M3 = list(name = "Mean-age centring", proxy = "individual mean age", dis = "age-independent", app = "not accounted for",
            text = "Centring-based. Each individual's age is split into its mean age (among-individual term) and the deviation from that mean (\u0394age, within-individual ageing), so among-individual differences such as selective disappearance do not bias the within-individual ageing terms (van de Pol & Wright 2009; Fay et al. 2022)."),
  M4 = list(name = "ALR interaction", proxy = "ALR (age at last record)", dis = "age-dependent (includes the age-independent ALR effect)", app = "not accounted for",
            text = "ALR-based. ALR interacts with the ageing terms, so the difference between longer- and shorter-lived individuals can change with age."),
  M5 = list(name = "Centring interaction", proxy = "individual mean age", dis = "age-dependent (includes the age-independent mean-age effect)", app = "not accounted for",
            text = "Centring-based. Individual mean age interacts with the within-individual ageing terms, so individuals sampled at older ages (usually the longer-lived) can age differently."),
  M6 = list(name = "LS interaction (positive control)", proxy = "known lifespan (LS)", dis = "age-dependent (includes the age-independent LS effect)", app = "not accounted for",
            text = "Known lifespan interacts with the ageing terms. Only available when lifespan is truly known; a benchmark for Models 4 and 5."),
  M7 = list(name = "Additive ALR + AFR", proxy = "ALR and AFR (age at first observation)", dis = "age-independent", app = "age-independent",
            text = "ALR-based. ALR and AFR enter as main effects: individuals sampled until later, or first observed later, may differ by constant amounts at all ages."),
  M8 = list(name = "Interactive ALR + AFR", proxy = "ALR and AFR (age at first observation)", dis = "age-dependent", app = "age-dependent",
            text = "ALR-based. ALR and AFR both interact with the ageing terms."),
  M9 = list(name = "Interactive ALR + additive AFR", proxy = "ALR and AFR (age at first observation)", dis = "age-dependent", app = "age-independent",
            text = "ALR-based. ALR interacts with the ageing terms; AFR enters as a main effect only."),
  M10 = list(name = "Additive ALR + interactive AFR", proxy = "ALR and AFR (age at first observation)", dis = "age-independent", app = "age-dependent",
             text = "ALR-based. AFR interacts with the ageing terms; ALR enters as a main effect only.")
)

# Internal term names in plain language, e.g. "cv_diet:f1" -> "diet x age"
readable_terms <- function(x, age_function = "Quadratic", cov_labels = character(0)) {
  age_lab <- switch(age_function, Logarithmic = "log(age)", "Asymptotic exponential" = "exp(\u2212z_age)", "age")
  s <- x
  for (cv in names(cov_labels)) {
    s <- gsub(paste0("\\b", gsub(".", "\\.", cv, fixed = TRUE), "\\b"), cov_labels[[cv]], s, perl = TRUE)
  }
  s <- gsub("\\bmean_f1\\b", paste0("mean(", age_lab, ")"), s, perl = TRUE)
  s <- gsub("\\bmean_f2\\b", "mean(age\u00b2)", s, perl = TRUE)
  s <- gsub("\\bmean_f3\\b", "mean(age\u00b3)", s, perl = TRUE)
  s <- gsub("\\bdelta_f1\\b", paste0("\u0394", age_lab), s, perl = TRUE)
  s <- gsub("\\bdelta_f2\\b", "\u0394age\u00b2", s, perl = TRUE)
  s <- gsub("\\bdelta_f3\\b", "\u0394age\u00b3", s, perl = TRUE)
  s <- gsub("\\bf1\\b", age_lab, s, perl = TRUE)
  s <- gsub("\\bf2\\b", "age\u00b2", s, perl = TRUE)
  s <- gsub("\\bf3\\b", "age\u00b3", s, perl = TRUE)
  s <- gsub("\\b(ALR|AFR|LS)2\\b", "\\1\u00b2", s, perl = TRUE)
  s <- gsub("\\b(ALR|AFR|LS)3\\b", "\\1\u00b3", s, perl = TRUE)
  s <- gsub(" * ", " \u00d7 ", s, fixed = TRUE)
  gsub(":", " \u00d7 ", s, fixed = TRUE)
}

model_info_content <- function(m, s, meta = NULL) {
  mm <- MODEL_MEANING[[m]]
  fn <- s$age_function %||% "Quadratic"
  fam <- s$family %||% "gaussian"
  covars <- meta$covars %||% character(0)
  cov_labels <- meta$cov_labels %||% character(0)
  cov_age <- intersect(meta$cov_age %||% character(0), covars)
  cov_pairs <- meta$cov_pairs %||% character(0)
  extra <- clean_extra_terms(s$extra, covars)
  rhs <- function(x) if (length(x)) paste(x, collapse = " + ") else "1"
  r_formula <- NULL
  if (identical(fn, A3_NONLINEAR) && identical(fam, "gaussian")) {
    if (m %in% c("M3", "M5")) {
      spec <- "Not available with the non-linear exponential function (the centring models need a linear-in-parameters age function)."
    } else {
      parts <- nonlinear_extra_parts(NONLINEAR_PARTS[[m]], extra[[m]] %||% character(0), covars)
      a <- union(c(covars, cov_pairs), parts$a)
      b <- union(cov_age, parts$b)
      spec <- paste0("trait = a \u00b7 exp(b \u00b7 z_age); level a ~ ", readable_terms(rhs(a), "Linear", cov_labels),
                     "; rate b ~ ", readable_terms(rhs(b), "Linear", cov_labels))
    }
  } else {
    fn2 <- if (identical(fn, A3_NONLINEAR)) "Linear" else fn
    b <- switch(fn2, Quadratic = c("f1", "f2"), Cubic = c("f1", "f2", "f3"), "f1")
    among <- s$among %||% "linear"
    fs <- apply_extra_terms(model_formula_strings(b, covars, among, cov_age, cov_pairs), b, extra, among)
    spec <- paste("trait ~", readable_terms(fs[[m]], fn2, cov_labels))
    r_formula <- paste("trait ~", fs[[m]], "+", random_term_string(b[[1]], s$random_slope, isTRUE(meta$has_group), meta$random_terms %||% character(0), has_group2 = isTRUE(meta$has_group2)))
  }
  list(title = paste0(model_label(m), " \u00b7 ", mm$name), meaning = mm$text,
       attributes = c(`Lifespan proxy` = mm$proxy, `Selective disappearance` = mm$dis, `Selective appearance` = mm$app),
       spec = spec, random = if (is.null(meta)) "(1 | ID)" else random_display(meta, s$random_slope),
       family = paste0(family_label(fam), if (fam %in% ZI_FAMILIES) paste0("; zero inflation ", s$zi %||% "~1") else ""),
       ageing = fn, extra = extra[[m]] %||% character(0), r_formula = r_formula)
}

# ---------------------------------------------------------------------------
# Exact model equation for the per-model "i": the fixed effects fully expanded (A * B written as A + B + A:B),
# each term with its own coefficient, the random effects and error distribution, and the R syntax as fitted.
# ---------------------------------------------------------------------------
html_esc <- function(x) htmltools::htmlEscape(as.character(x))

# One component of a model term in mathematical notation (HTML). i = individual, j = sampling occasion.
math_part <- function(p, age_function = "Quadratic", meta = NULL, level = NA_character_) {
  age <- switch(age_function, Logarithmic = "ln(age)", "Asymptotic exponential" = "exp(\u2212age)", "age")
  pw <- function(k) if (k > 1) paste0("<sup>", k, "</sup>") else ""
  if (grepl("^f[1-3]$", p)) {
    k <- as.integer(substring(p, 2))
    return(paste0(age, "<sub>ij</sub>", pw(k)))
  }
  if (grepl("^mean_f[1-3]$", p)) {
    k <- as.integer(substring(p, 7))
    return(paste0("mean(", age, pw(k), ")<sub>i</sub>"))
  }
  if (grepl("^delta_f[1-3]$", p)) {
    k <- as.integer(substring(p, 8))
    return(if (k > 1) paste0("\u0394(", age, pw(k), ")<sub>ij</sub>") else paste0("\u0394", age, "<sub>ij</sub>"))
  }
  if (grepl("^(ALR|AFR|LS)[23]?$", p)) {
    v <- sub("[23]$", "", p)
    k <- if (grepl("[23]$", p)) as.integer(substring(p, nchar(p))) else 1L
    return(paste0(v, "<sub>i</sub>", pw(k)))
  }
  lab <- html_esc((meta$cov_labels %||% character(0))[p] %|NA|% display_term(p))
  if (!is.na(level)) paste0("[", lab, " = ", html_esc(level), "]<sub>ij</sub>") else paste0(lab, "<sub>ij</sub>")
}
`%|NA|%` <- function(a, b) if (length(a) != 1 || is.na(a)) b else a

# Terms as R expands them (terms() order, as in the coefficient table), with one row per coefficient:
# categorical covariates get one indicator per non-reference level.
expand_model_terms <- function(rhs, factor_levels = list()) {
  if (!nzchar(trimws(rhs %||% ""))) return(list())
  tl <- attr(stats::terms(stats::as.formula(paste("~", rhs))), "term.labels")
  out <- list()
  for (t in tl) {
    parts <- strsplit(t, ":", fixed = TRUE)[[1]]
    choices <- lapply(parts, function(p) {
      lv <- factor_levels[[p]]
      if (length(lv) >= 2) as.character(lv[-1]) else NA_character_
    })
    combos <- expand.grid(choices, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    for (r in seq_len(nrow(combos))) {
      lev <- vapply(seq_along(parts), function(q) as.character(combos[r, q]), character(1))
      out[[length(out) + 1]] <- list(term = t, parts = parts, levels = lev,
                                     r_name = paste(ifelse(is.na(lev), parts, paste0(parts, lev)), collapse = ":"))
    }
  }
  out
}

model_equation <- function(m, s, meta = NULL, factor_levels = list(), dat = NULL) {
  fn <- s$age_function %||% "Quadratic"
  fam <- s$family %||% "gaussian"
  covars <- meta$covars %||% character(0)
  cov_age <- intersect(meta$cov_age %||% character(0), covars)
  cov_pairs <- meta$cov_pairs %||% character(0)
  extra <- clean_extra_terms(s$extra, covars)
  rs_req <- normalise_slope(s$random_slope)
  rs <- rs_req
  rs_note <- ""
  if (identical(rs_req, "auto")) {
    adv <- if (!is.null(dat)) tryCatch(random_slope_advice(individual_data_support(dat)), error = function(e) NULL) else NULL
    rs <- adv$recommended %||% "none"
    rs_note <- paste0("Automatic random effects: with the current data this is the ", slope_text(rs), " structure.")
  }
  sb <- function(x, k) paste0(x, "<sub>", k, "</sub>")
  bta <- function(k) sb("\u03b2", k)
  nonlinear <- identical(fn, A3_NONLINEAR) && identical(fam, "gaussian")
  fn_math <- if (identical(fn, A3_NONLINEAR)) "Linear" else fn
  slope_var <- math_part("f1", fn_math, meta)
  # ---- random effects (shared by both model types)
  re_terms <- character(0); re_dist <- character(0)
  if (isTRUE(meta$has_group)) {
    re_terms <- c(re_terms, sb("g", "k"))
    re_dist <- c(re_dist, paste0(sb("g", "k"), " ~ N(0, \u03c3<sup>2</sup><sub>group</sub>)", if (isTRUE(meta$nested)) ": individuals i are nested in groups k" else ": groups k crossed with individuals"))
  }
  if (isTRUE(meta$has_group2)) {
    re_terms <- c(re_terms, sb("h", "l"))
    re_dist <- c(re_dist, paste0(sb("h", "l"), " ~ N(0, \u03c3<sup>2</sup><sub>group2</sub>)", if (isTRUE(meta$nested)) ": groups k are nested in top-level groups l" else ": top-level groups l crossed with groups k and individuals"))
  }
  for (rt in meta$random_terms %||% character(0)) {
    lab <- html_esc(display_term(rt))
    re_terms <- c(re_terms, paste0("v<sub>", lab, "</sub>"))
    re_dist <- c(re_dist, paste0("v<sub>", lab, "</sub> ~ N(0, \u03c3<sup>2</sup><sub>", lab, "</sub>): one intercept per level of ", lab))
  }
  # ---- error distribution and link
  resp <- switch(fam, gaussian = sb("y", "ij"), poisson = , zip = paste0("log(", sb("\u03bb", "ij"), ")"),
                 binomial = , betabinomial = paste0("logit(", sb("p", "ij"), ")"), paste0("log(", sb("\u03bc", "ij"), ")"))
  fam_dist <- switch(fam,
    gaussian = paste0(sb("\u03b5", "ij"), " ~ N(0, \u03c3<sup>2</sup>)"),
    poisson = paste0(sb("y", "ij"), " ~ Poisson(", sb("\u03bb", "ij"), ")"),
    nbinom2 = paste0(sb("y", "ij"), " ~ NegBin(", sb("\u03bc", "ij"), ", \u03b8), Var(y) = \u03bc + \u03bc<sup>2</sup>/\u03b8"),
    nbinom1 = paste0(sb("y", "ij"), " ~ NegBin(", sb("\u03bc", "ij"), ", \u03c6), Var(y) = \u03bc(1 + \u03c6)"),
    zip = paste0(sb("y", "ij"), " = 0 with probability ", sb("\u03c0", "ij"), ", otherwise ", sb("y", "ij"), " ~ Poisson(", sb("\u03bb", "ij"), ")"),
    zinb = paste0(sb("y", "ij"), " = 0 with probability ", sb("\u03c0", "ij"), ", otherwise ", sb("y", "ij"), " ~ NegBin(", sb("\u03bc", "ij"), ", \u03b8)"),
    zinb1 = paste0(sb("y", "ij"), " = 0 with probability ", sb("\u03c0", "ij"), ", otherwise ", sb("y", "ij"), " ~ NegBin(", sb("\u03bc", "ij"), ", \u03c6), Var(y) = \u03bc(1 + \u03c6)"),
    binomial = paste0(sb("y", "ij"), " ~ Binomial(", sb("n", "ij"), ", ", sb("p", "ij"), "), with ", sb("n", "ij"), " = 1 for binary data and the trials column otherwise"),
    betabinomial = paste0(sb("y", "ij"), " ~ Beta-binomial(", sb("n", "ij"), ", ", sb("p", "ij"), ", \u03c6): binomial with extra variation between observations"))
  zi_line <- NULL; zi_coef <- NULL
  if (fam %in% ZI_FAMILIES) {
    zv <- tryCatch(all.vars(stats::as.formula(s$zi %||% "~1")), error = function(e) character(0))
    zt <- expand_model_terms(paste(zv, collapse = " + "), factor_levels)
    zi_line <- paste0("logit(", sb("\u03c0", "ij"), ") = ", sb("\u03b3", 0),
                      paste0(vapply(seq_along(zt), function(k) paste0(" + ", sb("\u03b3", k), " \u00b7 ",
                             paste(vapply(seq_along(zt[[k]]$parts), function(q) math_part(zt[[k]]$parts[[q]], fn_math, meta, zt[[k]]$levels[[q]]), character(1)), collapse = " \u00d7 ")),
                             character(1)), collapse = ""))
    zi_coef <- data.frame(Coefficient = c(sb("\u03b3", 0), vapply(seq_along(zt), function(k) sb("\u03b3", k), character(1))),
                          Multiplies = c("1 (zero-inflation intercept)", vapply(zt, function(z) paste(vapply(seq_along(z$parts), function(q) math_part(z$parts[[q]], fn_math, meta, z$levels[[q]]), character(1)), collapse = " \u00d7 "), character(1))),
                          R_term = c("zi~(Intercept)", vapply(zt, function(z) paste0("zi~", z$r_name), character(1))), stringsAsFactors = FALSE)
  }
  if (nonlinear) {
    if (m %in% c("M3", "M5")) return(list(available = FALSE, message = "Not available with the non-linear exponential function (the centring models need a linear-in-parameters age function)."))
    parts <- nonlinear_extra_parts(NONLINEAR_PARTS[[m]], extra[[m]] %||% character(0), covars)
    a_rhs <- paste(c("1", union(c(covars, cov_pairs), parts$a)), collapse = " + ")
    b_rhs <- paste(c("1", union(cov_age, parts$b)), collapse = " + ")
    ta <- expand_model_terms(a_rhs, factor_levels)
    tb <- expand_model_terms(b_rhs, factor_levels)
    lin <- function(sym, tt) paste0(sb(sym, 0), paste0(vapply(seq_along(tt), function(k) paste0(" + ", sb(sym, k), " \u00b7 ",
                    paste(vapply(seq_along(tt[[k]]$parts), function(q) math_part(tt[[k]]$parts[[q]], "Linear", meta, tt[[k]]$levels[[q]]), character(1)), collapse = " \u00d7 ")),
                    character(1)), collapse = ""))
    ua <- c(sb("u", "a,i"), re_terms)
    ub <- if (rs %in% c("correlated", "uncorrelated")) sb("u", "b,i") else character(0)
    eq <- c(paste0(sb("y", "ij"), " = ", sb("a", "i"), " \u00b7 exp(", sb("b", "i"), " \u00b7 age<sub>ij</sub>) + ", sb("\u03b5", "ij")),
            paste0(sb("a", "i"), " = ", lin("\u03b1", ta), paste0(" + ", ua, collapse = "")),
            paste0(sb("b", "i"), " = ", lin("\u03b2", tb), if (length(ub)) paste0(" + ", ub) else ""))
    dist <- c(paste0(sb("\u03b5", "ij"), " ~ N(0, \u03c3<sup>2</sup>)"),
              switch(rs,
                correlated = paste0("(", sb("u", "a,i"), ", ", sb("u", "b,i"), ") ~ MVN(0, \u03a3): SDs \u03c3<sub>a</sub>, \u03c3<sub>b</sub> and correlation \u03c1"),
                uncorrelated = c(paste0(sb("u", "a,i"), " ~ N(0, \u03c3<sup>2</sup><sub>a</sub>)"), paste0(sb("u", "b,i"), " ~ N(0, \u03c3<sup>2</sup><sub>b</sub>), independent of ", sb("u", "a,i"))),
                paste0(sb("u", "a,i"), " ~ N(0, \u03c3<sup>2</sup><sub>a</sub>)")),
              re_dist)
    coef <- data.frame(Coefficient = c(sb("\u03b1", 0), vapply(seq_along(ta), function(k) sb("\u03b1", k), character(1)),
                                       sb("\u03b2", 0), vapply(seq_along(tb), function(k) sb("\u03b2", k), character(1))),
                       Multiplies = c("level a: intercept", vapply(ta, function(z) paste("level a:", paste(vapply(seq_along(z$parts), function(q) math_part(z$parts[[q]], "Linear", meta, z$levels[[q]]), character(1)), collapse = " \u00d7 ")), character(1)),
                                      "rate b: intercept", vapply(tb, function(z) paste("rate b:", paste(vapply(seq_along(z$parts), function(q) math_part(z$parts[[q]], "Linear", meta, z$levels[[q]]), character(1)), collapse = " \u00d7 ")), character(1))),
                       R_term = c("a.(Intercept)", vapply(ta, function(z) paste0("a.", z$r_name), character(1)), "b.(Intercept)", vapply(tb, function(z) paste0("b.", z$r_name), character(1))),
                       stringsAsFactors = FALSE)
    rand <- switch(rs, correlated = "pdSymm(a + b ~ 1)", uncorrelated = "pdDiag(a + b ~ 1)", "pdSymm(a ~ 1)")
    rhs_or_1 <- function(x) {
      y <- sub("^1 \\+ ", "", x)
      if (!nzchar(y)) "1" else y
    }
    r_call <- paste0("nlme::nlme(trait ~ a * exp(b * f1), data = dat,\n           fixed = list(a ~ ", rhs_or_1(a_rhs), ", b ~ ", rhs_or_1(b_rhs), "),\n           random = ",
                     if (isTRUE(meta$has_group2)) paste0("list(group2 = pdSymm(a ~ 1), group = pdSymm(a ~ 1), id = ", rand, ")") else
                       if (isTRUE(meta$has_group)) paste0("list(group = pdSymm(a ~ 1), id = ", rand, ")") else rand,
                     ", groups = ", if (isTRUE(meta$has_group2)) "~ group2/group/id" else if (isTRUE(meta$has_group)) "~ group/id" else "~ id",
                     ", method = \"ML\", start = <starting values>)")
    return(list(available = TRUE, equation = eq, distributions = c(dist, if (nzchar(rs_note)) rs_note), coefficients = coef,
                basis = "age on its original scale; a (level) and b (rate) are estimated, so the curve is not linear in its coefficients",
                r_compact = paste0("trait ~ a * exp(b * f1); a ~ ", rhs_or_1(a_rhs), "; b ~ ", rhs_or_1(b_rhs)), r_expanded = NULL, r_call = r_call))
  }
  b <- switch(fn_math, Quadratic = c("f1", "f2"), Cubic = c("f1", "f2", "f3"), "f1")
  basis_note <- paste0(switch(fn_math,
                              Logarithmic = "f1 = log(age)",
                              "Asymptotic exponential" = "f1 = exp(\u2212age)",
                              Quadratic = "f1 = age, f2 = age\u00b2",
                              Cubic = "f1 = age, f2 = age\u00b2, f3 = age\u00b3",
                              "f1 = age"),
                       if (isTRUE(s$standardise)) ", centred and scaled across the sampled ages" else "")
  among <- s$among %||% "linear"
  fixed <- apply_extra_terms(model_formula_strings(b, covars, among, cov_age, cov_pairs), b, extra, among)[[m]]
  tt <- expand_model_terms(fixed, factor_levels)
  random_str <- random_term_string(b[[1]], rs, isTRUE(meta$has_group), meta$random_terms %||% character(0), has_group2 = isTRUE(meta$has_group2))
  id_re <- switch(rs, correlated = , uncorrelated = c(sb("u", "0i"), paste0(sb("u", "1i"), " \u00b7 ", slope_var)), sb("u", "0i"))
  terms_math <- vapply(tt, function(z) paste(vapply(seq_along(z$parts), function(q) math_part(z$parts[[q]], fn_math, meta, z$levels[[q]]), character(1)), collapse = " \u00d7 "), character(1))
  eq <- paste0(resp, " = ", bta(0), paste0(vapply(seq_along(tt), function(k) paste0(" + ", bta(k), " \u00b7 ", terms_math[[k]]), character(1)), collapse = ""),
               paste0(" + ", c(id_re, re_terms), collapse = ""), if (identical(fam, "gaussian")) paste0(" + ", sb("\u03b5", "ij")) else "")
  id_dist <- switch(rs,
    correlated = paste0("(", sb("u", "0i"), ", ", sb("u", "1i"), ") ~ MVN(0, \u03a3): SDs \u03c3<sub>u0</sub>, \u03c3<sub>u1</sub> and correlation \u03c1"),
    uncorrelated = c(paste0(sb("u", "0i"), " ~ N(0, \u03c3<sup>2</sup><sub>u0</sub>)"), paste0(sb("u", "1i"), " ~ N(0, \u03c3<sup>2</sup><sub>u1</sub>), independent of ", sb("u", "0i"))),
    paste0(sb("u", "0i"), " ~ N(0, \u03c3<sup>2</sup><sub>u0</sub>)"))
  coef <- data.frame(Coefficient = c(bta(0), vapply(seq_along(tt), function(k) bta(k), character(1))),
                     Multiplies = c("1 (intercept: all standardised variables at their centre, categorical covariates at their reference level)", terms_math),
                     R_term = c("(Intercept)", vapply(tt, function(z) z$r_name, character(1))), stringsAsFactors = FALSE)
  if (!is.null(zi_coef)) coef <- rbind(coef, zi_coef)
  tl <- vapply(tt, function(z) z$term, character(1))
  r_expanded <- paste0("trait ~ 1", paste0(" + ", unique(tl), collapse = ""), " + ", random_str)
  fam_call <- switch(fam, poisson = "poisson()", zip = "poisson()", nbinom1 = "glmmTMB::nbinom1()", zinb1 = "glmmTMB::nbinom1()",
                     binomial = "binomial()", betabinomial = "glmmTMB::betabinomial()", "glmmTMB::nbinom2()")
  r_call <- if (identical(fam, "gaussian")) {
    paste0("lme4::lmer(", r_expanded, ",\n           data = dat, REML = FALSE)")
  } else {
    paste0("glmmTMB::glmmTMB(", r_expanded, ",\n                 data = dat, family = ", fam_call, ", ziformula = ",
           if (fam %in% ZI_FAMILIES) (s$zi %||% "~1") else "~0",
           if (fam %in% BINOMIAL_FAMILIES && isTRUE(meta$has_trials)) ", weights = .trials" else "", ", REML = FALSE)")
  }
  list(available = TRUE, equation = c(eq, zi_line), distributions = c(if (identical(fam, "gaussian")) fam_dist else paste(fam_dist, "with the linear predictor above"),
                                                                     id_dist, re_dist, if (nzchar(rs_note)) rs_note),
       coefficients = coef, basis = basis_note, r_compact = paste("trait ~", fixed, "+", random_str), r_expanded = r_expanded, r_call = r_call)
}

# ---------------------------------------------------------------------------
# Section R scripts: the app's own functions that a section needs, found through their calls, written out
# so that a script runs without disappR
# ---------------------------------------------------------------------------
app_code_closure <- function(entry, env = environment(standardise_data)) {
  seen <- character(0)
  todo <- entry
  while (length(todo)) {
    f <- todo[[1]]
    todo <- todo[-1]
    if (f %in% seen || !exists(f, envir = env, inherits = FALSE)) next
    seen <- c(seen, f)
    obj <- get(f, envir = env, inherits = FALSE)
    if (is.function(obj) && requireNamespace("codetools", quietly = TRUE)) {
      g <- tryCatch(codetools::findGlobals(obj, merge = TRUE), error = function(e) character(0))
      todo <- c(todo, setdiff(g, seen))
    }
  }
  seen
}

app_code_definitions <- function(nms, env = environment(standardise_data)) {
  has_pkg <- c(HAS_GLMMTMB = "glmmTMB", HAS_LMERTEST = "lmerTest", HAS_DHARMA = "DHARMa")
  lab <- function(n) if (identical(make.names(n), n)) n else paste0("`", n, "`")
  objs <- lapply(nms, function(n) get(n, envir = env, inherits = FALSE))
  is_fun <- vapply(objs, is.function, logical(1))
  small <- vapply(objs, function(o) is.function(o) || as.numeric(utils::object.size(o)) < 2e5, logical(1))
  one <- function(n, o) {
    if (n %in% names(has_pkg)) return(sprintf("%s <- requireNamespace(\"%s\", quietly = TRUE)", n, has_pkg[[n]]))
    paste0(lab(n), " <- ", paste(deparse(o, width.cutoff = 500L), collapse = "\n"))
  }
  keep <- small
  c(unlist(Map(one, nms[keep & !is_fun], objs[keep & !is_fun]), use.names = FALSE),
    unlist(Map(one, nms[keep & is_fun], objs[keep & is_fun]), use.names = FALSE))
}

prepare_model_data <- function(dat, age_function = "Quadratic", covars = character(0), standardise = TRUE,
                               has_group = FALSE, random_terms = character(0), has_group2 = FALSE) {
  d <- dat[is.finite(dat$trait) & is.finite(dat$age), , drop = FALSE]
  if (nrow(d) < 2) return(NULL)
  im <- individual_metrics(dat)
  mp <- match(d$id, im$id)
  afr_i <- ifelse(is.finite(im$entry), im$entry, im$first_recorded)
  raw <- list(ALR = im$alr[mp], LS = im$lifespan[mp], AFR = afr_i[mp])
  first <- !duplicated(d$id)
  pp <- list()
  for (v in names(raw)) {
    x <- raw[[v]][first]
    x <- x[is.finite(x)]
    ctr <- if (isTRUE(standardise) && length(x)) mean(x) else 0
    s <- if (length(x) > 1) stats::sd(x) else NA_real_
    scl <- if (isTRUE(standardise) && is.finite(s) && s > 0) s else 1
    pp[[v]] <- c(centre = ctr, scale = scl)
    d[[paste0(v, "_raw")]] <- raw[[v]]
    d[[v]] <- (raw[[v]] - ctr) / scl
    d[[paste0(v, "2")]] <- d[[v]]^2
    d[[paste0(v, "3")]] <- d[[v]]^3
  }
  ap <- make_age_params(d$age, age_function, standardise)
  b <- age_basis(d$age, ap)
  d <- cbind(d, b)
  for (nm in names(b)) {
    mm <- stats::ave(d[[nm]], d$id, FUN = mean)
    d[[paste0("mean_", nm)]] <- mm
    d[[paste0("delta_", nm)]] <- d[[nm]] - mm
  }
  for (cv in covars) {
    if (cv %in% names(d) && !is.numeric(d[[cv]])) d[[cv]] <- factor(d[[cv]])
  }
  for (rt in random_terms) {
    if (rt %in% names(d)) d[[rt]] <- factor(d[[rt]])
  }
  d$id <- factor(d$id)
  if (isTRUE(has_group)) d$group <- factor(d$group)
  if (isTRUE(has_group2)) d$group2 <- factor(d$group2)
  list(data = d, age_params = ap, proxy_params = pp, basis = names(b))
}

# ok: the model can be fitted at all; default: selected by default (Model 6 is left
# unselected when LS is missing for some individuals, because the common-row AIC
# comparison would otherwise drop those individuals from every model).
model_availability <- function(dat, meta) {
  im <- individual_metrics(dat)
  im <- im[im$n_trait > 0, , drop = FALSE]
  ok <- stats::setNames(rep(TRUE, length(MODEL_IDS)), MODEL_IDS)
  why <- stats::setNames(rep("", length(MODEL_IDS)), MODEL_IDS)
  if (!isTRUE(stats::sd(im$alr, na.rm = TRUE) > 0)) {
    ok[c("M2", "M4", "M7", "M8")] <- FALSE
    why[c("M2", "M4", "M7", "M8")] <- "ALR does not vary among individuals."
  }
  if (!isTRUE(meta$has_life) || isTRUE(meta$life_auto)) {
    ok["M6"] <- FALSE
    why["M6"] <- "Needs a known lifespan column (automatic LS = last record is not a positive control)."
  } else if (!isTRUE(stats::sd(im$lifespan, na.rm = TRUE) > 0)) {
    ok["M6"] <- FALSE
    why["M6"] <- "LS does not vary among individuals."
  }
  afr <- ifelse(is.finite(im$entry), im$entry, im$first_recorded)
  if (!isTRUE(stats::sd(im$alr, na.rm = TRUE) > 0)) {
    ok[c("M9", "M10")] <- FALSE
    why[c("M9", "M10")] <- "ALR does not vary among individuals."
  }
  if (!isTRUE(stats::sd(afr, na.rm = TRUE) > 0)) {
    ok[c("M7", "M8", "M9", "M10")] <- FALSE
    why[c("M7", "M8", "M9", "M10")] <- "AFR does not vary among individuals."
  }
  default <- ok
  afr_f <- afr[is.finite(afr)]
  if (ok[["M7"]] && length(afr_f)) {
    mode_afr <- as.numeric(names(sort(table(afr_f), decreasing = TRUE))[[1]])
    n_diff <- sum(abs(afr_f - mode_afr) > 1e-8)
    if (n_diff < max(5, 0.05 * length(afr_f))) {
      default[c("M7", "M8", "M9", "M10")] <- FALSE
      why[c("M7", "M8", "M9", "M10")] <- sprintf("Only %d individuals have an AFR different from the most common value: Models 7-10 are left unselected.", n_diff)
    }
  }
  if (ok[["M6"]] && any(!is.finite(im$lifespan))) {
    default[["M6"]] <- FALSE
    why[["M6"]] <- sprintf("LS missing for %d individuals: selecting Model 6 drops them from every model.", sum(!is.finite(im$lifespan)))
  }
  list(ok = ok, default = default, why = why)
}

fit_one_model <- function(formula_str, random_str, data, family = "gaussian", zi_str = "~1") {
  ff <- stats::as.formula(paste("trait ~", formula_str, "+", random_str))
  msgs <- character(0)
  fit <- tryCatch(
    withCallingHandlers({
      if (identical(family, "gaussian")) {
        ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        if (HAS_LMERTEST) lmerTest::lmer(ff, data = data, REML = FALSE, control = ctrl)
        else lme4::lmer(ff, data = data, REML = FALSE, control = ctrl)
      } else {
        fam <- switch(family, poisson = stats::poisson(), zip = stats::poisson(), nbinom1 = glmmTMB::nbinom1(),
                      zinb1 = glmmTMB::nbinom1(), binomial = stats::binomial(), betabinomial = glmmTMB::betabinomial(),
                      glmmTMB::nbinom2())
        zi <- if (family %in% ZI_FAMILIES) stats::as.formula(zi_str) else ~0
        ctrl <- tryCatch(glmmTMB::glmmTMBControl(rank_check = "adjust"), error = function(e) glmmTMB::glmmTMBControl())
        # Binomial families: proportions are weighted by the number of trials (binary data have one trial).
        # The trials column exists (all NA) even when no trials column is mapped, so weights are used only when every
        # analysed row has a positive number of trials; otherwise binary 0/1 data would lose every row.
        trials <- if (family %in% BINOMIAL_FAMILIES && ".trials" %in% names(data)) as.numeric(data$.trials) else NULL
        wts <- if (length(trials) && all(is.finite(trials) & trials > 0)) trials else NULL
        if (is.null(wts)) {
          glmmTMB::glmmTMB(ff, data = data, family = fam, ziformula = zi, REML = FALSE, control = ctrl)
        } else {
          glmmTMB::glmmTMB(ff, data = data, family = fam, ziformula = zi, REML = FALSE, control = ctrl, weights = wts)
        }
      }
    }, warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      txt <- trimws(conditionMessage(m))
      if (grepl("rank deficient|dropping", txt)) msgs <<- c(msgs, txt)
      invokeRestart("muffleMessage")
    }),
    error = function(e) {
      em <- conditionMessage(e)
      hint <- if (grepl("number of observations", em) && grepl("random effects", em)) {
        " (too few records per individual for this random-effect structure: choose a random intercept or use individuals with more records)"
      } else if (grepl("contrasts", em)) {
        " (a factor covariate has only one level in the analysed rows)"
      } else ""
      msgs <<- c(msgs, paste0("Error: ", em, hint))
      NULL
    })
  notes <- unique(msgs)
  if (!is.null(fit)) {
    if (methods::is(fit, "merMod")) {
      if (isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) FALSE))) notes <- c(notes, "singular fit (a random-effect variance is ~0)")
    } else if (inherits(fit, "glmmTMB")) {
      if (isFALSE(fit$sdr$pdHess)) notes <- c(notes, "non-positive-definite Hessian")
      small <- tryCatch({
        vc <- summary(fit)$varcor$cond
        sds <- unlist(lapply(vc, function(m) attr(m, "stddev")))
        any(!is.finite(sds) | sds < 1e-3)
      }, error = function(e) FALSE)
      if (isTRUE(small)) notes <- c(notes, "random-effect SD \u2248 0 (boundary)")
    }
    if (!is.finite(tryCatch(stats::AIC(fit), error = function(e) NA_real_))) {
      notes <- c(notes, "AIC not finite")
    }
  }
  notes <- unique(notes)
  list(fit = fit, notes = notes, validity = classify_fit(fit, notes))
}

# Valid: no warnings or notes. Caution: boundary (singular) fits, convergence-gradient warnings, rank
# deficiency or other notes. Failed: no fit, a non-positive-definite or singular Hessian, or a non-finite AIC;
# such fits are excluded from automated "lowest AIC" interpretation unless the user overrides it.
classify_fit <- function(fit, notes = character(0)) {
  if (is.null(fit)) return("Failed")
  invalid <- "Hessian is numerically singular|non-positive-definite Hessian|degenerate +Hessian|unable to evaluate scaled gradient|Hessian is not positive definite|invalid Hessian|AIC not finite"
  if (any(grepl(invalid, notes, ignore.case = TRUE))) return("Failed")
  if (length(notes)) "Caution" else "Valid"
}

coef_table <- function(fit, label) {
  if (inherits(fit, "nlme")) {
    tt <- tryCatch(summary(fit)$tTable, error = function(e) NULL)
    if (is.null(tt)) return(NULL)
    return(data.frame(Model = label, Term = display_term(rownames(tt)), Raw_term = rownames(tt), Estimate = tt[, "Value"],
                      SE = tt[, "Std.Error"], Statistic = tt[, "t-value"], P_value = tt[, "p-value"],
                      stringsAsFactors = FALSE, row.names = NULL))
  }
  cf <- tryCatch({
    s <- summary(fit)$coefficients
    if (is.list(s) && !is.matrix(s)) s <- s$cond
    s
  }, error = function(e) NULL)
  if (is.null(cf) || !is.matrix(cf) || !"Estimate" %in% colnames(cf)) return(NULL)
  cn <- colnames(cf)
  stat_col <- grep("value$", cn, value = TRUE)
  p_col <- grep("^Pr\\(", cn, value = TRUE)
  est <- cf[, "Estimate"]
  se <- cf[, "Std. Error"]
  stat <- if (length(stat_col)) cf[, stat_col[[1]]] else est / se
  p <- if (length(p_col)) cf[, p_col[[1]]] else 2 * stats::pnorm(abs(stat), lower.tail = FALSE)
  data.frame(Model = label, Term = display_term(rownames(cf)), Raw_term = rownames(cf), Estimate = est, SE = se, Statistic = stat,
             P_value = p, stringsAsFactors = FALSE, row.names = NULL)
}

# Scaling constants behind each coefficient. A term's scale factor is the product of the scales of its
# standardised components (age term powers, ALR/LS/AFR powers); dividing the estimate (and SE) by it gives
# the coefficient per original unit of the centred variables, e.g. per year of ALR, or per (age - mean age)^2.
# Numeric covariates are not standardised; factor contrasts have no scale.
scaling_constants <- function(res) {
  ap <- res$age_params
  pp <- res$proxy_params
  age_var <- switch(ap$fun, Logarithmic = "log(age)", "Asymptotic exponential" = "exp(\u2212z_age)", "age")
  age_centre <- if (!isTRUE(ap$standardise)) 0 else switch(ap$fun, Logarithmic = ap$mean_log, "Asymptotic exponential" = ap$mean_exp, ap$mean_age)
  age_scale <- if (!isTRUE(ap$standardise)) 1 else switch(ap$fun, Logarithmic = ap$sd_log, "Asymptotic exponential" = ap$sd_exp, ap$sd_age)
  out <- data.frame(Variable = age_var, Centre = age_centre, Scale = age_scale, stringsAsFactors = FALSE)
  for (v in c("ALR", "LS", "AFR")) {
    if (!is.null(pp[[v]])) out <- rbind(out, data.frame(Variable = v, Centre = pp[[v]][["centre"]], Scale = pp[[v]][["scale"]], stringsAsFactors = FALSE))
  }
  out
}

term_scaling <- function(raw_terms, res) {
  sc <- scaling_constants(res)
  age_scale <- sc$Scale[[1]]
  age_var <- sc$Variable[[1]]
  poly <- res$age_params$fun %in% c("Linear", "Quadratic", "Cubic")
  one <- function(raw) {
    if (identical(raw, "(Intercept)")) return(c(scale = NA_real_, unit = "value at the centre of every scaled variable"))
    s <- 1
    units <- character(0)
    if (grepl("^[ab]\\.", raw)) {
      rate <- startsWith(raw, "b.")
      raw <- sub("^[ab]\\.", "", raw)
      if (rate) {
        s <- s * age_scale
        units <- c(units, "rate per age unit")
      }
      if (identical(raw, "(Intercept)")) return(c(scale = if (rate) s else NA_real_, unit = if (rate) "rate per age unit (in the exponent)" else "level a at the centre of every scaled variable"))
    }
    parts <- strsplit(raw, ":", fixed = TRUE)[[1]]
    for (p in parts) {
      base <- sub("^(mean_|delta_)", "", p)
      if (grepl("^f[1-3]$", base)) {
        pw <- if (poly) as.integer(substring(base, 2)) else 1L
        s <- s * age_scale^pw
        units <- c(units, paste0(if (grepl("^mean_", p)) "individual mean " else if (grepl("^delta_", p)) "within-individual deviation of " else "",
                                 if (pw > 1) paste0("(", age_var, " \u2212 centre)^", pw) else paste(age_var, "unit")))
      } else if (grepl("^(ALR|LS|AFR)[23]?$", p)) {
        v <- sub("[23]$", "", p)
        pw <- if (grepl("[23]$", p)) as.integer(substring(p, nchar(p))) else 1L
        s <- s * sc$Scale[sc$Variable == v][[1]]
        if (pw > 1) s <- s * sc$Scale[sc$Variable == v][[1]]^(pw - 1)
        units <- c(units, if (pw > 1) paste0("(", v, " \u2212 centre)^", pw) else paste(v, "unit"))
      } else {
        units <- c(units, paste0(display_term(p), " (not scaled)"))
      }
    }
    c(scale = s, unit = paste(units, collapse = " \u00d7 "))
  }
  m <- vapply(raw_terms, one, character(2))
  data.frame(Scale_factor = suppressWarnings(as.numeric(m["scale", ])), Per = m["unit", ], stringsAsFactors = FALSE, row.names = NULL)
}

# Random-effect variances, standard deviations and correlations of a fitted lme4 or glmmTMB model,
# plus the residual SD (Gaussian) or the family dispersion parameter (negative binomial).
variance_components <- function(fit, label) {
  tryCatch({
    if (inherits(fit, "nlme")) {
      vm <- unclass(nlme::VarCorr(fit))
      rn <- rownames(vm)
      grp <- "id"
      rows <- list()
      prev_term <- NULL
      for (i in seq_len(nrow(vm))) {
        if (grepl(" =$", rn[i])) {
          grp <- sub(" =$", "", trimws(rn[i]))
          next
        }
        sdv <- suppressWarnings(as.numeric(vm[i, "StdDev"]))
        if (!is.finite(sdv)) next
        g <- if (identical(trimws(rn[i]), "Residual")) "Residual" else grp
        rows[[length(rows) + 1]] <- data.frame(Model = label, Group = g, Term = trimws(rn[i]), Type = "SD", Estimate = sdv,
                                               Variance = suppressWarnings(as.numeric(vm[i, "Variance"])), stringsAsFactors = FALSE)
        if ("Corr" %in% colnames(vm)) {
          cr <- suppressWarnings(as.numeric(vm[i, "Corr"]))
          if (is.finite(cr) && !is.null(prev_term)) {
            rows[[length(rows) + 1]] <- data.frame(Model = label, Group = g, Term = paste(prev_term, "~", trimws(rn[i])), Type = "Correlation",
                                                   Estimate = cr, Variance = NA_real_, stringsAsFactors = FALSE)
          }
        }
        prev_term <- trimws(rn[i])
      }
      return(if (length(rows)) do.call(rbind, rows) else NULL)
    }
    if (methods::is(fit, "merMod")) {
      vc <- as.data.frame(lme4::VarCorr(fit))
      is_sd <- is.na(vc$var2)
      data.frame(Model = label, Group = vc$grp,
                 Term = ifelse(is.na(vc$var1), "Residual", ifelse(is_sd, vc$var1, paste(vc$var1, "~", vc$var2))),
                 Type = ifelse(is_sd, "SD", "Correlation"), Estimate = vc$sdcor,
                 Variance = ifelse(is_sd, vc$vcov, NA_real_), stringsAsFactors = FALSE)
    } else {
      vc <- summary(fit)$varcor$cond
      rows <- list()
      for (g in names(vc)) {
        sds <- attr(vc[[g]], "stddev")
        rows[[length(rows) + 1]] <- data.frame(Model = label, Group = g, Term = names(sds), Type = "SD",
                                               Estimate = as.numeric(sds), Variance = as.numeric(sds)^2,
                                               stringsAsFactors = FALSE)
        cr <- attr(vc[[g]], "correlation")
        if (!is.null(cr) && length(sds) > 1) {
          ij <- which(upper.tri(cr), arr.ind = TRUE)
          rows[[length(rows) + 1]] <- data.frame(Model = label, Group = g,
                                                 Term = paste(rownames(cr)[ij[, 1]], "~", colnames(cr)[ij[, 2]]),
                                                 Type = "Correlation", Estimate = as.numeric(cr[ij]), Variance = NA_real_,
                                                 stringsAsFactors = FALSE)
        }
      }
      disp <- tryCatch({
        fam <- stats::family(fit)$family
        if (grepl("nbinom", fam)) stats::sigma(fit) else NA_real_
      }, error = function(e) NA_real_)
      if (length(disp) == 1 && is.finite(disp)) {
        rows[[length(rows) + 1]] <- data.frame(Model = label, Group = "Family", Term = "Dispersion parameter",
                                               Type = "Dispersion", Estimate = disp, Variance = NA_real_,
                                               stringsAsFactors = FALSE)
      }
      if (length(rows)) do.call(rbind, rows) else NULL
    }
  }, error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Non-linear exponential mixed models (Gaussian): trait = a * exp(b * z_age), fitted with nlme by ML
# ---------------------------------------------------------------------------
# The level a and the rate b have their own linear predictors, so selective disappearance can act on the
# level (age-independent: a ~ ALR) or on the rate (age-dependent: b ~ ALR). Individuals have random levels
# (and optionally random rates). Age is always standardised (z_age) for numerical stability. For count
# families a * exp(b * age) is linear in age on the log link, so the Linear function is fitted instead.
NONLINEAR_PARTS <- list(
  M1 = list(a = character(0), b = character(0)),
  M2 = list(a = "ALR", b = character(0)),
  M4 = list(a = "ALR", b = "ALR"),
  M6 = list(a = "LS", b = "LS"),
  M7 = list(a = c("ALR", "AFR"), b = character(0)),
  M8 = list(a = c("ALR", "AFR"), b = c("ALR", "AFR")),
  M9 = list(a = c("ALR", "AFR"), b = "ALR"),
  M10 = list(a = c("ALR", "AFR"), b = "AFR")
)
# Adds extra-term keys to the level (a) and rate (b) parts of a non-linear model. Terms without age act on the
# level; terms that interact with age act on the level and the rate (as ALR x age does). Mean age is not used.
nonlinear_extra_parts <- function(p, keys, covars = character(0)) {
  for (k in keys %||% character(0)) {
    if (k %in% c("ALR", "AFR", "LS")) {
      p$a <- union(p$a, k)
    } else if (k %in% c("ALR_x_age", "AFR_x_age", "LS_x_age")) {
      v <- sub("_x_age$", "", k)
      p$a <- union(p$a, v)
      p$b <- union(p$b, v)
    } else if (grepl(":age$", k) && sub(":age$", "", k) %in% covars) {
      p$b <- union(p$b, sub(":age$", "", k))
    } else {
      for (cp in parse_built_term(k) %||% list()) {
        rest <- setdiff(cp, c("age", "mean_age"))
        if (!length(rest)) next
        lab <- attr(stats::terms(stats::as.formula(paste("~", paste(rest, collapse = " * ")))), "term.labels")
        p$a <- union(p$a, lab)
        if ("age" %in% cp) p$b <- union(p$b, lab)
      }
    }
  }
  p
}

fit_nonlinear_suite <- function(dat, meta, models = MODEL_IDS, random_slope = FALSE, include_invalid = FALSE,
                                extra = NULL, progress = NULL) {
  fail <- function(msg, status = list()) list(ok = FALSE, message = msg, fits = list(), aic = data.frame(),
                                              coefficients = data.frame(), status = status, lrt = data.frame())
  if (!requireNamespace("nlme", quietly = TRUE)) return(fail("The nlme package is needed for the non-linear exponential function."))
  covars <- meta$covars %||% character(0)
  cov_age <- intersect(meta$cov_age %||% character(0), covars)
  cov_pairs <- (meta$cov_pairs %||% character(0))[vapply(strsplit(meta$cov_pairs %||% character(0), ":", fixed = TRUE),
                                                          function(p) all(p %in% covars), logical(1))]
  extra <- clean_extra_terms(extra, covars)
  prep <- prepare_model_data(dat, "Linear", covars, TRUE, meta$has_group, character(0), has_group2 = isTRUE(meta$has_group2))
  if (is.null(prep)) return(fail("No usable trait data."))
  d <- prep$data
  status <- list()
  av <- model_availability(dat, meta)
  for (m in models) {
    if (m %in% c("M3", "M5")) {
      status[[m]] <- "Unavailable: the mean-age centring models are not defined for the non-linear exponential function."
    } else if (m %in% names(av$ok) && !av$ok[[m]]) {
      status[[m]] <- paste("Unavailable:", av$why[[m]])
    }
  }
  models <- setdiff(intersect(models, names(NONLINEAR_PARTS)), names(status))
  if (!length(models)) return(fail("None of the selected models can be fitted with the non-linear exponential function.", status = status))
  parts <- NONLINEAR_PARTS[models]
  if (length(extra)) {
    for (m in intersect(names(extra), models)) parts[[m]] <- nonlinear_extra_parts(parts[[m]], extra[[m]], covars)
  }
  for (m in models) {
    parts[[m]]$a <- union(c(covars, cov_pairs), parts[[m]]$a)
    parts[[m]]$b <- union(cov_age, parts[[m]]$b)
  }
  # interaction terms such as cv_a:cv_b are not columns: use their components
  vars <- unique(c(unlist(strsplit(unlist(lapply(parts, unlist)), ":", fixed = TRUE)), "f1", "id", if (isTRUE(meta$has_group)) "group", if (isTRUE(meta$has_group2)) "group2"))
  ok <- rep(TRUE, nrow(d))
  drop_by <- character(0)
  for (v in vars) {
    x <- d[[v]]
    bad <- if (is.null(x)) rep(TRUE, nrow(d)) else if (is.numeric(x)) !is.finite(x) else is.na(x)
    if (any(bad & ok)) drop_by <- c(drop_by, sprintf("%s: %d", display_term(v), sum(bad & ok)))
    ok <- ok & !bad
  }
  dd <- droplevels(d[ok, , drop = FALSE])
  if (nrow(dd) < 20 || length(unique(dd$id)) < 6) {
    return(fail("Too few complete observations shared by the selected models (need >= 20 rows from >= 6 individuals).", status = status))
  }
  if (length(unique(dd$age)) < 3) return(fail("The non-linear exponential function needs at least 3 distinct ages.", status = status))
  novar <- covars[!vapply(covars, function(cv) { x <- dd[[cv]]; if (is.numeric(x)) isTRUE(stats::sd(x) > 0) else length(unique(x)) >= 2 }, logical(1))]
  drop_term <- function(x) x[!vapply(strsplit(x, ":", fixed = TRUE), function(p) any(p %in% novar), logical(1))]
  for (m in models) {
    parts[[m]]$a <- drop_term(parts[[m]]$a)
    parts[[m]]$b <- drop_term(parts[[m]]$b)
  }
  if (length(novar)) drop_by <- c(drop_by, sprintf("covariate %s omitted (no variation in the analysed rows)", paste(unname(meta$cov_labels[novar]), collapse = ", ")))
  if (length(meta$random_terms)) drop_by <- c(drop_by, "additional crossed random intercepts are not used by the non-linear exponential models")
  rs_request <- normalise_slope(random_slope)
  advice <- random_slope_advice(individual_data_support(dd))
  rs <- if (identical(rs_request, "auto")) advice$recommended else rs_request
  random_note <- if (identical(rs_request, "auto")) paste0("Automatic random effects: ", slope_text(rs), ". ", advice$text) else ""
  id_pd <- switch(rs, correlated = nlme::pdSymm(a + b ~ 1), uncorrelated = nlme::pdDiag(a + b ~ 1), nlme::pdSymm(a ~ 1))
  rand <- if (isTRUE(meta$has_group2)) {
    list(group2 = nlme::pdSymm(a ~ 1), group = nlme::pdSymm(a ~ 1), id = id_pd)
  } else if (isTRUE(meta$has_group)) list(group = nlme::pdSymm(a ~ 1), id = id_pd) else id_pd
  grp <- if (isTRUE(meta$has_group2)) ~ group2 / group / id else if (isTRUE(meta$has_group)) ~ group / id else ~ id
  rstr <- paste0(if (isTRUE(meta$has_group2)) "random level a | group2 + " else "", if (isTRUE(meta$has_group)) "random level a | group + " else "",
                 switch(rs, correlated = "correlated random level a and rate b | id", uncorrelated = "uncorrelated random level a and rate b | id", "random level a | id"))
  pos <- dd$trait > 0
  st0 <- if (mean(pos) > 0.8 && sum(pos) >= 10) {
    c0 <- stats::coef(stats::lm(log(dd$trait[pos]) ~ dd$f1[pos]))
    c(exp(c0[[1]]), c0[[2]])
  } else c(mean(dd$trait), 0)
  rhs <- function(x) if (length(x)) paste(x, collapse = " + ") else "1"
  ncols <- function(x) ncol(stats::model.matrix(stats::as.formula(paste("~", rhs(x))), dd))
  fits <- list(); aic_rows <- list(); coef_rows <- list(); vc_rows <- list(); validity <- list()
  pseudo <- list(); shown <- list()
  for (i in seq_along(models)) {
    m <- models[[i]]
    if (is.function(progress)) progress(i, length(models), model_label(m))
    pa <- parts[[m]]$a
    pb <- parts[[m]]$b
    pseudo[[m]] <- paste(make.names(c(paste0("a_", c("int", pa)), paste0("b_", c("int", pb)))), collapse = " + ")
    shown[[m]] <- paste0("a\u00b7exp(b\u00b7z_age); a ~ ", rhs(pa), "; b ~ ", rhs(pb))
    msgs <- character(0)
    fx <- list(stats::as.formula(paste("a ~", rhs(pa))), stats::as.formula(paste("b ~", rhs(pb))))
    start <- c(st0[1], rep(0, ncols(pa) - 1), st0[2], rep(0, ncols(pb) - 1))
    fit <- tryCatch(withCallingHandlers(
      nlme::nlme(trait ~ a * exp(b * f1), data = dd, fixed = fx, random = rand, groups = grp, start = start, method = "ML",
                 control = nlme::nlmeControl(maxIter = 100, pnlsMaxIter = 20, msMaxIter = 200, returnObject = TRUE)),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        msgs <<- c(msgs, paste("Error:", conditionMessage(e)))
        NULL
      })
    if (!is.null(fit)) {
      if (is.character(fit$apVar)) msgs <- c(msgs, "non-positive-definite approximate variance-covariance matrix (invalid Hessian)")
      if (!is.finite(tryCatch(stats::AIC(fit), error = function(e) NA_real_))) msgs <- c(msgs, "AIC not finite")
    }
    msgs <- unique(msgs)
    validity[[m]] <- classify_fit(fit, msgs)
    if (is.null(fit)) {
      status[[m]] <- paste("Failed:", paste(msgs, collapse = " | "))
      next
    }
    fits[[m]] <- fit
    status[[m]] <- if (length(msgs)) paste0(validity[[m]], ": ", paste(msgs, collapse = " | ")) else "Valid"
    ll <- stats::logLik(fit)
    aic_rows[[m]] <- data.frame(Model = model_label(m), AIC = tryCatch(stats::AIC(fit), error = function(e) NA_real_), df = attr(ll, "df"),
                                logLik = as.numeric(ll), N = nrow(dd), Fit = validity[[m]], stringsAsFactors = FALSE)
    ct <- coef_table(fit, model_label(m))
    if (!is.null(ct)) coef_rows[[m]] <- ct
    vc <- variance_components(fit, model_label(m))
    if (!is.null(vc)) vc_rows[[m]] <- vc
  }
  aic <- if (length(aic_rows)) do.call(rbind, aic_rows) else data.frame()
  n_excluded <- 0L
  if (nrow(aic)) {
    aic$Eligible <- is.finite(aic$AIC) & (aic$Fit != "Failed" | isTRUE(include_invalid))
    n_excluded <- sum(!aic$Eligible)
    aic$Delta_AIC <- NA_real_
    aic$Weight <- NA_real_
    if (any(aic$Eligible)) {
      best <- min(aic$AIC[aic$Eligible])
      aic$Delta_AIC[aic$Eligible] <- aic$AIC[aic$Eligible] - best
      w <- exp(-aic$Delta_AIC[aic$Eligible] / 2)
      aic$Weight[aic$Eligible] <- w / sum(w)
    }
    aic <- aic[order(!aic$Eligible, aic$AIC), , drop = FALSE]
    rownames(aic) <- NULL
  }
  n_eligible <- if (nrow(aic)) sum(aic$Eligible) else 0L
  list(ok = n_eligible > 0,
       message = if (n_eligible > 0) "Non-linear exponential models fitted by maximum likelihood (nlme) on the same complete-case rows." else
         if (nrow(aic)) "Every non-linear fit has an invalid Hessian: simplify the random effects or tick 'Include fits with invalid Hessians' to inspect them." else
         "No non-linear model could be fitted (see the fitting status); non-linear fits need good starting values and positive trait values.",
       fits = fits, aic = aic, coefficients = if (length(coef_rows)) do.call(rbind, coef_rows) else data.frame(),
       status = status, data = dd, n_dropped = sum(!ok), drop_by = drop_by, formulas = shown, random = rstr,
       basis = "f1", age_params = prep$age_params, proxy_params = prep$proxy_params, family = "gaussian", zi = "~1",
       covars = setdiff(covars, novar), random_terms = character(0),
       varcomp = if (length(vc_rows)) do.call(rbind, vc_rows) else data.frame(),
       age_function = A3_NONLINEAR, random_slope = rs, standardise = TRUE, nonlinear = TRUE, parts = parts,
       validity = validity, n_excluded = n_excluded, include_invalid = isTRUE(include_invalid),
       random_request = rs_request, random_structure = rs, random_note = random_note, random_advice = advice,
       among = "linear", cov_age = setdiff(cov_age, novar), extra = extra,
       cov_pairs = cov_pairs[!vapply(strsplit(cov_pairs, ":", fixed = TRUE), function(p) any(p %in% novar), logical(1))],
       lrt = lrt_table(fits, validity, include_invalid, pseudo))
}

fit_model_suite <- function(dat, meta, models = MODEL_IDS, age_function = "Quadratic", family = "gaussian",
                            random_slope = FALSE, standardise = TRUE, zi_str = "~1", progress = NULL,
                            among = "linear", include_invalid = FALSE, extra = NULL) {
  fail <- function(msg, status = list()) list(ok = FALSE, message = msg, fits = list(), aic = data.frame(),
                                              coefficients = data.frame(), status = status, lrt = data.frame())
  nonlinear_note <- NULL
  if (identical(age_function, A3_NONLINEAR)) {
    if (identical(family, "gaussian")) {
      return(fit_nonlinear_suite(dat, meta, models, random_slope, include_invalid, extra, progress))
    }
    age_function <- "Linear"
    nonlinear_note <- "a\u00b7exp(b\u00b7age) is linear in age on the log or logit link of the non-Gaussian families, so it was fitted as the Linear function"
  }
  if (!identical(family, "gaussian") && !HAS_GLMMTMB) {
    return(fail("The selected family needs the glmmTMB package: install.packages('glmmTMB'), or choose Gaussian."))
  }
  if (family %in% COUNT_FAMILIES) {
    v <- dat$trait[is.finite(dat$trait)]
    if (any(v < 0) || any(abs(v - round(v)) > 1e-8)) {
      return(fail(paste("Count families need non-negative integer trait values.",
                        if (identical(meta$dup_action %||% "keep", "mean")) "Averaging duplicate ID x age records has made the trait non-integer: choose 'keep duplicates' on the Data tab to use a count family." else "")))
    }
  }
  if (family %in% BINOMIAL_FAMILIES) {
    v <- dat$trait[is.finite(dat$trait)]
    if (any(v < 0) || any(v > 1)) {
      return(fail("Binomial families need the trait to be 0/1 (binary) or a proportion between 0 and 1. For counts of successes, divide by the number of trials and choose that number under 'Weights' below the error family."))
    }
    binary <- all(v %in% c(0, 1))
    if (!binary && !isTRUE(meta$has_trials)) {
      return(fail("These are proportions, so the binomial families need the number of trials: choose the column with the number of trials (for example clutch size) under 'Weights' below the error family, or choose Gaussian."))
    }
    if (identical(family, "betabinomial") && binary && !isTRUE(meta$has_trials)) {
      return(fail("The beta-binomial family needs proportions with a trials column: with binary 0/1 data there is no extra-binomial variation to estimate, so use the binomial family."))
    }
  }
  covars <- meta$covars %||% character(0)
  rterms <- meta$random_terms %||% character(0)
  # Fewer than two AFR values: the appearance models cannot be told apart from Models 1-4.
  afr_levels <- {
    im_afr <- individual_metrics(dat)
    av_afr <- ifelse(is.finite(im_afr$entry), im_afr$entry, im_afr$first_recorded)
    length(unique(av_afr[is.finite(av_afr)]))
  }
  prep <- prepare_model_data(dat, age_function, covars, standardise, meta$has_group, rterms, has_group2 = isTRUE(meta$has_group2))
  if (is.null(prep)) return(fail("No usable trait data."))
  d <- prep$data
  b <- prep$basis
  cov_age <- intersect(meta$cov_age %||% character(0), covars)
  cov_pairs <- meta$cov_pairs %||% character(0)
  extra <- clean_extra_terms(extra, covars)
  fs <- apply_extra_terms(model_formula_strings(b, covars, among, cov_age, cov_pairs), b, extra, among)
  models <- intersect(models, names(fs))
  status <- list()
  av <- model_availability(dat, meta)
  for (m in models) {
    if (!av$ok[[m]]) status[[m]] <- paste("Unavailable:", av$why[[m]])
  }
  models <- setdiff(models, names(status))
  if (!length(models)) return(fail("None of the selected models can be fitted with these data.", status = status))

  vars <- unique(unlist(lapply(fs[models], function(s) all.vars(stats::as.formula(paste("trait ~", s))))))
  vars <- c(vars, "id", if (isTRUE(meta$has_group)) "group", if (isTRUE(meta$has_group2)) "group2", rterms,
            if (family %in% BINOMIAL_FAMILIES && isTRUE(meta$has_trials)) ".trials")
  ok <- rep(TRUE, nrow(d))
  drop_by <- character(0)
  for (v in vars) {
    x <- d[[v]]
    bad <- if (is.null(x)) rep(TRUE, nrow(d)) else if (identical(v, ".trials")) (!is.finite(x) | x <= 0) else if (is.numeric(x)) !is.finite(x) else is.na(x)
    if (any(bad & ok)) drop_by <- c(drop_by, sprintf("%s: %d", if (identical(v, ".trials")) "rows without a positive number of trials" else display_term(sub("_raw$", "", v)), sum(bad & ok)))
    ok <- ok & !bad
  }
  dd <- droplevels(d[ok, , drop = FALSE])
  if (nrow(dd) < 20 || length(unique(dd$id)) < 6) {
    return(fail("Too few complete observations shared by the selected models (need >= 20 rows from >= 6 individuals).", status = status))
  }
  n_age <- length(unique(dd$age))
  if (n_age < length(b) + 1) {
    return(fail(sprintf("The %s function needs at least %d distinct ages, but the analysed rows have %d. Choose a simpler ageing function.",
                        age_function, length(b) + 1, n_age), status = status))
  }
  # covariates without variation in the analysed rows would stop every model (e.g. contrasts error)
  keep_cv <- covars[vapply(covars, function(cv) {
    x <- dd[[cv]]
    if (is.numeric(x)) isTRUE(stats::sd(x) > 0) else length(unique(x)) >= 2
  }, logical(1))]
  for (cv in setdiff(covars, keep_cv)) {
    drop_by <- c(drop_by, sprintf("covariate %s omitted (no variation in the analysed rows)", meta$cov_labels[[cv]]))
  }
  if (!identical(keep_cv, covars)) {
    covars <- keep_cv
    cov_age <- intersect(cov_age, covars)
    extra <- clean_extra_terms(extra, covars)
    fs <- apply_extra_terms(model_formula_strings(b, covars, among, cov_age, cov_pairs), b, extra, among)
  }
  if (afr_levels < 2) {
    drop_by <- c(drop_by, "AFR takes fewer than two values, so selective appearance cannot be estimated: Models 7-10 are unavailable (their AFR terms would be dropped as rank deficient, repeating Models 1-4).")
  }
  if (!is.null(nonlinear_note)) drop_by <- c(drop_by, nonlinear_note)
  # zero-inflation predictors must exist and vary in the analysed rows
  zi_vars <- if (family %in% ZI_FAMILIES) tryCatch(all.vars(stats::as.formula(zi_str)), error = function(e) character(0)) else character(0)
  zi_keep <- zi_vars[vapply(zi_vars, function(v) {
    x <- dd[[v]]
    !is.null(x) && (if (is.numeric(x)) isTRUE(stats::sd(x) > 0) else length(unique(x)) >= 2)
  }, logical(1))]
  if (!identical(zi_keep, zi_vars)) {
    drop_by <- c(drop_by, sprintf("zero-inflation term %s omitted (missing or without variation)", paste(display_term(setdiff(zi_vars, zi_keep)), collapse = ", ")))
  }
  if (family %in% ZI_FAMILIES) zi_str <- if (length(zi_keep)) paste("~", paste(zi_keep, collapse = " + ")) else "~1"
  # random-effect structure for individuals (automatic choice from the data support)
  rs_request <- normalise_slope(random_slope)
  advice <- random_slope_advice(individual_data_support(dd))
  rs <- if (identical(rs_request, "auto")) advice$recommended else rs_request
  random_note <- if (identical(rs_request, "auto")) {
    paste0("Automatic random effects: ", slope_text(rs), ". ", advice$text)
  } else if (!identical(rs, "none") && identical(advice$level, "none")) {
    paste("Random slopes were requested although the data barely support them:", advice$facts)
  } else ""
  # additional random intercepts need at least two levels in the analysed rows
  keep_re <- rterms[vapply(rterms, function(rt) length(unique(dd[[rt]])) >= 2, logical(1))]
  for (rt in setdiff(rterms, keep_re)) {
    drop_by <- c(drop_by, sprintf("random intercept for %s omitted (fewer than 2 levels)", meta$random_labels[[rt]]))
  }
  rstr <- random_term_string(b[[1]], rs, meta$has_group, keep_re, has_group2 = isTRUE(meta$has_group2))
  fits <- list(); aic_rows <- list(); coef_rows <- list(); vc_rows <- list(); validity <- list()
  for (i in seq_along(models)) {
    m <- models[[i]]
    if (is.function(progress)) progress(i, length(models), model_label(m))
    res <- fit_one_model(fs[[m]], rstr, dd, family, zi_str)
    validity[[m]] <- res$validity
    if (is.null(res$fit)) {
      status[[m]] <- paste("Failed:", paste(res$notes, collapse = " | "))
      next
    }
    fits[[m]] <- res$fit
    status[[m]] <- if (length(res$notes)) paste0(res$validity, ": ", paste(res$notes, collapse = " | ")) else "Valid"
    ll <- stats::logLik(res$fit)
    aic_rows[[m]] <- data.frame(Model = model_label(m), AIC = tryCatch(stats::AIC(res$fit), error = function(e) NA_real_), df = attr(ll, "df"),
                                logLik = as.numeric(ll), N = stats::nobs(res$fit), Fit = res$validity, stringsAsFactors = FALSE)
    ct <- coef_table(res$fit, model_label(m))
    if (!is.null(ct)) coef_rows[[m]] <- ct
    vc <- variance_components(res$fit, model_label(m))
    if (!is.null(vc)) vc_rows[[m]] <- vc
  }
  aic <- if (length(aic_rows)) do.call(rbind, aic_rows) else data.frame()
  n_excluded <- 0L
  if (nrow(aic)) {
    # fits with an invalid Hessian (or non-finite AIC) are excluded from ranking unless the user overrides it
    aic$Eligible <- is.finite(aic$AIC) & (aic$Fit != "Failed" | isTRUE(include_invalid))
    n_excluded <- sum(!aic$Eligible)
    aic$Delta_AIC <- NA_real_
    aic$Weight <- NA_real_
    if (any(aic$Eligible)) {
      best <- min(aic$AIC[aic$Eligible])
      aic$Delta_AIC[aic$Eligible] <- aic$AIC[aic$Eligible] - best
      w <- exp(-aic$Delta_AIC[aic$Eligible] / 2)
      aic$Weight[aic$Eligible] <- w / sum(w)
    }
    aic <- aic[order(!aic$Eligible, aic$AIC), , drop = FALSE]
    rownames(aic) <- NULL
  }
  n_eligible <- if (nrow(aic)) sum(aic$Eligible) else 0L
  list(ok = n_eligible > 0,
       message = if (n_eligible > 0) {
         paste0("Models fitted by maximum likelihood on the same complete-case rows.",
                if (n_excluded > 0) sprintf(" %d fit(s) with an invalid Hessian were excluded from the ranking.", n_excluded) else "")
       } else if (nrow(aic)) {
         "Every fitted model has an invalid (non-positive-definite or singular) Hessian: simplify the random effects, zero-inflation or fixed effects, or tick 'Include fits with invalid Hessians' to inspect them."
       } else "No model could be fitted (see the fitting status).",
       validity = validity, n_excluded = n_excluded, include_invalid = isTRUE(include_invalid),
       random_request = rs_request, random_structure = rs, random_note = random_note, random_advice = advice,
       among = if (identical(among, "same") && length(b) > 1) "same" else "linear", cov_age = cov_age, extra = extra, nonlinear = FALSE,
       cov_pairs = cov_pairs[vapply(strsplit(cov_pairs, ":", fixed = TRUE), function(p) all(p %in% covars), logical(1))],
       fits = fits, aic = aic, coefficients = if (length(coef_rows)) do.call(rbind, coef_rows) else data.frame(),
       status = status, data = dd, n_dropped = sum(!ok), drop_by = drop_by, formulas = fs[models], random = rstr,
       basis = b, age_params = prep$age_params, proxy_params = prep$proxy_params, family = family, zi = zi_str,
       covars = covars, random_terms = keep_re,
       varcomp = if (length(vc_rows)) do.call(rbind, vc_rows) else data.frame(),
       age_function = age_function, random_slope = rs, standardise = isTRUE(standardise),
       lrt = lrt_table(fits, validity, include_invalid, fs))
}

lrt_table <- function(fits, validity = list(), include_invalid = FALSE, formulas = NULL) {
  pairs <- list(
    c("M1", "M2", "Age-independent selective disappearance (ALR main effect)"),
    c("M2", "M4", "Age-dependent selective disappearance (ALR \u00d7 ageing terms)"),
    c("M3", "M5", "Age-dependent selective disappearance (mean age \u00d7 \u0394age terms)"),
    c("M2", "M7", "Age-independent selective appearance (AFR main effect, given ALR)"),
    c("M4", "M8", "Age-dependent selective appearance (AFR \u00d7 ageing terms, given ALR \u00d7 ageing)"),
    c("M7", "M8", "Age-dependent disappearance + appearance vs additive"),
    c("M7", "M9", "Age-dependent disappearance (ALR \u00d7 age), given additive AFR"),
    c("M7", "M10", "Age-dependent appearance (AFR \u00d7 age), given additive ALR"),
    c("M4", "M9", "Age-independent appearance (AFR main effect), given ALR \u00d7 age"),
    c("M9", "M8", "Age-dependent appearance (AFR \u00d7 age), given ALR \u00d7 age"),
    c("M10", "M8", "Age-dependent disappearance (ALR \u00d7 age), given AFR \u00d7 age")
  )
  rows <- lapply(pairs, function(pr) {
    a <- fits[[pr[1]]]; b <- fits[[pr[2]]]
    if (is.null(a) || is.null(b)) return(NULL)
    va <- validity[[pr[1]]] %||% "Valid"
    vb <- validity[[pr[2]]] %||% "Valid"
    if (!isTRUE(include_invalid) && (identical(va, "Failed") || identical(vb, "Failed"))) return(NULL)
    if (!is.null(formulas) && !is.null(formulas[[pr[1]]]) && !is.null(formulas[[pr[2]]]) &&
        !formula_nested(formulas[[pr[1]]], formulas[[pr[2]]])) return(NULL)
    la <- stats::logLik(a); lb <- stats::logLik(b)
    ddf <- attr(lb, "df") - attr(la, "df")
    if (ddf <= 0) return(NULL)
    chi <- max(0, 2 * (as.numeric(lb) - as.numeric(la)))
    data.frame(Comparison = paste(model_label(pr[1]), "vs", model_label(pr[2])), Tests = pr[3],
               Chisq = chi, df = ddf, P_value = stats::pchisq(chi, ddf, lower.tail = FALSE),
               Fits = if (identical(va, "Valid") && identical(vb, "Valid")) "Valid" else paste(va, "/", vb), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Population-level trajectory from a fitted model (see the "i" help for the Models tab):
#   * random effects excluded (re.form = NA);
#   * numeric covariates held at their mean across individuals;
#   * factor covariates marginalised over the observed combinations of levels, weighted by the number of
#     individuals with each combination (unless `by` names a factor, which is then held at each level);
#   * ALR, LS, AFR and individual mean-age terms held at their individual-level means (or at `hold` values,
#     given on the model's standardised scale), with polynomial proxy terms computed from the held value;
#   * response scale (count models include the zero-inflation probability).
predict_population_curve <- function(fit, res, ages, hold = NULL, by = NULL) {
  d <- res$data
  ind <- d[!duplicated(d$id), , drop = FALSE]
  covs <- intersect(res$covars %||% character(0), names(d))
  fac <- covs[!vapply(d[covs], is.numeric, logical(1))]
  num <- setdiff(covs, fac)
  by <- if (length(by) == 1 && isTRUE(by %in% fac)) by else NULL
  other <- setdiff(fac, by)
  if (length(other)) {
    key <- do.call(paste, c(lapply(ind[other], as.character), sep = "\r"))
    wt <- as.numeric(table(key)[key[!duplicated(key)]])
    base <- ind[!duplicated(key), other, drop = FALSE]
  } else {
    base <- data.frame(.one = 1)
    wt <- 1
  }
  if (!is.null(by)) {
    lv <- levels(d[[by]])
    nb0 <- nrow(base)
    base <- base[rep(seq_len(nb0), times = length(lv)), , drop = FALSE]
    wt <- rep(wt, times = length(lv))
    base[[by]] <- factor(rep(lv, each = nb0), levels = lv)
    base$.level <- rep(lv, each = nb0)
  } else {
    base$.level <- "All"
  }
  nb <- nrow(base)
  nd <- base[rep(seq_len(nb), times = length(ages)), , drop = FALSE]
  nd$.w <- rep(wt, times = length(ages))
  nd$age <- rep(ages, each = nb)
  for (v in num) nd[[v]] <- mean(ind[[v]], na.rm = TRUE)
  bb <- age_basis(nd$age, res$age_params)
  for (nm in names(bb)) nd[[nm]] <- bb[[nm]]
  for (v in c("ALR", "LS", "AFR")) {
    val <- if (!is.null(hold[[v]])) hold[[v]] else if (v %in% names(ind)) mean(ind[[v]], na.rm = TRUE) else 0
    if (!is.finite(val)) val <- 0
    nd[[v]] <- val
    nd[[paste0(v, "2")]] <- val^2
    nd[[paste0(v, "3")]] <- val^3
  }
  for (nm in res$basis) {
    hk <- paste0("mean_", nm)
    mv <- if (!is.null(hold[[hk]])) hold[[hk]] else mean(ind[[hk]], na.rm = TRUE)
    nd[[hk]] <- mv
    nd[[paste0("delta_", nm)]] <- nd[[nm]] - mv
  }
  nd$id <- d$id[[1]]
  if ("group" %in% names(d)) nd$group <- d$group[[1]]
  if ("group2" %in% names(d)) nd$group2 <- d$group2[[1]]
  for (rt in intersect(res$random_terms %||% character(0), names(d))) nd[[rt]] <- d[[rt]][[1]]
  pr <- tryCatch({
    if (inherits(fit, "nlme")) as.numeric(stats::predict(fit, newdata = nd, level = 0))
    else as.numeric(stats::predict(fit, newdata = nd, re.form = NA, type = "response", allow.new.levels = TRUE))
  }, error = function(e) rep(NA_real_, nrow(nd)))
  ok <- is.finite(pr)
  if (!any(ok)) return(data.frame())
  grp <- paste(nd$.level[ok], format(nd$age[ok], digits = 15), sep = "\r")
  first <- !duplicated(grp)
  num_w <- tapply(pr[ok] * nd$.w[ok], grp, sum)
  den_w <- tapply(nd$.w[ok], grp, sum)
  out <- data.frame(level = nd$.level[ok][first], age = nd$age[ok][first],
                    fitted = as.numeric(num_w[grp[first]] / den_w[grp[first]]), stringsAsFactors = FALSE)
  out <- out[order(out$level, out$age), , drop = FALSE]
  rownames(out) <- NULL
  if (is.null(by)) out$level <- NULL
  out
}

# ---------------------------------------------------------------------------
# Plain-language interpretation of the selective disappearance / appearance terms of one model
# ---------------------------------------------------------------------------
# Uses Wald p-values of the terms, the nested likelihood-ratio test where available, and predicted
# population-level contrasts between individuals at the 10th and 90th percentiles of the proxy at a
# young, middle and old age (10th, 50th, 90th percentiles of the analysed ages).
interpret_model_terms <- function(res, m) {
  fit <- res$fits[[m]]
  if (is.null(fit)) return(character(0))
  ct <- res$coefficients[res$coefficients$Model == model_label(m), , drop = FALSE]
  if (!nrow(ct)) return("No coefficient table for this model.")
  # which selection terms does this model contain (standard models or extra terms)?
  term_parts <- strsplit(ct$Raw_term, ":", fixed = TRUE)
  has_term <- function(v) any(vapply(term_parts, function(pp) any(pp %in% c(v, paste0(v, 2:3), paste0("a.", v), paste0("b.", v))), logical(1)))
  has_mean <- any(grepl("^mean_f[1-3]$", ct$Raw_term) | grepl("(^|:)mean_f[1-3](:|$)", ct$Raw_term))
  if (!has_term("ALR") && !has_term("LS") && !has_term("AFR") && !has_mean) {
    return(paste0(model_label(m), " has no selective disappearance or appearance term: its trajectory describes the sampled records and is biased if individuals that disappear early differ from the others."))
  }
  count <- !identical(res$family, "gaussian")
  d <- res$data
  ind <- d[!duplicated(d$id), , drop = FALSE]
  ages <- as.numeric(stats::quantile(d$age, c(0.1, 0.5, 0.9), names = FALSE, type = 7))
  ages <- unique(ages)
  lrt_p <- function(cmp) {
    l <- res$lrt
    if (is.data.frame(l) && nrow(l) && cmp %in% l$Comparison) l$P_value[l$Comparison == cmp][[1]] else NA_real_
  }
  fmt_diff <- function(hi, lo) {
    if (count) {
      r <- hi / lo
      if (is.finite(r)) sprintf("\u00d7%s (%s)", format_num(r), if (r >= 1) sprintf("%.0f%% higher", 100 * (r - 1)) else sprintf("%.0f%% lower", 100 * (1 - r))) else "NA"
    } else {
      sprintf("%s%s", if (hi - lo >= 0) "+" else "\u2212", format_num(abs(hi - lo)))
    }
  }
  contrast <- function(label, var, hold_lo, hold_hi, lo_txt, hi_txt, main_terms, inter_terms, p_main_lrt, p_inter_lrt, lrt_names) {
    pm <- ct$P_value[ct$Raw_term %in% main_terms]
    pi <- ct$P_value[ct$Raw_term %in% inter_terms]
    p_main <- if (length(pm)) min(pm, na.rm = TRUE) else NA_real_
    p_int <- if (length(pi)) min(pi, na.rm = TRUE) else NA_real_
    lo <- predict_population_curve(fit, res, ages, hold = hold_lo)
    hi <- predict_population_curve(fit, res, ages, hold = hold_hi)
    if (!nrow(lo) || !nrow(hi)) return(paste0(label, ": predictions could not be computed."))
    dif <- if (count) hi$fitted / lo$fitted else hi$fitted - lo$fitted
    at <- paste(vapply(seq_along(ages), function(i) paste0(fmt_diff(hi$fitted[i], lo$fitted[i]), " at age ", format_num(ages[i])), character(1)), collapse = ", ")
    eff <- if (count) log(dif) else dif
    ok <- all(is.finite(eff))
    shape <- if (!ok || length(eff) < 2) "" else if (any(eff > 0) && any(eff < 0)) {
      "the direction of the difference reverses across ages"
    } else if (abs(eff[length(eff)]) > 1.5 * abs(eff[1]) && abs(eff[length(eff)]) > 0) {
      "the difference is larger at older ages"
    } else if (abs(eff[1]) > 1.5 * abs(eff[length(eff)]) && abs(eff[1]) > 0) {
      "the difference is larger at younger ages"
    } else {
      "the difference is similar in size across ages"
    }
    sig_int <- (is.finite(p_inter_lrt) && p_inter_lrt < 0.05) || (!is.finite(p_inter_lrt) && is.finite(p_int) && p_int < 0.05)
    sig_main <- (is.finite(p_main_lrt) && p_main_lrt < 0.05) || (is.finite(p_main) && p_main < 0.05)
    evidence <- paste0(
      if (length(main_terms)) sprintf("main-effect term(s) Wald p = %s", format_p(p_main)) else NULL,
      if (length(inter_terms)) sprintf("%sinteraction term(s) smallest Wald p = %s", if (length(main_terms)) "; " else "", format_p(p_int)) else NULL,
      if (is.finite(p_main_lrt)) sprintf("; LRT %s p = %s", lrt_names[[1]], format_p(p_main_lrt)) else "",
      if (is.finite(p_inter_lrt)) sprintf("; LRT %s p = %s", lrt_names[[2]], format_p(p_inter_lrt)) else "")
    what <- paste0("Predicted ", if (count) "ratio" else "difference", " between ", hi_txt, " and ", lo_txt, ": ", at, "; ", shape, ".")
    verdict <- if (length(inter_terms) && sig_int) {
      paste0(label, ": consistent with an age-dependent pattern \u2014 the gap between ", hi_txt, " and ", lo_txt, " changes with age (", evidence, ").")
    } else if (sig_main) {
      paste0(label, ": consistent with an age-independent pattern \u2014 ", hi_txt, " differ from ", lo_txt, " by a similar amount at all ages",
             if (length(inter_terms)) " (the interaction is not clearly supported)" else "", " (", evidence, ").")
    } else {
      paste0(label, ": no clear evidence for this term in this model (", evidence, "). Absence of evidence is not evidence of absence, especially with few long-lived individuals.")
    }
    c(verdict, what)
  }
  pct <- function(x) as.numeric(stats::quantile(x[is.finite(x)], c(0.1, 0.9), names = FALSE, type = 7))
  proxy_part <- function(v, label_hi, label_lo, lrt_main, lrt_inter, lrt_names) {
    raw <- ind[[paste0(v, "_raw")]]
    if (is.null(raw) || sum(is.finite(raw)) < 3) return(character(0))
    q <- pct(raw)
    cs <- res$proxy_params[[v]]
    z <- (q - cs[["centre"]]) / cs[["scale"]]
    main <- intersect(c(v, paste0(v, 2:3), paste0("a.", v)), ct$Raw_term)
    inter <- ct$Raw_term[vapply(strsplit(ct$Raw_term, ":", fixed = TRUE), function(pp) {
      length(pp) > 1 && any(pp %in% c(v, paste0(v, 2:3))) && any(grepl("^f[1-3]$", pp))
    }, logical(1)) | ct$Raw_term == paste0("b.", v)]
    hl <- stats::setNames(list(z[1]), v)
    hh <- stats::setNames(list(z[2]), v)
    contrast(paste(v, "term"), v, hl, hh,
             paste0(label_lo, " (", v, " = ", format_num(q[1]), ")"), paste0(label_hi, " (", v, " = ", format_num(q[2]), ")"),
             main, inter, lrt_main, lrt_inter, lrt_names)
  }
  out <- character(0)
  if (has_term("ALR")) {
    out <- c(out, proxy_part("ALR", "individuals with a late last record (longer-lived)", "individuals with an early last record (shorter-lived)",
                             if (m %in% c("M2", "M4", "M7", "M8", "M9", "M10")) lrt_p("Model 1 vs Model 2") else NA_real_,
                             if (m %in% c("M4", "M8")) lrt_p("Model 2 vs Model 4") else if (identical(m, "M9")) lrt_p("Model 7 vs Model 9") else NA_real_,
                             c("Model 1 vs 2", if (identical(m, "M9")) "Model 7 vs 9" else "Model 2 vs 4")))
  }
  if (has_term("LS")) {
    out <- c(out, proxy_part("LS", "longer-lived individuals", "shorter-lived individuals", NA_real_, NA_real_, c("", "")))
  }
  if (has_term("AFR")) {
    out <- c(out, proxy_part("AFR", "individuals first observed later", "individuals first observed earlier",
                             if (m %in% c("M7", "M8", "M9", "M10")) lrt_p("Model 2 vs Model 7") else NA_real_,
                             if (identical(m, "M8")) lrt_p("Model 4 vs Model 8") else if (identical(m, "M10")) lrt_p("Model 7 vs Model 10") else NA_real_,
                             c("Model 2 vs 7", if (identical(m, "M10")) "Model 7 vs 10" else "Model 4 vs 8")))
  }
  if (has_mean && !is.null(ind[["mean_f1"]])) {
    mf <- ind[["mean_f1"]]
    q <- pct(mf)
    near <- function(target) which.min(abs(mf - target))
    hold_for <- function(target) {
      j <- near(target)
      stats::setNames(lapply(res$basis, function(nm) ind[[paste0("mean_", nm)]][[j]]), paste0("mean_", res$basis))
    }
    main <- ct$Raw_term[grepl("^mean_f[1-3]$", ct$Raw_term)]
    inter <- ct$Raw_term[grepl("mean_f[1-3]", ct$Raw_term) & grepl("delta_f[1-3]", ct$Raw_term)]
    ap <- res$age_params
    to_age <- function(x) if (isTRUE(ap$standardise) && ap$fun %in% c("Linear", "Quadratic", "Cubic")) x * ap$sd_age + ap$mean_age else x
    out <- c(out, contrast("Mean-age term", "mean_f1", hold_for(q[1]), hold_for(q[2]),
                           paste0("individuals whose records centre on younger ages (mean age \u2248 ", format_num(to_age(q[1])), ")"),
                           paste0("individuals whose records centre on older ages (mean age \u2248 ", format_num(to_age(q[2])), ")"),
                           main, inter, NA_real_, if (identical(m, "M5")) lrt_p("Model 3 vs Model 5") else NA_real_, c("", "Model 3 vs 5")))
  }
  v <- res$validity[[m]] %||% "Valid"
  c(out,
    if (!identical(v, "Valid")) paste0("This fit is classified '", v, "': treat these estimates with caution.") else NULL,
    "Caveats: these are associations conditional on the model. ALR and mean age are proxies for lifespan that also depend on sampling; p-values are approximate; conclusions should agree with the visual diagnosis plots, the residual checks and alternative random-effect structures and error families.")
}

# ---------------------------------------------------------------------------
# performance package checks (optional). Every check, including the formatting of its result, runs inside its
# own error handler and time limit, so a check that fails or runs too long is reported instead of stopping the app.
# ---------------------------------------------------------------------------
with_time_limit <- function(expr, seconds) {
  setTimeLimit(elapsed = seconds, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
  expr
}

run_performance <- function(fit, family, time_limit = 60) {
  if (!requireNamespace("performance", quietly = TRUE)) {
    return(list(ok = FALSE, message = "Install the performance package (install.packages('performance')) for these checks."))
  }
  tryCatch(run_performance_checks(fit, family, time_limit),
           error = function(e) list(ok = FALSE, message = paste("The performance checks could not be run for this model:", conditionMessage(e))))
}

run_performance_checks <- function(fit, family, time_limit = 60) {
  num1 <- function(x) {
    x <- suppressWarnings(as.numeric(unlist(x)))
    x <- x[is.finite(x)]
    if (length(x)) x[[1]] else NA_real_
  }
  element <- function(x, name) if (is.list(x) && !is.null(x[[name]])) x[[name]] else NULL
  explain <- function(e) {
    msg <- conditionMessage(e)
    if (grepl("time limit", msg, ignore.case = TRUE)) paste0("stopped after ", time_limit, " s (too slow for this model)") else paste("not available for this model:", msg)
  }
  rows <- list()
  run_check <- function(label, note, code) {
    res <- tryCatch(with_time_limit(suppressWarnings(suppressMessages(code)), time_limit), error = explain)
    if (!is.character(res) || length(res) != 1 || is.na(res) || !nzchar(res)) res <- "not estimable"
    rows[[length(rows) + 1]] <<- data.frame(Check = label, Result = res, Note = note, stringsAsFactors = FALSE)
  }
  run_check("R\u00b2 (marginal / conditional)",
            "Nakagawa R\u00b2: variance explained by fixed effects (marginal) and by fixed plus random effects (conditional).", {
    r2 <- performance::r2(fit)
    vals <- suppressWarnings(as.numeric(unlist(r2)))
    nm <- names(unlist(r2))
    ok <- is.finite(vals)
    if (any(ok) && length(nm) == length(vals)) paste(paste0(sub("^R2_", "", nm[ok]), " = ", sprintf("%.3f", vals[ok])), collapse = "; ") else "not estimable"
  })
  run_check("ICC (adjusted)", "Share of variance due to the random effects.", {
    v <- num1(element(performance::icc(fit), "ICC_adjusted"))
    if (is.finite(v)) format_num(v) else "not estimable"
  })
  run_check("Singular fit", "Singular: a random-effect variance is estimated at (or near) zero; simplify the random effects.", {
    if (isTRUE(performance::check_singularity(fit))) "yes" else "no"
  })
  run_check("Convergence", "", {
    if (isTRUE(performance::check_convergence(fit))) "converged" else "possible convergence problem"
  })
  run_check("Collinearity (largest VIF)", "Polynomial and interaction terms are collinear by construction; high VIFs for them are expected and not a fault.", {
    coll <- performance::check_collinearity(fit)
    vif <- if (is.data.frame(coll) && "VIF" %in% names(coll)) suppressWarnings(as.numeric(coll$VIF)) else numeric(0)
    if (any(is.finite(vif))) {
      j <- which(vif == max(vif[is.finite(vif)]))[[1]]
      sprintf("%s (%s)", format_num(vif[[j]]), if ("Term" %in% names(coll)) as.character(coll$Term[[j]]) else "")
    } else "not estimable"
  })
  if (!identical(family, "gaussian")) {
    run_check("Overdispersion", if (family %in% BINOMIAL_FAMILIES) "Ratios well above 1 suggest a beta-binomial family (proportions) or a missing random effect." else "Ratios well above 1 suggest a negative binomial family.", {
      od <- performance::check_overdispersion(fit)
      ratio <- num1(element(od, "dispersion_ratio"))
      if (is.finite(ratio)) sprintf("dispersion ratio %s, p = %s", format_num(ratio), format_p(num1(element(od, "p_value")))) else "not estimable"
    })
    if (family %in% COUNT_FAMILIES) run_check("Zero inflation", "Ratios below 1 mean the model predicts fewer zeros than observed.", {
      zi <- performance::check_zeroinflation(fit)
      obs <- num1(element(zi, "observed.zeros"))
      if (is.finite(obs)) sprintf("observed %s vs predicted %s zeros (ratio %s)", format_num(obs), format_num(num1(element(zi, "predicted.zeros"))),
                                  format_num(num1(element(zi, "ratio")))) else "not estimable"
    })
  }
  cm <- NULL
  plot_message <- ""
  if (requireNamespace("see", quietly = TRUE)) {
    # zero-inflated and other glmmTMB count models: skip the panels that need many simulations or refits
    checks <- if (inherits(fit, "glmmTMB")) c("qq", "reqq", "homogeneity", "linearity", "vif") else "all"
    cm <- tryCatch(with_time_limit(suppressWarnings(suppressMessages(performance::check_model(fit, check = checks, verbose = FALSE))), time_limit),
                   error = function(e) e)
    if (inherits(cm, "error")) {
      plot_message <- paste("performance::check_model() panels are not available for this model:", explain(cm))
      cm <- NULL
    }
  } else {
    plot_message <- "Install the see package (install.packages('see')) to draw performance::check_model() panels."
  }
  list(ok = TRUE, table = do.call(rbind, rows), check_model = cm, plot_message = plot_message)
}

# ---------------------------------------------------------------------------
# A4-A7: disappearance diagnostics on person-occasion data
# ---------------------------------------------------------------------------
# One row per individual x sampled age (with a trait value). "Disappears" = dies before the next sampling
# occasion when lifespan is known; otherwise the record is the individual's last one. Records of censored
# individuals at their last record, and last records at the oldest sampled age of the study (their fate is
# unknown), are excluded. Trait is standardised within each age (z); change = trait minus the individual's
# record one sampling step earlier.
disappearance_data <- function(dat, meta, log_scale = FALSE) {
  ia <- id_age_means(dat)
  if (nrow(ia) < 10) return(NULL)
  if (isTRUE(log_scale)) ia$trait <- log1p(pmax(ia$trait, 0))
  step <- infer_age_step(ia$age, ia$id)
  if (!is.finite(step) || step <= 0) step <- 1
  im <- individual_metrics(dat)
  mi <- match(ia$id, im$id)
  life <- im$lifespan[mi]
  last <- stats::ave(ia$age, ia$id, FUN = max)
  cens <- ia$id %in% (meta$censored_ids %||% character(0))
  known <- isTRUE(meta$has_life) && !isTRUE(meta$life_auto) && sum(is.finite(im$lifespan)) >= 5
  tol <- 1e-8
  if (known) {
    event <- ifelse(is.finite(life), as.numeric(life < ia$age + step - tol), NA_real_)
    event[!is.finite(life) & ia$age < last] <- 0
  } else {
    event <- as.numeric(abs(ia$age - last) < tol)
    event[event == 1 & abs(last - max(ia$age)) < tol] <- NA_real_
  }
  event[cens & abs(ia$age - last) < tol] <- NA_real_
  ia$event <- event
  ia <- ia[order(ia$id, ia$age), , drop = FALSE]
  n <- nrow(ia)
  prev_same <- c(FALSE, ia$id[-1] == ia$id[-n] & abs(diff(ia$age) - step) <= 0.01 * step)
  ia$change <- NA_real_
  ia$change[prev_same] <- ia$trait[prev_same] - ia$trait[which(prev_same) - 1]
  zs <- stats::ave(ia$trait, ia$age, FUN = function(v) {
    s <- stats::sd(v)
    if (length(v) >= 3 && is.finite(s) && s > 0) (v - mean(v)) / s else rep(NA_real_, length(v))
  })
  ia$trait_z <- zs
  ia$step <- step
  attr(ia, "lifespan_known") <- known
  ia
}

# Discrete-time hazard: logit P(disappear before next occasion) = age-specific baseline + beta * trait_z
# (+ trait_z x age, + change). Age classes are the sampled ages (or 10 quantile classes if there are more
# than 30); classes with no events or only events carry no information and are removed.
disappearance_models <- function(ia) {
  if (is.null(ia)) return(list(ok = FALSE, message = "Too few records."))
  x <- ia[is.finite(ia$event) & is.finite(ia$trait_z), , drop = FALSE]
  if (nrow(x) < 30 || sum(x$event) < 5 || sum(1 - x$event) < 5) {
    return(list(ok = FALSE, message = sprintf("Too few records or disappearances (%d records, %d disappearances; need 30 and 5).", nrow(x), sum(x$event))))
  }
  ua <- sort(unique(x$age))
  x$age_class <- if (length(ua) <= 30) factor(x$age) else cut(x$age, unique(stats::quantile(x$age, seq(0, 1, 0.1))), include.lowest = TRUE)
  inf <- tapply(x$event, x$age_class, function(e) any(e == 1) && any(e == 0))
  x <- droplevels(x[x$age_class %in% names(inf)[inf %in% TRUE], , drop = FALSE])
  if (nrow(x) < 30 || nlevels(x$age_class) < 1) return(list(ok = FALSE, message = "No age class has both survivors and disappearances."))
  x$age_c <- (x$age - mean(x$age)) / (if (stats::sd(x$age) > 0) stats::sd(x$age) else 1)
  base_rhs <- if (nlevels(x$age_class) > 1) "age_class" else "1"
  glm_safe <- function(f, data) tryCatch(suppressWarnings(stats::glm(stats::as.formula(f), data = data, family = stats::binomial())), error = function(e) NULL)
  row_for <- function(fit, term, label, n, ev) {
    if (is.null(fit)) return(NULL)
    cf <- summary(fit)$coefficients
    if (!term %in% rownames(cf)) return(NULL)
    est <- cf[term, 1]; se <- cf[term, 2]
    data.frame(Model = label, Term = term, Odds_ratio = exp(est), Lower_95 = exp(est - 1.96 * se), Upper_95 = exp(est + 1.96 * se),
               P = cf[term, 4], Records = n, Disappearances = ev, stringsAsFactors = FALSE)
  }
  m1 <- glm_safe(paste("event ~", base_rhs, "+ trait_z"), x)
  m2 <- glm_safe(paste("event ~", base_rhs, "+ trait_z + trait_z:age_c"), x)
  xc <- x[is.finite(x$change), , drop = FALSE]
  m3 <- NULL
  if (nrow(xc) >= 30 && sum(xc$event) >= 5 && stats::sd(xc$change) > 0) {
    xc$change_z <- (xc$change - mean(xc$change)) / stats::sd(xc$change)
    xc <- droplevels(xc)
    m3 <- glm_safe(paste("event ~", if (nlevels(xc$age_class) > 1) "age_class" else "1", "+ trait_z + change_z"), xc)
  }
  tab <- rbind(row_for(m1, "trait_z", "Trait now", nrow(x), sum(x$event)),
               row_for(m2, "trait_z:age_c", "Trait effect changing with age (per SD of age)", nrow(x), sum(x$event)),
               row_for(m3, "trait_z", "Trait now, adjusting for change", nrow(xc), sum(xc$event)),
               row_for(m3, "change_z", "Change since the previous occasion", nrow(xc), sum(xc$event)))
  p_age <- if (!is.null(m1) && !is.null(m2)) stats::pchisq(max(0, m1$deviance - m2$deviance), 1, lower.tail = FALSE) else NA_real_
  ages <- as.numeric(stats::quantile(x$age, c(0.1, 0.5, 0.9), names = FALSE))
  or_at <- if (!is.null(m2)) {
    cf <- stats::coef(m2)
    vc <- stats::vcov(m2)
    ac <- (ages - mean(x$age)) / (if (stats::sd(x$age) > 0) stats::sd(x$age) else 1)
    data.frame(Age = ages, Odds_ratio = exp(cf[["trait_z"]] + cf[["trait_z:age_c"]] * ac),
               Lower_95 = exp(cf[["trait_z"]] + cf[["trait_z:age_c"]] * ac - 1.96 * sqrt(vc["trait_z", "trait_z"] + ac^2 * vc["trait_z:age_c", "trait_z:age_c"] + 2 * ac * vc["trait_z", "trait_z:age_c"])),
               Upper_95 = exp(cf[["trait_z"]] + cf[["trait_z:age_c"]] * ac + 1.96 * sqrt(vc["trait_z", "trait_z"] + ac^2 * vc["trait_z:age_c", "trait_z:age_c"] + 2 * ac * vc["trait_z", "trait_z:age_c"])))
  } else data.frame()
  # empirical disappearance rates by within-age trait quintile, in three age bands
  x$band <- cut(x$age, unique(stats::quantile(x$age, c(0, 1/3, 2/3, 1))), include.lowest = TRUE)
  x$tq <- cut(x$trait_z, unique(stats::quantile(x$trait_z, seq(0, 1, 0.2))), include.lowest = TRUE)
  emp <- stats::aggregate(cbind(event, trait_z) ~ band + tq, data = x, FUN = mean)
  emp$n <- stats::aggregate(event ~ band + tq, data = x, FUN = length)$event
  list(ok = !is.null(tab) && nrow(tab) > 0, table = tab, p_age_dependence = p_age, or_at_ages = or_at, empirical = emp,
       lifespan_known = isTRUE(attr(ia, "lifespan_known")), message = "")
}

# A6: selection differentials at each age, in within-age SD units. S = mean(survivors) - mean(all at that age);
# contrast = mean(survivors) - mean(disappearing). Ages need >= 3 survivors and >= 3 disappearances.
selection_differentials <- function(ia) {
  if (is.null(ia)) return(data.frame())
  x <- ia[is.finite(ia$event) & is.finite(ia$trait), , drop = FALSE]
  rows <- lapply(split(x, x$age), function(g) {
    s <- g$trait[g$event == 0]; d <- g$trait[g$event == 1]
    if (length(s) < 3 || length(d) < 3) return(NULL)
    sdv <- stats::sd(g$trait)
    if (!is.finite(sdv) || sdv <= 0) return(NULL)
    p <- length(d) / nrow(g)
    diff_sd <- (mean(s) - mean(d)) / sdv
    se_diff <- sqrt(stats::var(s) / length(s) + stats::var(d) / length(d)) / sdv
    data.frame(Age = g$age[[1]], N = nrow(g), Disappearing = length(d), Proportion_disappearing = p,
               Differential = p * diff_sd, Differential_SE = p * se_diff, Contrast = diff_sd, Contrast_SE = se_diff,
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$Age), , drop = FALSE]
}

# A7: life table of disappearance. Known lifespan: at risk at age a = individuals first observed at or before a
# with LS >= a (censored individuals stop contributing after their last record); events = LS in [a, a + step).
# Unknown lifespan: ALR replaces LS and events at the oldest sampled age are not counted.
life_table <- function(dat, meta) {
  im <- individual_metrics(dat)
  im <- im[is.finite(im$first_recorded), , drop = FALSE]
  if (nrow(im) < 5) return(list(ok = FALSE, message = "Too few individuals."))
  step <- infer_age_step(dat$age, dat$id)
  if (!is.finite(step) || step <= 0) step <- 1
  known <- isTRUE(meta$has_life) && !isTRUE(meta$life_auto) && sum(is.finite(im$lifespan)) >= 5
  cens <- im$id %in% (meta$censored_ids %||% character(0))
  start <- ifelse(is.finite(im$entry), pmin(im$entry, im$first_recorded), im$first_recorded)
  end <- if (known) ifelse(is.finite(im$lifespan), im$lifespan, im$last_recorded) else im$last_recorded
  dead <- if (known) is.finite(im$lifespan) & !cens else !cens & im$last_recorded < max(dat$age) - 1e-8
  ages <- seq(min(start), max(end), by = step)
  rows <- lapply(ages, function(a) {
    at_risk <- sum(start <= a + 1e-8 & end >= a - 1e-8)
    events <- sum(dead & end >= a - 1e-8 & end < a + step - 1e-8)
    if (at_risk == 0) return(NULL)
    q <- events / at_risk
    z <- 1.96
    centre <- (q + z^2 / (2 * at_risk)) / (1 + z^2 / at_risk)
    half <- z * sqrt(q * (1 - q) / at_risk + z^2 / (4 * at_risk^2)) / (1 + z^2 / at_risk)
    data.frame(Age = a, At_risk = at_risk, Disappearances = events, Hazard = q, Lower_95 = max(0, centre - half), Upper_95 = min(1, centre + half))
  })
  lt <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(lt) || !nrow(lt) || sum(lt$Disappearances) == 0) return(list(ok = FALSE, message = "No disappearances to estimate a hazard."))
  overall <- sum(lt$Disappearances) / sum(lt$At_risk)
  lt_fit <- lt[lt$At_risk > 0, , drop = FALSE]
  p_const <- tryCatch({
    m0 <- suppressWarnings(stats::glm(cbind(Disappearances, At_risk - Disappearances) ~ 1, data = lt_fit, family = stats::binomial()))
    m1 <- suppressWarnings(stats::glm(cbind(Disappearances, At_risk - Disappearances) ~ factor(Age), data = lt_fit, family = stats::binomial()))
    stats::pchisq(max(0, m0$deviance - m1$deviance), df = max(1, m0$df.residual - m1$df.residual), lower.tail = FALSE)
  }, error = function(e) NA_real_)
  list(ok = TRUE, table = lt, overall = overall, p_constant = p_const, lifespan_known = known, step = step)
}

# A5: trait against occasions before death (known LS) or before the last record, by lifespan tercile.
terminal_data <- function(dat, meta, n_min = 3) {
  ia <- id_age_means(dat)
  if (nrow(ia) < 10) return(data.frame())
  # deviation from the mean of all individuals at the same age, so that the population age trend does not
  # masquerade as a terminal change (occasions further before death are also younger ages)
  ia$resid <- ia$trait - stats::ave(ia$trait, ia$age, FUN = mean)
  step <- infer_age_step(ia$age, ia$id)
  if (!is.finite(step) || step <= 0) step <- 1
  im <- individual_metrics(dat)
  known <- isTRUE(meta$has_life) && !isTRUE(meta$life_auto) && sum(is.finite(im$lifespan)) >= 5
  cens <- im$id %in% (meta$censored_ids %||% character(0))
  endv <- if (known) im$lifespan else im$last_recorded
  endv[cens] <- NA_real_
  mi <- match(ia$id, im$id)
  ia$end <- endv[mi]
  ia <- ia[is.finite(ia$end), , drop = FALSE]
  if (nrow(ia) < 10) return(data.frame())
  ia$before <- floor((ia$end - ia$age) / step + 1e-8)
  ia <- ia[ia$before >= 0, , drop = FALSE]
  if (nrow(ia) < 3) return(data.frame())
  ends <- endv[is.finite(endv)]
  br <- unique(stats::quantile(ends, c(0, 1/3, 2/3, 1), names = FALSE))
  if (length(br) < 2) return(data.frame())
  ia$end_bin <- cut(ia$end, br, include.lowest = TRUE, labels = make.unique(paste0(format_num(br[-length(br)]), "\u2013", format_num(br[-1]))))
  z <- stats::aggregate(cbind(trait, resid) ~ end_bin + before, data = ia, FUN = mean)
  z$n <- stats::aggregate(trait ~ end_bin + before, data = ia, FUN = length)$trait
  z <- z[z$n >= n_min, , drop = FALSE]
  attr(z, "lifespan_known") <- known
  z
}

# ---------------------------------------------------------------------------
# Internal consistency of the model evidence
# ---------------------------------------------------------------------------
consistency_notes <- function(r, deviation = NULL, a2_trend = NULL) {
  out <- character(0)
  if (is.null(r) || !isTRUE(r$ok)) return(out)
  aic <- r$aic[r$aic$Eligible, , drop = FALSE]
  if (!nrow(aic)) return(out)
  best <- aic$Model[[1]]
  if (is.data.frame(deviation) && nrow(deviation)) {
    dm <- deviation[grepl("^Model", deviation$Method), , drop = FALSE]
    if (nrow(dm) >= 2 && best %in% dm$Method) {
      dbest <- dm$Method[which.min(dm$Mean_abs_D)]
      gap <- dm$Mean_abs_D[dm$Method == best] - min(dm$Mean_abs_D)
      if (!identical(dbest, best) && gap > max(2, 0.25 * min(dm$Mean_abs_D))) {
        out <- c(out, sprintf("AIC favours %s, but %s recovers the simulated population trajectory best (mean |D| %.1f%% vs %.1f%%). AIC rewards fit to the sampled records (including random effects, zero inflation and the error family), whereas deviation measures how well the population-level curve matches the simulated typical individual. They can disagree when late ages are sparsely sampled (e.g. missing when old or few occasions), when a model fits the observed range well but extrapolates poorly, when holding ALR or AFR at their means describes a different 'average' individual than the simulation, or when selective appearance and disappearance are confounded. Report both and favour conclusions that agree across AIC, likelihood-ratio tests and the visual diagnostics.",
                              best, dbest, dm$Mean_abs_D[dm$Method == best], min(dm$Mean_abs_D)))
      }
    }
  }
  l <- r$lrt
  if (is.data.frame(l) && nrow(l)) {
    p24 <- if ("Model 2 vs Model 4" %in% l$Comparison) l$P_value[l$Comparison == "Model 2 vs Model 4"][[1]] else NA_real_
    if (is.finite(p24)) {
      if (best %in% c("Model 4", "Model 8", "Model 9") && p24 >= 0.05) {
        out <- c(out, sprintf("%s has the lowest AIC, but the likelihood-ratio test of the ALR \u00d7 age terms (Model 2 vs 4) is not significant (p = %s): AIC penalises extra parameters less than a 5%% test, so the evidence for age-dependent selective disappearance is weak.", best, format_p(p24)))
      }
      if (best %in% c("Model 1", "Model 2") && p24 < 0.05) {
        out <- c(out, sprintf("%s has the lowest AIC, but the likelihood-ratio test of the ALR \u00d7 age terms is significant (p = %s): the two criteria disagree, so report both.", best, format_p(p24)))
      }
    }
  }
  best_id <- sub("^Model ", "M", best)
  if (identical(r$validity[[best_id]], "Caution")) {
    out <- c(out, sprintf("The lowest-AIC model (%s) is classified 'Caution' (%s): check whether a simpler random-effect structure or family gives the same conclusion.",
                          best, sub("^Caution: ", "", r$status[[best_id]] %||% "")))
  }
  inter_best <- best %in% c("Model 4", "Model 5", "Model 6", "Model 8", "Model 9", "Model 10")
  if (inter_best && identical(r$random_structure %||% "none", "none") && identical(r$random_advice$level %||% "none", "good")) {
    out <- c(out, sprintf("%s (with age-dependent terms) is favoured while only random intercepts are fitted, although the data support random slopes: individual differences in ageing rates can make interaction models look better without selective disappearance. Refit with a random slope and compare.", best))
  }
  if (is.data.frame(a2_trend) && nrow(a2_trend) == 1 && is.finite(a2_trend$P[[1]])) {
    inter <- best %in% c("Model 4", "Model 5", "Model 6", "Model 8", "Model 9", "Model 10")
    if (inter && a2_trend$P[[1]] >= 0.05) {
      out <- c(out, sprintf("The best model has age-dependent terms, but the lifespan\u2013trait coefficient does not change clearly across age bins (p = %s): check the visual diagnosis plots and the data support at old ages.", format_p(a2_trend$P[[1]])))
    } else if (!inter && best %in% c("Model 1", "Model 2", "Model 3", "Model 7") && a2_trend$P[[1]] < 0.05) {
      out <- c(out, sprintf("The best model has no age-dependent terms, but the lifespan\u2013trait coefficient changes across age bins (p = %s): the visual diagnostic and the model comparison disagree.", format_p(a2_trend$P[[1]])))
    }
  }
  out
}

# Columns usable for subsetting: text columns with 2-50 levels or numeric columns with 2-10 distinct values.
subset_candidates <- function(df) {
  if (is.null(df) || !ncol(df)) return(character(0))
  names(df)[vapply(df, function(x) {
    v <- x[!is.na(x)]
    if (!length(v)) return(FALSE)
    nu <- length(unique(v))
    num <- suppressWarnings(as.numeric(as.character(v)))
    if (mean(is.finite(num)) >= 0.95) nu >= 2 && nu <= 10 else nu >= 2 && nu <= 50
  }, logical(1))]
}

prediction_ages <- function(age) {
  u <- sort(unique(age[is.finite(age)]))
  if (length(u) <= 30) u else seq(min(u), max(u), length.out = 80)
}

# Relativised deviation from the simulated truth (Methods Eq. 11), restricted to ages
# with at least min_n observed individuals and a non-negligible true value.
deviation_from_truth <- function(curves, truth, obs_n, min_n = 10) {
  if (!nrow(curves) || !nrow(obs_n)) return(data.frame())
  ok_ages <- obs_n$age[obs_n$n >= min_n]
  if (!length(ok_ages)) return(data.frame())
  z <- curves[curves$age >= min(ok_ages) & curves$age <= max(ok_ages), , drop = FALSE]
  if (!nrow(z)) return(data.frame())
  z$true <- toy_true_curve(truth, z$age)
  z <- z[abs(z$true) > 0.05 * max(abs(z$true)), , drop = FALSE]
  if (!nrow(z)) return(data.frame())
  z$D <- 100 * (z$fitted - z$true) / z$true
  out <- do.call(rbind, lapply(split(z, z$Method), function(x) {
    x <- x[order(x$age), , drop = FALSE]
    data.frame(Method = x$Method[[1]], Mean_abs_D = mean(abs(x$D)), Max_abs_D = max(abs(x$D)),
               D_at_oldest_age = x$D[nrow(x)], Oldest_age = x$age[nrow(x)], stringsAsFactors = FALSE)
  }))
  out[order(out$Mean_abs_D), , drop = FALSE]
}

run_dharma <- function(fit, family) {
  if (!HAS_DHARMA) return(list(ok = FALSE, message = "Install the DHARMa package for simulation-based residual checks."))
  sim <- tryCatch(DHARMa::simulateResiduals(fit, n = 250, plot = FALSE), error = function(e) e)
  if (inherits(sim, "error")) return(list(ok = FALSE, message = paste("DHARMa failed:", conditionMessage(sim))))
  p_of <- function(expr) tryCatch(expr$p.value, error = function(e) NA_real_)
  tab <- data.frame(
    Test = c("Uniformity (KS)", "Dispersion", "Zero inflation"),
    P_value = c(p_of(DHARMa::testUniformity(sim, plot = FALSE)),
                p_of(DHARMa::testDispersion(sim, plot = FALSE)),
                if (identical(family, "gaussian")) NA_real_ else p_of(DHARMa::testZeroInflation(sim, plot = FALSE))),
    stringsAsFactors = FALSE)
  list(ok = TRUE, sim = sim, table = tab)
}

# ---------------------------------------------------------------------------
# Reproducible R code for the fitted comparison
# ---------------------------------------------------------------------------
model_r_code <- function(res, meta, source_label = "your_data.csv") {
  if (isTRUE(res$nonlinear)) return(model_r_code_nonlinear(res, meta, source_label))
  map <- meta$map
  q <- function(x) paste0('"', x, '"')
  ap <- res$age_params
  pp <- res$proxy_params
  num <- function(x) format(x, digits = 10)
  lines <- c(
    "# Code generated by disappR 0.9.4 - reproduces the app's model comparison",
    if (identical(res$family, "gaussian")) "library(lme4)" else "library(glmmTMB)",
    paste0("dat <- read.csv(", q(source_label), ", stringsAsFactors = FALSE, check.names = FALSE, colClasses = c(", q(map$id), " = \"character\"))"),
    if (!is.null(meta$subset)) sprintf("dat <- dat[!is.na(dat[[%s]]) & as.character(dat[[%s]]) %%in%% c(%s), ]  # subset used in the app",
                                       q(meta$subset$var), q(meta$subset$var), paste(q(meta$subset$levels), collapse = ", ")) else NULL,
    paste0("dat$id <- trimws(as.character(dat[[", q(map$id), "]]))"),
    if (isTRUE(meta$has_group)) paste0("dat$group <- as.character(dat[[", q(map$group), "]])") else NULL,
    if (isTRUE(meta$has_group2)) paste0("dat$group2 <- as.character(dat[[", q(map$group2), "]])") else NULL,
    if (isTRUE(meta$has_group2) && isTRUE(meta$nested)) "dat$group <- ifelse(is.na(dat$group) | is.na(dat$group2), dat$group, paste(dat$group2, dat$group, sep = \"/\"))  # groups nested in top-level groups" else NULL,
    if (isTRUE(meta$has_group) && isTRUE(meta$nested)) "dat$id <- ifelse(is.na(dat$group), dat$id, paste(dat$group, dat$id, sep = \"/\"))  # nested: (1 | group) + (1 | group/ID)" else NULL,
    paste0("dat$age <- as.numeric(dat[[", q(map$age), "]])"),
    if (is.finite(meta$age_round %||% NA_real_)) sprintf("dat$age <- round(dat$age / %s) * %s  # ages rounded as in the app", num(meta$age_round), num(meta$age_round)) else NULL,
    paste0("dat$trait <- as.numeric(dat[[", q(map$trait), "]])"),
    if (res$family %in% BINOMIAL_FAMILIES && isTRUE(meta$has_trials)) paste0("dat$.trials <- as.numeric(dat[[", q(meta$trials_col), "]])  # number of binomial trials (prior weights)") else NULL,
    "dat <- dat[!is.na(dat$id) & is.finite(dat$age), ]",
    if (isTRUE(meta$alr_mapped)) paste0("dat$ALR_raw <- as.numeric(dat[[", q(map$alr), "]])") else "dat$ALR_raw <- ave(dat$age, dat$id, FUN = max)  # age at last record",
    if (isTRUE(meta$entry_mapped)) paste0("dat$AFR_raw <- as.numeric(dat[[", q(map$entry), "]])") else "dat$AFR_raw <- ave(dat$age, dat$id, FUN = min)  # age at first record",
    if (isTRUE(meta$has_life) && !isTRUE(meta$life_auto)) paste0("dat$LS_raw <- as.numeric(dat[[", q(map$life), "]])") else "dat$LS_raw <- NA_real_",
    if (isTRUE(meta$has_censor)) sprintf("dat$LS_raw[dat$id %%in%% dat$id[as.character(dat[[%s]]) %%in%% %s]] <- NA  # censored individuals: lifespan unknown",
                                         q(map$censor), q(map$censor_value)) else NULL,
    if (identical(meta$dup_action, "mean")) "# NOTE: the app averaged duplicate ID x age records; do the same here before fitting" else NULL,
    "dat <- dat[is.finite(dat$trait), ]",
    if (res$family %in% BINOMIAL_FAMILIES && isTRUE(meta$has_trials)) "dat <- dat[is.finite(dat$.trials) & dat$.trials > 0, ]  # rows without a positive number of trials are dropped, as in the app" else NULL,
    "# individual-level centring/scaling constants used by the app",
    sprintf("dat$ALR <- (dat$ALR_raw - %s) / %s", num(pp$ALR[["centre"]]), num(pp$ALR[["scale"]])),
    sprintf("dat$AFR <- (dat$AFR_raw - %s) / %s", num(pp$AFR[["centre"]]), num(pp$AFR[["scale"]])),
    sprintf("dat$LS  <- (dat$LS_raw - %s) / %s", num(pp$LS[["centre"]]), num(pp$LS[["scale"]])),
    if (identical(res$among, "same")) "for (v in c(\"ALR\", \"AFR\", \"LS\")) { dat[[paste0(v, 2)]] <- dat[[v]]^2; dat[[paste0(v, 3)]] <- dat[[v]]^3 }  # polynomial among-individual terms" else NULL
  )
  z <- sprintf("(dat$age - %s) / %s", num(ap$mean_age), num(ap$sd_age))
  basis_lines <- switch(ap$fun,
    Linear = sprintf("dat$f1 <- %s", if (ap$standardise) z else "dat$age"),
    Quadratic = c(sprintf("dat$f1 <- %s", if (ap$standardise) z else "dat$age"), "dat$f2 <- dat$f1^2"),
    Cubic = c(sprintf("dat$f1 <- %s", if (ap$standardise) z else "dat$age"), "dat$f2 <- dat$f1^2", "dat$f3 <- dat$f1^3"),
    Logarithmic = {
      lx <- if (ap$all_positive) "log(dat$age)" else sprintf("log(dat$age - %s + 1)", num(ap$min_age))
      sprintf("dat$f1 <- %s", if (ap$standardise) sprintf("(%s - %s) / %s", lx, num(ap$mean_log), num(ap$sd_log)) else lx)
    },
    "Asymptotic exponential" = sprintf("dat$f1 <- %s", if (ap$standardise) sprintf("(exp(-%s) - %s) / %s", z, num(ap$mean_exp), num(ap$sd_exp)) else sprintf("exp(-%s)", z))
  )
  lines <- c(lines, basis_lines)
  for (nm in res$basis) {
    lines <- c(lines, sprintf("dat$mean_%s <- ave(dat$%s, dat$id, FUN = mean); dat$delta_%s <- dat$%s - dat$mean_%s", nm, nm, nm, nm, nm))
  }
  for (cv in meta$covars) {
    lines <- c(lines, if (identical(meta$cov_types[[cv]] %||% "", "continuous")) {
      sprintf("dat$%s <- suppressWarnings(as.numeric(dat[[%s]]))  # continuous covariate", cv, q(meta$cov_labels[[cv]]))
    } else {
      sprintf("dat$%s <- factor(dat[[%s]])  # categorical covariate", cv, q(meta$cov_labels[[cv]]))
    })
  }
  if (isTRUE(meta$has_group)) lines <- c(lines, "dat$group <- factor(dat$group)")
  if (isTRUE(meta$has_group2)) lines <- c(lines, "dat$group2 <- factor(dat$group2)")
  for (rt in res$random_terms %||% character(0)) {
    lines <- c(lines, sprintf("dat$%s <- factor(dat[[%s]])", rt, q(meta$random_labels[[rt]])))
  }
  lines <- c(lines, "dat$id <- factor(dat$id)")
  vars <- unique(unlist(lapply(res$formulas, function(s) all.vars(stats::as.formula(paste("trait ~", s))))))
  vars <- c(vars, "id", if (isTRUE(meta$has_group)) "group", if (isTRUE(meta$has_group2)) "group2", res$random_terms %||% character(0))
  lines <- c(lines, paste0("dat <- dat[complete.cases(dat[, c(", paste(q(vars), collapse = ", "), ")]), ]  # common rows for AIC"))
  call_for <- function(m) {
    f <- paste("trait ~", res$formulas[[m]], "+", res$random)
    if (identical(res$family, "gaussian")) {
      sprintf("%s <- lmer(%s, data = dat, REML = FALSE, control = lmerControl(optimizer = \"bobyqa\"))", tolower(m), f)
    } else {
      fam <- switch(res$family, poisson = "poisson()", zip = "poisson()", nbinom1 = "nbinom1()", zinb1 = "nbinom1()",
                    binomial = "binomial()", betabinomial = "betabinomial()", "nbinom2()")
      zi <- if (res$family %in% ZI_FAMILIES) res$zi else "~0"
      wt <- if (res$family %in% BINOMIAL_FAMILIES && isTRUE(meta$has_trials)) ", weights = .trials" else ""
      sprintf("%s <- glmmTMB(%s, family = %s, ziformula = %s, data = dat%s)", tolower(m), f, fam, zi, wt)
    }
  }
  fitted <- names(res$fits)
  lines <- c(lines, vapply(fitted, call_for, character(1)),
             paste0("AIC(", paste(tolower(fitted), collapse = ", "), ")"),
             if (length(fitted)) plot_code_lines(res, meta) else NULL)
  paste(lines, collapse = "\n")
}

# R code for the population-level trajectory figure: the fitted models' predictions (computed as in the app),
# observed means and the decomposition. Appended to the exported model script.
plot_code_lines <- function(res, meta) {
  ap <- res$age_params
  num <- function(x) format(x, digits = 10)
  q <- function(x) paste0('"', x, '"')
  zf <- sprintf("(age - %s) / %s", num(ap$mean_age), num(ap$sd_age))
  f1 <- if (isTRUE(ap$standardise)) zf else "age"
  basis <- switch(ap$fun,
    Linear = sprintf("data.frame(f1 = %s)", f1),
    Quadratic = sprintf("{ f1 <- %s; data.frame(f1 = f1, f2 = f1^2) }", f1),
    Cubic = sprintf("{ f1 <- %s; data.frame(f1 = f1, f2 = f1^2, f3 = f1^3) }", f1),
    Logarithmic = {
      lx <- if (isTRUE(ap$all_positive)) "log(age)" else sprintf("log(age - %s + 1)", num(ap$min_age))
      sprintf("data.frame(f1 = %s)", if (isTRUE(ap$standardise)) sprintf("(%s - %s) / %s", lx, num(ap$mean_log), num(ap$sd_log)) else lx)
    },
    "Asymptotic exponential" = sprintf("data.frame(f1 = %s)", if (isTRUE(ap$standardise)) sprintf("(exp(-%s) - %s) / %s", zf, num(ap$mean_exp), num(ap$sd_exp)) else sprintf("exp(-%s)", zf)))
  fitted <- names(res$fits)
  covs <- res$covars %||% character(0)
  rterms <- res$random_terms %||% character(0)
  c("",
    "# ---------------------------------------------------------------------------",
    "# Figure: population-level trajectories of the fitted models, observed means and the decomposition",
    "# Predictions as in the app: random effects excluded; numeric covariates at their mean; factor covariates",
    "# averaged over the observed combinations of levels (weighted by individuals); ALR, AFR, LS and mean-age",
    "# terms at their individual-level means; response scale.",
    "# ---------------------------------------------------------------------------",
    "library(ggplot2)",
    paste0("make_basis <- function(age) ", basis, "  # the same age transformation as above"),
    paste0("fits <- list(", paste(sprintf("%s = %s", q(model_label(fitted)), tolower(fitted)), collapse = ", "), ")"),
    paste0("covariates <- c(", paste(q(covs), collapse = ", "), ")"),
    "population_curve <- function(fit, dat, ages) {",
    "  ind <- dat[!duplicated(dat$id), , drop = FALSE]",
    "  fac <- covariates[!vapply(ind[covariates], is.numeric, logical(1))]",
    "  numv <- setdiff(covariates, fac)",
    "  if (length(fac)) {",
    "    key <- do.call(paste, c(lapply(ind[fac], as.character), sep = \"\\r\"))",
    "    base <- ind[!duplicated(key), fac, drop = FALSE]",
    "    w <- as.numeric(table(key)[key[!duplicated(key)]])",
    "  } else {",
    "    base <- data.frame(.one = 1)",
    "    w <- 1",
    "  }",
    "  nd <- base[rep(seq_len(nrow(base)), times = length(ages)), , drop = FALSE]",
    "  nd$.w <- rep(w, times = length(ages))",
    "  nd$age <- rep(ages, each = nrow(base))",
    "  for (v in numv) nd[[v]] <- mean(ind[[v]], na.rm = TRUE)",
    "  b <- make_basis(nd$age)",
    "  for (nm in names(b)) {",
    "    nd[[nm]] <- b[[nm]]",
    "    mv <- if (paste0(\"mean_\", nm) %in% names(ind)) mean(ind[[paste0(\"mean_\", nm)]], na.rm = TRUE) else 0",
    "    nd[[paste0(\"mean_\", nm)]] <- mv",
    "    nd[[paste0(\"delta_\", nm)]] <- nd[[nm]] - mv",
    "  }",
    "  for (v in c(\"ALR\", \"AFR\", \"LS\")) {",
    "    val <- if (v %in% names(ind)) mean(ind[[v]], na.rm = TRUE) else 0",
    "    if (!is.finite(val)) val <- 0",
    "    nd[[v]] <- val; nd[[paste0(v, 2)]] <- val^2; nd[[paste0(v, 3)]] <- val^3",
    "  }",
    "  nd$id <- dat$id[[1]]",
    if (isTRUE(meta$has_group)) "  nd$group <- dat$group[[1]]" else NULL,
    if (isTRUE(meta$has_group2)) "  nd$group2 <- dat$group2[[1]]" else NULL,
    if (length(rterms)) paste0("  for (rt in c(", paste(q(rterms), collapse = ", "), ")) nd[[rt]] <- dat[[rt]][[1]]") else NULL,
    "  pr <- if (inherits(fit, \"nlme\")) {",
    "    predict(fit, newdata = nd, level = 0)",
    "  } else if (inherits(fit, \"glmmTMB\")) {",
    "    predict(fit, newdata = nd, re.form = NA, type = \"response\", allow.new.levels = TRUE)",
    "  } else {",
    "    predict(fit, newdata = nd, re.form = NA)",
    "  }",
    "  num_w <- tapply(as.numeric(pr) * nd$.w, nd$age, sum)",
    "  den_w <- tapply(nd$.w, nd$age, sum)",
    "  data.frame(age = as.numeric(names(num_w)), fitted = as.numeric(num_w / den_w))",
    "}",
    "ages <- sort(unique(dat$age))",
    "if (length(ages) > 30) ages <- seq(min(ages), max(ages), length.out = 80)",
    "curves <- do.call(rbind, lapply(names(fits), function(m) cbind(population_curve(fits[[m]], dat, ages), Model = m)))",
    "# observed means: mean over individuals of each individual's mean at that age",
    "ia <- aggregate(trait ~ id + age, data = dat, FUN = mean)",
    "ia$id <- as.character(ia$id)",
    "ia <- ia[order(ia$id, ia$age), ]",
    "observed <- aggregate(trait ~ age, data = ia, FUN = mean)",
    "observed$n <- as.numeric(table(ia$age)[as.character(observed$age)])",
    "# decomposition (Rebke et al. 2010): the app's own functions, copied here so that the figure matches the app",
    app_code_definitions(app_code_closure("decomposition_trajectory")),
    "decomp <- decomposition_trajectory(dat)",
    "ggplot() +",
    "  geom_point(data = observed, aes(age, trait, size = n), colour = \"grey40\", alpha = 0.5) +",
    "  geom_line(data = curves, aes(age, fitted, colour = Model), linewidth = 1.1) +",
    "  geom_line(data = decomp, aes(age, fitted), colour = \"purple\", linetype = 2, linewidth = 1.1) +",
    "  scale_size_area(max_size = 4, guide = \"none\") +",
    "  labs(x = \"Age\", y = \"Trait\", colour = NULL, subtitle = \"Points: observed means; dashed purple: decomposition\") +",
    "  theme_minimal(base_size = 13)")
}

model_r_code_nonlinear <- function(res, meta, source_label = "your_data.csv") {
  # build a copy explicitly: utils::modifyList() merges nested lists, so it would not empty fits or formulas
  r0 <- res
  r0$nonlinear <- FALSE
  r0$formulas <- list()
  r0$fits <- list()
  r0$basis <- "f1"
  r0$family <- "gaussian"
  r0$random_terms <- character(0)
  r0$among <- "linear"
  base <- model_r_code(r0, meta, source_label)
  base <- sub("library\\(lme4\\)", "library(nlme)", base)
  base <- sub("\nAIC\\(\\)$", "", base)
  lines <- strsplit(base, "\n", fixed = TRUE)[[1]]
  lines <- lines[!grepl("^dat <- dat\\[complete.cases", lines)]
  rand <- switch(res$random_structure %||% "none", correlated = "pdSymm(a + b ~ 1)", uncorrelated = "pdDiag(a + b ~ 1)", "pdSymm(a ~ 1)")
  grp <- if (isTRUE(meta$has_group2)) "~ group2 / group / id" else if (isTRUE(meta$has_group)) "~ group / id" else "~ id"
  rand_arg <- if (isTRUE(meta$has_group2)) {
    sprintf("list(group2 = pdSymm(a ~ 1), group = pdSymm(a ~ 1), id = %s)", rand)
  } else if (isTRUE(meta$has_group)) sprintf("list(group = pdSymm(a ~ 1), id = %s)", rand) else rand
  calls <- vapply(names(res$fits), function(m) {
    pa <- res$parts[[m]]$a
    pb <- res$parts[[m]]$b
    rhs <- function(x) if (length(x)) paste(x, collapse = " + ") else "1"
    sprintf("%s <- nlme(trait ~ a * exp(b * f1), data = dat, fixed = list(a ~ %s, b ~ %s), random = %s, groups = %s, start = c(%s), method = \"ML\", control = nlmeControl(returnObject = TRUE))",
            tolower(m), rhs(pa), rhs(pb), rand_arg, grp,
            paste(format(nlme::fixef(res$fits[[m]]), digits = 6), collapse = ", "))
  }, character(1))
  paste(c(lines, "# non-linear exponential models: trait = a * exp(b * z_age); starting values are the app's estimates",
          calls, paste0("AIC(", paste(tolower(names(res$fits)), collapse = ", "), ")"),
          if (length(res$fits)) plot_code_lines(res, meta) else NULL), collapse = "\n")
}

# Agreement between age at first record and age at last record. A strong positive association means
# individuals recorded late are also recorded late in life, so entry and exit are not independent and the
# observation window, rather than lifespan alone, drives who is seen at which age.
afr_alr_agreement <- function(dat, meta) {
  im <- individual_metrics(dat)
  im <- im[im$n_trait > 0, , drop = FALSE]
  afr <- ifelse(is.finite(im$entry), im$entry, im$first_recorded)
  alr <- im$alr
  ok <- is.finite(afr) & is.finite(alr)
  if (sum(ok) < 5 || stats::sd(afr[ok]) == 0 || stats::sd(alr[ok]) == 0) {
    return(list(ok = FALSE, message = if (sum(ok) < 5) "Fewer than 5 individuals with both AFR and ALR." else
      "AFR or ALR does not vary among individuals, so the two cannot be correlated."))
  }
  r <- stats::cor(afr[ok], alr[ok])
  rho <- suppressWarnings(stats::cor(afr[ok], alr[ok], method = "spearman"))
  span <- alr[ok] - afr[ok]
  list(ok = TRUE, n = sum(ok), r = r, rho = rho, mean_span = mean(span), sd_span = stats::sd(span),
       table = data.frame(
         Measure = c("Individuals", "Correlation AFR vs ALR (Pearson)", "Correlation AFR vs ALR (Spearman)",
                     "Observation window (ALR \u2212 AFR)"),
         Value = c(format(sum(ok), big.mark = ","), format_num(r), format_num(rho),
                   sprintf("mean %s, SD %s", format_num(mean(span)), format_num(stats::sd(span)))),
         stringsAsFactors = FALSE))
}

# ---------------------------------------------------------------------------
# B1: error-family check. AICs are comparable among count families (all are
# probability mass functions for the same counts) but NOT with a Gaussian density.
# ---------------------------------------------------------------------------
compare_families <- function(dat, meta, model = "M4", age_function = "Quadratic", random_slope = FALSE,
                             standardise = TRUE, zi_str = "~1", families = COUNT_FAMILIES, progress = NULL,
                             among = "linear", include_invalid = FALSE) {
  if (!HAS_GLMMTMB) return(list(ok = FALSE, message = "Install glmmTMB to compare count families."))
  v <- dat$trait[is.finite(dat$trait)]
  binom_set <- all(families %in% BINOMIAL_FAMILIES)
  if (binom_set) {
    if (any(v < 0) || any(v > 1)) return(list(ok = FALSE, message = "The trait is not binary or a proportion, so the binomial families do not apply."))
  } else if (any(v < 0) || any(abs(v - round(v)) > 1e-8)) {
    return(list(ok = FALSE, message = "The trait is not a non-negative integer count, so count families do not apply (use Gaussian)."))
  }
  rows <- list()
  for (i in seq_along(families)) {
    fam <- families[[i]]
    if (is.function(progress)) progress(i, length(families), family_label(fam))
    res <- fit_model_suite(dat, meta, models = model, age_function = age_function, family = fam,
                           random_slope = random_slope, standardise = standardise, zi_str = zi_str,
                           among = among, include_invalid = include_invalid)
    if (!isTRUE(res$ok) || !nrow(res$aic)) {
      rows[[fam]] <- data.frame(Family = family_label(fam), AIC = NA_real_, df = NA_real_, N = NA_real_, Fit = res$validity[[model]] %||% "Failed",
                                Status = paste(unlist(res$status), res$message, collapse = " "), stringsAsFactors = FALSE)
    } else {
      rows[[fam]] <- data.frame(Family = family_label(fam), AIC = res$aic$AIC[[1]], df = res$aic$df[[1]], N = res$aic$N[[1]], Fit = res$aic$Fit[[1]],
                                Status = unlist(res$status)[[1]], stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  if (any(is.finite(tab$AIC))) tab$Delta_AIC <- tab$AIC - min(tab$AIC, na.rm = TRUE) else tab$Delta_AIC <- NA_real_
  tab <- tab[order(tab$AIC, na.last = TRUE), c("Family", "AIC", "Delta_AIC", "df", "N", "Fit", "Status"), drop = FALSE]
  rownames(tab) <- NULL
  list(ok = TRUE, table = tab, model = model)
}

# ---------------------------------------------------------------------------
# B2: comparison of ageing functions for one fixed-effect structure (Model 1 on the B2 tab, any
# model in the within-model check), with the chosen family, covariates and random effects.
# All functions use the same rows, so AICs are comparable.
# ---------------------------------------------------------------------------
compare_population_functions <- function(dat, meta, family = "gaussian", random_slope = FALSE, standardise = TRUE,
                                         zi_str = "~1", progress = NULL, functions = AGE_FUNCTIONS, model = "M1",
                                         among = "linear", include_invalid = FALSE) {
  functions <- intersect(functions, MODEL_FUNCTIONS)
  rows <- list()
  curves <- list()
  for (i in seq_along(functions)) {
    fn <- functions[[i]]
    if (is.function(progress)) progress(i, length(functions), fn)
    res <- fit_model_suite(dat, meta, models = model, age_function = fn, family = family,
                           random_slope = random_slope, standardise = standardise, zi_str = zi_str,
                           among = among, include_invalid = include_invalid)
    if (!isTRUE(res$ok) || !nrow(res$aic)) {
      rows[[fn]] <- data.frame(Function = fn, AIC = NA_real_, df = NA_real_, N = NA_real_, Fit = res$validity[[model]] %||% "Failed",
                               Status = paste(unlist(res$status), res$message, collapse = " "), stringsAsFactors = FALSE)
      next
    }
    rows[[fn]] <- data.frame(Function = fn, AIC = res$aic$AIC[[1]], df = res$aic$df[[1]], N = res$aic$N[[1]], Fit = res$aic$Fit[[1]],
                             Status = unlist(res$status)[[1]], stringsAsFactors = FALSE)
    cv <- predict_population_curve(res$fits[[model]], res, prediction_ages(res$data$age))
    if (nrow(cv)) {
      cv$Function <- fn
      curves[[fn]] <- cv
    }
  }
  if (!length(rows)) return(list(table = data.frame(), curves = data.frame(), model = model, functions = functions))
  tab <- do.call(rbind, rows)
  tab$Delta_AIC <- if (any(is.finite(tab$AIC))) tab$AIC - min(tab$AIC, na.rm = TRUE) else NA_real_
  tab <- tab[order(tab$AIC, na.last = TRUE), c("Function", "AIC", "Delta_AIC", "df", "N", "Fit", "Status"), drop = FALSE]
  rownames(tab) <- NULL
  list(table = tab, curves = if (length(curves)) do.call(rbind, curves) else data.frame(), model = model, functions = functions)
}

# ---------------------------------------------------------------------------
# "i" help content: what each section does, how to use and read it, and cautions
# ---------------------------------------------------------------------------
info_entry <- function(title, what, how, caution) list(title = title, what = what, how = how, caution = caution)

INFO <- list(
  overview = info_entry("The workflow",
    "The tabs follow the manuscript's recommended workflow: diagnose selective disappearance and appearance visually (Visual diagnosis, including the disappearance diagnostics, and Individual and population trajectories), check sampling (Missingness and proxies), then fit and compare mixed models (Modelling).",
    "Work through the tabs in order. Learn on simulated data, where the answer is known, before analysing your own data. Use 'Save to summary' on any result you want in the exported report.",
    "Visual diagnostics suggest selection processes but cannot prove them. AIC measures relative support, not causation. Many tests are run across ages, proxies and models with no correction for multiple testing, so treat the diagnostics as exploratory: an association with ALR or AFR is a statistical pattern, compatible with selective disappearance but also with observation processes, cohort or environmental effects, and among-individual heterogeneity."),
  quickstart = info_entry("Quick start",
    "Loads a simulated teaching dataset (where the answer is known) or one of the bundled published datasets.",
    "Simulated data show textbook patterns with a known truth; the fly data are a real reanalysis with count data, covariates, nesting and censoring. Saved results are collected on the 'Summary and report' tab.",
    "Teaching simulations use deliberately strong effects."),
  data_source = info_entry("Data source and simulation",
    "Chooses simulated data (known truth), one of the bundled empirical examples, or your own CSV with one row per individual \u00d7 age. The empirical examples cover laboratory systems (fruit-fly and seed-beetle fecundity, leafcutting-bee locomotor activity) and wild populations (common tern immunity and navigation, painted turtle reproduction, great tit recruitment, eastern chipmunk reproduction, Soay sheep reproduction); each opens with the mapping, error family and models set close to the analysis reported in its paper. The simulator generates individual trajectories with a chosen ageing form (linear, quadratic, cubic, logarithmic or exponential) and a biologically motivated shape within that form, for example slow, fast or no senescence, an early or late peak followed by senescence, a decline followed by improvement, or a rapid or gradual decline (or growth) that levels off. The default shape of each form is the one used by earlier versions. Selective disappearance can be absent, age-independent (lifespan linked to the individual's level), age-dependent (lifespan linked to its rate of ageing) or both, with a positive or negative direction. With individual-specific AFR (age at first observation), selective appearance is set in the same way. Diet lowers the trait only; families add a nested random effect.",
    "Change one setting at a time, press 'Simulate & use this dataset', and compare the plots and model rankings with 'What you should see'. The mean lifespan sets how many occasions each individual is sampled. Use 'Subset the data' to analyse one group (for example one sex or treatment). For your own CSV, age must be numeric (not age classes or text), and missing values can be blank cells or NA.",
    "Large individual differences in ageing rates make interaction models look better even without selective disappearance unless random slopes are fitted. Very short lifespans make individual fits and interaction models unreliable. The exponential form is defined on age scaled by the SD of the expected ages, like the app's Exponential function, so it is recovered only approximately when missingness changes the sampled ages. The empirical examples are exploratory: the defaults broadly reproduce the published patterns, but they need not match a paper's numbers exactly, because published analyses differ in how the data were subset, in covariates and random effects that are not always fully reported, in software and estimation, and in decisions about which records to exclude. Treat any difference as a reason to inspect the settings rather than as a failure of either analysis."),
  afr_expression = info_entry("Counting missed occasions: from AFR or from AFE",
    "Sets where each individual's expected-occasion window opens. 'Age at first record (AFR)' is when the individual enters the data, which is how it is mapped on the Data tab. 'Age at first trait expression (AFE)' is the age from which the trait exists and could in principle have been measured, entered as one age that applies to every individual: body size is expressed from birth even if recording begins at first breeding.",
    "Use AFR when the trait genuinely begins when the individual enters the data (a first clutch, for example). Use AFE when individuals could have been measured earlier than they were, and give the age at which the trait starts.",
    "This setting changes the missingness figures only. It affects this tab's grid, missed-occasion counts and coverage; it does not enter any model, proxy or statistic anywhere in the app, all of which use AFR. Counting from AFE usually raises the estimated missingness, sometimes substantially, and changes its age profile at young ages, so the two definitions are not comparable: say which one you used when reporting."),
  afr_alr = info_entry("Agreement between AFR and ALR",
    "Correlation between each individual's age at first record and its age at last record, with the mean and SD of the observation window (ALR \u2212 AFR). The figure shows the dashed 1:1 line and the linear regression of ALR on AFR (solid black, with its 95% band).",
    "Under a fixed entry age, AFR does not vary and the correlation is undefined. A correlation near zero means entry and exit are independent, so ALR reflects lifespan rather than when an individual was first seen. A strong positive correlation means individuals that enter late also leave late, so the two proxies carry overlapping information and AFR terms (Models 7-10) and ALR terms (Models 2, 4) compete for the same variance.",
    "A high correlation can arise from biology (late starters live longer) or purely from the study design (a short study window forces late entrants to have late exits), and this diagnostic cannot separate the two. With few individuals the correlation is unstable. Mean age is affected by both ends of the window, so a variable AFR degrades it as a proxy of lifespan more than it degrades ALR."),
  among_order = info_entry("Among-individual terms: polynomial order",
    "Sets the functional form of the among-individual proxy (ALR, LS, AFR or mean age). 'Linear' enters the proxy as a single term. 'Same polynomial order as the ageing function' gives the proxy the same form as the age terms: with a quadratic ageing function the model gains LS and LS\u00b2 (or ALR and ALR\u00b2), and their interactions with each age term.",
    "Use the linear form unless there is a reason to expect a curved relationship between the proxy and the trait. Choose the matched order when the trait is not monotonic in lifespan, for example when intermediate-lifespan individuals have the highest trait values, or when a plot of the trait against ALR or lifespan is clearly U-shaped or humped. A linear ALR term cannot represent such a pattern and will report no selective disappearance even when it is strong.",
    "Matching the order adds parameters quickly (a quadratic ageing function doubles the proxy terms and their interactions), so it needs more individuals and a wider spread of lifespans. The extra terms are collinear with the linear ones by construction, so individual coefficients are hard to interpret: compare models by fit, and read the trajectories rather than the coefficients."),
  mapping = info_entry("Map columns",
    "Tells the app which columns hold ID, age and trait, plus optional ALR, lifespan, AFR (age at first observation), continuous and categorical fixed-effect covariates, their interactions, grouping, additional random intercepts, censoring and an optional age resolution. Uploaded columns keep their names.",
    "Continuous covariates (e.g. temperature) enter the models as linear effects; categorical covariates (e.g. treatment, sex, diet) enter as factors with one coefficient per level. Interactions can be a covariate \u00d7 age (the covariate changes the shape of the ageing trajectory, e.g. diet \u00d7 (age + age\u00b2)) or between two covariates (e.g. diet \u00d7 sex, added to every model with both main effects). Covariates and interactions are used in the mixed models on tab 5 (Modelling). Map lifespan (LS) only if it is truly known; use nesting when IDs repeat across groups; map a censoring column for individuals alive at the end of the study. With a 'next level up' column the random effects have three levels, e.g. individuals within fathers within families: (1 | family) + (1 | family:father) + (1 | family:father:ID). The CONTINUOUS and CATEGORICAL boxes warn straight away when a column looks as if it belongs in the other box.",
    "Values of a continuous covariate that are not numbers are treated as missing, which drops those rows from the models. Categorical covariates with many levels or rare combinations make interaction models hard to estimate. Additional random intercepts are crossed with ID unless nested. Age must be numeric (whole numbers or decimals); it cannot be categorical."),
  subset = info_entry("Subset the data",
    "Restricts every analysis to the rows whose value of one categorical variable (a text column with up to 50 levels, or a numeric column with up to 10 distinct values) is among the chosen levels.",
    "Pick the variable, then tick the levels to keep (for example females only, or one treatment). Leave the variable empty to use all rows. The subset applies to the integrity checks, all diagnostics, individual fits and models, and to the exported R code only through the data you supply.",
    "Subsetting reduces the number of individuals, so rare long-lived individuals may disappear from the analysis. Subsetting on a variable that is affected by survival (for example 'still alive') creates selection bias of its own. To compare groups within one analysis use panels in the visual diagnosis or covariates instead."),
  integrity = info_entry("Data integrity checks",
    "Flags problems that change results: duplicate records, IDs in several groups, missing or impossible lifespans, censoring, the sampling schedule (step, off-schedule rows, staggered or irregular schedules), few records per individual, AFR variation, AFR\u2013ALR correlation, missing-value codes (-99, -999, -9999), ages at or below zero, covariates without variation, with many levels or with missing values, empty columns, extreme trait values and the trait distribution.",
    "'Warning' rows need action (fix the data or the mapping); 'Note' rows are informative. The sampling step is the most common interval between an individual's consecutive records, so staggered cohorts and floating-point ages are handled.",
    "The checks cannot find every recording error. The suggested family is only a starting point: confirm it with the error-family check."),
  distribution = info_entry("Distribution of a variable",
    "Histogram (numbers) or bar chart (categories) of any column, per row or per individual.",
    "Look for skew, zeros, outliers and unexpected categories before choosing an error family or covariates.",
    "'Per individual' averages the raw ID column, so individuals nested in groups with repeated IDs are merged here."),
  visual_settings = info_entry("Settings for the lifespan-group figures",
    "The variable that groups individuals (ALR, mean age, LS or AFR), the number of bins, panels and the trait scale, applied to both lifespan-group figures. Settings that affect a single figure sit just above that figure.",
    "The grouping variable forms the bins of the trajectory plot and the x-axis of the lifespan plot; the number of bins sets the number of groups in the trajectory plot and of age bins in the lifespan plot. Choose ALR, mean age or LS to diagnose selective disappearance and AFR to diagnose selective appearance. Every bin \u00d7 age point in the trajectory plot and every age bin in the lifespan plot needs at least 3 individuals, so sparse bins are not drawn. With no more distinct values than bins, each value is its own bin, and the same bin boundaries are used in every panel.",
    "Few bins hide patterns, many bins give noisy points. Facets split the sample, so each panel has fewer individuals. A time-varying covariate is summarised by each individual's most common value. On the log scale, differences are ratios."),
  a1 = info_entry("Trajectories within bins",
    "Mean trait at each age for individuals grouped into bins of ALR, LS, mean age or AFR (each individual counted once), optionally in separate panels per covariate level. Points need at least 3 individuals. The bin boundaries (equal width or quantiles), \u00b11 SE and the panels of this figure (any categorical variable, or as in the settings above) are set above the figure; the bin differences below use the same panels.",
    "Overlapping bins: no selective disappearance. Parallel but offset bins: age-independent selection. Bins that diverge or converge with age: age-dependent selection. Use AFR bins for selective appearance. Compare panels to see whether the pattern differs between treatments.",
    "At late ages only long-lived bins remain. ALR-based bins shift when final occasions are missed. If the bin means are non-monotonic (for example highest in the middle bins), the lifespan effect is not linear: Models 2 and 4 use a linear ALR term and will find nothing, so add ALR\u00b2 with 'Same polynomial order' or as an extra term."),
  a1_diff = info_entry("Difference between consecutive bins at each age",
    "Signed difference in mean trait between consecutive bins (higher bin minus the next lower bin) at each age, with least-squares trend lines of difference against age: one line per bin pair, or a single line pooling all pairs (chosen above the figure), each with its 95% confidence band (shaded).",
    "Differences near zero with flat lines: no selection. Constant non-zero differences (flat lines away from zero): age-independent selection. Lines that rise or fall: age-dependent selection. The single pooled line answers 'do neighbouring lifespan bins become more (or less) different with age on average?'; separate lines show whether particular bins drive the pattern, for example only the longest-lived bin diverging at late ages. The saved table reports each line's slope per unit age.",
    "A point appears only where both bins have at least 3 individuals; differences at late ages rest on few individuals. The pooled line does not weight points by sample size, so treat its slope as descriptive. The confidence bands come from the same unweighted least-squares fits (each plotted difference counts as one point), so they are descriptive too; a line through only two points has no band."),
  a2 = info_entry("Trait against lifespan within age bins",
    "Each individual's mean trait within an age bin, plotted against its lifespan (or ALR / AFR), with a straight-line fit and its 95% confidence band per age bin; age bins need at least 3 individuals. The table lists the regression coefficient (slope) in each consecutive age bin, its 95% CI, the correlation r, and the change from the previous bin, plus a weighted trend of the coefficient across age.",
    "Coefficients near zero in every bin: no selection. Similar non-zero coefficients in every bin: age-independent selection. Coefficients that change across consecutive age bins: age-dependent selection. Untick 'Show individual points' above the figure to see only the fitted lines, which makes differences in slope between age bins easier to judge.",
    "Older age bins contain only long-lived individuals, so their lifespan range is narrow and coefficients are imprecise (wide CIs). When count means and variances fall with age, coefficients can change through scale alone: compare with the log scale or r. The age-specific coefficients share individuals, so the trend across ages treats correlated estimates as independent."),
  a5_terminal = info_entry("Trait before death (terminal trajectories)",
    "Mean trait (absolute values) against occasions before death (known lifespan; 0 = the last occasion before death) or before the last record, for individuals in the shortest-, middle- and longest-lived thirds.",
    "Aligning individuals on time before death separates changes linked to approaching death from chronological ageing. A drop (terminal decline) or rise (terminal investment) in the last occasions means dying individuals leave the population with atypical values. Because short-lived individuals reach their final occasions at young ages and long-lived ones at old ages, the same terminal change removes atypical individuals at different ages: if the terminal change has the same size in every lifespan group, selective disappearance through terminal effects is roughly age-independent in strength; if it is larger in one lifespan group, its strength depends on age, which is age-dependent selective disappearance. The note below the figure compares the last occasion with earlier ones after removing the average age trend, because occasions further before death are also younger ages.",
    "Needs individuals that died during the study (censored individuals are excluded). With unknown lifespan the last record may not be close to death. Occasions far before death come only from long-lived individuals. A terminal change on its own is enough to make the interaction models (4, 5) win on the Modelling tab even when the trait does not change with age at all, so use this figure before reading such a win as ageing-related selective disappearance. Time to death = ALR \u2212 age, so a terminal effect is the same model written differently: a change proportional to time to death is an additive ALR term (the age-independent signature), while an effect confined to the last occasions needs ALR \u00d7 age terms (the age-dependent signature). In simulations with no selection at all, a drop over the final two occasions made Models 4 and 5 beat Model 1 by roughly 1450 AIC. Fitting age \u00d7 time-to-death alongside age \u00d7 ALR separates the two, except in the strictly linear case where they are algebraically identical; they also answer different questions, because a trajectory anchored to the date of death is undefined beyond it."),
  a6_selection = info_entry("Selection differentials by age",
    "At each age, the difference in mean trait between the individuals that survive to their next occasion and all individuals present (the selection differential), or between survivors and those that disappear, in within-age SD units, with 95% intervals and a weighted trend across age.",
    "Positive values: survivors have higher trait values than the population they came from, so the population mean rises through selective disappearance; negative values: the reverse. Differentials that are similar at all ages indicate age-independent selective disappearance; a trend across age indicates age-dependent selective disappearance. Summing differentials over ages approximates how far cross-sectional means drift from the individual trajectory.",
    "Ages need at least three survivors and three disappearing individuals. Disappearance is defined as for the other disappearance diagnostics. Differentials describe selection on the observed trait, not its cause."),
  a7_hazard = info_entry("Disappearance hazard",
    "At each age, the probability of disappearing before the next occasion (with 95% Wilson intervals, point size = individuals at risk), the overall, age-independent hazard (dashed line: all disappearances / all individual-occasions at risk), and a likelihood-ratio test of a constant hazard against age-specific hazards.",
    "A flat hazard close to the dashed line means the chance of dying does not change with age; a rising hazard indicates actuarial senescence. The hazard shows how quickly each age class is thinned: selective disappearance can only bias late ages strongly if many individuals disappear before them, so read the other diagnostics with the hazard in mind.",
    "With unknown lifespan, disappearances include missed detections and emigration, and individuals alive at the study end must be marked as censored. Hazards at old ages rest on few individuals (wide intervals). The hazard is fitted with ordinary binomial GLMs on individual \u00d7 occasion records: the discrete-time likelihood is standard, but unmodelled differences in frailty between individuals are not accounted for."),
  sampling_guidance = info_entry("Sampling summary and proxy guidance",
    "Summarises detection and proxy quality to guide the choice between ALR (Models 2, 4) and mean age (Models 3, 5). Expected occasions end at the known lifespan, or at the last record if lifespan is unknown.",
    "Read the guidance alongside the proxy-agreement plots.",
    "The missingness grid assumes regular, scheduled sampling occasions: it is most reliable when every individual is due to be sampled at a common interval. With irregular ages (for example ages calculated from dates) expected occasions are approximate, so missingness percentages, patterns and drivers should be read qualitatively; round ages to a sampling resolution on the Data tab. Without known lifespan, occasions missed after the last record are invisible, so detection is overestimated. If individuals genuinely follow different schedules (for example some sampled every year and others every third year), the common step is the shortest one and the less frequently sampled individuals appear here as heavily missing: that is study design, not failed detection."),
  heatmap = info_entry("Sampling grid",
    "One row per individual and one column per expected occasion: observed, missed (expected but not recorded), the final expected occasion (death if LS is known, otherwise ALR), and white cells where no observations are expected (before AFR, after death or ALR). Each individual's schedule starts from its own first record (or the common start age), so staggered cohorts are not given spurious missed occasions.",
    "Scattered red cells: random gaps. Red concentrated at young or old ages: age-dependent sampling. Order rows by ALR, AFR, mean trait, ID or at random to see whether gaps line up with lifespan, entry age or trait values.",
    "Up to 500 individuals (set with 'Individuals to show in the grid') and 60,000 cells are drawn, evenly spaced in the chosen order. With irregular ages the occasions are approximate; round ages on the Data tab. Very fine or irregular ages are coarsened so that the grid stays below 300,000 cells, which is reported on this tab."),
  missing_by_age = info_entry("Missingness and sample size by age",
    "Individuals observed at each age (bars) and the percentage of expected occasions missed (line), with a table.",
    "Missingness rising with age suggests 'missing when old'; high early missingness suggests 'missing when young'.",
    "Percentages at the oldest ages rest on few individuals."),
  missing_vs_var = info_entry("Missingness against other variables",
    "Missed-occasion rate against an individual-level variable, and logistic regressions of missingness on age, LS, AFR, the previous trait value and condition.",
    "Clear associations suggest that missingness is not completely at random.",
    "Associations in observed data cannot prove MCAR, MAR or MNAR, and running many tests inflates false positives. Missingness that depends on the trait itself (the 'Prior observed trait' row) is the most serious case for this app: individuals whose trait is low are recorded less often, so their records end earlier and their trait values are the low ones. That produces the same evidence as age-dependent selective disappearance - the interaction and centring models (3, 4, 5) can beat Model 1 decisively when there is no selection at all - and no model on this data can separate the two. Treat a win for those models as conditional on this row, and check the trait-before-death figure on the Visual diagnosis tab as well. Each individual contributes many rows to these regressions, so the p-values are optimistic, and failing to detect an association is weak evidence for MCAR."),
  proxy_agreement = info_entry("ALR, mean age and lifespan",
    "Pairwise agreement between ALR, mean sampled age and known lifespan (one point per individual, with r, a dotted 1:1 line and a solid black linear regression line).",
    "High r(ALR, LS) supports ALR-based models (2, 4); high r(mean age, LS) supports centring models (3, 5). Under complete sampling ALR and mean age are nearly interchangeable. The regression line shows how one measure scales with the other: departures from the 1:1 line reveal systematic under- or overestimation, for example ALR falling increasingly short of lifespan when final occasions are missed.",
    "Missed final occasions make ALR underestimate lifespan; intermittent gaps shift mean age. With unknown lifespan only ALR vs mean age is shown."),
  a3_settings = info_entry("Individual fit settings",
    "Chooses the ageing function fitted to each individual, how the population curve is reconstructed, and which individuals are drawn. Individual fits are a visual diagnostic of how individuals age.",
    "Compare functions in the table. 'Use this ageing function for the mixed models' makes the selected function the ageing function on the Modelling tab (otherwise the Modelling tab starts with the function with the lowest mean \u0394AICc across individuals). 'Exponential (a \u00b7 exp(b \u00b7 age))' is fitted by non-linear least squares and is the only function that is non-linear in its parameters. Draw the longest- or shortest-lived individuals to see how fits change with the number of records.",
    "Individuals need at least 2 records at different ages for a linear fit, 3 for quadratic, logarithmic and exponential fits, and 4 for a cubic fit. With exactly as many records as parameters the curve passes through every point, so it carries no residual information (no AICc) and is left out of the function comparison. Non-linear fits that do not converge are dropped."),
  a3_fits = info_entry("Individual trajectories and data support",
    "Each panel shows one individual's records and fitted function, with free axes. The support metrics above summarise how much information there is to estimate individual ageing slopes: individuals, total trait observations, median and interquartile range of observations per individual, the percentage of individuals with at least 2, 3 and 4 time points (distinct ages), the number of distinct ages and the median age span per individual.",
    "Look for systematic misfit, for example curvature that a straight line misses. A straight line can be fitted to 2 time points, but it passes through both exactly; slopes that can be separated from residual noise need at least 3 time points, and a quadratic fit that is not exact needs 4 or more.",
    "With few records the fit is nearly exact and says little about the true shape. If few individuals have \u2265 3 time points, individual slopes (individual fits) and random slopes (Models tab) are poorly estimated. A quadratic fitted to 3 records, or a cubic to 4, interpolates them exactly (no residual degrees of freedom): the curve is drawn but says nothing about how well that function describes the individual."),
  a3_mean = info_entry("Mean-coefficient and mean-function trajectories",
    "Two reconstructions from the individual fits f(age; \u03b8\u1d62), i = 1\u2026N. Mean of coefficients: f(age; \u03b8\u0304) with \u03b8\u0304 = (1/N) \u03a3\u1d62 \u03b8\u1d62. Mean of functions: f\u0304(age) = (1/N) \u03a3\u1d62 f(age; \u03b8\u1d62). When f is linear in its parameters, f(age; \u03b8) = \u03a3\u2096 \u03b8\u2096 g\u2096(age), so f\u0304(age) = \u03a3\u2096 \u03b8\u0304\u2096 g\u2096(age) = f(age; \u03b8\u0304) and the two curves are identical (all polynomial, logarithmic and exponential-basis functions here). When f is non-linear in \u03b8 (e.g. a\u00b7exp(b\u00b7age)), they differ by Jensen's inequality: f(age; \u03b8\u0304) is the trajectory of the average individual, whereas f\u0304(age) is the average trajectory of the population. Bands: 95% intervals from SD\u1d62[f(age; \u03b8\u1d62)]/\u221aN (functions) and the delta method with Cov(\u03b8\u1d62)/N (coefficients).",
    "Unlike observed means, averaged within-individual fits are not shifted by who survives. Without selective disappearance the band should cover the true trajectory. For non-linear functions choose the estimand: the average individual (mean of coefficients, comparable to a GLMM's fixed-effect curve on the link scale) or the population average (mean of functions, comparable to observed means in the absence of selection).",
    "These reconstructions are unreliable under high missingness, few time steps, few records per individual or a small number of individuals, and should then be read as illustrative only. Both reconstructions can depart from the true latent trajectory for reasons that are not selective disappearance, and the departure is usually largest at the youngest and oldest ages. Each individual's fitted function is evaluated at ages the individual itself never supplied (before it entered, after it died), so the average includes extrapolation; with few records per individual the fits are close to saturated, which makes that extrapolation unstable; and for count traits the individual fits are on the response scale while the simulated or modelled truth is a curve on the link scale, so the average of individual curves estimates the mean of exponentials rather than the exponential of the mean. Under missingness these effects grow. Always read these curves alongside the population-level model predictions on the Modelling tab: if the model predictions track the truth and these curves do not, the discrepancy is in the individual fits, not in the models. Individuals with few records give very imprecise coefficients, so the curves can miss the truth at young and old ages even without bias: read the bands. The mean of functions is sensitive to individuals with extreme non-linear fits. Curves are drawn only across ages observed in at least 10 individuals. Short-lived individuals with too few records cannot be fitted, which reintroduces survivor bias, and robust averages such as medians are biased under age-dependent selection."),
  a3_compare = info_entry("Individual-level function comparison",
    "Fit of each function to individuals, including the non-linear exponential: mean R\u00b2, adjusted R\u00b2, and mean \u0394AICc over individuals for which every function is estimable.",
    "The lowest mean \u0394AICc (and highest adjusted R\u00b2) indicates the best individual-level function. Treat the ranking as suggestive, not prescriptive: it describes which function is estimable from these individuals, not which function generated the data, and a difference of a few AICc units is not evidence for one shape over another.",
    "This whole comparison is unreliable when missingness is high, when individuals have few time steps or few records, or when the sample of individuals is small: with little data per individual the simplest function wins almost automatically, and the mean \u0394AICc is driven by a handful of well-sampled individuals. N_common can be small, and functions with more parameters are fitted to fewer individuals. With few occasions per individual the simplest function usually wins on AICc even when the data were generated by a more complex one: a quadratic needs at least six records for AICc to be defined at all, so in short-lived systems it is often undefined or penalised out, and a linear fit ranks first. This says which function is estimable from these individuals, not which function generated the data. The exponential model a \u00b7 exp(b \u00b7 age) is fitted here on raw age but on standardised age in the population models, so its rate b is not comparable between the two."),
  a3_coefs = info_entry("Individual coefficients",
    "Mean and SD of the individual coefficients on the raw age scale.",
    "The SD combines true among-individual differences and estimation error.",
    "Coefficients are not comparable across functions."),
  b2 = info_entry("Population-level ageing function",
    "Evaluated at the population level rather than per individual: fits Model 1 as a mixed model with each selected ageing function (including the non-linear exponential for Gaussian traits), using the error family, covariates and random effects chosen on the Modelling tab, and compares AIC and fitted curves. Points are observed means.",
    "Compare the result with the individual fits above: the individual comparison shows which shape describes individuals, this one which shape describes the population trajectory of the sampled records. Tick only the functions you want to compare to reduce clutter.",
    "Model 1 ignores selective disappearance, so the preferred function can change once proxies are added. Cubic curves extrapolate poorly. With count families the exponential model a \u00b7 exp(b \u00b7 age) equals the Linear function on the log scale and is left out. 'Asymptotic exponential' is an intercept plus \u03b2 \u00b7 exp(\u2212standardised age): it is linear in its coefficients and flattens towards a plateau, with the rate fixed by the age scale. It is a different model from 'Exponential (a \u00b7 exp(b \u00b7 age))', whose rate b is estimated, and its shape depends on the age distribution of the data analysed."),
  model_settings = info_entry("Model settings",
    "Error family, ageing function, the individual-level random effects, how among-individual terms enter, standardisation, whether invalid fits may be ranked, and which models to compare. Each model row has a tick box, an 'i' that explains what the model tests (ALR- or centring-based; selective disappearance and appearance, age-dependent or age-independent) and gives its exact specification with the current covariates, interactions and random effects, and a menu to add terms to that model only (for example ALR additively to Model 5 or 6, or a covariate \u00d7 age). The wrench next to each menu builds a term from up to three chosen terms (age, ALR, AFR, LS, mean age or a covariate), each joined to the previous one by + (additive) or \u00d7 (interaction with main effects), for example ALR \u00d7 AFR \u00d7 age; built terms appear in that model's menu and can be removed there. Binomial families show a 'Weights' menu for the number of trials behind each proportion.",
    "Random effects: '(1 | ID)' random intercept; uncorrelated '(1 | ID) + (0 + age | ID)'; correlated '(1 + age | ID)'; 'Automatic' chooses from the data support shown below the menu (supported when at least 30 individuals and half of all individuals have \u2265 3 distinct ages). Among-individual terms (ALR, LS, AFR, mean age) enter linearly by default even with quadratic or cubic ageing; 'Same polynomial order' adds their squares (and cubes). The non-linear exponential fits trait = a \u00b7 exp(b \u00b7 z_age) with nlme (Gaussian traits); for counts it equals the Linear function on the log link. Press 'Fit models' after any change.",
    "AIC compares models only within the same family and random-effect structure. Random slopes need several records per individual. Terms added to one model change what it tests, so nested likelihood-ratio tests are only reported for genuinely nested pairs. Including invalid fits in the ranking is for inspection only. Built interactions with age multiply the ageing terms (a three-way interaction with a quadratic adds many coefficients), so they can over-fit, fail to converge or be hard to interpret: compare them with simpler models. The zero-inflated negative binomial comes in two variance forms, nbinom2 (variance \u03bc + \u03bc\u00b2/\u03b8) and nbinom1 (variance \u03bc(1 + \u03c6)); compare them with the error-family check."),
  model_support = info_entry("Model support and fit validity",
    "\u0394AIC, Akaike weights, nested likelihood-ratio tests and fitting status for Models 1\u201310 (and any extra terms), all on the same rows. Each fit is classified: Valid (no warnings); Caution (boundary or singular fit, convergence-gradient warnings, rank-deficient terms dropped, or other notes); Failed (no fit, a non-positive-definite or singular Hessian, or a non-finite AIC).",
    "\u0394AIC < 2 means similar support. Model 2 vs 4 tests age-dependent selective disappearance through ALR, 3 vs 5 through mean age; 2 vs 7 and 4 vs 8 test selective appearance; Models 9 and 10 let disappearance and appearance differ in age dependence (9: ALR \u00d7 age + AFR; 10: ALR + AFR \u00d7 age). Likelihood-ratio tests are only reported for pairs whose terms are genuinely nested. Failed fits are excluded from the ranking, the Akaike weights, the tests, the prediction plot and all automated interpretation unless 'Include fits with invalid Hessians' is ticked.",
    "Missingness that depends on the trait itself produces the same ranking as selective disappearance (Models 3 and 5 can win with no true selection at all), so read this ranking together with the missingness drivers on the 'Missingness and proxies' tab before concluding. Terminal declines shortly before death do the same: check the trait-before-death figure on the disappearance-diagnostics tab. Integer data treated as Gaussian: check that the trait is continuous rather than a count or a 0/1 outcome, to avoid a misfit family. Without random slopes, interaction models can be favoured even without selective disappearance when individuals differ in ageing rates. A Caution fit can still be informative, but compare it with a simpler random-effect structure. Rows missing any model variable are dropped from all models. Likelihood-ratio p-values are approximate for mixed models (boundary variances, singular fits, zero inflation), and the ranking shows relative support among the models fitted, not that the best one generated the data."),
  consistency = info_entry("Internal consistency of the evidence",
    "Automatic cross-checks between the model ranking and other evidence: for simulations, whether the lowest-AIC model also recovers the true trajectory best; whether likelihood-ratio tests agree with AIC; whether the age dependence in the best model agrees with the trend of the lifespan\u2013trait coefficient across age bins; whether the best fit is only classified 'Caution'; and whether interaction models win while random slopes are not fitted although the data support them.",
    "Agreement across criteria strengthens a conclusion. When they disagree, report the disagreement and prefer the more conservative reading. AIC measures how well a model describes the sampled records; recovering the latent population trajectory is a different goal, and with sparse late ages, few occasions, zero inflation or confounded appearance and disappearance the two can point to different models.",
    "The checks are heuristic and cannot detect every inconsistency; the absence of a flag is not proof that the evidence is coherent."),
  varcomp = info_entry("Random-effect variances",
    "Estimated random-effect SDs (and correlations) for each fitted model.",
    "With 'Random slope' ticked, the individual level shows an SD for the age term (f1) and its correlation with the intercept, which confirms that random slopes were fitted.",
    "SDs near zero indicate boundary (singular) fits."),
  predictions = info_entry("Population-level trajectories",
    "Each model's predicted population-level trajectory, calculated by: (1) excluding random effects (re.form = NA), so the curve describes a typical individual rather than any sampled individual; (2) holding numeric covariates at their mean across individuals; (3) marginalising over factor covariates \u2014 predictions are made for every observed combination of factor levels and averaged, weighted by the number of individuals with that combination (or shown separately for each level of the factor chosen in 'Show predictions by'); (4) holding ALR, LS, AFR and individual mean-age terms at their individual-level means (polynomial among-individual terms use the square or cube of that mean value); (5) reporting the response scale, so count models include the zero-inflation probability. Observed means, the decomposition, the reconstruction from individual fits (mean of coefficients and mean of individual functions, using the function chosen on the 'Individual and population trajectories' tab) and, for simulations, the truth can be overlaid. Choose which models to draw to reduce clutter.",
    "Models that account for selective disappearance should track the truth in simulations; observed means and the decomposition are biased under age-dependent selection. With covariate \u00d7 age interactions, show predictions by that covariate to see how its levels age differently. Tick 'Individual records (jittered)' to add every record's raw trait value as a small transparent point, jittered slightly along age only (trait values are not changed). The decomposition (Rebke et al. 2010) chains the mean within-individual change between successive sampling occasions, using only the individuals sampled at both (survivor-restricted); on an irregular schedule the records are first placed on a common grid of occasions, and a caution appears below the figure.",
    "Curves at sparsely sampled late ages are uncertain. For counts the curve averages predictions across factor combinations on the response scale, which is not the same as a prediction at average covariate values. Failed fits are not drawn unless included; Caution fits are dashed. When selective disappearance is age-independent (lifespan linked to an individual's level, not to its rate of ageing), every model returns almost the same curve even when their AIC values differ by hundreds of units. This is expected, not a fault: the proxy is held at its mean, so additive and interaction proxy terms contribute nothing at that value, and the age slope is unbiased in all of the models. The AIC difference comes from explaining among-individual variation in level, which improves fit without changing the trajectory. Curves separate when disappearance is age-dependent, because then the uncorrected models estimate a biased slope."),
  accuracy = info_entry("Accuracy against the simulated truth",
    "Relativised deviation D = 100 \u00d7 (estimate \u2212 truth) / truth (Methods Eq. 11) for each model, observed means and the decomposition.",
    "Smaller mean |D| is better.",
    "Simulations only; restricted to ages with at least 10 observed individuals."),
  coefficients = info_entry("Coefficients, scaling and interpretation",
    "Fixed-effect estimates, standard errors and tests for the selected model, with the scaling constant behind each coefficient, the estimate per original unit, the random-effect variances, and a plain-language reading of the selective disappearance or appearance terms. Scaling: with standardisation, age = (age \u2212 centre) / SD and ALR, LS and AFR are standardised across individuals; a term's scale factor is the product of the SDs of its components (for example SD\u00b2 for age\u00b2, SD_age \u00d7 SD_ALR for age \u00d7 ALR). Estimate per original unit = estimate / scale factor: the coefficient per unit of the centred original variables (for example per year of ALR, or per (age \u2212 centre)\u00b2). Numeric covariates are not standardised. Count models are on the log scale.",
    "The interpretation compares predicted trajectories for individuals at the 10th and 90th percentiles of ALR (or LS, AFR, mean age) at a young, middle and old age. If only the additive term is supported, longer- and shorter-lived individuals differ by a similar amount at all ages (age-independent selective disappearance). If the interaction is supported, the gap changes with age (age-dependent), for example longer-lived individuals having higher trait values mainly at older ages.",
    "Standardisation is recommended for numerical stability; the per-unit column translates coefficients back to biological units. P-values are approximate (Wald; likelihood-ratio tests where nested models were fitted). The interpretation is an association conditional on the model, not proof of selection: check it against the visual diagnosis plots, residual diagnostics and alternative random-effect structures. Standardised coefficients are expressed on the age scale of the data analysed, so they change if the sample or its age range changes."),
  definitions = info_entry("Model definitions",
    "Fixed effects and the question addressed by each ticked model, with the chosen covariates and interactions.",
    "Models 3 and 5 separate within- and among-individual age effects; Model 6 needs known lifespan; Models 7\u201310 add AFR (age at first observation). Click the 'i' next to a model in the settings for its exact specification, including terms added to that model.",
    "Covariates enter every model; interactions chosen on the Data tab are added to every model; terms added in a model's menu apply only to that model."),
  family_check = info_entry("Error-family check",
    "Refits one model with each count family and compares AIC.",
    "Choose the lowest AIC among count families, then check residuals.",
    "Not comparable with a Gaussian AIC; zero-inflated models can be slow or fail to converge."),
  dharma = info_entry("Residual checks (DHARMa)",
    "Simulation-based residual tests for uniformity, dispersion and zero inflation, with diagnostic plots.",
    "Small p-values indicate that the error distribution or model structure fits poorly.",
    "In large datasets trivial deviations become significant: judge the plots."),
  performance = info_entry("Model diagnostics (performance)",
    "Checks from the performance package for the selected model: Nakagawa R\u00b2, intraclass correlation, singularity, convergence, collinearity and, for counts, overdispersion and zero inflation, plus the performance::check_model() diagnostic panels (needs the see package).",
    "Use these alongside the DHARMa residual checks: singular fits suggest simplifying the random effects; overdispersion suggests a negative binomial family; a zero-inflation ratio well below 1 suggests a zero-inflated family.",
    "Not every check supports every model type (for example some zero-inflated glmmTMB models): unavailable checks are reported rather than stopping the app. Polynomial and interaction terms are collinear by construction. check_model() can be slow on large datasets."),
  rcode = info_entry("Reproducible R code",
    "A script that reproduces the model comparison with lme4 or glmmTMB.",
    "Run it on your data file to reproduce or extend the analysis.",
    "It assumes the same column names and settings as the app."),
  saved = info_entry("Saved results",
    "The results you saved from other tabs. They make up the exported HTML report (with plots) or text report.",
    "Remove items you no longer need before exporting.",
    "Saved items are snapshots: they do not update when the data or settings change."),
  overview_auto = info_entry("Automatic overview",
    "An automatic summary of the current data and settings.",
    "Use it as a checklist of where to look, not as a conclusion. Save it to include it in the report.",
    "The wording is deliberately cautious: automated rules cannot judge biological plausibility, sampling design or model adequacy. Model rankings depend on the error family, the random-effect structure, the ageing function, proxy quality, sampling and fit validity; confirm with the plots and diagnostics before reporting.")
)

INFO$a2_table <- INFO$a2
INFO$a2_table$title <- "Regression coefficient across consecutive age bins"

# Collapsible box with the R code that reproduces a section (the minus sign hides it)
section_code_box <- function(section) {
  fluidRow(box(width = 12, title = "R code for this section", status = "info", solidHeader = TRUE, collapsible = TRUE,
    p(class = "small-note", "A script that reproduces this section with the current settings. The app's own helper functions are included, so it runs on its own (packages: ggplot2; lme4, and glmmTMB or nlme where models are fitted). If you uploaded a file, put it in R's working directory."),
    div(class = "code-toolbar",
        tags$button(type = "button", class = "btn btn-default btn-sm", onclick = sprintf("disapprCopy('code_%s_text')", section), icon("copy"), " Copy"),
        downloadButton(paste0("download_code_", section), "Download .R", class = "btn-default btn-sm")),
    div(class = "section-code", verbatimTextOutput(paste0("code_", section, "_text")))))
}

# Video walkthroughs hosted on the project's OSF page (https://osf.io/kevnm/); streamed, so they need internet access.
HELP_VIDEOS <- list(
  part1 = list(title = "How to use disappR, part 1", guid = "u6mxk", minutes = "4:46",
               text = "The first half of a screen-recorded walkthrough of the app."),
  part2 = list(title = "How to use disappR, part 2", guid = "5krbc", minutes = "5:02",
               text = "The second half of the walkthrough, continuing from part 1.")
)

help_video_modal <- function(v) {
  modalDialog(
    title = v$title, size = "l", easyClose = TRUE,
    tags$video(src = paste0("https://osf.io/download/", v$guid, "/"), controls = NA, autoplay = NA, preload = "auto",
               playsinline = NA, style = "width:100%; max-height:70vh; background:#000; border-radius:8px;",
               "Your browser cannot play this video. ",
               tags$a(href = paste0("https://osf.io/", v$guid, "/"), target = "_blank", "Open it on OSF instead.")),
    p(class = "small-note", style = "margin-top:8px;", "Streamed from OSF, so it needs an internet connection and may take a few seconds to start. ",
      tags$a(href = paste0("https://osf.io/", v$guid, "/"), target = "_blank", "Open on OSF"), "."),
    footer = modalButton("Close")
  )
}

info_title <- function(title, key) {
  tags$span(title, actionLink(paste0("info_", key), label = "", icon = icon("circle-info"),
                              class = "info-link", title = "What is this, and how do I read it?"))
}

info_body <- function(e) {
  tagList(
    h5(strong("What it does")), p(e$what),
    h5(strong("How to use and interpret it")), p(e$how),
    h5(strong("Cautions")), p(e$caution)
  )
}

save_button <- function(id, label = "Save to summary") {
  actionButton(id, label, icon = icon("bookmark"), class = "btn-default btn-sm save-btn")
}

# ---------------------------------------------------------------------------
# Saved results: rendering and export
# ---------------------------------------------------------------------------
draw_saved_plot <- function(pl) {
  if (is.function(pl)) pl() else print(pl)
  invisible(NULL)
}

html_escape <- function(x) htmltools::htmlEscape(as.character(x))

table_to_html <- function(df, digits = 4) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return("")
  cell <- function(v) if (is.numeric(v)) format(signif(v, digits), trim = TRUE) else as.character(v)
  head <- paste0("<tr>", paste0("<th>", html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(df)), function(i) {
    vals <- vapply(df, function(col) cell(col[i]), character(1))
    paste0("<tr>", paste0("<td>", html_escape(vals), "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<table>", head, paste(body, collapse = ""), "</table>")
}

plot_to_data_uri <- function(pl, width = 8, height = 5) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  f <- tempfile(fileext = ".png")
  opened <- tryCatch({
    grDevices::png(f, width = width, height = height, units = "in", res = 110)
    TRUE
  }, error = function(e) FALSE)
  if (!opened) return(NULL)
  dev <- grDevices::dev.cur()
  ok <- tryCatch({
    draw_saved_plot(pl)
    TRUE
  }, error = function(e) FALSE)
  grDevices::dev.off(dev)
  if (!ok || !file.exists(f)) return(NULL)
  paste0("data:image/png;base64,", jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size)))
}

html_report <- function(entries, heading = "disappR report") {
  css <- paste0("body{font-family:Helvetica,Arial,sans-serif;max-width:980px;margin:30px auto;color:#3d342d;line-height:1.5;padding:0 16px}",
                "h1{color:#5a4033}h2{color:#7b5a45;border-bottom:1px solid #e5d6c6;padding-bottom:4px;margin-top:34px}",
                ".meta{color:#8a7564;font-size:13px}table{border-collapse:collapse;margin:10px 0;font-size:13px}",
                "th,td{border:1px solid #e0d2c2;padding:4px 8px;text-align:left}th{background:#f4e8da}",
                "img{max-width:100%;margin:8px 0}pre{background:#2e2620;color:#f3e9dc;padding:12px;border-radius:8px;overflow:auto;font-size:12px}")
  parts <- c("<!DOCTYPE html><html><head><meta charset='utf-8'><title>disappR report</title><style>", css, "</style></head><body>",
             paste0("<h1>", html_escape(heading), "</h1>"),
             paste0("<p class='meta'>Generated ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M")), " with disappR 0.7.0</p>"))
  if (!length(entries)) parts <- c(parts, "<p>No results were saved.</p>")
  for (e in entries) {
    parts <- c(parts, paste0("<h2>", html_escape(paste(e$section, "\u00b7", e$title)), "</h2>"),
               paste0("<p class='meta'>Saved ", html_escape(e$time), " \u00b7 ", html_escape(e$data), "</p>"))
    if (length(e$text)) parts <- c(parts, "<ul>", paste0("<li>", html_escape(e$text), "</li>"), "</ul>")
    tabs <- e$tables
    for (j in seq_along(tabs)) {
      cap <- names(tabs)[j]
      if (!is.null(cap) && !is.na(cap) && nzchar(cap)) parts <- c(parts, paste0("<p><b>", html_escape(cap), "</b></p>"))
      parts <- c(parts, table_to_html(tabs[[j]]))
    }
    for (pl in e$plots) {
      uri <- plot_to_data_uri(pl)
      if (!is.null(uri)) parts <- c(parts, paste0("<img src='", uri, "'/>"))
    }
    if (!is.null(e$code)) parts <- c(parts, paste0("<pre>", html_escape(e$code), "</pre>"))
  }
  c(parts, "</body></html>")
}

text_report <- function(entries, heading = "disappR report") {
  lines <- c(heading, paste("Generated:", format(Sys.time())), "")
  if (!length(entries)) lines <- c(lines, "No results were saved.")
  for (e in entries) {
    lines <- c(lines, paste0("== ", e$section, " \u00b7 ", e$title, " =="), paste("Saved:", e$time, "|", e$data))
    if (length(e$text)) lines <- c(lines, paste("-", e$text))
    tabs <- e$tables
    for (j in seq_along(tabs)) {
      cap <- names(tabs)[j]
      if (!is.null(cap) && !is.na(cap) && nzchar(cap)) lines <- c(lines, paste0("[", cap, "]"))
      lines <- c(lines, utils::capture.output(print(as.data.frame(tabs[[j]]), row.names = FALSE)))
    }
    if (length(e$plots)) lines <- c(lines, sprintf("(%d plot(s) included in the HTML report)", length(e$plots)))
    if (!is.null(e$code)) lines <- c(lines, e$code)
    lines <- c(lines, "")
  }
  lines
}
