# Figure 7 of Sanghvi et al. 2026 (Methods in Ecology and Evolution) and its model comparisons, refitted with random
# slopes: fruit flies with a random slope and curvature (linear and quadratic age) per daughter, great tits with a
# random slope per female. Fixed effects, data preparation, prediction grid and plotting follow the authors' code
# ('sel dis code.Rmd', OSF kevnm, sections 'Empirical data' and 'Model predictions'). Random slopes are on standardised
# age and uncorrelated with the intercept (glmmTMB diag()), the structure that converges for both datasets.
#   Rscript dev/figure7_random_slopes.R <folder with 'Fruit fly data.csv' and 'Great tits data.csv'> <output folder> [slopes|intercept]
# 'intercept' reproduces the published figure (random intercepts only), for comparison.
suppressPackageStartupMessages({ library(glmmTMB); library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
dat_dir <- args[1]; out_dir <- args[2]; re_mode <- if (length(args) >= 3) args[3] else "slopes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) { p <- file.path(dat_dir, f); if (!file.exists(p)) stop("Missing ", p); read.csv(p, stringsAsFactors = TRUE) }
log_lines <- character(0); say <- function(...) { x <- sprintf(...); cat(x, "\n"); log_lines <<- c(log_lines, x) }

# ---------------------------------------------------------------- fruit flies (as in the authors' code) ----------
dataAMNAT <- rd("Fruit fly data.csv")
dataAMNAT$lifespan <- as.numeric(dataAMNAT$LS); dataAMNAT$F1_Age <- as.numeric(dataAMNAT$F1_Age); dataAMNAT$ALR <- as.numeric(dataAMNAT$ALR)
dataAMNAT <- dataAMNAT %>% ungroup() %>% mutate(F1_Age_scaled = as.numeric(scale(F1_Age)), ALR_scaled = as.numeric(scale(ALR))) %>%
  group_by(F1_ID) %>% mutate(Mean_age = mean(F1_Age), Delta_age = F1_Age - Mean_age, Age2 = F1_Age^2, Mean_age2 = mean(Age2)) %>% ungroup()
fly_age_mu <- mean(dataAMNAT$F1_Age); fly_age_sd <- stats::sd(dataAMNAT$F1_Age)
dataAMNAT$a <- (dataAMNAT$F1_Age - fly_age_mu) / fly_age_sd
fly_re <- if (re_mode == "slopes") "(1 | F0_ID) + diag(1 + a + I(a^2) | F0_ID:F1_ID)" else "(1 | F0_ID/F1_ID)"
fly_fixed <- c(
  Model1_Naive     = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + F1_Age + I(F1_Age^2)",
  Model2_AddALR    = "F2_Count ~ Rep + ALR + Paternal_age + Paternal_sperm_age + F1_Age + I(F1_Age^2)",
  Model3_AddCenter = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + Mean_age + Delta_age + I(Age2 - Mean_age2)",
  Model4_IntALR    = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + ALR*F1_Age + ALR*I(F1_Age^2)",
  Model5_IntCenter = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + (Mean_age)*(Delta_age) + I(Age2-Mean_age2)*(Mean_age)")
fly_ctl <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
# Rescaled versions of the same models (ages / 10, squared terms / 100, mean age centred): an exact reparametrisation
# with the same likelihood, which converges where the raw scale does not.
mm_fly <- mean(dataAMNAT$Mean_age)
add_scaled <- function(x) dplyr::mutate(x, s_age = F1_Age / 10, ALR_s = ALR / 10, mc = (Mean_age - mm_fly) / 10, da = Delta_age / 10, dq = (Age2 - Mean_age2) / 100)
dataAMNAT <- add_scaled(dataAMNAT)
fly_scaled <- c(
  Model1_Naive     = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + s_age + I(s_age^2)",
  Model2_AddALR    = "F2_Count ~ Rep + ALR_s + Paternal_age + Paternal_sperm_age + s_age + I(s_age^2)",
  Model3_AddCenter = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + mc + da + dq",
  Model4_IntALR    = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + ALR_s*s_age + ALR_s*I(s_age^2)",
  Model5_IntCenter = "F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + mc*da + dq*mc")
# Each model is fitted in the authors' raw form and in the rescaled form, each with the authors' optimiser (BFGS) and
# glmmTMB's default (nlminb). All four are the same model, so the fit with the highest likelihood and a valid Hessian
# is the maximum-likelihood fit and is kept. With random slopes, the raw-scale fits of the mean-centring models 3 and 5
# stop far from the maximum (about 2,500 log-likelihood units below it) with both optimisers.
best_fit <- function(forms, ...) {
  fits <- list()
  for (k in names(forms)) {
    fits[[paste(k, "BFGS")]] <- tryCatch(glmmTMB(forms[[k]], ..., control = fly_ctl), error = function(e) NULL)
    fits[[paste(k, "nlminb")]] <- tryCatch(glmmTMB(forms[[k]], ...), error = function(e) NULL)
  }
  fits <- Filter(Negate(is.null), fits)
  ll <- vapply(fits, function(m) { v <- as.numeric(logLik(m)); if (isTRUE(m$sdr$pdHess) && is.finite(v)) v else -Inf }, numeric(1))
  m <- fits[[which.max(ll)]]; attr(m, "optimiser") <- names(fits)[which.max(ll)]
  attr(m, "loglik_by_optimiser") <- ll; m
}
fly_fit <- lapply(names(fly_fixed), function(nm) best_fit(list(raw = as.formula(paste(fly_fixed[[nm]], "+", fly_re)), rescaled = as.formula(paste(fly_scaled[[nm]], "+", fly_re))),
                                                          family = nbinom2(), ziformula = ~1, data = dataAMNAT))
names(fly_fit) <- names(fly_fixed)
for (nm in names(fly_fit)) say("Fly %s: kept %s (log-likelihood %s)", nm, attr(fly_fit[[nm]], "optimiser"),
                               paste(sprintf("%s %.2f", names(attr(fly_fit[[nm]], "loglik_by_optimiser")), attr(fly_fit[[nm]], "loglik_by_optimiser")), collapse = ", "))
# Model 4 rescaled (ALR and age standardised) for coefficients and p values, as in the authors' code
fly_m4_scaled <- best_fit(list(scaled = as.formula(paste("F2_Count ~ Rep + Paternal_age + Paternal_sperm_age + ALR_scaled * poly(F1_Age, 2) +", fly_re))),
                          family = nbinom2(), ziformula = ~1, data = dataAMNAT)

# ---------------------------------------------------------------- great tits ----------------------------------------
dataPROCB <- rd("Great tits data.csv")
dataPROCBLD <- dataPROCB %>% filter(LD > -1) %>% mutate(f_min_age = as.numeric(f_min_age), f_ALR = as.numeric(f_ALR)) %>%
  group_by(female) %>% mutate(mean_age_LD = mean(f_min_age, na.rm = TRUE), delta_age_LD = f_min_age - mean_age_LD,
                              age_sq = f_min_age^2, mean_age_sq_LD = mean(age_sq, na.rm = TRUE)) %>% ungroup()
tit_age_mu <- mean(dataPROCBLD$f_min_age); tit_age_sd <- stats::sd(dataPROCBLD$f_min_age)
dataPROCBLD$a <- (dataPROCBLD$f_min_age - tit_age_mu) / tit_age_sd
tit_re <- if (re_mode == "slopes") "diag(1 + a | female) + (1|year) + (1|area)" else "(1|female) + (1|year) + (1|area)"
tit_env <- "+ loc_density + pred + f_status + YR_FL"
tit_fixed <- c(
  laying_Model1_Naive     = "LD ~ f_min_age + I(f_min_age^2)",
  laying_Model2_AddALR    = "LD ~ f_min_age + I(f_min_age^2) + f_ALR",
  laying_Model3_AddCenter = "LD ~ mean_age_LD + delta_age_LD + I(age_sq - mean_age_sq_LD)",
  laying_Model4_IntALR    = "LD ~ f_min_age * f_ALR + I(f_min_age^2) * f_ALR",
  laying_Model5_IntCenter = "LD ~ mean_age_LD * delta_age_LD + mean_age_LD * I(age_sq - mean_age_sq_LD)")
tit_fit <- lapply(tit_fixed, function(f) glmmTMB(as.formula(paste(f, tit_env, "+", tit_re)), dataPROCBLD, family = gaussian))

# ---------------------------------------------------------------- model comparison tables ---------------------------
aic_table <- function(fits, label) {
  a <- vapply(fits, AIC, numeric(1)); i4 <- grep("Model4", names(fits))
  data.frame(Dataset = label, Model = sub(".*Model([0-9]).*", "Model \\1", names(fits)), AIC = round(a, 2), dAIC_vs_Model4 = round(a - a[i4], 2),
             dAIC_vs_best = round(a - min(a), 2), valid_Hessian = vapply(fits, function(m) isTRUE(m$sdr$pdHess), logical(1)), row.names = NULL)
}
tab <- rbind(aic_table(fly_fit, "Fruit fly fecundity"), aic_table(tit_fit, "Great tit laying date"))
say("Random effects: fly %s; great tit %s", fly_re, tit_re)
print(tab, row.names = FALSE)
utils::write.csv(tab, file.path(out_dir, sprintf("figure7_model_comparison_%s.csv", re_mode)), row.names = FALSE)
coef_table <- function(m, label) {
  s <- summary(m)$coefficients$cond
  data.frame(Dataset = label, Term = rownames(s), Estimate = signif(s[, 1], 4), SE = signif(s[, 2], 4), z = round(s[, 3], 2), p = signif(s[, 4], 3), row.names = NULL)
}
coefs <- rbind(coef_table(fly_m4_scaled, "Fruit fly, Model 4 (ALR and age standardised)"), coef_table(tit_fit$laying_Model4_IntALR, "Great tit, Model 4"))
utils::write.csv(coefs, file.path(out_dir, sprintf("figure7_model4_coefficients_%s.csv", re_mode)), row.names = FALSE)
vc <- function(m) { v <- VarCorr(m)$cond; paste(vapply(names(v), function(g) paste0(g, ": SD ", paste(signif(attr(v[[g]], "stddev"), 3), collapse = ", ")), character(1)), collapse = "; ") }
say("Fly Model 4 random-effect SDs: %s", vc(fly_fit$Model4_IntALR)); say("Great tit Model 4 random-effect SDs: %s", vc(tit_fit$laying_Model4_IntALR))

# ---------------------------------------------------------------- predictions (authors' grid) -----------------------
cross_ages <- function(covars, ages, age_name) { out <- covars[rep(seq_len(nrow(covars)), times = length(ages)), , drop = FALSE]; out[[age_name]] <- rep(ages, each = nrow(covars)); rownames(out) <- NULL; as.data.frame(out) }
grid_predict <- function(model, newdata, type = "response") as.numeric(predict(model, newdata = newdata, re.form = NA, type = type, allow.new.levels = TRUE))
fly_ind <- dataAMNAT %>% distinct(F1_ID, .keep_all = TRUE)
fly_ages <- sort(unique(dataAMNAT$F1_Age))
fly_pred_frame <- cross_ages(fly_ind %>% select(Rep, Paternal_age, Paternal_sperm_age), fly_ages, "F1_Age") %>%
  mutate(Age2 = F1_Age^2, ALR = mean(fly_ind$ALR, na.rm = TRUE), Mean_age = mean(fly_ind$Mean_age, na.rm = TRUE), Mean_age2 = mean(fly_ind$Mean_age2, na.rm = TRUE),
         Delta_age = F1_Age - Mean_age, a = (F1_Age - fly_age_mu) / fly_age_sd, F0_ID = dataAMNAT$F0_ID[1], F1_ID = dataAMNAT$F1_ID[1]) %>% add_scaled()
for (k in 1:5) fly_pred_frame[[paste0("Model_", k)]] <- grid_predict(fly_fit[[k]], fly_pred_frame)
fly_grid_long <- fly_pred_frame %>% group_by(F1_Age) %>% summarise(across(Model_1:Model_5, mean), .groups = "drop") %>%
  pivot_longer(-F1_Age, names_to = "Model", values_to = "Predicted_Count")
tit_ind <- dataPROCBLD %>% distinct(female, .keep_all = TRUE)
tit_ages <- sort(unique(dataPROCBLD$f_min_age))
tit_pred_frame <- cross_ages(dataPROCBLD %>% add_count(female, name = "n_rec") %>% mutate(w = 1 / n_rec) %>% select(loc_density, pred, f_status, YR_FL, w), tit_ages, "f_min_age") %>%
  mutate(age_sq = f_min_age^2, f_ALR = mean(tit_ind$f_ALR, na.rm = TRUE), mean_age_LD = mean(tit_ind$mean_age_LD, na.rm = TRUE),
         mean_age_sq_LD = mean(tit_ind$mean_age_sq_LD, na.rm = TRUE), delta_age_LD = f_min_age - mean_age_LD, a = (f_min_age - tit_age_mu) / tit_age_sd,
         female = dataPROCBLD$female[1], year = dataPROCBLD$year[1], area = dataPROCBLD$area[1])
for (k in 1:5) tit_pred_frame[[paste0("Model_", k)]] <- grid_predict(tit_fit[[k]], tit_pred_frame)
tit_grid_long <- tit_pred_frame %>% group_by(f_min_age) %>% summarise(across(Model_1:Model_5, ~ stats::weighted.mean(.x, w)), .groups = "drop") %>%
  pivot_longer(-f_min_age, names_to = "Model", values_to = "Predicted_Val")
utils::write.csv(fly_grid_long, file.path(out_dir, sprintf("figure7A_predictions_%s.csv", re_mode)), row.names = FALSE)
utils::write.csv(tit_grid_long, file.path(out_dir, sprintf("figure7C_predictions_%s.csv", re_mode)), row.names = FALSE)

# ---------------------------------------------------------------- plots (authors' styling) --------------------------
model_cols <- c(Model_1 = "darkgrey", Model_2 = "hotpink1", Model_3 = "darkgreen", Model_4 = "orange", Model_5 = "steelblue2")
model_labs <- c(Model_1 = "Model 1", Model_2 = "Model 2", Model_3 = "Model 3", Model_4 = "Model 4", Model_5 = "Model 5")
model_ltys <- c(Model_1 = "dashed", Model_2 = "dashed", Model_3 = "dashed", Model_4 = "dashed", Model_5 = "dashed")
leg_theme <- theme(legend.position = "top", legend.title = element_text(size = 11), legend.text = element_text(size = 10), legend.key.size = unit(0.7, "lines"),
                   legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -8, 0), legend.box.spacing = unit(2, "pt"),
                   legend.spacing.y = unit(0, "pt"), axis.title = element_text(size = 16))
p1 <- ggplot(fly_grid_long, aes(x = F1_Age, y = Predicted_Count, colour = Model)) + geom_line(aes(linetype = Model), linewidth = 1) + geom_point(size = 1.8) +
  scale_colour_manual(values = model_cols, labels = model_labs, breaks = names(model_cols)) + scale_linetype_manual(values = model_ltys, guide = "none") +
  labs(x = "Age (days)", y = "Predicted fecundity", colour = "Model") + theme_bw(base_size = 13) + leg_theme +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) + scale_x_continuous(breaks = seq(4, 65, 14))
n_min <- 3
fly_brks <- seq(min(dataAMNAT$ALR, na.rm = TRUE), max(dataAMNAT$ALR, na.rm = TRUE), length.out = 5)
dataAMNAT_complete <- dataAMNAT %>% mutate(lifespan_bin = cut(ALR, breaks = fly_brks, include.lowest = TRUE, labels = paste0(round(head(fly_brks, -1)), "–", round(fly_brks[-1]))))
fly_bin_means <- dataAMNAT_complete %>% group_by(lifespan_bin, F1_ID, F1_Age) %>% summarise(y = mean(F2_Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(lifespan_bin, F1_Age) %>% summarise(mean_y = mean(y), n = n(), .groups = "drop") %>% filter(n >= n_min)
p2 <- ggplot(dataAMNAT_complete, aes(x = F1_Age, y = F2_Count, colour = lifespan_bin)) +
  geom_line(data = fly_bin_means, aes(x = F1_Age, y = mean_y, colour = lifespan_bin), inherit.aes = FALSE, linewidth = 1) +
  geom_point(data = fly_bin_means, aes(x = F1_Age, y = mean_y, colour = lifespan_bin), inherit.aes = FALSE, size = 2.4) +
  labs(x = "Age (days)", y = "Fecundity", colour = "ALR bin") + theme_bw(base_size = 13) + leg_theme +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) + scale_x_continuous(breaks = seq(4, 65, 14))
p3 <- ggplot(tit_grid_long, aes(x = f_min_age, y = Predicted_Val, colour = Model)) + geom_line(aes(linetype = Model), linewidth = 1) + geom_point(size = 1.8) +
  scale_colour_manual(values = model_cols, labels = model_labs, breaks = names(model_cols)) + scale_linetype_manual(values = model_ltys, guide = "none") +
  labs(x = "Age (years)", y = "Predicted laying date", colour = "Model") + theme_bw(base_size = 13) + leg_theme + theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(1, 9, 1))
tit_brks <- seq(0, ceiling(max(dataPROCBLD$f_ALR, na.rm = TRUE) / 2) * 2, by = 2)
dataPROCB_complete <- dataPROCBLD %>% mutate(f_ALR_bin = droplevels(cut(f_ALR, breaks = tit_brks, include.lowest = TRUE, labels = paste0(head(tit_brks, -1) + 1, "–", tit_brks[-1]))))
tit_bin_means <- dataPROCB_complete %>% group_by(f_ALR_bin, female, f_min_age) %>% summarise(y = mean(LD, na.rm = TRUE), .groups = "drop") %>%
  group_by(f_ALR_bin, f_min_age) %>% summarise(mean_y = mean(y), n = n(), .groups = "drop") %>% filter(n >= n_min)
p4 <- ggplot(dataPROCB_complete, aes(x = f_min_age, y = LD, colour = f_ALR_bin)) +
  geom_line(data = tit_bin_means, aes(x = f_min_age, y = mean_y, colour = f_ALR_bin), inherit.aes = FALSE, linewidth = 1) +
  geom_point(data = tit_bin_means, aes(x = f_min_age, y = mean_y, colour = f_ALR_bin), inherit.aes = FALSE, size = 2.4) +
  labs(x = "Age (years)", y = "Laying date", colour = "ALR bin") + theme_bw(base_size = 13) + leg_theme +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) + scale_x_continuous(breaks = seq(1, 9, 1))
composite_plot <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "A") & theme(plot.margin = margin(2, 4, 2, 2))
ggsave(file.path(out_dir, sprintf("figure7_%s.png", re_mode)), composite_plot, width = 12, height = 8.5, dpi = 300)
ggsave(file.path(out_dir, sprintf("figure7_%s.pdf", re_mode)), composite_plot, width = 12, height = 8.5)
writeLines(log_lines, file.path(out_dir, sprintf("figure7_notes_%s.txt", re_mode)))
cat("Saved to", out_dir, "\n")
