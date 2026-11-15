library(MASS)
library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)

# Random slope and curvature terms in the quadratic design: models 1-6 fitted with (i) a random intercept only,
# (ii) a random intercept and random slope, (iii) a random intercept, slope and curvature (all uncorrelated, on
# standardised age). 9 of the 27 scenarios, complete sampling, 100 individuals and 25 replicates per scenario.
# Fixed effects are entered on a rescaled but equivalent scale (age/10, age^2/100, ALR/10, ...) for convergence;
# this does not change the AICs.

run_random_structure_scenario <- function(cor_lifespan_intercept, cor_lifespan_slope, cor_lifespan_quad_slope, replicate,
                                          n_individuals = 100,
                                          mean_lifespan = 25, sd_lifespan = 0.2 * 25,
                                          mean_intercept_fert = 600, sd_intercept_fert = 0.2 * 600,
                                          mean_slope_fert = 9, sd_slope_fert = 0.2 * abs(9),
                                          mean_quad_slope_fert = -0.3, sd_quad_slope_fert = 0.2 * abs(-0.3)) {

  mu <- c(mean_lifespan, mean_intercept_fert, mean_slope_fert, mean_quad_slope_fert)
  sds <- c(sd_lifespan, sd_intercept_fert, sd_slope_fert, sd_quad_slope_fert)
  R <- matrix(c(
    1,                       cor_lifespan_intercept, cor_lifespan_slope, cor_lifespan_quad_slope,
    cor_lifespan_intercept,  1,                      0,                  0,
    cor_lifespan_slope,      0,                      1,                  0,
    cor_lifespan_quad_slope, 0,                      0,                  1
  ), nrow = 4, byrow = TRUE)
  Sigma <- diag(sds) %*% R %*% diag(sds)

  ind <- as_tibble(MASS::mvrnorm(n = n_individuals, mu = mu, Sigma = Sigma), .name_repair = "minimal") %>%
    setNames(c("lifespan_raw", "intercept", "slope", "quad_slope")) %>%
    mutate(id = 1:n_individuals, lifespan = pmax(1, round(lifespan_raw)))

  d <- expand_grid(id = 1:n_individuals, age = 1:max(ind$lifespan)) %>%
    left_join(ind, by = "id") %>%
    filter(age <= lifespan) %>%
    mutate(fertility = pmax(0, intercept + slope * age + quad_slope * age^2) + rnorm(n(), 0, 5)) %>%
    group_by(id) %>%
    mutate(ALR = max(age), mean_age = mean(age), mean_age2 = mean(age^2)) %>%
    ungroup() %>%
    mutate(id = factor(id),
           a = age / 10, a2 = age^2 / 100, alr = ALR / 10, ls = lifespan / 10,
           ma = mean_age / 10, da = (age - mean_age) / 10, da2 = (age^2 - mean_age2) / 100,
           age_s = as.numeric(scale(age)))

  fixed <- c(
    "Model 1" = "fertility ~ a + a2",
    "Model 2" = "fertility ~ a + a2 + alr",
    "Model 3" = "fertility ~ ma + da + da2",
    "Model 4" = "fertility ~ a * alr + a2 * alr",
    "Model 5" = "fertility ~ ma * da + da2 * ma",
    "Model 6" = "fertility ~ a * ls + a2 * ls")
  random <- c(
    "Random intercept"                   = "(1 | id)",
    "+ random slope"            = "(1 | id) + (0 + age_s | id)",
    "+ random slope and shape"  = "(1 | id) + (0 + age_s | id) + (0 + I(age_s^2) | id)")
  ctrl <- lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore")

  out <- expand_grid(model = names(fixed), random_structure = names(random))
  out$AIC <- mapply(function(m, r) {
    fit <- tryCatch(lmer(as.formula(paste(fixed[[m]], "+", random[[r]])), data = d, REML = FALSE, control = ctrl), error = function(e) NULL)
    if (is.null(fit)) NA_real_ else AIC(fit)
  }, out$model, out$random_structure)
  out %>% mutate(cor_lifespan_intercept = cor_lifespan_intercept, cor_lifespan_slope = cor_lifespan_slope,
                 cor_lifespan_quad_slope = cor_lifespan_quad_slope, replicate = replicate)
}

nine_scenarios_rs <- tibble::tribble(
  ~cor_lifespan_intercept, ~cor_lifespan_slope, ~cor_lifespan_quad_slope,
   0,    0,    0,
   0.5,  0,    0,
  -0.5,  0,    0,
   0,    0.5,  0,
   0,   -0.5,  0,
   0,    0,    0.5,
   0,    0,   -0.5,
   0.5,  0.5,  0.5,
  -0.5, -0.5, -0.5)

n_replicates_rs <- 25
grid_rs <- tidyr::crossing(nine_scenarios_rs, replicate = 1:n_replicates_rs)

n_cores_rs <- if (.Platform$OS.type == "windows") 1 else max(1, parallel::detectCores() - 1)
set.seed(42)
RNGkind("L'Ecuyer-CMRG")
results_rs <- bind_rows(parallel::mclapply(seq_len(nrow(grid_rs)), function(k) {
  g <- grid_rs[k, ]
  run_random_structure_scenario(g$cor_lifespan_intercept, g$cor_lifespan_slope, g$cor_lifespan_quad_slope, g$replicate)
}, mc.cores = n_cores_rs, mc.set.seed = TRUE))

sym <- function(x) ifelse(x > 0, "+", ifelse(x < 0, "-", "0"))
results_rs <- results_rs %>%
  mutate(scenario = paste(sym(cor_lifespan_intercept), sym(cor_lifespan_slope), sym(cor_lifespan_quad_slope)),
         scenario = factor(scenario, levels = with(nine_scenarios_rs, paste(sym(cor_lifespan_intercept), sym(cor_lifespan_slope), sym(cor_lifespan_quad_slope)))),
         random_structure = factor(random_structure, levels = c("Random intercept", "+ random slope", "+ random slope and shape"))) %>%
  group_by(scenario, random_structure, replicate) %>%
  mutate(dAIC = AIC - min(AIC, na.rm = TRUE)) %>%
  ungroup()

aic_summary_rs <- results_rs %>%
  group_by(scenario, random_structure, model) %>%
  summarise(mean_dAIC = mean(dAIC, na.rm = TRUE), se_dAIC = sd(dAIC, na.rm = TRUE) / sqrt(sum(!is.na(dAIC))), .groups = "drop")

interaction_support_rs <- results_rs %>%
  select(scenario, random_structure, replicate, model, AIC) %>%
  pivot_wider(names_from = model, values_from = AIC) %>%
  group_by(scenario, random_structure) %>%
  summarise(`Model 4 better than Model 2 by >2 AIC (%)` = round(100 * mean(`Model 2` - `Model 4` > 2, na.rm = TRUE)),
            `Model 5 better than Model 3 by >2 AIC (%)` = round(100 * mean(`Model 3` - `Model 5` > 2, na.rm = TRUE)),
            .groups = "drop")
print(interaction_support_rs, n = Inf)
cat("Model fits that failed:", sum(is.na(results_rs$AIC)), "of", nrow(results_rs), "\n")

model_cols_rs <- c("Model 1" = "darkgrey", "Model 2" = "hotpink1", "Model 3" = "darkgreen",
                   "Model 4" = "orange", "Model 5" = "steelblue2", "Model 6" = "purple")

plot_random_structure <- ggplot(aic_summary_rs, aes(x = model, y = mean_dAIC, colour = model)) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_dAIC - se_dAIC, ymax = mean_dAIC + se_dAIC), width = 0.3) +
  facet_grid(random_structure ~ scenario, scales = "free_y") +
  scale_colour_manual(values = model_cols_rs) +
  labs(x = NULL, y = expression(Delta * "AIC (mean " %+-% " SE)"), colour = "Model") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top",
        strip.text.y = element_text(size = 9))
print(plot_random_structure)
