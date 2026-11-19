# Does the age x ALR interaction need random slopes? (`tests/scripts/interaction_false_positives.R`)

Design of Sanghvi et al. 2026 (complete-sampling simulation; step 1's paper simulation): 250 individuals, lifespan
25 (SD 5), quadratic trajectories whose intercept, slope and shape vary among individuals (CV 0.2) and are correlated
with lifespan as set (0 or 0.5); 25 replicates per scenario; models fitted by ML in the app's engine. "Declared" means
the interaction model is more than 2 AIC better than its additive counterpart (Model 4 vs Model 2; Model 5 vs Model 3
gives the same percentages).

| Selective disappearance | Random effects | Interaction declared | median ΔAIC (M2 - M4) |
|---|---|---|---|
| none | random intercept | **88%** (false positive) | 24.2 |
| none | + random slope | **60%** (false positive) | 3.8 |
| none | + random slope and shape | **4%** (false positive) | -3.0 |
| intercept only (age-independent) | random intercept / + slope / + slope and shape | 88% / 60% / 4% (false positive) | 24.2 / 3.8 / -3.0 |
| slope (age-dependent) | random intercept / + slope / + slope and shape | 100% / 100% / 100% (power) | 537 / 44 / 62 |
| shape (age-dependent) | random intercept / + slope / + slope and shape | 100% / 100% / 100% (power) | 514 / 417 / 53 |

The intercept-only scenario reuses the seeds of the scenario without selective disappearance, and ALR absorbs an
intercept-lifespan correlation in both models, so its ΔAICs match to within 0.12.

**Reading.** When individuals differ in their ageing trajectories (as in every simulated scenario of the paper), an
interaction model fitted with random intercepts only absorbs that heterogeneity and is preferred 88% of the time even
with no selective disappearance; a random slope alone halves the excess but leaves the heterogeneity in shape, so the
false-positive rate is still 60%; with random slope and shape terms it falls to 4%, close to what is expected for a
2-AIC rule on two extra parameters (about 5%), while true age-dependent selective disappearance is still detected in
every dataset. The random-effect structure should match the ageing function: for a quadratic function, random slopes
for both age terms.

**The two empirical examples** (direct glmmTMB / lme4 refits with uncorrelated random slopes on scaled age):

- Fruit fly (zero-inflated negative binomial): random intercept, M2 - M4 = 24.0 AIC; with random slope and shape
  terms (which lower the AIC of both models by 12-15), M2 - M4 = 21.0. The age-dependent selective disappearance is
  robust to the random-effect structure. (With a random slope only, Model 4 did not converge.)
- Great tit laying date (Gaussian): random intercept, M2 - M4 = 9.4; with a random slope (SD 3.2, lowering the AIC by
  about 10), M2 - M4 = 6.4. Robust but weaker. A random shape term cannot be estimated: half the females have one
  record and 11% have three or more.

## Are the interaction terms themselves 'significant'? (`interaction_terms_without_slopes.txt`)

With no selective disappearance and random intercepts only, the interaction terms of Model 4 are jointly significant
(Wald p < 0.05) in 88% of datasets and at least one is individually significant in 92% (median |z| 3.0 and 3.9): a
false positive for the effect, not only a better fit. Their signs are random (linear term positive in 52%, quadratic in
32%): unmodelled individual slopes make the standard errors far too small. The population trajectory at the mean ALR is
unaffected (mean deviation from the truth 1.03% for Model 4, 1.05% for Model 2); the error lies in the ALR-specific
trajectories, i.e. in the claim of age-dependent selective disappearance. With random slope and shape terms: 4% jointly
significant, 8% for either term, median |z| 0.6 and 0.5.

## The methods paper's nine scenarios (`dev/random_slope_shape_simulation_chunk.R`)

The Rmd chunk written for the methods paper (its own simulation code, not the app's engine): 9 of the 27 scenarios, 100
individuals, 25 replicates, models 1-6 with a random intercept, + random slope, + random slope and shape (uncorrelated,
standardised age); all 4,050 fits converged. Interaction declared (Model 4 > Model 2 by >2 AIC) without age-dependent
selective disappearance: no selection 76% / 68% / 0%; intercept + 92% / 56% / 0%; intercept - 88% / 48% / 4%. With
age-dependent selection (slope or shape linked to lifespan): 92-100% with every random structure.
(`random_slope_shape_nine_scenarios.txt`, `random_slope_shape_nine_scenarios.png`)

## The GAMM step (`tests/scripts/gamm_proxy_false_positives.R`)

Same simulation with no selective disappearance (250 individuals), step 7's proxy comparison: the age-dependent ALR
form, s(age, by = ALR), was preferred in 12/12 datasets with a random intercept, 7/12 with a random intercept and
slope, and 0/8 with a random smooth per individual. The app's GAMM guide therefore asks for a random smooth per
individual before an age-dependent proxy effect is interpreted.
