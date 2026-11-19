# Validation of the GAMM step (step 7, disappR 0.25.1)

Every result here comes from the real engines: mgcv (or gamm4) for the GAMMs and lme4 for the GLMMs. The 0.23 layer
was validated on a Python re-implementation of a GAMM. That port modelled age-dependent selective disappearance
differently (`s(age, by = ALR)`) from the R code (`ti(age, ALR)`), so its results did not validate the app.

| File | Script | What it checks |
|---|---|---|
| `replication_results.txt` | `tests/scripts/gamm_replication.R` | Two published GAMM analyses, from the authors' data and model code, against the app's engine |
| `gamm_vs_glmm_results.txt`, `gamm_vs_glmm_replicates.csv` | `tests/scripts/gamm_vs_glmm.R` | GAMMs against the parametric models of step 5 on 240 simulated datasets with a known answer |
| `gamm_vs_glmm_empirical.txt` | the same, `DISAPPR_GVG_REPS=0` | The same comparison on the bundled Gaussian examples |
| `stress_results.txt` | `tests/scripts/gamm_stress.R` | Every bundled example, every out-of-sample dataset in `tests/stress_data`, and 31 awkward datasets |
| `out_of_sample_results.txt` | `tests/scripts/gamm_out_of_sample.R` | The app's simulator (56 datasets with a known true trajectory) and the 13 parametric-study examples through every GAMM option |

## Replication (bundled examples `tsui_mass` and `brown_haematocrit`)

### Tsui et al. 2026, *Evol. Lett.*: Seychelles warbler body mass

Run on the authors' data and `model1.R`.

| Check | Paper | Authors' code, refitted | disappR |
|---|---|---|---|
| Records / birds | — | 2,560 / 1,172 | 2,559 / 1,172 |
| Age smooth | ageing in body mass | edf 3.05, F = 4.46, p = 0.002 | edf 3.03, F = 4.50, p = 0.002 |
| Lifespan, per year | linear term | 0.0252 (SE 0.0055) | 0.0244 (SE 0.0055) |
| Sex curves | G / GS / S chosen by AIC | S best; GS +4.5, G +4.3 | S best, G +5.7 (mgcv: choose the bam engine); in fast mode, the default for these 1,172 birds, G best (GS and S cannot be fitted) |
| Onset of decline | 4.7 years (limits 1–8.5) | — | final decline from 4.9 years (95% limits 3.5–14.8) |
| Age × lifespan | reported as significant for body mass | `ti` p = 0.056 | supported: `s(age, by = LS)` best, linear lifespan +3.0 (marginal likelihood, fast mode); before 0.25.1 the mgcv comparison by conditional AIC found no support (`ti` +8.5) |

- **Records:** one record differs because its observer is blank. disappR treats a blank as missing; the authors'
  `factor()` kept it as a level.
- **Onset of decline:** with pointwise intervals, decline is significant from 8.6 to 12.2 years. With simultaneous
  intervals no period is significant.
- **Age × lifespan:** the only result that does not match the paper. Both the authors' formula refitted here and the
  app find it weak or absent. It may come from a different specification in the authors' AIC scripts, which we did
  not run.

### Brown et al. 2021, *Ecol. Evol.*: Seychelles warbler haematocrit

| Check | Paper (gamm4) | Authors' specification refitted | disappR |
|---|---|---|---|
| Records / birds | 1,379 / 730 | 1,379 / 730 | 1,379 / 730 |
| Age smooth | edf 7.805, F = 16.73, p < .001 | gamm4 edf 7.81, F = 17.5; mgcv edf 7.81, F = 15.5 | edf 7.81, F = 15.5 |
| Peak | about 1.5 years | — | 1.62 years (95% limits 1.44–1.88) |
| Selective disappearance | none detected | — | ALR p = 0.45: none |

## GAMMs against GLMMs (simulations)

The design has four true ageing shapes, three kinds of selective disappearance and 20 replicates each. Each approach
fits its version of Models 1, 2 and 4 by ML. "Age-dependent declared" means the age-dependent model is more than
2 AIC better than the additive one.

| Approach | Falsely age-dependent (no true age-dependence) | Power (true age-dependence) | Mean shape error |
|---|---|---|---|
| GAMM | 1.3% | 100% | 0.34 |
| GLMM, function chosen by AIC | 30% | 100% | 0.37 |
| GLMM, quadratic | 75% | 100% | 0.41 |

The GAMM is the most reliable of the three. Up to 0.25.0 its false-positive rate was 11%: the proxy comparison used
mgcv's AIC, which is conditional on the random effects. Since 0.25.1 it uses the marginal likelihood (as step 5's
mixed models and gamm4 do), and the rate is 1.3% (3 of 240 datasets) with unchanged power. The conditional AIC also
missed age-independent selective disappearance in simulated data (see Out-of-sample tests): the individual random
intercepts absorb an individual-level term such as ALR almost for free.

## Performance

These figures are for the great tit example: 7,126 records from 4,848 birds.

- **mgcv:** every individual's random effect is a dense column. A `bam` fit took 202 s and 5.6 GB, and `summary()`
  with random-effect tests another 151 s.
- **gamm4:** lme4's sparse random effects fit the same model in 97 s and 1.4 GB. The app uses gamm4 automatically
  above 1,500 individuals when the model allows it.
- **Summaries:** `summary()` now skips mgcv's random-effect tests. They are slow, and they failed on some ML fits.
  Variance components are reported instead.

All timings are from a 4-core container with OpenBLAS. R's reference BLAS, the default on Windows, is several
times slower for these models.

## Test results on the final code

| Tier | Result |
|---|---|
| Unit tests, GAMM (`tests/testthat/test-gamm.R`) | all pass (73 expectations) |
| App server test (`tests/scripts/app_test.R`) | all outputs render; GAMM fit, three comparisons, six saves; methods paragraph; changing a GAMM setting no longer alters the fingerprint of a step-5 result |
| Stress and break tests (`stress_results.txt`) | 144 fits succeed, 40 return an intended plain message, no error escapes the engine |
| Replication (`replication_results.txt`) | runs without errors; results above |
| Full unit suite | passes except `test-combinations.R:142` (a beta-binomial cubic fit with an invalid Hessian in step 5) |
| Every example through steps 1-7 (`tests/scripts/app_examples_test.R`) | all 27 examples and the 7 other published traits of the multi-trait files, with each trait's published settings |
| Out-of-sample tests (`out_of_sample_results.txt`) | 56 simulated datasets and 13 parametric-study examples under every GAMM option: no errors, no broken invariants |

The `test-combinations.R:142` failure also occurs on unmodified 0.24.9 in this environment (glmmTMB 1.1.15.2,
OpenBLAS), so it predates the GAMM step.

The break tests found and fixed these engine bugs:

- `bam` failed on small data with random effects (an mgcv quirk). Such fits now use `gam`.
- Extremely large traits gave an mgcv error. They now get a message to rescale.
- Data with a single individual gave a contrasts error. They now get a plain message.
- `gamm4` reports deviance explained as `numeric(0)`, which broke the reading and the comparison tables.
- A `chol` error escaped mgcv's random-effect tests in `summary()`.

## Out-of-sample tests (`out_of_sample_results.txt`, `tests/scripts/gamm_out_of_sample.R`)

Data the GAMM step was not built or tuned on, run as a user would, through every post-fit function (tables, model
output, guide, checks, curves, derivative, peak, residuals, cautions, R code) with invariants checked (finite curves,
intervals containing the fit, R code that parses). Fits run under the app's time limits.

**A. The app's simulator, 56 datasets with a known true trajectory**: every trait type (Gaussian, Poisson, negative
binomial, zero-inflated, the paper's simulation) and age function; selective disappearance none, age-independent,
age-dependent or both, in both directions; selective appearance; four missingness patterns; grouped data with a diet
covariate. All 56 fit, with no errors.

| Check | Result |
|---|---|
| Age-independent selective disappearance: ALR term in the simulated direction, p < 0.05 | 23 of 23 datasets |
| Age-dependent selective disappearance, median mean \|D\| from the true curve | linear ALR 9.8%; `ti(age, ALR)` 2.7%; best GLMM (Model 4) 1.4% |
| Proxy comparison names the simulated kind of selective disappearance | 15 of 18 (marginal-likelihood AIC) |

The tests found that the proxy comparison, then by mgcv's conditional AIC, could miss strong age-independent selective
disappearance: on one dataset the model without ALR was 1.4 AIC worse while the marginal log-likelihood differed by 100.
The individual random intercepts absorb an individual-level term almost for free. Since 0.25.1 the comparison uses the
marginal likelihood (`gamm_marginal_aic()`, identical to lme4's AIC for the same model), and the guide asks for a refit
when the comparison supports another proxy form than the one fitted. A linear ALR term does not correct age-dependent
selective disappearance (9.8% from the truth); the age-dependent form does (2.7%). The three misses are borderline:
two age-independent datasets called age-dependent (for one, step 5 also chose Model 4) and one dataset without
selection called age-independent (for it, step 5 also chose Model 4).

**B. The 13 examples analysed with parametric models (steps 1-6)**, each through 16-22 GAMM settings: ALR linear,
smooth, `s(age, by = ALR)` and `ti(age, ALR)`; no proxy; known lifespan; AFR; random slopes; no individual effect;
terminal record; observation-level random effect; k = 5; the bam and gamm4 engines; no covariates or random terms;
numeric covariates as smooths; structures GS, GI, S and I on a categorical covariate; and the shape, proxy and group
comparisons. No errors and no broken invariants; no fit beyond 300 s and no comparison beyond 1,800 s (the slowest,
the great tit proxy comparison on 7,341 records, took 12 minutes under load).

| Agreement with step 5 (Model 2 against the default GAMM) | Result |
|---|---|
| ALR effect in the same direction | 13 of 13 examples; never opposite and significant |
| Correlation of the population curves | median 0.97 |

The two low correlations are not GAMM failures: for the seed-beetle fecundity data (r = 0.19) the quadratic GLMM on the
log scale explodes at the oldest ages (beyond 10^8 eggs a day), and for the chipmunk data (r = -0.09) both curves rise
but the smooth (edf 4.8) has a different shape from the quadratic, as the shape comparison also prefers.

## Known limitations

- **ML comparisons with many random-effect levels are slow.** Proxy forms, group structures and polynomials are all
  compared by ML, and ML in mgcv cannot use `bam`'s discretisation. The haematocrit example (730 birds, 747 breed
  groups) took around 30 minutes for the eight proxy models on this container. In the app, a comparison running
  longer than 300 seconds is stopped with a message.
- **Age-dependent proxy terms need mgcv.** With more than 1,500 individuals they are left out of the proxy
  comparison, because gamm4 cannot fit `ti()`.
- **False positives.** The rule that ΔAIC > 2 declares age-dependent selective disappearance gives about 1-2% false
  positives in the simulations.

## Fast comparisons (gamm4) against mgcv

Above 300 random-effect levels, comparisons are fitted with gamm4. gamm4 uses lme4's sparse random effects and the
marginal likelihood, and its age-dependent form is `s(age, by = proxy)`. On the haematocrit example (1,491 levels)
the proxy comparison went from about 30 minutes to 10 s, polynomials to 1 s and sex groups to 6 s.

The same simulations in fast mode (`fast/`, forced with `DISAPPR_GAMM_FAST_LEVELS=0`):

| Approach | Falsely age-dependent | Power |
|---|---|---|
| GAMM, mgcv `ti()` (conditional AIC, up to 0.25.0) | 11% | 100% |
| GAMM, mgcv `ti()` (marginal likelihood, 0.25.1) | 1.3% | 100% |
| GAMM, gamm4 `s(age, by = ALR)` (marginal likelihood) | 1.9% | 100% |

The simulated age-dependence is linear in ALR, which matches the `s(age, by = ALR)` form, so this design favours it.

Group comparison, 10 replicates per scenario:

- Both engines detect groups that truly age differently 10 times out of 10.
- When groups differ only in level, gamm4 never declares a difference in shape; mgcv does so once in 10.

Fast mode can change conclusions on real data. On the Tsui body-mass data:

| Comparison | mgcv | gamm4 (fast mode) |
|---|---|---|
| Lifespan | linear best; `ti` +8.5 AIC (conditional AIC, up to 0.25.0) | `s(age, by = LS)` best; linear +3.0 AIC |
| Sex | S best; G +5.7 | G best; GI +4.0 |

The gamm4 lifespan result agrees with the paper's reported age × lifespan effect.

## Replication: superb fairy-wren and Asian elephant (`replication_results_2.txt`, `tests/scripts/gamm_replication_2.R`)

### Cooper et al. 2021, *Am. Nat.*: superb fairy-wren

Age-smooth edf. Significance stars: \*\* P < .01.

| Trait | Paper (Table 1) | Authors' model refitted (gamm4) | disappR engine |
|---|---|---|---|
| Female independent offspring | 4.540\*\* | 4.21, p < .001 | 4.86, p < .001 |
| Male total offspring | 3.567\*\* | 3.07, p < .001 | 4.51, p < .001 |
| Male extra-group offspring | 3.731\*\* | 3.13, p < .001 | 4.77, p < .001 |
| Male within-group offspring | 2.019, n.s. | 1.75, p = 0.26 | 2.66, p = 0.09 |
| Female survival | 1.000\*\* | 2.80, p < .001 | 3.39, p < .001 |

**Reproduction.** Every conclusion matches the paper: total and extra-group siring and female output change with age,
and within-group siring does not. Peaks from the disappR curves are 5.2 years (male total) and 6.0 years (extra-group).
The paper's break points were about 4.3 years (extra-group) and 5.7 years (within-group).

**Why the edfs differ.** The app's edfs are somewhat higher because it has no observation-level random effect for
counts. Choosing the negative binomial family is the app's alternative.

**Survival.** This does not replicate. The paper found a straight-line decline (edf 1.000); both the refitted model and
the app find a curve. The deposit has no survival column, so survival was derived from age and longevity. That coding
is the likely cause, and it needs the authors' code to resolve.

### Lalande et al. 2022, *J. Evol. Biol.*: Asian elephant body mass

| Check | Paper | Authors' model refitted (gamm4, ML) | disappR engine |
|---|---|---|---|
| By-sex smooth vs common smooth | better by 65.7 AICc | better by 67.3 AIC | group comparison: GI best, G 19.2 AIC worse |
| Female delta-age smooth edf | 6.899 | 7.07 | — |
| Male delta-age smooth edf | 6.096 | 5.49 | — |

The disappR group comparison ran in fast mode (gamm4), with GI best and I +1.3. Its fit uses a smooth of age per sex
with ALR (k = 15). The rise to a peak is at 55 years in females and 60 in males.

## Replication: the other Tsui et al. 2026 traits (`replication_results_tsui.txt`, `tests/scripts/gamm_replication_tsui.R`)

The deposit has data for 12 of the paper's 13 GAMMs (the telomere file is missing). Body mass is above; the other 11
are compared here: the authors' model G formula refitted with mgcv::gam (REML) against the disappR engine with the
example's settings and each trait's own covariates.

| Model | Trait (family) | Records | Age smooth edf, authors / disappR | p, authors / disappR |
|---|---|---|---|---|
| 2 | Wing length (Gaussian) | 2,623 / 2,622 | 3.77 / 3.78 | <0.001 / <0.001 |
| 3 | Tarsus (Gaussian) | 2,640 / 2,639 | 1.00 / 1.00 | 0.115 / 0.119 |
| 4 | Fat score (Poisson) | 1,453 / 1,453 | 1.00 / 1.00 | 0.251 / 0.250 |
| 5 | Oxidative damage, ROMs (Gaussian) | 639 / 639 | 1.99 / 1.99 | 0.241 / 0.240 |
| 6 | Antioxidant capacity, OXY (Gaussian) | 638 / 638 | 1.00 / 1.00 | 0.413 / 0.413 |
| 7 | Haematocrit (Gaussian) | 1,744 / 1,744 | 1.01 / 1.00 | 0.141 / 0.140 |
| 8 | White blood cells (Gaussian) | 1,744 / 1,744 | 4.48 / 4.48 | <0.001 / <0.001 |
| 9 | Malaria infection (binomial) | 3,197 / 3,197 | 5.47 / 5.35 | <0.001 / <0.001 |
| 11 | Provisioning rate (Gaussian) | 1,004 / 1,004 | 1.55 / 1.55 | 0.009 / 0.009 |
| 12 | Annual reproductive success (Poisson) | 5,176 / 5,176 | 9.66 / 9.61 | <0.001 / <0.001 |
| 13 | Annual survival (binomial) | 5,174 / 5,098 | 2.28 / 2.25 | 0.051 / 0.051 |

Every conclusion about the age smooth and the lifespan, ALR and AFR terms agrees. The one-record differences are the
bird with a blank observer, which disappR treats as missing. For survival, the example leaves out the 76 records whose
previous-year status is 'dead' (none survived: complete separation), which brings the age smooth closer to the published one. Differences in edf come from the engine (bam for large
data) and, for reproductive success, from using the ages of first and last reproduction where the authors used the
calendar years (the version with ages, refitted with the authors' code, gives edf 9.65).

## Replication: Sanghvi et al. 2021, *Evol. Ecol.*: seed beetles (`replication_results_beetles.txt`, `tests/scripts/gamm_replication_beetles.R`)

The paper compared group-specific age curves with one shared curve (beetle and container random effects, adult
lifespan as a fixed term), by AIC.

| Trait | Published | mgcv refit (ML) | disappR group comparison |
|---|---|---|---|
| Fecundity of mated females, by density | density-specific best, dAIC 27.5 | Gaussian: density-specific better by 26.7 (24.2 with lifespan) | GS best; one curve (G) 28.3 worse |
| Flying ability of virgins, by sex | sex-specific best, dAIC 21.8 (ordinal family) | Gaussian: sex-specific better by 36.0 | GI best; G 26.2 worse |
| Flying propensity of virgins, by sex | sex-specific best, dAIC 118 | negative binomial: sex-specific better by 94.2 | GI best; G 72.0 worse |

All three conclusions replicate. The fecundity result holds with a Gaussian family; with a Poisson or negative binomial
family one shared curve fits as well, so the published model was presumably Gaussian, and the example uses it.
Flying ability is ordinal in the paper (mgcv's ocat family, not offered by disappR) and continuous here. The flight
comparisons ran in fast mode (gamm4; GS and S not fitted) because the two flight examples have 377 beetles.
