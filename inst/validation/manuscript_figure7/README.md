# Figure 7 of Sanghvi et al. 2026 with random slopes (`dev/figure7_random_slopes.R`)

Refit of the two empirical analyses with random slopes, following the authors' code (OSF kevnm, 'Empirical data' and
'Model predictions'): fixed effects, data preparation, prediction grid (proxies at their individual-level means,
covariates marginalised) and plotting unchanged. Only the random effects differ:

- fruit fly: `(1 | F0_ID) + diag(1 + a + I(a^2) | F0_ID:F1_ID)` - random intercept, slope and curvature per daughter
  (uncorrelated; a = standardised age), random intercept per father;
- great tit: `diag(1 + a | female) + (1 | year) + (1 | area)` - random intercept and slope per female. A curvature term
  cannot be estimated (half the females have one laying-date record; 11% have three or more).

Each fly model was fitted on the raw scale and in an exactly equivalent rescaled form, each with BFGS (the authors'
optimiser) and nlminb, and the highest-likelihood fit with a valid Hessian kept: with random slopes, the raw-scale fits
of Models 3 and 5 stop about 2,500 log-likelihood units below the maximum. All kept fits have valid Hessians.
`*_intercept.*` files reproduce the published analysis (random intercepts only) for comparison.

| ΔAIC vs Model 4 | Model 1 | Model 2 | Model 3 | Model 5 |
|---|---|---|---|---|
| Fly, random slope and curvature | 36.05 | 20.99 | 30.51 | 3.50 |
| Fly, published (random intercept) | 38.33 | 24.05 | 27.27 | 2.66 |
| Great tit, random slope | 8.66 | 6.44 | 50.98 | 43.91 |
| Great tit, published (random intercept) | 11.77 | 9.43 | 53.29 | 48.62 |

Model 4 is best in both datasets with either random structure. Fly Model 4 (standardised ALR and age): ALR × age
z = 4.78, p = 1.7e-6; ALR × age² z = 3.24, p = 0.001. Great tit Model 4: ALR × age p = 0.18; ALR × age² p = 0.004.
Panels B and D (raw means by ALR bin) do not depend on the model and are unchanged.
