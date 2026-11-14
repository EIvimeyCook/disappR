# Full test sweep, disappR 0.25.1 (GAMM branch), 14 Nov 2026

## Every bundled example through steps 1-7 in the app server (tests/scripts/app_examples_test.R)

```
fly                      records  1901  step 5 Model 4   GAMM bam, 1901 rows     shape Cubic age      saved 3   115 s
bichet                   records  1097  step 5 Model 4   GAMM gam, 1023 rows     shape Linear age     saved 3    53 s
bichet_marmot            records   173  step 5 Model 2   GAMM gam, 169 rows      shape Smooth of age  saved 3    30 s
moullec_swift            records  2084  step 5 Model 8   GAMM bam, 2069 rows     shape Smooth of age  saved 3    59 s
pasztor_apollo           records  3888  step 5 Model 4   GAMM bam, 3150 rows     shape Smooth of age  saved 3    56 s
wynn                     records  2689  step 5 Model 3   GAMM gam, 2689 rows     shape No age effect  saved 3    70 s
sanghvi_female           records  6001  step 5 Model 5   GAMM bam, 5997 rows     shape Smooth of age  saved 3   294 s
allain                   records   487  step 5 Model 6   GAMM gam, 487 rows      shape Smooth of age  saved 3    48 s
bouwhuis                 records  7341  step 5 Model 4   GAMM gamm4, 7126 rows   shape Cubic age      saved 3   785 s
warner                   records  2205  step 5 Model 4   GAMM bam, 1839 rows     shape Cubic age      saved 3    42 s
mckennaell_breeding      records  3173  step 5 Model 5   GAMM bam, 3173 rows     shape Quadratic age  saved 3   159 s
mckennaell_weight        records  2317  step 5 Model 3   GAMM bam, 2317 rows     shape Cubic age      saved 3    48 s
szejnersigal_activity    records   808  step 5 Model 6   GAMM gam, 808 rows      shape Smooth of age  saved 3    35 s
tsui_mass                records  3705  step 5 Model 6   GAMM bam, 2559 rows     shape Quadratic age  saved 3    65 s
brown_haematocrit        records  1383  step 5 Model 4   GAMM bam, 1379 rows     shape Smooth of age  saved 3    57 s
cooper_wren_females      records  1885  step 5 Model 4   GAMM gamm4, 1361 rows   shape Cubic age      saved 3    78 s
cooper_wren_males        records  2310  step 5 Model 2   GAMM gamm4, 2310 rows   shape Quadratic age  saved 3   201 s
lalande_elephants        records  3886  step 5 Model 4   GAMM bam, 3886 rows     shape Cubic age      saved 3    85 s
```

One output error, fly b2_plot, was a test artefact (the Modelling inputs were not yet set for the first example); fixed in the test, confirmed working.

## Test tiers (tests/scripts/run_all.R and re-runs after fixes)

| Tier | Result |
|---|---|
| 0 Unit tests | 2,104 pass; 1 known failure (test-combinations.R:142, beta-binomial cubic fit with an invalid Hessian; also fails on unmodified 0.24.9) |
| 1 Smoke | 305 checks, no errors (after fixing two script bugs present since 0.24.9) |
| 3 Robustness | passed |
| 4 App server | passed |
| 5 Stress (out-of-sample data) | passed |
| 6 Empirical and simulations | 35 pass, 3 fail: simulated-scenario ranking expectations (S1, S3, S5a); unmodified 0.24.9 fails 4 of the same checks |
| 7 Sweep of 0.20.9-0.20.19 changes | passed (after fixing two script bugs present since 0.24.9) |
| 8 Complex scenarios | passed |
| 9 Trait-specific ALR/AFR snapshot | passed (no errors) |
| 10 GAMMs against GLMMs | passed |
| 11 GAMM stress and break tests | passed |
| 12 GAMM replication | passed |

## Problems found and fixed in this sweep

- Lalande elephants froze the hazard (life-table) panel: continuous ages gave ~235 age classes and an hours-long frailty
  GLMM. Now at most ~40 classes; the frailty model is skipped on very large tables.
- Fly (negative binomial) GAMM comparisons took hours in mgcv ML. Now gamm4 with theta fixed at the REML estimate (~20 s).
- An older installed disappR could be used by the new app ("unused argument" errors): build stamp check added (0.25.1).
