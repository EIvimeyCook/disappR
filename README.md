# disappR

> **Status:** version 0.21.5. Run `Rscript tests/scripts/run_all.R` after any change to re-run the full test suite.
> for an analysis.

disappR is a Shiny app for estimating ageing from longitudinal data on individuals followed over time.

Ecologists want to know the average within-individual ageing pattern, but they only have population-level data on
individuals followed over time, and those data are usually incomplete. When an individual's phenotype is associated
with its entry into or removal from the sample (selective appearance and selective disappearance), standard analyses
can return biased ageing patterns.

disappR diagnoses these problems in your data, visually first and then through sampling and missingness, and fits
and compares the models that address them. It is meant as a first port of call for empiricists analysing
longitudinal data, especially on ageing. It covers the mixed models and diagnostics behind a publishable analysis,
with guidance at each step.

## Launch the app

You need R 4.1 or later. Download and unzip disappR, then run in R:

```r
# 1. Install the packages the app uses (once)
install.packages(c("shiny", "shinydashboard", "ggplot2", "lme4",               # required
                   "glmmTMB", "lmerTest", "nlme", "DHARMa", "performance",     # optional: every feature
                   "see", "codetools", "callr"))

# 2. Install disappR from the unzipped folder (once)
install.packages("path/to/disappR", repos = NULL, type = "source")

# 3. Launch the app (each session)
disappR::run_app()
```

To run without installing disappR, replace steps 2 and 3 with `shiny::runApp("path/to/disappR/inst/app")`.

The app opens in a browser window or the RStudio viewer. Without the optional packages it still runs, with fewer
analyses; the *Start here* page lists any that are missing.

| Optional package | Enables |
|---|---|
| glmmTMB | Count, proportion and binary traits |
| lmerTest | p-values for Gaussian models |
| nlme | The exponential ageing function |
| DHARMa | Simulation-based residual checks |
| performance, see | Model diagnostics and their plots |
| codetools | Exported R code |
| callr | Step 4 runs in a separate R process, so a failure there cannot stop the app |

## What the app does

The sidebar takes you through six steps: load the data, diagnose, model, and summarise.

### Load and check data (step 1)

* Use simulated data with a known answer, one of 13 datasets from 12 published studies, or your own CSV. Each published
  dataset opens with the settings that most closely reproduce the analysis in its paper.
* Your data need one row per individual per occasion, with columns for the individual, its age and the trait.
  Lifespan and age at first record are optional but make more models available.
* Map continuous and categorical covariates and their interactions, groups, and censoring.
* Check that ages, lifespans, age at first record and censoring are read correctly. An individual must have one
  lifespan, one age at last record and one age at first record: if its rows disagree, loading stops and names it.

### Diagnose what could bias the ageing pattern (steps 2 to 4)

Before any model is fitted, the app shows whether selection, sampling or the shape of ageing could bias the result.

| Source of bias | Level | What the app shows | Step |
|---|---|---|---|
| Selective disappearance and appearance | Population | Trait trajectories of lifespan or age-at-first-record groups and their differences across age; the trait against lifespan within age bins; the trait before death; selection differentials by age; the disappearance hazard | 2 · Visual diagnosis |
| Incomplete sampling | Data | The sampling grid; missingness by age and against the trait and other variables; agreement between lifespan proxies, to choose one for the models | 3 · Missingness and proxies |
| A misspecified ageing shape | Individual and population | A function fitted to each individual; the population curve rebuilt from the individual fits; ageing functions compared | 4 · Trajectories |

### Model the population ageing pattern (step 5)

* Fit ten mixed models, from no correction to models for selective disappearance and appearance acting on the level
  of the trait or on how it changes with age.
* Choose the error family (continuous traits, counts including overdispersed and zero-inflated counts, proportions,
  binary traits), random intercepts and slopes, the ageing function, and covariates and their interactions with age.
* Compare models by AIC and likelihood-ratio tests, and plot their predicted age trajectories.
* Check the fitted models: residual diagnostics, the ageing function, sensitivity to random slopes, effect sizes, and
  a null-model test of the lifespan effect.

### Summarise and report (step 6)

* Results you save are gathered into an automatic, cautious summary of the evidence for selective disappearance and
  appearance.
* Export an HTML or text report, the analysis data, and R code that reproduces the analysis outside the app.

### Learn with simulated data

* Generate datasets whose true ageing trajectory is known, controlling data structure, missingness, the organism's
  biology, sample size, and the type and strength of selection.
* See what selective disappearance looks like before meeting it in real data, and compare each model's estimate
  with the truth.

## Scripted use

The analyses also run from R without the app:

```r
library(disappR)
x   <- disappr_prepare(my_data, disappr_mapping(id = "ID", age = "age", trait = "mass", life = "lifespan"))
fit <- disappr_fit(x, models = c("M1", "M2", "M4"), age_function = "Quadratic")
disappr_compare(fit)
```

`vignette("disappR", package = "disappR")` lists the functions.

## Documentation

Every panel in the app has an ⓘ button explaining what it shows and how to read it. The vignettes cover the
workflow, the models, validation and replication, and development and testing: `browseVignettes("disappR")`.

## Testing

* `Rscript tests/scripts/run_all.R` runs every test tier.
* `Rscript tests/scripts/release_gate.R` also runs R CMD check.
* `python3 tools/static_audit/run_all.py` runs static checks that need no R installation.

## Citation

`citation("disappR")` gives both references:

* Sanghvi, K., Ivimey-Cook, E. 2026. disappR: a shiny app to model ageing and selective [dis]appearance.
* Sanghvi, K., Ivimey-Cook, E.R., Bouwhuis, S., Sepil, I. and van de Pol, M. 2026. A comparison of methods to assess
  selective disappearance and quantify ageing. EcoEvoRxiv.

## Maintainer and licence

E. R. Ivimey-Cook (E.Ivimey-Cook@uea.ac.uk). MIT licence.
