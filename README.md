<p align="center">
  <img src="inst/app/www/disappR_logo.png" width = "200"/>
</p>

<div align="center">
 <h1>disappR</h1>
</div>

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE.md)
<!-- badges: end -->

disappR is an R package with an interactive Shiny app for **diagnosing and modelling
selective disappearance and selective appearance** in longitudinal ageing data.

Ecologists want to know the average within-individual ageing pattern, but they only
have population-level data on individuals followed over time, and those data are
usually incomplete. When an individual's phenotype is associated with its entry into
or removal from the sample — selective appearance and selective disappearance —
standard analyses can return biased ageing patterns. disappR diagnoses these
problems in your data, visually first and then through sampling and missingness,
then fits and compares the mixed models that account for them. It is meant as a
first port of call for empiricists analysing longitudinal data, especially on
ageing, covering the diagnostics and models behind a publishable analysis with
guidance at every step.

## Features

- **Visual-first diagnosis.** Trait trajectories of lifespan or age-at-first-record
  groups and their differences across age, the trait against lifespan within age
  bins, the trait before death, selection differentials by age, and the
  disappearance hazard — before any model is fitted.
- **Sampling and missingness checks.** The expected-occasion grid, missingness by
  age and by other variables, and agreement between lifespan proxies, to choose one
  for the models.
- **Individual and population trajectories.** A function fitted to each individual
  (linear to a non-linear exponential), the population curve rebuilt from the
  individual fits, and ageing functions compared at the population level.
- **Ten comparative mixed models.** From a negative control to models for
  age-independent and age-dependent selective disappearance and appearance, with
  Gaussian, Poisson, negative binomial, zero-inflated, binomial and beta-binomial
  families, random slopes, fixed covariates and their interactions with age, and
  nested (up to three levels) or crossed random effects.
- **Teaching simulations with a known answer.** Choose the ageing form, the
  organism's biology, the sampling scheme, and the type and strength of selection,
  then see how each model's estimate compares with the truth.
- **13 published datasets to practise on**, from 12 studies in laboratory and wild
  populations, each pre-set to the analysis reported in its paper.
- **Reproducible output.** Residual diagnostics, a null-model test of the lifespan
  effect, downloadable R code for every analysis, and HTML or text reports of saved
  results.
- **Runs entirely on your machine.** Your data are read locally and never uploaded.

## Installation

disappR is not on CRAN. Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("EIvimeyCook/disappR")
```

Only `shiny`, `shinydashboard`, `ggplot2` and `lme4` are required. Everything
else is optional and enables specific features:

| Optional package | Enables |
| :---------------- | :------ |
| `glmmTMB` | Count, proportion and binary traits |
| `lmerTest` | p-values for Gaussian models |
| `nlme` | The non-linear exponential ageing function |
| `DHARMa` | Simulation-based residual checks |
| `performance`, `see` | Model diagnostics and their plots |
| `codetools` | Exported R code |
| `callr` | Runs a step in a separate R process, so a failure there cannot stop the app |

```r
install.packages(c("glmmTMB", "lmerTest", "nlme", "DHARMa", "performance", "see", "codetools", "callr"))
```

Without them the app still runs, with fewer analyses; the *Start here* page lists
any that are missing.

## Usage

The package exports one function that launches the app:

```r
library(disappR)
run_app()
```

This opens the app in your browser or the RStudio viewer, arguments are passed to
`shiny::runApp()`. The left-hand sidebar takes you through six steps, plus a
teaching-simulation mode:

1. **Data:** load simulated data with a known answer, one of the 13 published
   datasets, or your own CSV. Map the individual, age and trait columns, plus
   lifespan, age at first record, covariates, groups and censoring; the app checks
   that these are read correctly and names any individual whose rows disagree.
2. **Visual diagnosis:** look for the signatures of selective disappearance and
   appearance before fitting anything.
3. **Missingness and proxies:** check how complete sampling is, and which lifespan
   proxy — age at last record or mean age — best matches known lifespan.
4. **Trajectories:** fit each individual's ageing shape, rebuild the population
   curve from the individual fits, and compare ageing functions.
5. **Modelling:** fit and compare the ten models, with the error family, random
   effects and ageing function you choose; check the fits and see the coefficients,
   predictions and diagnostics.
6. **Summary and report:** gather every result you saved into a cautious automatic
   summary, and export it as an HTML or text report together with the data and R
   code.

### Diagnosing bias before modelling

| Source of bias | Level | What the app shows | Step |
| :-------------- | :---- | :------------------ | :--- |
| Selective disappearance and appearance | Population | Trait trajectories of lifespan or age-at-first-record groups and their differences across age; the trait against lifespan within age bins; the trait before death; selection differentials by age; the disappearance hazard | 2 · Visual diagnosis |
| Incomplete sampling | Data | The sampling grid; missingness by age and against the trait and other variables; agreement between lifespan proxies | 3 · Missingness and proxies |
| A misspecified ageing shape | Individual and population | A function fitted to each individual; the population curve rebuilt from the individual fits; ageing functions compared | 4 · Trajectories |

### The ten models

Fitted with `Rscript`-free, one-click comparison in the app, from a negative
control (Model 1, no lifespan or entry term) through models with an additive or
interactive age-at-last-record term (age-independent and age-dependent selective
disappearance), to models that add an age-at-first-record term as well (selective
appearance), in every combination of additive and interactive. Every model shares
common complete-case rows, the same error family, and the same random-effect
structure, so their AICs are directly comparable.

### Scripted use

The same analyses run from R without the app:

```r
library(disappR)
x   <- disappr_prepare(my_data, disappr_mapping(id = "ID", age = "age", trait = "mass", life = "lifespan"))
fit <- disappr_fit(x, models = c("M1", "M2", "M4"), age_function = "Quadratic")
disappr_compare(fit)
```

`disappr_mapping()`, `disappr_prepare()`, `disappr_fit()`, `disappr_compare()`,
`disappr_predict()`, `disappr_interpret()`, `disappr_diagnose()`, `disappr_code()`,
`disappr_provenance()`, `disappr_read()`, `disappr_example()` and
`disappr_simulate()` are all exported; `vignette("disappR", package = "disappR")`
lists what each one does.

## Getting help inside the app

Every panel has an ⓘ button explaining what it shows and how to read it — click it
for the specifics of any figure, table or setting. The vignettes cover the broader
picture:

```r
browseVignettes("disappR")
```

- **disappR** — an overview of the package and its workflow.
- **workflow** — a walkthrough of the six steps, tab by tab.
- **models-and-methods** — the ten models, their fixed effects, and what each
  comparison tests.
- **validation-and-replication** — how the app's output was checked against
  independent implementations and the published analyses it reproduces.
- **development-and-testing** — the test suite and how to run it.

## Example data

13 published datasets are bundled, each pre-mapped to the model reported in its
paper, spanning laboratory systems and wild populations:

- Sanghvi et al. 2025, *American Naturalist* — *Drosophila melanogaster* daily fecundity
- Sanghvi et al. 2022, *Evolution* — seed beetle female fecundity
- Bichet et al. 2022, *Journal of Animal Ecology* — common tern immune parameters
- Bichet et al. 2022, *Ecology and Evolution* — Alpine marmot immune parameters
- Wynn et al. 2025, *Journal of Animal Ecology* — common tern navigational efficiency
- Moullec et al. 2023, *Frontiers in Ecology and Evolution* — Alpine swift reproduction
- Pásztor et al. 2022, *Ecology and Evolution* — clouded Apollo butterfly body size
- Allain et al. 2023, *Oikos* — eastern chipmunk reproduction
- Bouwhuis et al. 2009, *Proc. R. Soc. B* — great tit recruit production
- Warner et al. 2016, *PNAS* — painted turtle reproduction
- McKenna-Ell et al. 2023, *Biology Letters* — Soay sheep breeding probability and offspring survival
- McKenna-Ell et al. 2023, *Biology Letters* — Soay sheep offspring birth weight
- Szejner-Sigal et al. 2025, *Proc. R. Soc. B* — alfalfa leafcutting bee locomotor activity

Sources and licences are listed in
[`inst/app/data/PROVENANCE.md`](inst/app/data/PROVENANCE.md).

## Data format

One row per individual per occasion. Required columns are the individual's ID, its
age and the trait. Optional columns — known lifespan, age at first record,
covariates, a grouping variable and censoring — unlock more models and figures but
are not required to get started; age at last record and age at first record are
computed automatically when not mapped.

## Interpretation caveats

Associations between missingness and observed variables cannot prove data are
missing completely at random, at random, or not at random. AIC measures relative
support, not causation. The app is a first port of call, not a substitute for
checking your own model's assumptions and residuals — the vignettes and the ⓘ
panels say where each of those caveats applies.

## Bug reports and contributions

Please file issues and feature requests at
<https://github.com/EIvimeyCook/disappR/issues>. Pull requests are welcome.

## Related tools

- [**DCQC**](https://github.com/SORTEE/DCQC) — data and code quality control for
  ecology and evolutionary biology
- [**READMEBuilder**](https://github.com/EIvimeyCook/READMEBuilder) — document a
  finished analysis project for archiving

## Citation

If disappR helps with your work, please cite it:

> Sanghvi, K., & Ivimey-Cook, E. R. (2026). *disappR: a shiny app to model ageing
> and selective [dis]appearance.* R package.
> <https://github.com/EIvimeyCook/disappR>

and the methods paper it implements:

> Sanghvi, K., Ivimey-Cook, E. R., Bouwhuis, S., Sepil, I., & van de Pol, M. (2026).
> A comparison of methods to assess selective disappearance and quantify ageing.
> *EcoEvoRxiv*.

`citation("disappR")` gives both references. A machine-readable
[`CITATION.cff`](CITATION.cff) is included, so GitHub's "Cite this repository"
button gives formatted APA and BibTeX.

## Contact

E. R. Ivimey-Cook — <e.ivimeycook@gmail.com>

## License

Released under the [MIT License](LICENSE.md).
