<p align="center">
  <img src="inst/app/www/logo.png" width="220"/>
</p>

<div align="center">
 <h1>disappR</h1>
 <p><b>Diagnose and model selective disappearance and appearance in longitudinal ageing data</b></p>
</div>

<!-- badges: start -->
<p align="center">
  <a href="https://doi.org/10.17605/OSF.IO/KEVNM"><img src="https://img.shields.io/badge/OSF-10.17605%2FOSF.IO%2FKEVNM-blue" alt="OSF DOI"></a>
  <a href="LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green" alt="License: MIT"></a>
  <img src="https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg" alt="R >= 4.1">
</p>
<!-- badges: end -->

disappR is an R package providing a Shiny app for **diagnosing and modelling
selective disappearance and selective appearance** in longitudinal ageing data.
It follows a visual-first workflow: look for the signatures of selection in the
data, check sampling and the shape of ageing, then choose the error family,
random effects and ageing function and compare a set of mixed models.

Longitudinal ageing trajectories are biased when individuals that die (or enter
the study) early differ from the others in their trait values, or in how the
trait changes with age. Population-level curves then mix within-individual
ageing with changes in who is still being sampled. disappR makes these biases
visible before any model is fitted, and then fits and compares the models that
account for them, with simulations where the right answer is known and
published datasets to practise on.

## Video walkthrough

Two short screen recordings show the app in use. They are also available from
the **Help and videos** tab inside the app.

<p align="center">
  <a href="https://osf.io/u6mxk/"><img src="https://img.shields.io/badge/%E2%96%B6%20Part%201-How%20to%20use%20disappR%20(4%3A46)-A8644E?style=for-the-badge" alt="Watch part 1"></a>
  &nbsp;
  <a href="https://osf.io/5krbc/"><img src="https://img.shields.io/badge/%E2%96%B6%20Part%202-How%20to%20use%20disappR%20(5%3A02)-A8644E?style=for-the-badge" alt="Watch part 2"></a>
</p>

## Features

- **Visual diagnosis first.** Trait trajectories within lifespan (ALR, mean age,
  LS) and AFR bins, differences between bins across age, the trait against
  lifespan within age bins, the trait before death, selection differentials by
  age and the disappearance hazard.
- **Sampling and missingness checks.** The expected-occasion grid, missingness
  by age and by variable, and agreement between ALR, mean age and lifespan, to
  guide the choice of lifespan proxy.
- **Individual and population trajectories.** Parametric fits (linear to
  exponential) for each individual, population curves rebuilt from individual
  fits, and a population-level comparison of ageing functions.
- **Ten comparative mixed models.** Gaussian, Poisson, negative binomial,
  zero-inflated (nbinom1 and nbinom2) and binomial families, random slopes,
  fixed covariates and their interactions, a term builder for extra terms, and
  nested (up to three levels) or crossed random effects.
- **Teaching simulations with known answers.** Choose the ageing form, the type
  of selection and the sampling scheme, and see how each model scores against the
  true trajectory.
- **Published data to practise on.** Ten empirical datasets from laboratory
  and wild populations, each pre-set to the analysis reported in its paper.
- **Reproducible output.** Residual diagnostics, downloadable R code for every
  section, and HTML or text reports of saved results.

## Installation

disappR is not on CRAN. Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("EIvimeyCook/disappR")
```

Count families, p-values and residual checks use optional packages:

```r
install.packages(c("glmmTMB", "lmerTest", "DHARMa"))
```

## Usage

The package exports one function:

```r
library(disappR)
run_app()
```

Arguments are passed to `shiny::runApp()`, e.g. `run_app(port = 4000)`.

## Workflow

| Tab | What it does |
| :-- | :----------- |
| **1 · Data** | Simulated, published or uploaded data; column mapping, covariates and interactions, nested grouping and extra random intercepts, censoring, AFR, subsetting and **data-integrity checks** |
| **2 · Visual diagnosis** | Trait trajectories within ALR / mean age / LS / AFR bins, differences between bins, the trait against lifespan within age bins, and disappearance diagnostics |
| **3 · Missingness and proxies** | Expected-occasion grid, missingness by age and by variable, drivers of missingness, and ALR vs mean age as lifespan proxies |
| **4 · Individual and population trajectories** | Individual parametric fits with data-support metrics, mean-of-coefficients and mean-of-functions trajectories, and ageing-function comparison |
| **5 · Modelling** | Models 1–10 with Gaussian and count families, random slopes, predictions, coefficients in original units, interpretation, DHARMa diagnostics and downloadable R code |
| **6 · Summary and report** | Results saved from every tab, an automatic cautious overview, and HTML and text reports |
| **Help and videos** | The video walkthroughs and the citation |

Constant differences between lifespan groups across age suggest
age-independent selective disappearance; differences that change with age
suggest age-dependent selective disappearance.

## Models

With `b(age)` the chosen ageing basis (e.g. age + age²), `mean(·)` individual
means and `Δ` within-individual deviations:

| Model | Fixed effects | Tests |
| :---- | :------------ | :---- |
| 1 | b(age) | negative control |
| 2 | b(age) + ALR | age-independent selective disappearance |
| 3 | mean(age) + Δb(age) | within/among separation (Fay et al. 2022 for quadratics) |
| 4 | b(age) × ALR | age-dependent selective disappearance (ALR proxy) |
| 5 | mean(age) × Δb(age) | age-dependent selective disappearance (mean-age proxy) |
| 6 | b(age) × LS | positive control with known lifespan |
| 7 | b(age) + ALR + AFR | age-independent disappearance and appearance |
| 8 | b(age) × ALR + b(age) × AFR | age-dependent disappearance and appearance |
| 9 | b(age) × ALR + AFR | age-dependent disappearance, age-independent appearance |
| 10 | b(age) + ALR + b(age) × AFR | age-independent disappearance, age-dependent appearance |

Mapped covariates enter every model additively. Random effects are `(1 | id)`,
plus `(1 | group)` if a grouping column is mapped; with nesting this equals
`(1 | group/id)`. All models use common complete-case rows and maximum
likelihood. Age and proxies are standardised by default, which is a linear
reparameterisation and leaves AIC unchanged.

## Example data

Ten published datasets are bundled, each pre-mapped to the model reported in
its paper:

- Sanghvi et al. 2025, *American Naturalist*: *Drosophila melanogaster* daily fecundity
- Sanghvi et al. 2022, *Evolution*: seed beetle female fecundity
- Bichet et al. 2022, *Journal of Animal Ecology*: common tern immune parameters
- Wynn et al. 2025, *Journal of Animal Ecology*: common tern navigational efficiency
- Allain et al. 2023, *Oikos*: eastern chipmunk reproduction
- Bouwhuis et al. 2009, *Proc. R. Soc. B*: great tit recruit production
- Warner et al. 2016, *PNAS*: painted turtle reproduction
- McKenna-Ell et al. 2023, *Biology Letters*: Soay sheep breeding probability and offspring survival (binomial)
- McKenna-Ell et al. 2023, *Biology Letters*: Soay sheep offspring birth weight
- Szejner-Sigal et al. 2025, *Proc. R. Soc. B*: alfalfa leafcutting bee locomotor activity

Sources and licences are listed in
[`inst/app/data/PROVENANCE.md`](inst/app/data/PROVENANCE.md).

## Data format

One row per individual × age. Required columns are individual ID, age and trait.
Optional columns are known lifespan (LS); age at last record (ALR) and age at
first record (AFR), both computed automatically if absent; condition;
covariates; and a grouping variable. Count families need non-negative integer
traits.

## Interpretation caveats

Associations between missingness and observed variables cannot prove MCAR, MAR
or MNAR. AIC measures relative support, not causation. The decomposition is not
recommended for recovering latent trajectories when selective disappearance is
age-dependent, missingness is substantial, or few individuals are sampled at
consecutive ages.

## Bug reports and contributions

Please file issues and feature requests at
<https://github.com/EIvimeyCook/disappR/issues>. Pull requests are welcome.

## Related tools

- [**shinyDigitise**](https://github.com/EIvimeyCook/shinyDigitise) — extract data
  from published figures
- [**metRscreen**](https://github.com/EIvimeyCook/metRscreen) — title and abstract
  screening for meta-analyses
- [**READMEBuilder**](https://github.com/EIvimeyCook/READMEBuilder) — README
  generation for reproducible research projects

## Citation

If disappR helps with your work, please cite the software and the paper:

> Sanghvi, K., & Ivimey-Cook, E. R. (2026). *disappR: a shiny app to model
> ageing and selective [dis]appearance.* R package.
> <https://github.com/EIvimeyCook/disappR>

> Sanghvi, K., Ivimey-Cook, E. R., Bouwhuis, S., Sepil, I., & van de Pol, M.
> (2026). A comparison of methods to assess selective disappearance and
> quantify ageing. *EcoEvoRxiv*.

Data, code and supplementary material for the paper are on OSF:
<https://doi.org/10.17605/OSF.IO/KEVNM>.

In R, `citation("disappR")` gives the same references. A machine-readable
[`CITATION.cff`](CITATION.cff) is included, so GitHub's "Cite this repository"
button gives formatted APA and BibTeX.

## Contact

Edward R. Ivimey-Cook — <e.ivimeycook@gmail.com> —
[ORCID 0000-0003-4910-0443](https://orcid.org/0000-0003-4910-0443)

## License

Released under the [MIT License](LICENSE.md).
