# disappR 0.9.6

* Start page: a short Purpose section above the workflow sets out why the app exists: selective disappearance and
  appearance are pervasive, the methods that handle them are recent and were built for simple cases, some are biased
  under missing data or under age-dependent selection (which is rarely tested), and there is little guidance on which
  model to fit.
* Decomposition on any sampling schedule. It always uses two successive sampling occasions and only the individuals
  sampled at both (survivor-restricted). On a common schedule the result is unchanged (identical to 0.9.5 for every
  bundled dataset except the chipmunk data, whose ages are irregular). On an irregular schedule the records are now
  placed on a common grid of occasions, spaced by the median interval between an individual's successive records,
  before the chain is built, and a caution below the population-trajectory figure says the result is approximate.
  Schedules such as "age 1, then every third age" are now decomposed through every occasion instead of stopping after
  the first two.
* The exported model script uses the app's own decomposition functions, so its figure matches the app.

# disappR 0.9.5

* Three new empirical examples, each pre-set to the published model: Soay sheep breeding probability and offspring
  first-winter survival (binomial) and Soay sheep offspring birth weight (Gaussian) from McKenna-Ell et al. 2023,
  Biology Letters, with age, age at last observation, bred as a yearling and early-life recruitment (each interacting
  with age; lamb capture age, sex and twin status for weight) and random intercepts for female, year and cohort; and
  leafcutting-bee locomotor activity from Szejner-Sigal et al. 2025, Proc. R. Soc. B, with a quadratic age function,
  lifespan added to Model 1 and a random intercept per bee, opening on females as in the paper.
* Examples can open with a data subset and with extra model terms; both are reset when another example is chosen.
* An example now keeps the ageing function of its published model when the Modelling tab opens; before, the tab
  replaced it with the function that fitted individuals best (the tab's own default for simulated and uploaded data).
* Decomposition: consecutive records are linked when their interval rounds to one sampling step, instead of only
  when it equals the step exactly. Schedules with an unequal first interval (day 1, then weekly on days 7, 14, 21, as
  in the bee data) are now decomposed through every occasion rather than stopping after the first two ages. Regular
  schedules give exactly the same result as before, and a missed occasion still breaks the chain. The exported
  figure code uses the same rule.
* Sampling step: a smallest common interval that is close to but not a divisor of a clearly dominant interval (6 and
  7 days in the bee data) is treated as an unequal first interval, and the dominant interval becomes the step, so the
  sampling grid follows the weekly schedule. No other bundled or test dataset changes step.
* The CONTINUOUS-box note about few distinct whole-number values no longer fires for counts from zero (0, 1, 2, 3).

# disappR 0.9.4

* Binomial models now fit binary (0/1) traits. Without a mapped trials column, the internal trials column (all
  missing) was still passed to glmmTMB as prior weights, so every row was dropped and the fit stopped with
  "contrasts can be applied only to factors with 2 or more levels". Weights are now used only when every analysed
  row has a positive number of trials. The exported model script wrote the negative binomial family for binomial
  fits and omitted the weights; it now writes binomial() or betabinomial() with weights = .trials.
* The number of binomial trials ("Weights") is chosen on the Modelling tab below the error family and is shown only
  for the binomial and beta-binomial families; a note pops up when one of them is selected. Changing it no longer
  resets the family and the model selection.
* New error family: zero-inflated negative binomial with linear variance (nbinom1), next to the nbinom2 version, in
  the models, the error-family check, the model equations and the exported code.
* Three-level nesting: an optional "next level up" grouping column, e.g. individuals within fathers within families,
  (1 | family) + (1 | family:father) + (1 | family:father:ID), in every model (including the non-linear exponential
  models), the equations, the variance table, the figures' panels and the exported code, with an integrity row for
  group labels reused across top-level groups.
* Term builder: the wrench on each model row adds a term built from up to three chosen terms (age, ALR, AFR, LS,
  mean age or any covariate) joined by + or x, including three-way interactions; interactions with age carry an
  over-fitting caution.
* Simulator: each ageing form offers biologically motivated shapes (slow, fast or no senescence; early or late peak;
  U-shape; rapid or gradual decline or growth that levels off; late-life increase). The default shape of each form
  is unchanged, so earlier simulations are reproduced exactly.
* Data tab: CONTINUOUS and CATEGORICAL covariate boxes that warn immediately when a column looks as if it belongs in
  the other box; a note and a check that age must be numeric; the upload note states that blank cells and NA are
  both read as missing.
* Figures: 95% confidence bands on the bin-difference trend lines and on the lines of the trait-against-lifespan
  figure; solid black regression lines on the ALR, mean age, LS and AFR agreement plots; a 'Panels by' menu
  (categorical variables) in the trajectory figure; optional jittered individual records in the population-level
  trajectory plot.
* The references to cite are shown in the sidebar.

# disappR 0.9.3

* Age at first trait expression (AFE) now lives only where it applies: the "Count missed occasions from" control
  on the Sampling and missingness tab, with AFR the default and a numeric AFE entry appearing only when AFE is
  chosen. Its "i" panel states that the setting changes the missingness figures on that tab and nothing else.
* The Data tab's AFR control is back to two options, automatic from each individual's first record (recommended)
  or a mapped column, and says plainly that AFR is what every model and statistic uses.
* The 0.9.2 option that set AFR equal to AFE is withdrawn, together with the constant-AFR path it used in
  standardise_data, so no entered age can reach the models. AFR, the ALR and mean-age proxies, the appearance
  models and every other statistic are unchanged by this release.

# disappR 0.9.2

* The Data tab now separates two ages explicitly. "Age at first trait expression (AFE)" is an optional single
  age entered as a number: the age from which the trait exists and could have been measured, even if recording
  began later. "Age at first record (AFR)" is when the individual enters the data and offers three options -
  automatic from the first record, the same as AFE (which requires AFE to be filled and makes AFR constant, so
  Models 7-10 become unavailable), or a column, whose dropdown appears only when that option is chosen. The
  column-based expression input added in 0.9.0 is withdrawn.
* The AFR-ALR agreement is now shown as a plot (with the 1:1 line and a fitted line) beside its table.
* The start page and the data-source panel no longer refer to "the fruit-fly data"; they refer to the published
  empirical datasets. The count-as-Gaussian note beneath the model table is shortened to the essential check.
* The predictions panel now explains why every model can return almost the same trajectory while their AIC values
  differ by hundreds of units: under age-independent selective disappearance the proxy is held at its mean, so
  proxy terms contribute nothing at that value and the age slope is unbiased in all models; the curves separate
  only when disappearance is age-dependent. Verified by simulation.

# disappR 0.9.1

* Fixes a parse error introduced in 0.9.0: the rewritten empirical-example notes were missing the comma between
  each example's own text and the shared closing paragraph inside paste(), so sourcing global.R failed. All seven
  notes are corrected.

# disappR 0.9.0

* Age at first observation is now separated from age at first trait expression. The Data tab maps entry into the
  dataset (AFR, used in Models 7-10) and, optionally, an earlier age at which the trait could first have been
  recorded. The distinction changes the missingness window: with entry only it runs from AFR to ALR, so earlier
  occasions are "not yet expressing"; with an expression age it opens there, so occasions between expression and
  entry count as missed. A new "i" panel explains the two and warns that the resulting missingness estimates are
  not comparable.
* New diagnostic on the missingness tab: the correlation between AFR and ALR, with the mean and SD of the
  observation window, and an "i" panel on how to read it (design-driven versus biological association, and what
  it implies for ALR and mean age as proxies).
* New "i" panel on the among-individual polynomial order, with a worked example (a quadratic ageing function
  gives ALR and ALR^2 plus their interactions with each age term) and the case for using it: a linear ALR term
  cannot represent a non-monotonic relation between lifespan and the trait, and will report no selective
  disappearance when the pattern is U-shaped or humped.
* Individual-fit panels: the function comparison and the mean-coefficient and mean-function reconstructions now
  state that they are unreliable under high missingness, few time steps, few records per individual or small
  samples, and that the dAICc ranking is suggestive rather than prescriptive. A warning now appears beside the
  random-slope advice whenever the data give weak support for individual slopes.
* The integer-data caution is shortened to the essential check: integer data treated as Gaussian, confirm the
  trait is continuous rather than a count or 0/1 outcome.
* The great tit example opens on laying date (LD), with the paper's other traits available from the trait
  dropdown. The seed-beetle male weight example has been removed; the female fecundity example remains.

# disappR 0.8.8

* The empirical-example notes are now qualitative. Each states what the paper did, what kind of model it fitted
  and what it generally found, then says that the app's defaults are set close to the published fit for that
  trait, can be changed, and generally reproduce the published pattern but may differ because of subsetting,
  model terms or covariates that were not fully specified, random-effect structures, or software and estimation.
  Specific published coefficients, AIC gaps and provenance wording have been removed from the notes; the file
  histories remain in data/PROVENANCE.md.
* The start page and the "Data source and simulation" panel now say that several empirical examples are
  available, from laboratory systems (fruit flies, seed beetles) and wild populations (terns, turtles, great
  tits, chipmunks), and the panel's cautions carry the same reproducibility caveat.
* New cautions on the individual-fit panels: the mean-coefficient and mean-function reconstructions can depart
  from the true latent trajectory for reasons other than selective disappearance - extrapolation beyond each
  individual's own observed ages, near-saturated fits when records are few, and, for counts, averaging on the
  response scale while the truth is a curve on the link scale - and these grow under missingness. Users are
  asked to read them alongside the population-level model predictions on the Modelling tab. The function
  comparison panel now notes that with few occasions per individual the simplest function usually wins on AICc
  (a quadratic needs at least six records for AICc to be defined), which reflects what is estimable rather than
  what generated the data.

# disappR 0.8.7

* The great tit example now ships the paper's full sample, 7,341 breeding attempts by 4,935 females, in place of
  the 7,126-row subset that omitted 215 attempts without a lay date. The replication improves: the
  individual-level peak is 2.80 years against the published 2.80 (was 2.85), the cross-sectional peak 3.48
  against 3.45, ALR +0.062 against +0.04, the post-peak slope -0.105 against -0.11, and clutch, brood and
  fledgling coefficients +0.499/-0.058, +0.594/-0.096 and +0.637/-0.109 against +0.52/-0.06, +0.62/-0.10 and
  +0.63/-0.11. The note records the counts, the lay-date missingness and the 361 females whose supplied ALR
  differs from their last recorded age.
* Built on 0.8.5, whose REPLICATION.md refits each example under its own default settings including the extra
  random intercepts. Carried over from the 0.8.6 branch: the repaired turtle sentence and a pointer from the
  empirical-example panel to the comparison table.

# disappR 0.8.5

* New REPLICATION.md compares, example by example, what each paper reported with what the app returns on the
  bundled data under the pre-filled settings, and says plainly which comparisons are quantitative, which are
  qualitative and which could not be verified offline (the fly and beetle-fecundity count families).
* Two figures quoted earlier were corrected by refitting with the examples' own default settings: the turtle
  hatchlings-alive comparison is weak (Model 4 leads Model 1 by 2.3 AIC, not 37.8, once the year random intercept
  is included) and the great tit recruit ranking under Gaussian settings gives Model 1 +6.1 rather than +13.
* Start-up checks: all three files balance, no empty arguments, every symbol used in ui.R is defined before the
  server starts, all eight example files exist with resolvable mappings, valid families, model ids and ageing
  functions, no duplicate or orphan input/output ids, and every info panel key referenced in the UI is defined.

# disappR 0.8.4

* Restores three details to the example notes: the tern stopover-threshold sensitivity (the age signature holds
  at 50 and 200 km but disappears with no filter, so the conclusion is conditional on it) and why 24 birds are
  lost; the evidence for harmonising the tern sex codes (nineteen individuals appear under both codings, and 1
  never coincides with f nor 2 with m); and the names of the restored columns, Wt_D0 in the beetle files and the
  nest microhabitat variables in the turtle file.

# disappR 0.8.3

* Merges the regenerated datasets and provenance record with the corrections from 0.8.1/0.8.2 that the previous
  build did not carry: the painted-turtle years and usable-record counts, the chipmunk sexes and the Model 6
  row restriction, the great tit sample description, and the app test's data source.
* The tern note now gives the published track and bird totals against the bundled ones, and states that the paper
  applied FLightR particle filtering before computing deflection, so the derived series reproduces the age
  signature but not every coefficient.
* The empirical-example panel points to data/PROVENANCE.md.

# disappR 0.8.2

* Data provenance audit of the bundled examples. Every file is now a verbatim copy of its archived original with
  all columns retained (earlier builds dropped a few unused columns from the beetle and turtle files, including
  the beetles emergence weight Wt_D0 and the turtle nest microhabitat variables, and renamed the turtle columns).
  The turtle example mapping uses the original column names again.
* The tern navigation example, the only derived dataset, now defines the autumn goal as the individuals mean
  wintering position rather than the tracks, which is what the paper describes: the among- and within-individual
  age effects move from -0.64 / +1.58 to -0.77 / +0.90 against the published -0.870 / +0.851. The file now also
  carries move_bearing, goal_bearing, goal_lat, goal_lon, step_km, the migration dates and sex_raw so that the
  derivation can be audited, and data/PROVENANCE.md documents every transformation.

# disappR 0.8.1

* Eighth empirical example: Bouwhuis et al. 2009, Proc. R. Soc. B - great tit (recruit production), the
  49-year Wytham dataset. Age at last reproduction is supplied in the data, so the example also demonstrates a
  mapped ALR that deliberately disagrees with the last recorded age.
* The example reproduces the published analysis closely: cross-sectional peak 3.43 years against 3.45, moving
  to 2.85 against 2.80 once ALR is fitted; ALR +0.06 against +0.04; post-peak slope -0.110 against -0.11;
  clutch size, brood size and fledgling numbers +0.48/-0.056, +0.60/-0.097 and +0.65/-0.112 against
  +0.52/-0.06, +0.62/-0.10 and +0.63/-0.11, with no ALR effect on any of the three.

# disappR 0.8.0

## Empirical examples
* The "Fruit-fly fecundity" data source becomes **Empirical examples (published data)**, a drop-down of seven
  published datasets listed as Surname et al., year, journal, species:
  - Sanghvi et al. 2025, American Naturalist - Drosophila melanogaster (daily fecundity)
  - Bichet et al. 2022, Journal of Animal Ecology - common tern (immune parameters)
  - Wynn et al. 2025, Journal of Animal Ecology - common tern (navigational efficiency)
  - Sanghvi et al. 2022, Evolution - seed beetle (male weight)
  - Sanghvi et al. 2022, Evolution - seed beetle (female fecundity)
  - Allain et al. 2023, Oikos - eastern chipmunk (reproduction)
  - Warner et al. 2016, PNAS - painted turtle (reproduction)
* Each example ships with the published data in inst/app/data, a pre-filled column mapping (covariates,
  interactions, nesting, extra random intercepts, censoring) and pre-filled model settings - error family,
  ageing function and the selected models - chosen to reproduce the published analysis. A note under the
  drop-down states what the paper found and what the app returns.
* Derived variables are included in the bundled files so that nothing has to be computed by hand: age in the
  Bichet data (year of sampling minus year of birth) and instantaneous deflection in the Wynn data (the angle
  between the movement bearing and the bearing to the goal, computed from the archived geolocator fixes, with
  steps under 100 km removed).
* New helpers: EXAMPLES, EXAMPLE_CHOICES, example_map() and load_example_file(). FLY_MAPPING and
  load_fly_example() are retained.

# disappR 0.7.6

* The two exponential functions are renamed so that they are not confused. The linear-in-coefficients basis
  (intercept + beta * exp(-age)) is now "Asymptotic exponential"; the model with an estimated rate is now
  "Exponential (a * exp(b * age))". Internal keys, the simulator, the exported R code and the tests follow the
  new names.
* Each model's "i" panel now states what the age terms are for the chosen function (for example
  "f1 = exp(-age), centred and scaled"), so models that share the formula trait ~ f1 + ... but differ in the
  age term are no longer ambiguous.
* Tests assert that the two exponential functions give different equations, coefficient names and fitting calls
  for every model where both are defined, that Models 3 and 5 remain unavailable for the estimated-rate model,
  and that each ageing function builds a different age term.

# disappR 0.7.5

* Caveats only, no change to any calculation. Added to the "i" panels: terminal effects as a rival explanation
  for the model ranking (trait before death); clustering of repeated records in the missingness and proxy-trend
  tests; hazards not accounting for frailty; saturated individual fits (quadratic on 3 records, cubic on 4);
  the two exponential parameterisations not being comparable; standardised coefficients depending on the sample's
  age range; likelihood-ratio p-values being approximate; model ranking as relative support; and multiple testing
  with no correction.

# disappR 0.7.4

* Trait-dependent missingness is now flagged as a rival explanation for the model ranking: the caveat in the
  "Missingness against other variables" help entry says that it produces the same evidence as age-dependent
  selective disappearance, and a warning appears under that plot whenever the prior-trait (or condition)
  association is significant.

# disappR 0.7.3

* AFR with fewer than two values is now a Warning in the integrity table (it names Models 7-10 and
  explains that they would repeat Models 1-4 because their AFR terms are dropped as rank deficient),
  and the same warning appears with the model results.
* Removed the calendar-year age warning added in 0.7.1.
* Binomial families: rows without a positive number of trials are dropped from the analysed rows and
  reported, instead of reaching glmmTMB as zero or negative weights.

# disappR 0.7.2

## Binomial families
* New error families: binomial (binary 0/1, or a proportion of successes) and beta-binomial for
  overdispersed proportions, both through glmmTMB with a logit link.
* New optional "Number of binomial trials" column on the Data tab. For proportion traits the trials
  are passed to glmmTMB as prior weights, reported in the integrity table, and rows without a
  positive number of trials are dropped from binomial models.
* suggest_family() now recognises binary traits (suggests binomial) and proportions between 0 and 1
  (suggests binomial and asks for a trials column); the Data tab reports the trait kind.
* Model equations, the exported R code, the error-family check and the performance checks all follow
  the binomial families (zero-inflation checks are limited to the count families).
* Guidance added: fitting binary or proportion traits as Gaussian can favour the interaction models
  even with no selective disappearance, because their variance changes with the predicted mean.

# disappR 0.7.1

## Fixes from the empirical stress test
* Fixed a stray empty argument in ui.R (`),,` before the sampling section code box) that stopped the
  app from starting.
* Binary (0/1) traits are recognised as binary rather than counts: the Data tab says so, Gaussian is
  suggested as a linear-probability approximation, and the text points to binomial() in the exported code.
* Underdispersed counts (clutch and litter sizes, breeding success) no longer get a Poisson suggestion:
  below a Pearson dispersion of 0.7 the suggestion is Gaussian, with the reason given.
* New integrity warning when ages look like calendar years (all ages between 1800 and 2100), the most
  common mapping mistake in datasets with no age column.
* The zero ratio in the trait-distribution line is capped (">999x") instead of printing huge numbers.
* The duplicate ID x age advice now covers IDs reused in different groups: map the grouping column and
  keep nesting instead of averaging.
* When a count family is refused because the trait is not integer, the message says so if averaging
  duplicates caused it.

## Guidance added
* Trajectories within bins: non-monotonic bin means mean the lifespan effect is not linear, so ALR^2 is
  needed (Models 2 and 4 use a linear ALR term and will find nothing).
* Model support: trait-dependent missingness and terminal declines both reproduce the ranking expected
  from selective disappearance, with pointers to the missingness drivers and the trait-before-death figure.
* Trait before death: a terminal change alone can make Models 4 and 5 win with no ageing at all.
* Sampling summary: individuals on genuinely different schedules appear as heavy missingness, which is
  study design rather than failed detection.

# disappR 0.7.0

## Stability, exact model equations and section scripts (latest round)
* Fixed the start-up crash "attempt to set an attribute on NULL": standardise_data() called setNames() on the
  (NULL) names of an empty covariate list whenever no covariates were mapped.
* Uploads: any error while reading a file, or a malformed result, is reported as a message before the data are used.
* Every observer runs inside an error guard, so an unexpected error is shown as a notification (and in the R
  console) instead of disconnecting the session.
* Each model's "i" shows the exact model: the fixed effects fully expanded as R expands them (A * B = A + B + A:B,
  one coefficient per term and per level of a categorical covariate), random effects and their distributions,
  the error family and zero-inflation part, and the R syntax as written, fully expanded and as fitted.
* Population-level trajectories: unticking every model no longer redraws all models (which made ticking one look
  like the others disappearing).
* Data, Visual diagnosis, Missingness and proxies, and Individual and population trajectories each end with a
  collapsible "R code for this section" box (copy or download) that reproduces the section with the current
  settings, including the app's own helper functions.

## Final polish
* Individual fits need fewer records: 2 for a linear fit, 3 for quadratic, logarithmic and exponential fits, 4 for a
  cubic fit (fits through every point carry no AICc and are left out of the function comparison).
* Workflow codes (A1-A7, B1-B5) removed from all titles, labels and help text.
* Removed the "Lifespan beyond the last record" plot.
* performance checks can no longer stop the app: each check and its formatting runs with its own error handler
  and time limit, glmmTMB models use lighter check_model() panels, and failures are reported as messages.
* The sampling grid shows up to 500 individuals (slider), with the plot height adjusted to the number shown.

## Layout and model specification
* Disappearance diagnostics merged into Visual diagnosis; one shared settings row for A1-A2 and figure-specific
  settings (bin boundaries and SE, trend lines, jittered points, survivors comparison) directly above each figure.
* A1 bin differences always compare consecutive bins; every bin x age point and A2 age bin needs at least 3
  individuals (the minimum-individuals setting was removed).
* Removed: the trait-now disappearance model plot, the A5 relative-to-age-mean option, the A7 hazard table and the
  within-model ageing-function comparison.
* Data: continuous and categorical fixed-effect covariates are mapped separately; interactions between a covariate
  and age or between two covariates are added to every model.
* Modelling: full-width settings with one row per model (tick box, "i" with the model's meaning and exact
  specification under the current settings, and a menu of terms added to that model only); model support below.
* Individual fits: "Use this ageing function for the mixed models" sets the Modelling tab's function; the
  population-level function comparison includes the non-linear exponential (Gaussian traits).
* Observed means in the population trajectory plots are drawn with a dark outline.

## Navigation and export
* Browser tab title; tabs named for what they do (Data, Visual diagnosis, Disappearance diagnostics, Missingness and
  proxies, Individual and population trajectories, Modelling, Summary and report); Start here describes each tab.
* Visual diagnosis uses one grouping variable (ALR, mean age, LS or AFR) and one number of bins for both figures.
* The population-level ageing-function comparison moved into 'Individual and population trajectories'.
* The modelling ageing function defaults to the function with the lowest mean ΔAICc across individual fits (with a
  note that simpler functions may be preferable).
* Data tab: distribution and preview before the integrity checks.
* Exported R code now also reproduces the subset, censoring and the trajectory figure (model predictions,
  observed means and the decomposition); fixed the non-linear script builder, which kept stale model lines.

## Modelling, diagnostics and interpretation
* New tab A4-A7 (disappearance diagnostics): A4 discrete-time hazard of disappearing before the next occasion
  from the current trait (optionally its change and a trait x age term); A5 trait before death by lifespan
  tercile (optionally relative to the mean at the same age); A6 selection differentials by age; A7 life table
  with the age-independent and age-specific hazard.
* Models 9 (ALR x age + AFR) and 10 (ALR + AFR x age), so appearance and disappearance can differ in age
  dependence; extra terms (ALR, AFR, LS, mean age, additively or x age, and covariates) can be added to chosen
  models; likelihood-ratio tests only for genuinely nested pairs.
* Random effects: intercept, uncorrelated or correlated intercept and slope, or automatic choice from the data
  support; polynomial among-individual terms (optional); covariate x age interactions chosen on the Data tab;
  zero inflation depending on age, ALR, AFR and covariates.
* Non-linear exponential ageing a * exp(b * age) fitted with nlme (Gaussian traits); the linear-in-parameters
  exponential basis is now exp(-z), and the simulated exponential form declines steeply early and slowly late.
* Fits classified Valid / Caution / Failed; Failed fits (e.g. non-positive-definite Hessian) are excluded from
  ranking, tests, prediction plots and automated interpretation unless included explicitly.
* Coefficients with scaling constants and estimates per original unit, random-effect variances, and a cautious
  plain-language interpretation of selection terms; performance::check_model() and related checks (fail-safe).
* Predictions: choose which models to draw, overlay the A3 reconstruction, and show curves by factor level.
* Internal consistency checks (AIC vs recovery of the simulated truth, AIC vs likelihood-ratio tests, best-fit
  validity, random-slope support, A2 trend) and a more cautious automatic summary.
* Data: optional subset by a categorical variable; AFR defined as age at first observation; lifespan margin
  removed; clearer A1-A2 colours; A1 difference trend as one pooled line or one per bin pair; A2 points optional;
  a prominent caveat when sampling is irregular; A3 data-support metrics; uploads up to 200 MB.

## Robustness (stress-tested on 22 out-of-sample datasets)
* Sampling step inferred from the most common interval between an individual's consecutive records, so
  staggered cohorts, floating-point ages and occasional extra records no longer shrink it; irregular ages use
  the median interval and are flagged. Expected occasions are anchored on each individual's first record.
  The grid is capped at 300,000 cells (coarsening reported) and the heatmap at 60,000 tiles. Results on
  regular schedules (including the fruit-fly data) are unchanged.
* Decomposition built from same-individual records one sampling step apart, following one schedule.
* Optional "Round ages to multiples of" for irregular ages (also applied in the exported R script).
* Uploads: automatic separator (comma, semicolon, tab, pipe) and decimal-comma detection, byte-order mark,
  common NA codes, trimmed white space, IDs kept as text ("007" != "7"), repaired blank or duplicated
  column names, latin1 fallback, clear messages for Excel workbooks, header-only files and dates in the age column.
* Integrity checks for irregular or offset schedules, -99/-999/-9999 codes, ages <= 0, extreme values,
  covariates without variation, with many levels or missing values, and empty columns.
* Models: covariates without variation are omitted with a note instead of stopping every model; ageing
  functions that need more distinct ages than the data have are refused with a message; covariate and
  random-effect names that collide after make.names are kept apart; glmmTMB adjusts rank-deficient designs;
  plain-language hints for unidentifiable random slopes. Individual fits use lm.fit, and the A3 function
  comparison uses at most 2,000 individuals.
* New `tests/stress_test.R` with the datasets in `tests/stress_data`.

## Latest revisions
* Simulator: five ageing forms (linear, quadratic, cubic, logarithmic, exponential, using the same bases as
  the models). Selective disappearance and selective appearance are each set as none / age-independent /
  age-dependent / both, with a positive or negative direction (replacing the intercept/slope/shape selectors).
  Individual differences in ageing rates (small/large) and mean lifespan (3-30 occasions); diet covariate that
  lowers the trait without affecting lifespan; AFR and ALR uncorrelated by construction.
* A1-A2: facet (panel) by any covariate or treatment with 2-8 levels; number of bins from 3 to the number of
  distinct ages (A1: values of the grouping variable) observed in at least the minimum number of individuals,
  counting single-record individuals; one bin per value when bins >= distinct values. A1 bin differences use
  age on the x-axis with a least-squares trend line per pair. A2 adds a table of the regression coefficient in
  each consecutive age bin (95% CI, r, change from the previous bin) and a weighted trend across age.
* A3: mean of coefficients and mean of individual functions (identical for functions linear in their
  parameters; formulas in the help), each with a 95% band; adds a non-linear individual function
  a * exp(b * age); curves restricted to ages observed in >= 10 individuals; minimum residual df fixed at 1.
* Sampling grid: observed / missed / death or ALR / not expected, ordered by ALR, AFR, mean trait, ID or at random.
* Models: random-effect variances, SDs, correlations and dispersion shown with the coefficients (and for all
  models); B2 function comparison with selectable functions and distinct colours; new within-model comparison
  of ageing functions for any of Models 1-8.
* "i" help for every section (what it does, how to read it, cautions) and "Save to summary" buttons; the
  Summary tab lists saved results and exports them as an HTML report (with plots) or a text report.

## Earlier revisions
* A2 now plots the trait against LS (or ALR / mean age / AFR) within age bins, with a linear fit per age bin.
* New A1 panel: signed difference in mean trait between proxy bins at each age (successive bins or all pairs).
  The A2 interpretation box was removed; bins are numbered in the A1 legend.
* Teaching data: AFR is fixed (everyone starts at age 1) or variable. Variable AFR is either random or
  selective appearance (positive or negative). MCAR samples ~75% of occasions (the paper simulation keeps the
  Methods' 0.5). Trait- and condition-dependent missingness were rescaled to ~25-30% missing.
* Data: age at first trait expression is either each individual's AFR or a common age for everyone, so that
  occasions between that age and the first record count as missing. Adds a censoring indicator (LS set to
  missing for censored individuals; fly data: `Censored = 0`), additional crossed random intercepts `+ (1 | X)`,
  and a distribution plot for any variable.
* Sampling: ALR vs mean age, ALR vs LS and mean age vs LS as separate labelled panels (ALR vs mean age only
  when LS is unknown).
* A3: free x and y scales per individual.

Rebuilt after an out-of-sample check on the fruit-fly fecundity data (Sanghvi et al. 2025, Am. Nat.).

## Why the rebuild: the fly data exposed a modelling problem
* 0.6.x fitted Gaussian LMMs only. On the fly counts (30% zeros, strongly overdispersed) a Gaussian
  Model 1-5 comparison favours **Model 5** (Model 4 is about 19 AIC worse). The manuscript's zero-inflated
  negative binomial favours **Model 4** (dAIC of Models 1, 2, 3, 5 = 38.4, 24.1, 28.0, 2.7). An independent
  re-implementation of the ZINB GLMM reproduced these values (38.3, 24.0, 27.2, 2.6). The error family
  changes the conclusion, so it is now a first-class choice.

## Modelling (B1-B5)
* Error families: Gaussian (lme4), Poisson, negative binomial (nbinom2, nbinom1), zero-inflated Poisson
  and zero-inflated NB (glmmTMB), with constant or age-dependent zero inflation.
* Fixed-effect covariates, and a higher-level random effect that can be nested (individuals identified by
  group + ID, i.e. `(1 | group/ID)`) or crossed.
* Age and proxies standardised before fitting (a linear reparameterisation: likelihood and AIC are
  unchanged, but raw age x ALR interactions no longer produce NaN/non-positive-definite Hessians).
* All models are fitted on common complete-case rows. Model 6 is left unselected by default when LS is
  missing for some individuals, and Models 7-8 when AFR is (nearly) constant, so Models 1-5 are not
  silently truncated.
* Nested likelihood-ratio tests (1 vs 2, 2 vs 4, 3 vs 5, 2 vs 7, 4 vs 8, 7 vs 8), Akaike weights,
  convergence/singularity/boundary notes.
* Population predictions follow the manuscript: random effects excluded, proxies held at individual-level
  means, factor covariates averaged over individuals, response scale.
* Error-family check (AIC across count families), DHARMa residual checks, and a generated R script that
  reproduces the comparison outside the app.
* Fitting runs on a button press, with a staleness guard when data or settings change.

## Diagnostics
* New data-integrity table: duplicate ID x age records, IDs linked to more than one group, mapped ALR
  vs last record, missing LS, LS earlier than the last record, LS - ALR gap, off-grid ages, records per
  individual, AFR variation, and a suggested family (Poisson GLM Pearson dispersion and zero excess).
* A1 bins are formed among individuals (not rows), with equal-width (as in the manuscript) or quantile
  boundaries and a minimum number of individuals per bin x age point.
* A2 (new): proxy-trait slope or correlation at each age with 95% CIs and an inverse-variance weighted
  test for change with age. It replaces the nlme "simple cross-check" models.
* A3: minimum residual df per individual (exact fits no longer report R2 = 1), a warning when fitted
  individuals are longer-lived than the population (survivor bias), and a choice of which individuals
  to draw (instead of "0 = all").
* Sampling: expected windows now end at the last occasion at or before LS (floor, not rounding), plus an
  optional lifespan margin. B3 guidance is based on detection probability and proxy correlations (the
  3/p - 1 threshold was removed). New panels: ALR / mean age vs LS, LS - ALR gap histogram, and
  coverage by age.
* Decomposition: starts at the first two consecutive ages at which individuals were sampled at both.
  The starting value is the mean trait at the first of these ages of those individuals. The 0.6.1
  all-individuals anchor and the anchor selector were removed.

## Teaching simulations
* Traits: continuous body mass (Gaussian), fecundity counts (Poisson, negative binomial, ZINB; log link),
  or the exact manuscript parameterisation.
* "Dramatic" strength (|rho| = 0.85 through a latent-quality factor) makes patterns obvious.
  "Moderate" uses |rho| = 0.5.
* Independent sign choices for lifespan correlations with intercept, slope and shape.
* Sampling options: complete, MCAR, MWO, MWY, trait-dependent and condition-dependent.
* Optional selective appearance, a selective treatment covariate, and 25 families (nested random effect).
* Each simulation states what the user should see. Model predictions and the decomposition are scored
  against the true trajectory (Methods Eq. 11).

## Empirical example
* The fruit-fly dataset is bundled, with mapping, covariates, nesting and family pre-filled.

## Other
* Removed: nlme dependency, the generic example CSV, the proxy-by-age-class plot.
* Tests: `tests/smoke_test.R` (all analysis functions, plus the fly replication checks) and
  `tests/app_test.R` (drives the real server with `shiny::testServer()`).

# disappR 0.6.1

## Fixes that stopped the app running
* `server.R`: `download_example` content handler contained `;` inside an argument list (parse error).
* Package layout: app files moved to `inst/app/` (logo to `inst/app/www/`); `run_app()` calls `shiny::runApp()`.

## Bugs
* Long model formulas were rebuilt via `deparse()`, which splits at ~60 characters (cubic Model 8 always failed).
* AFR variance check could return `NA` inside `if()`.
* Individual trajectories: "0 = all individuals" drew no curves; empty fits crashed the facet plot.
* Individual dAICc taken over different candidate sets per individual.
* Log basis choice fixed once per dataset.
* Missingness "After last observation" relied on an unmapped ALR column.
* Column auto-detection mangled upper-case names.
* Prediction transform parameters and within-individual centring corrected.
* Decomposition anchor changed to all individuals at the anchor age (replaced in 0.7.0).
* Font Awesome 4 icon names replaced.

## Improvements
* Model fitting status table; quantile bins for continuous proxies; corrected the claim that Models 4-6 are
  equivalent under complete sampling (Models 4 and 6 are; Model 5 is not in the same design space).
