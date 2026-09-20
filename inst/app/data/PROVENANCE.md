# Provenance of the bundled empirical examples

Every file below was checked cell by cell against the archived original.

| File | Rows vs original | Columns | Values | Derived columns |
|---|---|---|---|---|
| fly_fecundity.csv | unchanged | unchanged | unchanged | none (as supplied with the app) |
| bichet_2022_tern_immunity.csv | 1097 / 1097 | all kept; `year of sampling` renamed `year_sampling` | identical | `age` = year_sampling - year_of_birth |
| sanghvi_2022_beetle_male_weight.csv | 5354 / 5354 | all kept | identical | none; re-saved as UTF-8 |
| sanghvi_2022_beetle_female_fecundity.csv | 6001 / 6001 | all kept | identical | none; re-saved as UTF-8 |
| allain_2023_chipmunk_reproduction.csv | 487 / 487 | all kept | identical | none |
| warner_2016_turtle_reproduction.csv | 2205 / 2205 | all kept | identical | none |
| bouwhuis_2009_great_tit_recruitment.csv | 7126 / 7126 | all kept | identical (md5 match) | none |
| wynn_2025_tern_navigation.csv | **2689 / 126483** | subset + derived | identical for retained fields | `deflection`, `season`, `step_km`, `move_bearing`, `goal_bearing`, `goal_lat`, `goal_lon`, `sex` |

## The tern navigation file is the only transformed dataset

The response variable analysed by Wynn et al. (2025), instantaneous deflection from the goal, is not in the
archived file, which contains raw geolocator positions. It is computed here as:

1. each track is split into autumn and spring migration using the departure and arrival dates already in the file;
2. the goal is the individual's mean position between arrival at and departure from the wintering area (autumn)
   or the breeding colony at 53.5111 N, 8.1056 E (spring);
3. deflection is the absolute angle between the bearing to the next fix and the bearing to the goal;
4. steps shorter than 100 km are dropped, following the paper's stopover filter.

This keeps 2,689 fixes from 186 tracks and 78 of the 102 birds; tracks without wintering dates cannot be
classified and are lost. The published estimates are an among-individual effect of -0.870 degrees per year and
a within-individual effect of +0.851; this derivation gives -0.769 and +0.903. The sign pattern is unchanged at
50 km (-0.66 / +1.24) and 200 km (-0.55 / +1.58) thresholds, but disappears with no distance filter
(-0.06 / +3.61), so the stopover filter matters, exactly as the paper states.

`sex_raw` preserves the original coding. `sex` harmonises the numeric codes to letters: nineteen individuals
appear under both codings, and code 1 always coincides with m (nine individuals) and 2 with f (nine), never the
reverse, so the mapping is supported by the data rather than assumed.

## Added in 0.9.5

Supplied by Krish from the papers' data archives (Soay sheep: doi:10.5061/dryad.stqjq2c7s; bees:
doi:10.5061/dryad.j9kd51cq3) and bundled byte for byte under new file names.

| File | Rows | Columns | Values | Derived columns |
|---|---|---|---|---|
| mckennaell_2023_soay_breeding_survival.csv (`fecundityoffsurv.csv`) | 3173 records of 762 females (the paper's breeding-probability sample; 2573 non-blank offspring-survival records of 714 females) | unchanged | unchanged | none |
| mckennaell_2023_soay_offspring_weight.csv (`offspringwt.csv`) | 2317 lambs of 649 mothers (the paper's sample) | unchanged | unchanged | none |
| szejnersigal_2025_bee_activity.csv (`Activity.csv`) | 808 records of 184 bees (414 female, 394 male records); the paper's text reports 199 bees in the activity experiment | unchanged | unchanged | none |
