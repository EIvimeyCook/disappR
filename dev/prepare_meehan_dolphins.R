# Builds inst/app/data/meehan_2026_dolphin_group_size.csv from the authors' shared data:
# Meehan et al. 2026, PLOS Biology 10.1371/journal.pbio.3004022; data and code github.com/k-meeha/socialagingindolphins,
# archived at Zenodo 10.5281/zenodo.22675472 (CC-BY 4.0). Run from the package root:
#   Rscript dev/prepare_meehan_dolphins.R path/to/socialaging_main.csv
# Every sighting is kept (one per dolphin per day, as in the paper). The observation date is left out: the paper's
# random effect of date has 3,512 levels, too many for the app; year is kept. Ages are converted from days to years; ALS (age at
# the last sighting, computed by the authors as each dolphin's oldest sighting) is kept as the ALR column.
src <- commandArgs(TRUE)[1]
d <- utils::read.csv(src, stringsAsFactors = FALSE)
stopifnot(all(c("Observation.Date", "Year", "Sex", "ObsAge", "ALS", "Repro.Stat", "Associates", "anon_id") %in% names(d)))
out <- data.frame(dolphin = d$anon_id,
                  age_years = round(d$ObsAge / 365.25, 3),
                  sex = tolower(d$Sex),
                  repro_status = c(Calf = "calf", Newborn = "newborn", `Not Calf` = "independent")[d$Repro.Stat],
                  year = d$Year,
                  age_last_seen = round(d$ALS / 365.25, 3),
                  associates = d$Associates,
                  stringsAsFactors = FALSE)
stopifnot(!anyNA(out), nrow(out) == 64914L, length(unique(out$dolphin)) == 457L)
utils::write.csv(out, file.path("inst", "app", "data", "meehan_2026_dolphin_group_size.csv"), row.names = FALSE)
cat(nrow(out), "sightings of", length(unique(out$dolphin)), "dolphins written\n")
