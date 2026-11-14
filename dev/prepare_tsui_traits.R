# Prepares the other Seychelles warbler traits of Tsui et al. 2026 (Evol. Lett. 10: 506-520) as bundled examples.
# Source: Zenodo 10.5281/zenodo.21101620 (CC BY 4.0). Run with the folder of the unzipped deposit:
#   Rscript dev/prepare_tsui_traits.R "path/to/deposit"
# Columns are selected (no values are changed) except: catch time converted from hh:mm:ss to minutes after midnight
# (as the authors' scripts do with chron), lifespan also given in years (newlifespan / 365.25), and the malaria and
# provisioning files read with their ';' separator and ',' decimal mark. The telomere file is not in the deposit.
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[[1]] else stop("Give the folder of the deposit.")
out <- file.path("inst", "app", "data")
mins <- function(x) { p <- strsplit(as.character(x), ":", fixed = TRUE)
  vapply(p, function(v) if (length(v) >= 2 && all(nzchar(v))) 60 * as.numeric(v[1]) + as.numeric(v[2]) + (if (length(v) > 2) as.numeric(v[3]) / 60 else 0) else NA_real_, numeric(1)) }
yrs <- function(d) round(d / 365.25, 4)
w <- function(d, f) { utils::write.csv(d, file.path(out, f), row.names = FALSE); cat(f, nrow(d), "rows,", length(unique(d$BirdID)), "birds\n") }

oxy <- utils::read.csv(file.path(src, "oxy_28_5.csv"), stringsAsFactors = FALSE)
oxy$CatchTime_mins <- mins(oxy$CatchTime)
w(data.frame(BirdID = oxy$BirdID, age_year = oxy$age_year, sum_roms_z = oxy$sum_roms_z, sum_oxy_z = oxy$sum_oxy_z, SexEstimate = oxy$SexEstimate,
             newstat = oxy$newstat, avg_invert = oxy$avg_invert, CatchTime_mins = oxy$CatchTime_mins, birthyear = oxy$birthyear,
             newlifespan = oxy$newlifespan, lifespan_years = yrs(oxy$newlifespan)), "tsui_2026_warbler_oxidative.csv")

h <- utils::read.csv(file.path(src, "haematocrit_28_5.csv"), stringsAsFactors = FALSE)
h$CatchTime_mins <- mins(h$CatchTime)
w(data.frame(BirdID = h$BirdID, age = h$age, Cells_z = h$Cells_z, wbc_z = h$wbc_z, SexEstimate = h$SexEstimate, newstat = h$newstat,
             avg_invert = h$avg_invert, CatchTime_mins = h$CatchTime_mins, birthyear = h$birthyear, Observer = h$Observer,
             newlifespan = h$newlifespan, lifespan_years = yrs(h$newlifespan)), "tsui_2026_warbler_blood.csv")

m <- utils::read.csv(file.path(src, "malaria_28_5.csv"), stringsAsFactors = FALSE, sep = ";", dec = ",")
w(data.frame(BirdID = m$BirdID, age = m$age, TestResult = m$TestResult, SexEstimate = m$SexEstimate, newstat = m$newstat,
             avg_invert = as.numeric(m$avg_invert), birthyear = m$birthyear, newlifespan = m$newlifespan, lifespan_years = yrs(m$newlifespan)),
  "tsui_2026_warbler_malaria.csv")

p <- utils::read.csv(file.path(src, "provisioning_28_5.csv"), stringsAsFactors = FALSE, sep = ";", dec = ",")
w(data.frame(BirdID = p$BirdID, age = p$age, prate_z = p$prate_z, SexEstimate = p$SexEstimate, BroodSize = p$BroodSize, nr_helpers = p$nr_helpers,
             newstat = p$newstat, avg_invert = p$avg_invert, birthyear = p$birthyear, WatchType = p$WatchType, Observer = p$Observer,
             newlifespan = p$newlifespan, lifespan_years = yrs(p$newlifespan)), "tsui_2026_warbler_provisioning.csv")

r <- utils::read.csv(file.path(src, "repro_28_5.csv"), stringsAsFactors = FALSE)
w(data.frame(BirdID = r$BirdID, parent_age = r$parent_age, ars = r$ars, minrepro = r$minrepro, maxrepro = r$maxrepro,
             minreproage = r$minreproage, maxreproage = r$maxreproage, SexEstimate = r$sex, flno = r$flno, newstat = r$newstat,
             birthyear = r$birthyear, newlifespan = r$newlifespan, lifespan_years = yrs(r$newlifespan)), "tsui_2026_warbler_reproduction.csv")

s <- utils::read.csv(file.path(src, "survival_model.csv"), stringsAsFactors = FALSE)
s <- s[s$age > 0 & !is.na(s$SexEstimate), , drop = FALSE]                                  # as the authors' model13.R
w(data.frame(BirdID = s$BirdID, age = s$age, survival = s$survival, SexEstimate = s$SexEstimate, avg_bug = s$avg_bug, yearlytq = s$yearlytq,
             statbefore = s$statbefore, occasionyear = s$occasionyear, birthyear = s$birthyear), "tsui_2026_warbler_survival.csv")
