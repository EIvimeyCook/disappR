# Prepares the seed-beetle data of Sanghvi, Zajitschek, Iglesias-Carrasco and Head 2021 (Evolutionary Ecology 35:
# 367-385, 10.1007/s10682-021-10115-y) as bundled examples. Run with the authors' data workbook:
#   Rscript dev/prepare_sanghvi2021_beetles.R "path/to/EVEC_paper_final_data_Sanghvi_et_al_2021.xlsx"
# Uses the sheets 'Long format data flight' and 'Long format fecundity'. Columns are renamed and recoded for
# readability only: density H/L -> High/Low, sex F/M -> Female/Male, experiment 2/3 -> Virgin/Mated (the workbook's
# Notes: experiment 2 = flight assay only, 3 = mating and flight assay). No values are changed.
args <- commandArgs(trailingOnly = TRUE)
xl <- if (length(args)) args[[1]] else stop("Give the path of the workbook.")
out <- file.path("inst", "app", "data")
f <- as.data.frame(readxl::read_excel(xl, sheet = "Long format data flight"))
flight <- data.frame(id = f$id, container = f$batch, density = ifelse(f$density == "H", "High", "Low"),
                     sex = ifelse(f$sex == "F", "Female", "Male"), mating = ifelse(f$exp == 3, "Mated", "Virgin"),
                     dev_time = f$devtime, lifespan = f$lifespan, age = f$age, ability = f$ability, propensity = f$propensity)
utils::write.csv(flight, file.path(out, "sanghvi_2021_beetle_flight.csv"), row.names = FALSE)
cat("flight:", nrow(flight), "rows,", length(unique(flight$id)), "beetles\n")
e <- as.data.frame(readxl::read_excel(xl, sheet = "Long format fecundity"))
fec <- data.frame(id = e$ID, container = e$Batch, density = ifelse(e$Density == "H", "High", "Low"), dev_time = e[["Dev.time"]],
                  lifespan = e[["Adult.Lspan"]], day = e$Day, eggs = e$Eggs)
utils::write.csv(fec, file.path(out, "sanghvi_2021_beetle_fecundity.csv"), row.names = FALSE)
cat("fecundity:", nrow(fec), "rows,", length(unique(fec$id)), "females\n")
