# App-level audit of step 4's "All individuals' trajectories": every bundled example (and extra trait) is loaded through
# the real server (shiny::testServer); the figure and its note are rendered with every colour option, the ageing
# functions are cycled on a subset, and the save button is pressed. Run from the package root:
#   Rscript tests/scripts/all_trajectories_app_audit.R     (exits with status 1 on any problem)
suppressMessages(library(shiny))
app <- shinyAppDir(normalizePath(file.path("inst", "app")))
problems <- character(0); rows <- list()
render <- function(expr, label) {
  v <- tryCatch(force(expr), shiny.silent.error = function(e) structure("validation", class = "silent"), error = function(e) e)
  if (inherits(v, "error")) { problems <<- c(problems, paste0(label, ": ", conditionMessage(v))); return(FALSE) }
  TRUE
}
testServer(app, {
  session$setInputs(model_family = "gaussian", model_age_function = "Quadratic", random_structure = "none", among_order = "linear",
                    standardise = TRUE, include_invalid = FALSE, a3_scale = "auto", a3_min_records = 0, a3_n = 20, a3_order = "random",
                    a3_all_points = TRUE, a3_all_alpha = 0.35)
  for (k in names(EXAMPLES)) {
    ex <- EXAMPLES[[k]]
    for (tr in unique(c(ex$mapping$trait, names(ex$traits)))) {
      m <- ex$mapping; tp <- example_trait_preset(ex, tr)
      cov <- tp$covars %||% m$covars %||% character(0); fac <- tp$cov_factor %||% m$cov_factor %||% character(0)
      lab <- paste0(k, if (!identical(tr, m$trait)) paste0(":", tr) else "")
      t0 <- Sys.time()
      session$setInputs(data_source = "example", example_id = k)
      session$setInputs(col_id = m$id, col_age = m$age, col_trait = tr, col_alr = m$alr %||% "__AUTO_LAST__", col_entry = m$entry %||% "__AUTO_FIRST__",
                        col_life = m$life %||% "", col_condition = m$condition %||% "", col_cov_num = setdiff(cov, fac), col_cov_fac = intersect(fac, cov),
                        col_cov_int = c(if (length(m$cov_age)) paste(m$cov_age, "age", sep = "|||"), m$cov_int %||% character(0)),
                        col_group = m$group %||% "", col_group2 = m$group2 %||% "", group_nested = isTRUE(m$nested), col_random = m$random %||% character(0),
                        col_censor = m$censor %||% "", censor_value = m$censor_value %||% "", col_trials = m$trials %||% "",
                        exclude_inconsistent = identical(m$inconsistent, "exclude"), trait_ages = !isFALSE(m$trait_specific_ages),
                        start_mode = m$start_mode %||% "afr", start_age = m$start_age %||% NA, age_round = m$age_round %||% NA)
      session$setInputs(subset_var = ex$subset$var %||% "", subset_levels = ex$subset$levels %||% character(0), a3_function = "Quadratic")
      d <- tryCatch(dat(), error = function(e) NULL)
      if (is.null(d) || !nrow(d)) { problems <- c(problems, paste0(lab, ": data did not load")); next }
      ch <- tryCatch(a3_all_colour_choices(), error = function(e) { problems <<- c(problems, paste0(lab, ": colour choices: ", conditionMessage(e))); character(0) })
      render(output$a3_all_colour_ui, paste(lab, "colour menu"))
      for (col in ch) {
        session$setInputs(a3_all_colour = col)
        render(output$a3_all_plot, paste(lab, "figure, colour", if (nzchar(col)) col else "none"))
      }
      ok_note <- render(output$a3_all_note, paste(lab, "note"))
      note <- tryCatch(paste(unlist(output$a3_all_note), collapse = " "), error = function(e) "")
      # cycle every ageing function on a few examples with different data shapes
      if (k %in% c("fly", "bouwhuis", "tsui_malaria", "meehan_dolphins", "sanghvi2021_flight"))
        for (fn in c(AGE_FUNCTIONS, A3_NONLINEAR)) { session$setInputs(a3_function = fn); render(output$a3_all_plot, paste(lab, "figure,", fn)) }
      session$setInputs(a3_function = "Quadratic", save_a3_all = sample.int(1e6, 1))
      rows[[length(rows) + 1]] <- data.frame(Example = lab, Individuals = length(unique(d$id)), `Colour options` = length(ch),
                                             Seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                                             Note = substr(gsub("<[^>]+>|\\s+", " ", note), 1, 110), check.names = FALSE)
      cat(sprintf("%-32s %5d individuals  %d colour options  %3.0f s\n", lab, length(unique(d$id)), length(ch), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  n_saved <- sum(vapply(saved_results(), function(e) identical(e$title, "All individuals' trajectories"), logical(1)))
  cat(sprintf("\n%d examples rendered; %d saved results of the figure\n", length(rows), n_saved))
  if (n_saved < length(rows)) problems <<- c(problems, sprintf("only %d of %d saves were stored", n_saved, length(rows)))
})
cat(if (length(problems)) paste(c("Problems:", paste(" -", unique(problems))), collapse = "\n") else "No problems.", "\n")
if (length(problems)) quit(status = 1)
