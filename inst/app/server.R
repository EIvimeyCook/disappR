server <- function(input, output, session) {

  # Errors inside observers would otherwise end the session (the page turns grey and stops responding).
  # Every observer runs inside this guard: the error is reported on screen and in the R console instead.
  disappr_guard <- function(where, expr) {
    tryCatch(expr,
             shiny.silent.error = function(e) invisible(NULL),
             error = function(e) {
               msg <- paste0("Something went wrong (", where, "): ", conditionMessage(e))
               message("disappR: ", msg)
               try(showNotification(msg, type = "error", duration = 20), silent = TRUE)
               invisible(NULL)
             })
  }

  HEAT_COLOURS <- c("Observed" = "#7C8060", "Missed" = "#C0392B", "Death or ALR" = "#E67E22", "Not expected" = "#FFFFFF")

  # Help tab: video walkthroughs in a popup
  lapply(names(HELP_VIDEOS), function(k) {
    observeEvent(input[[paste0("play_", k)]], disappr_guard("play video", {
      showModal(help_video_modal(HELP_VIDEOS[[k]]))
    }))
  })

  # "i" help modals for every section
  lapply(names(INFO), function(key) {
    observeEvent(input[[paste0("info_", key)]], disappr_guard("input[[paste0('info_', key)]]", {
      if (is.function(session$sendModal)) {
        showModal(modalDialog(title = INFO[[key]]$title, info_body(INFO[[key]]), easyClose = TRUE, footer = modalButton("Close")))
      }
    }), ignoreInit = TRUE)
  })

  num_input <- function(x, default) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != 1 || !is.finite(v)) default else v
  }
  safe_get <- function(expr) tryCatch(expr, error = function(e) NULL)
  # Progress bars and notifications only when the session supports them (so the server
  # can also be exercised with shiny::testServer()).
  notify <- function(msg, type = "message", duration = 5) {
    if (is.function(session$sendNotification)) showNotification(msg, type = type, duration = duration)
    invisible(NULL)
  }
  run_with_progress <- function(msg, f) {
    if (is.function(session$sendProgress)) {
      withProgress(message = msg, value = 0, f(function(i, n, lab) incProgress(1 / n, detail = lab)))
    } else {
      f(NULL)
    }
  }
  metric_card <- function(label, value, sub = NULL, width = 3) {
    column(width, div(class = "metric-card", div(class = "metric-label", label), div(class = "metric-value", value),
                      if (!is.null(sub)) div(class = "metric-sub", sub)))
  }

  output$package_status <- renderUI({
    miss <- c(if (!HAS_GLMMTMB) "glmmTMB (count families)", if (!HAS_DHARMA) "DHARMa (residual checks)",
              if (!HAS_LMERTEST) "lmerTest (Gaussian p-values)")
    if (!length(miss)) return(p(class = "small-note", "All optional packages are installed."))
    div(class = "truth-card", strong("Optional packages not installed: "), paste(miss, collapse = "; "), ".")
  })

  observeEvent(input$go_toy, disappr_guard("input$go_toy", {
    updateRadioButtons(session, "data_source", selected = "toy")
    updateTabItems(session, "tabs", "data")
  }))
  observeEvent(input$go_fly, disappr_guard("input$go_fly", {
    updateRadioButtons(session, "data_source", selected = "example")
    updateTabItems(session, "tabs", "data")
  }))

  # ======================================================================
  # 1. Data
  # ======================================================================
  read_toy_inputs <- function() {
    individual <- identical(input$toy_afr_mode, "individual")
    sd_type <- input$toy_sd_type %||% "both"
    sa_type <- if (individual) (input$toy_sa_type %||% "none") else "none"
    list(
      trait = input$toy_trait %||% "mass",
      form = input$toy_form %||% "Quadratic",
      strength = input$toy_strength %||% "dramatic",
      sd_type = sd_type,
      sd_dir = if (identical(sd_type, "none")) 1 else num_input(input$toy_sd_dir, 1),
      rate_var = input$toy_rate_var %||% "low",
      mean_ls = min(30, max(3, round(num_input(input$toy_mean_ls, 20)))),
      missingness = input$toy_missingness %||% "complete",
      afr_mode = input$toy_afr_mode %||% "same",
      sa_type = sa_type,
      sa_dir = if (identical(sa_type, "none")) 1 else num_input(input$toy_sa_dir, 1),
      diet = isTRUE(input$toy_diet), groups = isTRUE(input$toy_groups),
      n_id = min(2000, max(30, round(num_input(input$toy_n, 300)))), seed = round(num_input(input$toy_seed, 1))
    )
  }
  toy_cfg <- reactiveVal(TOY_DEFAULTS)
  observeEvent(input$toy_commit, disappr_guard("input$toy_commit", {
    toy_cfg(read_toy_inputs())
    notify("Simulated dataset ready. All tabs now use it; refit models on the Models tab.", duration = 4)
  }))

  output$toy_status <- renderUI({
    active <- toy_cfg()
    pending <- read_toy_inputs()
    txt <- toy_truth_text(active)
    changed <- !identical(lapply(pending, as.character), lapply(utils::modifyList(TOY_DEFAULTS, active), as.character))
    tagList(
      div(class = "truth-card", strong("Active simulation. "), txt$sim, tags$ul(lapply(txt$expect, tags$li))),
      if (changed) div(class = "diagnosis-card", div(class = "diagnosis-detail", strong("Settings changed: "), "press 'Simulate & use this dataset' to apply them.")) else NULL
    )
  })

  raw_data <- reactive({
    src <- input$data_source %||% "toy"
    if (identical(src, "upload")) {
      validate(need(input$data_file, "Upload a CSV file on the Data tab to begin."))
      rd <- tryCatch(read_user_csv(input$data_file$datapath, input$data_file$name %||% ""),
                     error = function(e) list(data = NULL, note = paste("The file could not be read:", conditionMessage(e))))
      if (!is.list(rd)) rd <- list(data = NULL, note = "The file could not be read.")
      note <- if (is.character(rd$note) && length(rd$note) == 1 && !is.na(rd$note)) rd$note else "The file could not be read."
      x <- rd$data
      validate(need(is.data.frame(x), note))
      validate(need(ncol(x) >= 3 && nrow(x) >= 10, paste(note, "The app needs at least 3 columns and 10 rows: check the separator and header.")))
      attr(x, "read_note") <- note
      x
    } else if (identical(src, "example")) {
      ex <- current_example()
      x <- load_example_file(ex$file)
      validate(need(!is.null(x), paste0("The bundled dataset (data/", ex$file, ") was not found.")))
      x
    } else {
      simulate_toy_data(toy_cfg())
    }
  })

  current_example <- reactive({
    id <- input$example_id %||% names(EXAMPLES)[[1]]
    if (!id %in% names(EXAMPLES)) id <- names(EXAMPLES)[[1]]
    EXAMPLES[[id]]
  })
  output$example_note <- renderUI({
    ex <- current_example()
    div(class = "truth-card", strong("Empirical example. "), ex$note)
  })
  is_toy <- reactive(identical(input$data_source %||% "toy", "toy"))
  truth <- reactive(if (is_toy()) attr(raw_data(), "truth") else NULL)

  preset_map <- reactive({
    df <- raw_data()
    switch(input$data_source %||% "toy", toy = toy_mapping(df), example = current_example()$mapping, guess_mapping(df))
  })

  output$mapping_ui <- renderUI({
    df <- raw_data()
    pm <- preset_map()
    cols <- names(df)
    pick <- function(x, fallback = "") if (length(x) == 1 && !is.na(x) && x %in% cols) x else fallback
    pm_cov <- intersect(pm$covars %||% character(0), cols)
    pm_fac <- if (!is.null(pm$cov_factor)) intersect(pm$cov_factor, pm_cov) else pm_cov[!vapply(pm_cov, function(nm) is.numeric(maybe_numeric(df[[nm]])), logical(1))]
    pm_num <- setdiff(pm_cov, pm_fac)
    map_row <- fluidRow(
      column(4,
        selectInput("col_id", "Individual ID *", choices = cols, selected = pick(pm$id, cols[[1]])),
        selectInput("col_age", "Age / time *", choices = cols, selected = pick(pm$age, cols[[min(2, length(cols))]])),
        selectInput("col_trait", "Trait *", choices = cols, selected = pick(pm$trait, cols[[min(3, length(cols))]]))
      ),
      column(4,
        selectInput("col_alr", "ALR (age at last record)", choices = c("Automatic: last recorded age" = "__AUTO_LAST__", cols),
                    selected = pick(pm$alr, "__AUTO_LAST__")),
        selectInput("col_life", "Known lifespan (LS)", choices = c("(not available)" = "", "Automatic: last recorded age (proxy only)" = "__AUTO_LAST__", cols),
                    selected = if (identical(pm$life, "__AUTO_LAST__")) "__AUTO_LAST__" else pick(pm$life, ""))
      ),
      column(4,
        selectizeInput("col_cov_num", "Continuous fixed-effect covariates (linear effects)", choices = cols, selected = pm_num,
                       multiple = TRUE, options = list(placeholder = "e.g. temperature, body size")),
        selectizeInput("col_cov_fac", "Categorical fixed-effect covariates (factors)", choices = cols, selected = pm_fac,
                       multiple = TRUE, options = list(placeholder = "e.g. treatment, sex, diet")),
        uiOutput("cov_int_ui"),
        box_note("Fixed-effect covariates and their interactions are used in the mixed models on tab 5 (Modelling). Values of a continuous covariate that are not numbers are treated as missing."),
        selectInput("col_group", "Higher-level random effect", choices = c("(none)" = "", cols), selected = pick(pm$group, "")),
        checkboxInput("group_nested", "Individuals nested within this group", value = isTRUE(pm$nested)),
        selectizeInput("col_random", "Additional random intercepts: + (1 | X)", choices = cols,
                       selected = intersect(pm$random %||% character(0), cols), multiple = TRUE,
                       options = list(placeholder = "e.g. year, observer, replicate")),
        selectInput("col_condition", "Condition / state (optional)", choices = c("(none)" = "", cols), selected = pick(pm$condition, ""))
      )
    )
    read_note <- attr(df, "read_note")
    age_vals <- safe_numeric(df[[pick(pm$age, cols[[min(2, length(cols))]])]])
    min_age <- if (any(is.finite(age_vals))) min(age_vals[is.finite(age_vals)]) else 1
    tagList(
      if (!is.null(read_note)) p(class = "small-note", icon("file-csv"), " ", read_note) else NULL,
      map_row,
      conditionalPanel("input.data_source != 'toy'",
        fluidRow(column(4, numericInput("age_round", "Round ages to multiples of (optional)",
                                        value = if (isTRUE(is.finite(pm$age_round))) pm$age_round else NA, min = 0)),
                 column(8, p(class = "small-note", style = "margin-top:26px;",
                             "Leave blank for scheduled sampling. For irregular ages (e.g. from dates) enter the sampling resolution, such as 1 (year) or 7 (days), to place records on a common schedule.")))),
      tags$hr(),
      conditionalPanel("input.data_source != 'toy'",
        fluidRow(
          column(7,
            radioButtons("afr_mode", "Age at first record (AFR): when the individual enters the data",
                         choices = c("Automatic: each individual's first record (recommended)" = "auto",
                                     "Choose a column" = "column"),
                         selected = if (isTRUE(nzchar(pick(pm$entry, ""))) && !identical(pm$entry, "__AUTO_FIRST__")) "column" else "auto"),
            conditionalPanel("input.afr_mode == 'column'",
              selectInput("col_entry", "Column holding each individual's AFR", choices = cols,
                          selected = pick(pm$entry, cols[[1]]))),
            box_note("AFR is what every model and statistic uses: it is the AFR term in Models 7\u201310 and the start of each individual's observed window. An age at first trait expression can be set on the 'Sampling and missingness' tab, where it changes the missingness figures only.")
          ),
          column(5,
            selectInput("col_trials", "Number of binomial trials (optional)", choices = c("(none)" = "", cols), selected = pick(pm$trials, "")),
            selectInput("col_censor", "Censoring indicator (optional)", choices = c("(none)" = "", cols), selected = pick(pm$censor, "")),
            uiOutput("censor_value_ui")
          )
        )
      ),
      conditionalPanel("input.data_source == 'toy'",
        p(class = "small-note", "Simulated data: AFR (age at first observation), censoring and the start of sampling follow the simulation settings on the left."))
    )
  })
  try(outputOptions(output, "mapping_ui", suspendWhenHidden = FALSE), silent = TRUE)

  output$censor_value_ui <- renderUI({
    df <- raw_data()
    cc <- input$col_censor %||% ""
    if (!nzchar(cc) || !cc %in% names(df)) return(NULL)
    vals <- sort(unique(as.character(df[[cc]][!is.na(df[[cc]])])))
    if (!length(vals)) return(p(class = "small-note", "This column has no values."))
    if (length(vals) > 50) return(p(class = "small-note", "More than 50 distinct values: choose an indicator column (e.g. 0/1)."))
    pm <- isolate(preset_map())
    cur <- isolate(input$censor_value)
    sel <- if (isTRUE(cur %in% vals)) cur else if (isTRUE(pm$censor_value %in% vals)) pm$censor_value else vals[[1]]
    tagList(
      selectInput("censor_value", "Value marking a censored individual", choices = vals, selected = sel),
      box_note("LS is set to missing for censored individuals (e.g. escaped, or alive when the study ended).")
    )
  })
  try(outputOptions(output, "censor_value_ui", suspendWhenHidden = FALSE), silent = TRUE)

  # interactions: any covariate x age, or any pair of covariates ("a|||b"; "a|||age")
  output$cov_int_ui <- renderUI({
    cv <- unique(c(input$col_cov_num, input$col_cov_fac))
    if (!length(cv)) return(NULL)
    ch <- stats::setNames(paste(cv, "age", sep = "|||"), paste(cv, "\u00d7 age"))
    if (length(cv) > 1) {
      pr <- utils::combn(cv, 2)
      ch <- c(ch, stats::setNames(paste(pr[1, ], pr[2, ], sep = "|||"), paste(pr[1, ], "\u00d7", pr[2, ])))
    }
    pm <- isolate(preset_map())
    cur <- isolate(input$col_cov_int)
    pre <- c(if (length(pm$cov_age)) paste(pm$cov_age, "age", sep = "|||") else character(0), pm$cov_int %||% character(0))
    selectizeInput("col_cov_int", "Interactions (optional)", choices = ch, selected = intersect(if (!is.null(cur)) cur else pre, ch),
                   multiple = TRUE, options = list(placeholder = "e.g. diet \u00d7 age, diet \u00d7 sex"))
  })
  try(outputOptions(output, "cov_int_ui", suspendWhenHidden = FALSE), silent = TRUE)

  # ---- optional subset of the rows
  output$subset_ui <- renderUI({
    df <- safe_get(raw_data())
    if (is.null(df)) return(NULL)
    cand <- subset_candidates(df)
    cur <- isolate(input$subset_var)
    tagList(
      fluidRow(
        column(4, selectInput("subset_var", "Keep only rows where", choices = c("(no subsetting)" = "", cand),
                              selected = if (isTRUE(cur %in% cand)) cur else "")),
        column(8, uiOutput("subset_levels_ui"))
      ),
      uiOutput("subset_status")
    )
  })
  output$subset_levels_ui <- renderUI({
    df <- safe_get(raw_data())
    v <- input$subset_var %||% ""
    if (is.null(df) || !nzchar(v) || !v %in% names(df)) return(p(class = "small-note", style = "margin-top:26px;", "All rows are used."))
    lv <- sort(unique(as.character(df[[v]][!is.na(df[[v]])])))
    cur <- isolate(input$subset_levels)
    selectizeInput("subset_levels", paste(v, "is one of"), choices = lv, multiple = TRUE,
                   selected = if (length(intersect(cur, lv))) intersect(cur, lv) else lv[[1]])
  })
  output$subset_status <- renderUI({
    df <- safe_get(raw_data())
    du <- safe_get(raw_used())
    if (is.null(df) || is.null(du) || nrow(du) == nrow(df)) return(NULL)
    p(class = "small-note", sprintf("Subset: %s of %s rows are analysed (%s = %s).", format(nrow(du), big.mark = ","), format(nrow(df), big.mark = ","),
                                    input$subset_var, paste(input$subset_levels, collapse = ", ")))
  })
  try(outputOptions(output, "subset_ui", suspendWhenHidden = FALSE), silent = TRUE)
  try(outputOptions(output, "subset_levels_ui", suspendWhenHidden = FALSE), silent = TRUE)
  raw_used <- reactive({
    df <- raw_data()
    v <- input$subset_var %||% ""
    lv <- input$subset_levels %||% character(0)
    if (!nzchar(v) || !v %in% names(df) || !length(lv)) return(df)
    out <- df[!is.na(df[[v]]) & as.character(df[[v]]) %in% lv, , drop = FALSE]
    validate(need(nrow(out) >= 10, "The subset keeps fewer than 10 rows: choose more levels."))
    out
  })

  # ---- distribution of any variable
  output$dist_var_ui <- renderUI({
    cols <- names(raw_data())
    cur <- isolate(input$dist_var)
    tr <- preset_map()$trait
    sel <- if (isTRUE(cur %in% cols)) cur else if (isTRUE(tr %in% cols)) tr else cols[[1]]
    selectInput("dist_var", "Variable", choices = cols, selected = sel)
  })
  dist_values <- reactive({
    df <- raw_used()
    v <- input$dist_var
    validate(need(isTRUE(v %in% names(df)), "Choose a variable."))
    x <- df[[v]]
    if (identical(input$dist_unit, "individual")) {
      idc <- input$col_id
      validate(need(isTRUE(idc %in% names(df)), "Map the individual ID column first."))
      ids <- as.character(df[[idc]])
      ok <- !is.na(ids)
      xn <- safe_numeric(x)
      if (sum(is.finite(xn)) >= 0.8 * sum(!is.na(x) & nzchar(as.character(x)))) {
        x <- as.numeric(tapply(xn[ok], ids[ok], finite_mean))
      } else {
        x <- as.character(tapply(as.character(x)[ok], ids[ok], mode_or_na))
      }
    }
    x
  })
  output$dist_plot <- renderPlot(dist_plot_obj())
  dist_plot_obj <- reactive({
    x <- dist_values()
    v <- input$dist_var
    present <- !is.na(x) & nzchar(as.character(x))
    validate(need(any(present), "No non-missing values."))
    xn <- safe_numeric(x)
    unit <- if (identical(input$dist_unit, "individual")) "individuals" else "rows"
    if (sum(is.finite(xn)) >= 0.8 * sum(present)) {
      xn <- xn[is.finite(xn)]
      ux <- sort(unique(xn))
      sub <- sprintf("n = %d %s; mean %s, median %s, SD %s; %.0f%% zeros; %d missing",
                     length(xn), unit, format_num(mean(xn)), format_num(stats::median(xn)),
                     if (length(xn) > 1) format_num(stats::sd(xn)) else "NA", 100 * mean(xn == 0), sum(!present))
      p <- if (length(ux) <= 30 && all(abs(ux - round(ux)) < 1e-8)) {
        ggplot(data.frame(value = xn), aes(value)) + geom_bar(fill = warm_palette[[2]])
      } else {
        ggplot(data.frame(value = xn), aes(value)) + geom_histogram(bins = 30, fill = warm_palette[[2]], colour = "#FFFDF9")
      }
      p + labs(x = v, y = paste("Number of", unit), subtitle = sub) + theme_disappR(12)
    } else {
      tab <- sort(table(as.character(x[present])), decreasing = TRUE)
      n_lev <- length(tab)
      if (n_lev > 30) tab <- tab[1:30]
      df <- data.frame(level = factor(names(tab), levels = rev(names(tab))), n = as.numeric(tab))
      ggplot(df, aes(level, n)) + geom_col(fill = warm_palette[[2]], width = 0.7) + coord_flip() +
        labs(x = v, y = paste("Number of", unit),
             subtitle = paste0(n_lev, " distinct values", if (n_lev > 30) " (30 most common shown)" else "")) +
        theme_disappR(12)
    }
  })

  current_map <- reactive({
    req(input$col_id, input$col_age, input$col_trait)
    list(id = input$col_id, age = input$col_age, trait = input$col_trait,
         alr = input$col_alr %||% "__AUTO_LAST__", life = input$col_life %||% "",
         entry = if (identical(input$afr_mode %||% "auto", "column")) (input$col_entry %||% "__AUTO_FIRST__") else "__AUTO_FIRST__",
         condition = input$col_condition %||% "", trials = input$col_trials %||% "",
         covars = unique(c(input$col_cov_num, input$col_cov_fac)), cov_factor = input$col_cov_fac %||% character(0),
         cov_int = input$col_cov_int %||% character(0), group = input$col_group %||% "",
         nested = isTRUE(input$group_nested), random = input$col_random %||% character(0),
         censor = input$col_censor %||% "", censor_value = input$censor_value %||% "",
         cov_age = character(0),
         age_round = if (identical(input$data_source, "toy")) NA_real_ else num_input(input$age_round, NA_real_))
  })

  bundle <- reactive({
    df <- raw_used()
    mp <- current_map()
    validate(
      need(all(c(mp$id, mp$age, mp$trait) %in% names(df)), "Map the ID, age and trait columns."),
      need(length(unique(c(mp$id, mp$age, mp$trait))) == 3, "ID, age and trait must be three different columns.")
    )
    b <- standardise_data(df, mp, input$dup_action %||% "keep")
    sv <- input$subset_var %||% ""
    b$meta$subset <- if (nzchar(sv) && length(input$subset_levels)) list(var = sv, levels = input$subset_levels) else NULL
    d <- b$data
    validate(
      need(nrow(d) >= 10, paste0("Need at least 10 rows with an ID and a numeric age (found ", nrow(d),
                                 "). Check the age column: dates must be converted to ages, and decimal commas or text codes are not numbers.")),
      need(length(unique(d$id)) >= 3, "Need at least three individuals: check the ID column."),
      need(sum(is.finite(d$trait)) >= 10, paste0("Need at least 10 numeric trait values (found ", sum(is.finite(d$trait)),
                                                "): check the trait column for text, units or decimal commas."))
    )
    b
  })
  dat <- reactive(bundle()$data)
  meta <- reactive(bundle()$meta)
  imet <- reactive(individual_metrics(dat()))
  integrity <- reactive(data_integrity(dat(), meta()))
  life_known <- reactive(isTRUE(meta()$has_life) && !isTRUE(meta()$life_auto))

  data_sig <- reactive({
    d <- dat()
    m <- meta()
    paste(input$data_source, input$example_id %||% "", nrow(d), length(unique(d$id)), signif(sum(d$age), 12), signif(sum(d$trait, na.rm = TRUE), 12),
          signif(sum(d$alr, na.rm = TRUE), 12), signif(sum(d$life, na.rm = TRUE), 12), signif(sum(d$entry, na.rm = TRUE), 12),
          paste(m$covars, collapse = ","), paste(m$random_terms, collapse = ","), m$has_group, m$nested, m$life_auto,
          m$dup_action, m$n_censored, m$trials_col %||% "", input$afr_mode %||% "auto", paste(m$cov_age, collapse = ","), paste(m$cov_types, collapse = ","),
          paste(m$cov_pairs, collapse = ","), input$subset_var %||% "",
          paste(input$subset_levels %||% character(0), collapse = ","), sep = "|")
  })

  # Sensible defaults whenever the data change: family and model selection.
  observeEvent(safe_get(data_sig()), disappr_guard("data_sig()", {
    ex <- if (identical(input$data_source %||% "toy", "example")) current_example() else NULL
    fam <- if (is_toy()) {
      toy_family(toy_cfg()$trait)
    } else if (!is.null(ex$family)) {
      ex$family
    } else {
      integrity()$family$family
    }
    if (!identical(fam, "gaussian") && !HAS_GLMMTMB) fam <- "gaussian"
    updateSelectInput(session, "model_family", selected = fam)
    if (!is.null(ex$age_function)) updateSelectInput(session, "model_age_function", selected = ex$age_function)
    if (!is.null(ex$among)) updateRadioButtons(session, "among_order", selected = ex$among)
    av <- model_availability(dat(), meta())
    for (mid in MODEL_IDS) {
      on <- if (!is.null(ex$models)) (mid %in% ex$models && isTRUE(av$ok[[mid]])) else isTRUE(av$default[[mid]])
      updateCheckboxInput(session, paste0("use_", mid), value = on)
    }
  }))

  output$download_toy <- downloadHandler(
    filename = function() paste0("disappR_simulated_", toy_cfg()$trait, "_", toy_cfg()$missingness, ".csv"),
    content = function(file) utils::write.csv(simulate_toy_data(toy_cfg()), file, row.names = FALSE)
  )

  output$data_metrics <- renderUI({
    d <- dat()
    im <- imet()
    st <- infer_age_step(d$age, d$id)
    fam <- integrity()$family
    tagList(
      metric_card("Rows / individuals", paste0(format(nrow(d), big.mark = ","), " / ", format(nrow(im), big.mark = ",")),
                  paste0(sum(is.finite(d$trait)), " trait values")),
      metric_card("Ages", paste0(format_num(min(d$age)), "\u2013", format_num(max(d$age))),
                  paste0(length(unique(d$age)), " distinct; step ", if (is.finite(st)) format_num(st) else "NA")),
      metric_card("Records per individual", format_num(stats::median(im$n_trait)), "median"),
      metric_card("Trait type", fam$kind %||% (if (identical(fam$family, "gaussian")) "Continuous" else "Count"),
                  paste("Suggested:", family_label(fam$family)))
    )
  })

  output$integrity_table <- renderTable({
    integrity()$table
  }, striped = TRUE, spacing = "s")

  output$integrity_extra <- renderUI({
    bad <- integrity()$bad_ls_ids
    if (!length(bad)) return(NULL)
    p(class = "small-note", strong("IDs with LS < last recorded age: "), paste(utils::head(bad, 30), collapse = ", "),
      if (length(bad) > 30) " \u2026" else "")
  })

  output$data_preview <- renderTable(utils::head(raw_data(), 10), striped = TRUE, spacing = "xs")

  # ======================================================================
  # 2. A1-A2 visual diagnosis
  # ======================================================================
  output$visual_proxy_ui <- renderUI({
    ch <- proxy_choices(imet(), meta(), "all")
    labs <- c(ALR = "ALR (age at last record)", `Mean age` = "Mean age of the individual's records", LS = "LS (known lifespan)",
              AFR = "AFR (age at first observation): selective appearance")
    cur <- isolate(input$visual_proxy)
    tagList(
      selectInput("visual_proxy", "Group individuals by (both figures)", choices = stats::setNames(ch, unname(labs[ch])),
                  selected = if (isTRUE(cur %in% ch)) cur else "ALR"),
      if (!"AFR" %in% ch) p(class = "small-note", "AFR does not vary among individuals, so selective appearance cannot be diagnosed with these data.") else NULL
    )
  })
  visual_target <- reactive(if (identical(input$visual_proxy, "AFR")) "appearance" else "disappearance")
  try(outputOptions(output, "visual_proxy_ui", suspendWhenHidden = FALSE), silent = TRUE)

  # Facets: any column with 2-8 distinct values (mapped covariates first); one level per individual
  output$facet_ui <- renderUI({
    df <- safe_get(raw_data())
    m <- safe_get(meta())
    if (is.null(df) || is.null(m)) return(NULL)
    map <- m$map
    excl <- c(map$id, map$age, map$trait, map$alr, map$life, map$entry, map$censor)
    cand <- setdiff(names(df), excl)
    few <- cand[vapply(cand, function(nm) {
      v <- as.character(df[[nm]])
      u <- unique(v[!is.na(v) & nzchar(v)])
      length(u) >= 2 && length(u) <= 8
    }, logical(1))]
    ch <- unique(c(intersect(unname(m$cov_labels), few), few))
    cur <- isolate(input$facet_var)
    selectInput("facet_var", "Facet (separate panels) by", choices = c("(none)" = "", ch),
                selected = if (isTRUE(cur %in% ch)) cur else "")
  })
  try(outputOptions(output, "facet_ui", suspendWhenHidden = FALSE), silent = TRUE)
  facet_map <- reactive({
    v <- input$facet_var %||% ""
    if (!nzchar(v)) return(NULL)
    df <- raw_data()
    m <- meta()
    if (!v %in% names(df)) return(NULL)
    map <- m$map
    ids <- as.character(df[[map$id]])
    if (isTRUE(m$has_group) && isTRUE(m$nested)) {
      g <- as.character(df[[map$group]])
      ids <- ifelse(is.na(g), ids, paste(g, ids, sep = "/"))
    }
    lv <- as.character(df[[v]])
    ok <- !is.na(ids) & nzchar(ids) & !is.na(lv) & nzchar(lv)
    if (!any(ok)) return(NULL)
    per <- tapply(lv[ok], ids[ok], mode_or_na)
    stats::setNames(paste0(v, ": ", as.character(per)), names(per))
  })
  facet_layer <- function() if (!is.null(facet_map())) facet_wrap(~ facet) else NULL
  facet_text <- function() if (nzchar(input$facet_var %||% "")) paste0("; panels by ", input$facet_var) else ""

  # Bins range from 3 to the number of distinct values observed in >= the minimum number of individuals
  observe(disappr_guard("observer", {
    d <- safe_get(dat())
    im <- safe_get(imet())
    px <- safe_get(a1_px())
    if (is.null(d) || is.null(im) || is.null(px)) return(invisible(NULL))
    nmin <- 3
    step <- infer_age_step(d$age, d$id)
    mx1 <- max_bins_for(proxy_values(im, px), nmin, step)
    mx2 <- max_bins_for(id_age_means(d)$age, nmin, step)
    mx <- max(mx1, mx2)
    cur <- isolate(input$n_bins) %||% 4
    try(updateSliderInput(session, "n_bins", min = 3, max = mx, value = min(max(3, cur), mx)), silent = TRUE)
  }))

  visual_dat <- reactive({
    d <- dat()
    if (identical(input$trait_scale, "log1p")) {
      validate(need(all(d$trait[is.finite(d$trait)] > -1), "log(trait + 1) needs trait values above -1."))
      d$trait <- log1p(d$trait)
    }
    d
  })
  pick_proxy <- function(value, fallback_order) {
    ch <- proxy_choices(imet(), meta(), "all")
    validate(need(length(ch) > 0, "No usable grouping variable."))
    if (isTRUE(value %in% ch)) return(value)
    hit <- intersect(fallback_order, ch)
    if (length(hit)) hit[[1]] else ch[[1]]
  }
  vis_px <- reactive(pick_proxy(input$visual_proxy, c("ALR", "LS", "Mean age", "AFR")))
  a1_px <- reactive(vis_px())
  a2_px <- reactive(vis_px())
  trait_label <- reactive({
    lab <- current_map()$trait
    if (identical(input$trait_scale, "log1p")) paste0("log(", lab, " + 1)") else lab
  })
  a1_bins <- reactive({
    binned_trajectory(visual_dat(), proxy_values(imet(), a1_px()), input$n_bins %||% 4, input$bin_method %||% "equal",
                      3, facet = facet_map())
  })

  output$a1_plot <- renderPlot(a1_plot_obj())
  a1_plot_obj <- reactive({
    px <- a1_px()
    s <- a1_bins()
    validate(need(nrow(s) > 0, "Not enough variation in the proxy, or too few individuals per bin \u00d7 age point."))
    cols <- ordered_colours(length(levels(s$bin)))
    p <- ggplot(s, aes(age, mean, colour = bin, group = bin))
    if (isTRUE(input$show_se)) p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0, alpha = 0.6, na.rm = TRUE)
    p + geom_line(linewidth = 1.1) + geom_point(aes(size = n), alpha = 0.95) +
      scale_colour_manual(values = cols, drop = FALSE) + scale_size_area(max_size = 4.5, guide = "none") +
      labs(x = "Age", y = paste("Mean", trait_label()), colour = paste(px, "bin"),
           title = paste(trait_label(), "across age by", px, "bin"),
           subtitle = "Each line joins bin means of individual \u00d7 age means; point size = number of individuals") +
      facet_layer() + theme_disappR(13)
  })

  a2_data <- reactive({
    z <- trait_by_age_bins(visual_dat(), proxy_values(imet(), a2_px()), input$n_bins %||% 4, facet = facet_map())
    validate(need(nrow(z) > 0, "No individuals with both trait records and this variable."))
    n_ind <- stats::ave(rep(1, nrow(z)), z$facet, z$age_bin, FUN = length)
    z <- z[n_ind >= 3, , drop = FALSE]
    validate(need(nrow(z) > 0, "No age bin has at least 3 individuals."))
    z
  })
  a2_slopes <- reactive(a2_bin_slopes(a2_data()))
  a2_slope_display <- function(tab) {
    fmt_ci <- function(lo, hi) ifelse(is.finite(lo) & is.finite(hi),
                                      paste0(vapply(lo, format_num, character(1)), " to ", vapply(hi, format_num, character(1))), "NA")
    out <- data.frame(`Age bin` = tab$Age_bin, `Mean age` = signif(tab$Mean_age, 3), N = tab$N,
                      Coefficient = signif(tab$Coefficient, 3), SE = signif(tab$SE, 3), `95% CI` = fmt_ci(tab$Lower_95, tab$Upper_95),
                      r = round(tab$r, 2), P = vapply(tab$P, format_p, character(1)),
                      `Change from previous bin` = signif(tab$Change_from_previous, 3), check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(facet_map())) out <- cbind(Panel = tab$Facet, out, stringsAsFactors = FALSE)
    out
  }
  output$a2_slope_table <- renderTable({
    tab <- a2_slopes()
    validate(need(nrow(tab) > 0, "No age bin with enough individuals."))
    a2_slope_display(tab)
  }, striped = TRUE, spacing = "xs")
  a2_trend_lines <- reactive({
    tab <- a2_slopes()
    tr <- a2_slope_trend(tab)
    px <- a2_px()
    if (!nrow(tr)) return("Too few age bins with estimable coefficients (at least 3 needed) to test for a trend across age.")
    vapply(seq_len(nrow(tr)), function(i) {
      sprintf("%sThe coefficient of %s changes by %s per unit age across %d age bins (inverse-variance weighted trend, p = %s).",
              if (!is.null(facet_map())) paste0(tr$Facet[i], ": ") else "", px, format_num(tr$Change_per_age[i]),
              tr$N_bins[i], format_p(tr$P[i]))
    }, character(1))
  })
  output$a2_trend_note <- renderUI({
    txt <- safe_get(a2_trend_lines())
    if (is.null(txt)) return(NULL)
    div(class = "small-note", lapply(txt, p))
  })
  output$a2_plot <- renderPlot(a2_plot_obj())
  a2_plot_obj <- reactive({
    px <- a2_px()
    z <- a2_data()
    cols <- ordered_colours(length(levels(z$age_bin)))
    p <- ggplot(z, aes(proxy, trait, colour = age_bin))
    if (isTRUE(input$a2_points %||% TRUE)) {
      jw <- diff(range(z$proxy, na.rm = TRUE))
      jw <- if (is.finite(jw) && jw > 0) 0.01 * jw else 0
      p <- p + geom_point(alpha = 0.35, size = 1.4, position = position_jitter(width = jw, height = 0, seed = 1))
    }
    p +
      geom_smooth(aes(group = age_bin), method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.2, na.rm = TRUE) +
      scale_colour_manual(values = cols, drop = FALSE) +
      labs(x = px, y = paste("Mean", trait_label(), "within age bin"), colour = "Age",
           title = paste(trait_label(), "against", px, "within age bins"),
           subtitle = if (isTRUE(input$a2_points %||% TRUE)) "Points (jittered): each individual's mean within an age bin; lines: linear fit per age bin" else "Lines: linear fit of the trait on the proxy within each age bin (points hidden)") +
      facet_layer() + theme_disappR(13)
  })

  bin_diff_data <- reactive({
    s <- a1_bins()
    validate(need(nrow(s) > 0, "No bins to compare (see the trajectory plot)."))
    dz <- bin_differences(s, "successive")
    validate(need(nrow(dz) > 0, "No age at which two bins both meet the minimum number of individuals."))
    dz
  })
  bin_diff_trends <- reactive(bin_difference_trends(bin_diff_data()))
  # one least-squares line of difference against age across all bin pairs (per panel)
  bin_diff_pooled <- reactive({
    dz <- bin_diff_data()
    rows <- lapply(split(dz, dz$facet), function(x) {
      if (length(unique(x$age)) < 2) return(NULL)
      fit <- stats::lm(difference ~ age, data = x)
      cf <- suppressWarnings(summary(fit)$coefficients)
      if (!"age" %in% rownames(cf)) return(NULL)
      data.frame(facet = x$facet[[1]], x = min(x$age), xend = max(x$age),
                 y = cf[1, 1] + cf["age", 1] * min(x$age), yend = cf[1, 1] + cf["age", 1] * max(x$age),
                 Slope_per_age = cf["age", 1], P = if (nrow(x) > 2 && ncol(cf) >= 4) cf["age", 4] else NA_real_,
                 N_points = nrow(x), stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) do.call(rbind, rows) else data.frame()
  })
  output$bin_diff_plot <- renderPlot(bin_diff_obj())
  bin_diff_obj <- reactive({
    dz <- bin_diff_data()
    tr <- bin_diff_trends()
    ages <- sort(unique(dz$age))
    st <- if (length(ages) > 1) min(diff(ages)) else 1
    np <- length(levels(dz$pair))
    dz$age_plot <- dz$age + (as.integer(dz$pair) - (np + 1) / 2) * st * min(0.6 / np, 0.15)
    p <- ggplot(dz, aes(age_plot, difference, colour = pair)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
      geom_point(size = 3, alpha = 0.9)
    if (identical(input$diff_lines, "pooled")) {
      pl <- bin_diff_pooled()
      if (nrow(pl)) {
        p <- p + geom_segment(data = pl, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, colour = "black", linewidth = 1.2)
      }
    } else if (nrow(tr)) {
      seg <- data.frame(facet = tr$Facet, pair = factor(tr$Pair, levels = levels(dz$pair)), x = tr$Age_from, xend = tr$Age_to,
                        y = tr$Intercept + tr$Slope_per_age * tr$Age_from, yend = tr$Intercept + tr$Slope_per_age * tr$Age_to)
      p <- p + geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend, colour = pair), inherit.aes = FALSE, linewidth = 1)
    }
    p +
      scale_colour_manual(values = distinct_colours(np), drop = FALSE) +
      scale_x_continuous(breaks = if (length(ages) <= 15) ages else waiver()) +
      labs(x = "Age", y = paste("Difference in mean", trait_label(), "(higher \u2212 lower", a1_px(), "bin)"), colour = NULL,
           title = paste("Difference between", a1_px(), "bins at each age"),
           subtitle = if (identical(input$diff_lines, "pooled")) "Bin numbers as in the trajectory legend; black line: least-squares trend across all pairs. A rising or falling line indicates age-dependent selection." else "Bin numbers as in the trajectory legend; lines: least-squares trend of each difference against age. Trends that rise or fall indicate age-dependent selection.") +
      facet_layer() + theme_disappR(13)
  })

  output$toy_card_visual <- renderUI({
    if (!is_toy()) return(NULL)
    txt <- toy_truth_text(toy_cfg())
    div(class = "truth-card", strong("What you should see"), tags$ul(lapply(txt$expect, tags$li)))
  })

  # ======================================================================
  # 3. Sampling, missingness and proxies
  # ======================================================================
  # The missingness window may open at the age at first trait expression; nothing else in the app uses it.
  miss_start_age <- reactive({
    if (identical(input$miss_start %||% "afr", "afe")) num_input(input$afe_age, NA_real_) else NA_real_
  })
  grid_r <- reactive({
    afe <- miss_start_age()
    build_missing_grid(dat(), margin = 0,
                       start_mode = if (isTRUE(is.finite(afe))) "same" else "afr",
                       start_age = afe)
  })
  msum <- reactive(missingness_summary(dat(), grid_r(), 0.05, life_known()))

  output$sampling_metrics <- renderUI({
    ms <- msum()
    tagList(
      metric_card("Missing expected occasions", if (is.finite(ms$percent)) sprintf("%.1f%%", ms$percent) else "NA",
                  ms$type, width = 3),
      metric_card("Detection p", if (is.finite(ms$p)) sprintf("%.2f", ms$p) else "NA",
                  "Share of expected occasions with a trait record", width = 3),
      metric_card("r (mean age, ALR)", if (is.finite(ms$r_mean_alr)) sprintf("%.3f", ms$r_mean_alr) else "NA",
                  "Near 1 = the two proxies are interchangeable", width = 3),
      metric_card("r (ALR, LS) / r (mean age, LS)",
                  if (is.finite(ms$r_alr_ls)) sprintf("%.2f / %.2f", ms$r_alr_ls, ms$r_mean_ls) else "LS unknown",
                  "How well each proxy tracks known lifespan", width = 3)
    )
  })

  output$sampling_guidance <- renderUI({
    ms <- msum()
    tagList(
      div(class = "diagnosis-card", div(class = "diagnosis-title", "Proxy choice (ALR or mean age)"), div(class = "diagnosis-detail", ms$guidance)),
      div(class = "diagnosis-detail small-note", strong("Missing cells: "), ms$pattern)
    )
  })

  heat_n <- reactive(min(500, max(20, round(num_input(input$heat_n, 150)))))
  output$heatmap <- renderPlot(heatmap_obj(), height = function() {
    n <- safe_get(min(heat_n(), length(unique(grid_r()$id)))) %||% 150
    max(460, min(1400, round(2.2 * n) + 200))
  })
  heatmap_obj <- reactive({
    g <- grid_r()
    validate(need(nrow(g) > 0, "Could not build the sampling grid."))
    im <- imet()
    ids <- unique(g$id)
    mi <- match(ids, im$id)
    ord <- switch(input$heat_order %||% "alr",
                  alr = ids[order(im$alr[mi], im$first_recorded[mi])],
                  afr = ids[order(im$entry[mi], im$alr[mi])],
                  trait = ids[order(im$trait_mean[mi], im$alr[mi])],
                  id = sort(ids),
                  random = ids[order(stats::runif(length(ids)))])
    n_show <- heat_n()
    if (length(ord) > n_show) ord <- ord[unique(round(seq(1, length(ord), length.out = n_show)))]
    x <- grid_display(g, ord)
    validate(need(nrow(x) > 0, "Could not build the sampling grid."))
    n_shown <- attr(x, "n_shown") %||% length(ord)
    final_lab <- if (life_known()) "Death (known LS)" else "Last record (ALR)"
    labs_map <- c("Observed" = "Observed", "Missed" = "Missed (expected, not recorded)", "Death or ALR" = final_lab,
                  "Not expected" = "No observations expected (before AFR, after death/ALR)")
    x$id_f <- factor(x$id, levels = rev(ord))
    x$status <- factor(x$status, levels = names(HEAT_COLOURS))
    ggplot(x, aes(age, id_f, fill = status)) +
      geom_tile(width = g$step[[1]], colour = "#E9E1D6", linewidth = 0.15) +
      scale_fill_manual(values = HEAT_COLOURS, breaks = names(HEAT_COLOURS), labels = unname(labs_map[names(HEAT_COLOURS)]), drop = FALSE) +
      labs(x = "Age (scheduled occasions)", y = NULL, fill = NULL,
           subtitle = if (n_shown < length(ids)) paste("Showing", n_shown, "of", length(ids), "individuals (evenly spaced in the chosen order)") else NULL) +
      theme_disappR(12) +
      theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA), legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 2))
  })

  output$missing_by_age <- renderPlot(missing_by_age_obj())
  missing_by_age_obj <- reactive({
    g <- grid_r()
    validate(need(nrow(g) > 0, "No sampling grid."))
    cv <- coverage_by_age(g)
    sc <- max(cv$Expected) / 100
    ggplot(cv, aes(Age)) +
      geom_col(aes(y = Observed), fill = "#E5D2BE", width = 0.8 * g$step[[1]]) +
      geom_line(aes(y = (100 - Percent_observed) * sc), colour = warm_palette[[1]], linewidth = 1.1) +
      geom_point(aes(y = (100 - Percent_observed) * sc), colour = warm_palette[[1]], size = 2.2) +
      scale_y_continuous("Individuals observed (bars)", sec.axis = sec_axis(~ . / sc, name = "% of expected occasions missing (line)")) +
      labs(x = "Age") + theme_disappR(12)
  })

  output$coverage_table <- renderTable({
    cv <- coverage_by_age(grid_r())
    if (!nrow(cv)) return(NULL)
    if (nrow(cv) > 25) cv <- cv[unique(round(seq(1, nrow(cv), length.out = 25))), , drop = FALSE]
    cv$Age <- format_num(cv$Age)
    cv
  }, striped = TRUE, spacing = "xs", digits = 1)

  output$missing_vs_var <- renderPlot(missing_vs_var_obj())
  missing_vs_var_obj <- reactive({
    g <- grid_r()
    validate(need(nrow(g) > 0, "No sampling grid."))
    x <- missingness_by_variable(dat(), g, input$miss_var %||% "ALR")
    validate(need(nrow(x) > 0, "This variable is not available for these data."))
    if (identical(x$type[[1]], "categorical")) {
      ggplot(x, aes(stats::reorder(group, missing_rate), 100 * missing_rate)) + geom_col(fill = warm_palette[[1]], width = 0.65) +
        coord_flip() + labs(x = NULL, y = "% missing") + theme_disappR(12)
    } else {
      ggplot(x, aes(value, 100 * missing_rate)) + geom_line(colour = warm_palette[[1]], linewidth = 1) +
        geom_point(colour = warm_palette[[1]], size = 2.2) +
        labs(x = input$miss_var, y = "% missing", subtitle = if (identical(x$type[[1]], "binned")) "Quantile bins of individuals" else NULL) +
        theme_disappR(12)
    }
  })

  # A trait-dependent missingness signal mimics age-dependent selective disappearance, so it is flagged here
  output$missing_trait_note <- renderUI({
    d <- safe_get(msum())
    d <- if (is.null(d)) NULL else d$drivers
    if (is.null(d) || !nrow(d)) return(NULL)
    hit <- d[d$Predictor %in% c("Prior observed trait", "Condition") & d$Signal %in% "Associated", , drop = FALSE]
    if (!nrow(hit)) return(NULL)
    div(class = "diagnosis-card",
        div(class = "diagnosis-detail", strong("Missingness depends on the trait itself "),
            sprintf("(%s: %s, p = %s).", hit$Predictor[[1]], hit$Effect[[1]], format_p(hit$P_value[[1]])),
            " Occasions are missed according to the trait value, so individuals with low values are recorded less often and their records end earlier.",
            " This produces the same model ranking as age-dependent selective disappearance: Models 3, 4 and 5 can beat Model 1 decisively when no selection is present.",
            " A win for those models cannot be read as evidence of selective disappearance here; check the trait-before-death figure and report this association alongside the model comparison."))
  })

  output$drivers_table <- renderTable({
    d <- msum()$drivers
    if (!nrow(d)) return(data.frame(Note = "Too few missing cells for association tests."))
    d$P_value <- vapply(d$P_value, format_p, character(1))
    d
  }, striped = TRUE, spacing = "xs")

  proxy_im <- reactive({
    im <- imet()
    im[im$n_trait > 0, , drop = FALSE]
  })
  proxy_pair_plot <- function(x, y, xlab, ylab) {
    ok <- is.finite(x) & is.finite(y)
    validate(need(sum(ok) >= 4, "Not enough individuals with both values."))
    df <- data.frame(x = x[ok], y = y[ok])
    r <- cor_safe(df$x, df$y)
    ggplot(df, aes(x, y)) +
      geom_count(colour = warm_palette[[3]], alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey50") +
      annotate("text", x = -Inf, y = Inf, label = if (is.finite(r)) sprintf("r = %.3f", r) else "r = NA",
               hjust = -0.2, vjust = 1.5, colour = "#453A32") +
      scale_size_area(max_size = 6, guide = "none") +
      labs(x = xlab, y = ylab, subtitle = paste(ylab, "vs", xlab, "(one point per individual; dotted line: 1:1)")) +
      theme_disappR(12)
  }
  output$afr_alr_plot <- renderPlot({
    d <- safe_get(dat())
    if (is.null(d)) return(NULL)
    im <- individual_metrics(d)
    im <- im[im$n_trait > 0, , drop = FALSE]
    df <- data.frame(AFR = ifelse(is.finite(im$entry), im$entry, im$first_recorded), ALR = im$alr)
    df <- df[is.finite(df$AFR) & is.finite(df$ALR), , drop = FALSE]
    validate(need(nrow(df) >= 5, "Fewer than 5 individuals with both AFR and ALR."))
    validate(need(stats::sd(df$AFR) > 0, "AFR is the same for every individual, so it cannot be correlated with ALR."))
    r <- stats::cor(df$AFR, df$ALR)
    ggplot(df, aes(AFR, ALR)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
      geom_count(alpha = 0.55, colour = warm_palette[[2]]) +
      geom_smooth(method = "lm", formula = y ~ x, colour = warm_palette[[1]], fill = warm_palette[[1]], alpha = 0.18) +
      scale_size_area(max_size = 6, guide = "none") +
      labs(x = "Age at first record (AFR)", y = "Age at last record (ALR)",
           subtitle = sprintf("r = %s; points on the dashed line were recorded once", format_num(r))) +
      theme_disappR(12)
  })

  output$afr_alr_table <- renderTable({
    a <- safe_get(afr_alr_agreement(dat(), meta()))
    if (is.null(a) || !isTRUE(a$ok)) return(data.frame(Measure = "AFR vs ALR", Value = if (is.null(a)) "not available" else a$message))
    a$table
  }, striped = TRUE, spacing = "xs", width = "100%")

  output$proxy_plots_ui <- renderUI({
    if (life_known()) {
      fluidRow(column(4, plotOutput("proxy_alr_mean", height = 330)),
               column(4, plotOutput("proxy_alr_ls", height = 330)),
               column(4, plotOutput("proxy_mean_ls", height = 330)))
    } else {
      fluidRow(column(6, plotOutput("proxy_alr_mean", height = 340)),
               column(6, p(class = "small-note", style = "margin-top:20px;",
                           "Lifespan (LS) is not known for these data, so only ALR and mean age are compared. Map a known LS column on the Data tab to add ALR vs LS and mean age vs LS.")))
    }
  })
  proxy_alr_mean_obj <- reactive({
    im <- proxy_im()
    proxy_pair_plot(im$mean_age, im$alr, "Mean age", "ALR")
  })
  proxy_alr_ls_obj <- reactive({
    validate(need(life_known(), "LS unknown."))
    im <- proxy_im()
    proxy_pair_plot(im$lifespan, im$alr, "LS", "ALR")
  })
  proxy_mean_ls_obj <- reactive({
    validate(need(life_known(), "LS unknown."))
    im <- proxy_im()
    proxy_pair_plot(im$lifespan, im$mean_age, "LS", "Mean age")
  })
  output$proxy_alr_mean <- renderPlot(proxy_alr_mean_obj())
  output$proxy_alr_ls <- renderPlot(proxy_alr_ls_obj())
  output$proxy_mean_ls <- renderPlot(proxy_mean_ls_obj())


  # ======================================================================
  # 4. A3 individual parametric fits
  # ======================================================================
  # Minimum records per function: see min_records_for() (Linear 2; Quadratic, logarithmic, exponential 3; Cubic 4).
  a3_min_df <- function() 1
  a3_n_par <- function() {
    fn <- input$a3_function %||% "Quadratic"
    if (identical(fn, A3_NONLINEAR)) 2 else ncol(age_basis(1:3, make_age_params(1:10, fn, FALSE))) + 1
  }

  a3_ids <- reactive({
    im <- imet()
    im <- im[im$n_trait > 0, , drop = FALSE]
    k <- min(nrow(im), input$a3_n %||% 20)
    switch(input$a3_order %||% "random",
           longest = im$id[order(-im$alr)][seq_len(k)],
           shortest = im$id[order(im$alr)][seq_len(k)],
           im$id[order(stats::runif(nrow(im)))][seq_len(k)])
  })

  a3_fit <- reactive(fit_individual_function(dat(), input$a3_function %||% "Quadratic", a3_min_df(), a3_ids()))
  a3_compare <- reactive(compare_individual_functions(dat(), a3_min_df()))

  output$a3_metrics <- renderUI({
    z <- a3_fit()
    im <- imet()
    im <- im[im$n_trait > 0, , drop = FALSE]
    fitted <- im$id %in% z$fitted_ids
    alr_all <- mean(im$alr, na.rm = TRUE)
    alr_fit <- if (any(fitted)) mean(im$alr[fitted], na.rm = TRUE) else NA_real_
    sd_alr <- stats::sd(im$alr, na.rm = TRUE)
    biased <- is.finite(alr_fit) && isTRUE(sd_alr > 0) && (alr_fit - alr_all) / sd_alr > 0.25
    tagList(
      div(class = "metric-card", div(class = "metric-label", "Individuals fitted"),
          div(class = "metric-value", paste0(length(z$fitted_ids), " / ", z$n_total)),
          div(class = "metric-sub", paste("Need \u2265", min_records_for(input$a3_function %||% "Quadratic"), "records at different ages"))),
      div(class = "metric-card", div(class = "metric-label", "Mean ALR: fitted vs all"),
          div(class = "metric-value", if (is.finite(alr_fit)) paste0(format_num(alr_fit), " vs ", format_num(alr_all)) else "NA"),
          div(class = "metric-sub", if (biased) "Fitted individuals are longer-lived: the individual fits are survivor-biased here." else "No strong lifespan selection among fitted individuals.")),
      if (nrow(z$fit_stats)) div(class = "metric-card", div(class = "metric-label", "Median adjusted R\u00b2"),
          div(class = "metric-value", format_num(stats::median(z$fit_stats$Adjusted_R2, na.rm = TRUE))),
          div(class = "metric-sub", paste("Mean R\u00b2", format_num(mean(z$fit_stats$R2, na.rm = TRUE))))) else NULL
    )
  })

  output$a3_plot <- renderPlot({
    d <- dat()
    ids <- a3_ids()
    z <- a3_fit()
    raw <- d[d$id %in% ids & is.finite(d$trait), , drop = FALSE]
    validate(need(nrow(raw) > 0, "No trait records for the selected individuals."))
    ids <- ids[ids %in% raw$id]
    stat <- z$fit_stats
    lab <- stats::setNames(paste0(ids, "\n", "not estimable"), ids)
    if (nrow(stat)) {
      hit <- ids[ids %in% stat$id]
      s <- stat[match(hit, stat$id), , drop = FALSE]
      lab[hit] <- paste0(hit, "\nadj R\u00b2 ", ifelse(is.finite(s$Adjusted_R2), sprintf("%.2f", s$Adjusted_R2), "NA"))
    }
    raw$panel <- factor(lab[raw$id], levels = unique(lab[ids]))
    p <- ggplot(raw, aes(age, trait)) + geom_point(colour = "#8A7564", size = 1.3, alpha = 0.8)
    if (nrow(z$curves)) {
      cv <- z$curves[z$curves$id %in% ids, , drop = FALSE]
      cv$panel <- factor(lab[cv$id], levels = levels(raw$panel))
      p <- p + geom_line(data = cv, aes(age, fitted), colour = warm_palette[[1]], linewidth = 0.9)
    }
    p + facet_wrap(~ panel, ncol = 5, scales = "free") +
      labs(x = "Age", y = current_map()$trait, title = paste(input$a3_function, "fits to individual trajectories")) +
      theme_disappR(10) + theme(strip.text = element_text(size = 8, colour = "#5B493C"))
  }, height = function() 170 * ceiling(max(1, length(a3_ids())) / 5) + 60)

  output$a3_mean_plot <- renderPlot(a3_mean_obj())
  a3_mean_obj <- reactive({
    z <- a3_fit()
    validate(need(nrow(z$mean_curve) > 0, "No individual has enough records to fit this function."))
    mc <- z$mean_curve
    recon <- input$a3_recon %||% "both"
    lab_c <- "Mean of coefficients"
    lab_f <- "Mean of individual functions"
    lab_t <- "True (simulated)"
    obs <- observed_trajectory(dat())
    obs <- obs[obs$age >= min(mc$age) & obs$age <= max(mc$age), , drop = FALSE]
    lines <- data.frame(age = numeric(0), fitted = numeric(0), Method = character(0))
    bands <- data.frame(age = numeric(0), lo = numeric(0), hi = numeric(0), Method = character(0))
    if (recon %in% c("coef", "both")) {
      lines <- rbind(lines, data.frame(age = mc$age, fitted = mc$fitted, Method = lab_c))
      bands <- rbind(bands, data.frame(age = mc$age, lo = mc$lo, hi = mc$hi, Method = lab_c))
    }
    if (recon %in% c("fun", "both")) {
      lines <- rbind(lines, data.frame(age = mc$age, fitted = mc$fitted_fun, Method = lab_f))
      bands <- rbind(bands, data.frame(age = mc$age, lo = mc$lo_fun, hi = mc$hi_fun, Method = lab_f))
    }
    tc <- truth()
    if (!is.null(tc)) lines <- rbind(lines, data.frame(age = mc$age, fitted = toy_true_curve(tc, mc$age), Method = lab_t))
    lines <- lines[is.finite(lines$fitted), , drop = FALSE]
    bands <- bands[is.finite(bands$lo) & is.finite(bands$hi), , drop = FALSE]
    cols <- stats::setNames(c("#D55E00", "#0072B2", "black"), c(lab_c, lab_f, lab_t))
    note <- if (isTRUE(z$nonlinear)) {
      "Non-linear in its parameters: the two reconstructions differ (Jensen's inequality)."
    } else {
      "Linear in its parameters: the two reconstructions are identical (lines overlap)."
    }
    p <- ggplot()
    if (nrow(bands)) p <- p + geom_ribbon(data = bands, aes(x = age, ymin = lo, ymax = hi, fill = Method), alpha = 0.15)
    p <- p + geom_point(data = obs, aes(age, fitted, size = n), colour = "grey55", alpha = 0.7) +
      geom_line(data = lines, aes(age, fitted, colour = Method, linetype = Method), linewidth = 1.15) +
      scale_colour_manual(values = cols) + scale_fill_manual(values = cols, guide = "none") +
      scale_linetype_manual(values = stats::setNames(c("solid", "dotdash", "longdash"), c(lab_c, lab_f, lab_t))) +
      scale_size_area(max_size = 4, guide = "none") +
      labs(x = "Age", y = current_map()$trait, colour = NULL, linetype = NULL,
           subtitle = paste0(note, " Bands: 95% intervals from ", length(z$fitted_ids),
                             " individual fits; grey points: observed means; ages observed in \u2265 10 individuals")) +
      theme_disappR(12)
    if (nrow(obs)) {
      # keep the plot readable if an extreme non-linear fit inflates the mean of functions
      span <- diff(range(obs$fitted))
      if (!is.finite(span) || span <= 0) span <- max(abs(obs$fitted), 1)
      lim <- range(obs$fitted) + c(-3, 3) * span
      yv <- c(obs$fitted, lines$fitted, bands$lo, bands$hi)
      yv <- yv[is.finite(yv) & yv >= lim[1] & yv <= lim[2]]
      if (length(yv)) p <- p + coord_cartesian(ylim = range(yv))
    }
    p
  })

  output$a3_compare_table <- renderTable({
    x <- a3_compare()
    if (!nrow(x)) return(data.frame(Note = "No function could be fitted to any individual."))
    x$Mean_R2 <- round(x$Mean_R2, 3)
    x$Mean_adj_R2 <- round(x$Mean_adj_R2, 3)
    x$Median_adj_R2 <- round(x$Median_adj_R2, 3)
    x$Mean_dAICc <- round(x$Mean_dAICc, 2)
    x
  }, striped = TRUE, spacing = "xs")

  output$a3_coef_table <- renderTable({
    cf <- a3_fit()$coefs
    if (!nrow(cf)) return(data.frame(Note = "No coefficients."))
    nms <- setdiff(names(cf), "id")
    data.frame(Coefficient = nms,
               Mean = vapply(cf[nms], function(v) mean(v, na.rm = TRUE), numeric(1)),
               SD = vapply(cf[nms], function(v) stats::sd(v, na.rm = TRUE), numeric(1)),
               N = vapply(cf[nms], function(v) sum(is.finite(v)), integer(1)), check.names = FALSE)
  }, digits = 5, striped = TRUE, spacing = "xs")

  # ======================================================================
  # 6. Model settings shared by B2 and B4
  # ======================================================================
  model_settings <- reactive({
    fam <- input$model_family %||% "gaussian"
    zi_vars <- input$zi_vars %||% character(0)
    ex <- stats::setNames(lapply(MODEL_IDS, function(mid) as.character(input[[paste0("extra_", mid)]] %||% character(0))), MODEL_IDS)
    ex <- ex[vapply(ex, length, integer(1)) > 0]
    list(family = fam,
         zi = if (fam %in% c("zip", "zinb") && length(zi_vars)) paste("~", paste(zi_vars, collapse = " + ")) else "~1",
         age_function = input$model_age_function %||% "Quadratic",
         models = MODEL_IDS[vapply(MODEL_IDS, function(mid) isTRUE(input[[paste0("use_", mid)]]), logical(1))],
         random_slope = normalise_slope(input$random_structure %||% "none"),
         standardise = isTRUE(input$standardise),
         among = if (identical(input$among_order, "same")) "same" else "linear",
         include_invalid = isTRUE(input$include_invalid),
         extra = if (length(ex)) ex else NULL)
  })
  settings_sig <- function(s) {
    paste(s$family, s$zi, s$random_slope, s$standardise, s$among, s$include_invalid, sep = "||")
  }
  model_sig <- reactive({
    s <- model_settings()
    paste(data_sig(), settings_sig(s), s$age_function, paste(sort(s$models), collapse = ","),
          paste(names(s$extra), vapply(s$extra, function(x) paste(sort(x), collapse = "+"), character(1)), sep = "=", collapse = ";"), sep = "||")
  })

  fn_default <- reactiveValues(sig = NULL, best = NULL, n = NA, nonlinear_best = FALSE, manual = NULL)
  output$use_a3_function_ui <- renderUI({
    fn <- input$a3_function %||% "Quadratic"
    tagList(
      actionButton("use_a3_function", paste0("Use this ageing function (", fn, ") for the mixed models"),
                   class = "btn-default btn-block-space", icon = icon("arrow-right"), style = "white-space: normal; width: 100%;"),
      if (identical(isolate(fn_default$manual), fn)) p(class = "small-note", "The Modelling tab uses this function.") else NULL
    )
  })
  observeEvent(input$use_a3_function, disappr_guard("input$use_a3_function", {
    fn <- input$a3_function %||% "Quadratic"
    updateSelectInput(session, "model_age_function", selected = fn)
    fn_default$manual <- fn
    notify(paste0("The mixed models on the Modelling tab will use the ", fn, " ageing function."), duration = 4)
  }))
  observeEvent(list(input$tabs, safe_get(data_sig())), disappr_guard("list(input$tabs, safe_get(data_sig()))", {
    if (!identical(input$tabs, "models")) return(invisible(NULL))
    sig <- safe_get(data_sig())
    if (is.null(sig) || identical(fn_default$sig, sig)) return(invisible(NULL))
    fn_default$sig <- sig
    tab <- safe_get(a3_compare())
    lin <- if (is.data.frame(tab) && nrow(tab)) tab$Function[tab$Function %in% AGE_FUNCTIONS] else character(0)
    if (!length(lin)) {
      fn_default$best <- NULL
      return(invisible(NULL))
    }
    fn_default$best <- lin[[1]]
    fn_default$n <- tab$N_common[[1]]
    fn_default$nonlinear_best <- identical(tab$Function[[1]], A3_NONLINEAR)
    if (!is.null(fn_default$manual)) return(invisible(NULL))
    updateSelectInput(session, "model_age_function", selected = lin[[1]])
  }))
  output$function_default_note <- renderUI({
    b <- fn_default$best
    if (!is.null(fn_default$manual)) {
      return(p(class = "small-note", sprintf("Chosen on tab 4 (individual fits): %s.%s You can still change it here.", fn_default$manual,
                                             if (!is.null(b)) sprintf(" Lowest mean \u0394AICc across individuals: %s.", b) else "")))
    }
    if (is.null(b)) {
      return(p(class = "small-note", "Default: Quadratic. The default follows the individual fits once they can be compared (tab 4)."))
    }
    cur <- input$model_age_function %||% b
    p(class = "small-note",
      sprintf("Default: %s, the function with the lowest mean \u0394AICc across individuals (%s individuals with every function estimable)%s.",
              b, format(fn_default$n, big.mark = ","), if (isTRUE(fn_default$nonlinear_best)) "; the non-linear exponential fitted individuals even better and can be chosen above" else ""),
      " This is a starting point: a simpler function (e.g. Linear or Quadratic) may be preferable, and the population-level comparison (tab 4) can favour a different function.",
      if (!identical(cur, b)) strong(" You have chosen a different function.") else NULL)
  })

  output$zi_ui <- renderUI({
    m <- safe_get(meta())
    covs <- if (is.null(m) || !length(m$covars)) character(0) else stats::setNames(m$covars, unname(m$cov_labels[m$covars]))
    ch <- c("Age" = "f1", "ALR" = "ALR", "AFR (age at first observation)" = "AFR", covs)
    cur <- isolate(input$zi_vars)
    selectizeInput("zi_vars", "Zero inflation depends on", choices = ch, selected = intersect(cur %||% character(0), ch), multiple = TRUE,
                   options = list(placeholder = "constant (~1) if empty"))
  })
  try(outputOptions(output, "zi_ui", suspendWhenHidden = FALSE), silent = TRUE)
  output$random_support_note <- renderUI({
    d <- safe_get(dat())
    if (is.null(d)) return(NULL)
    sup <- individual_data_support(d)
    adv <- random_slope_advice(sup)
    col <- switch(adv$level, good = "#1B7837", limited = "#B35806", "#A50026")
    warn <- data_support_warning(sup)
    tagList(
      div(class = "small-note", style = paste0("border-left: 3px solid ", col, "; padding-left: 7px; margin-bottom: 8px;"), adv$text),
      if (nzchar(warn)) div(class = "small-note", style = "border-left: 3px solid #A50026; padding-left: 7px; margin-bottom: 8px;", strong(warn))
    )
  })
  extra_term_choices <- function(m) {
    covs <- if (is.null(m) || !length(m$covars)) character(0) else {
      labs <- unname(m$cov_labels[m$covars])
      c(stats::setNames(m$covars, paste0(labs, " (additive)")), stats::setNames(paste0(m$covars, ":age"), paste0(labs, " \u00d7 age")))
    }
    c(EXTRA_TERM_KEYS, covs)
  }
  extra_term_label <- function(x, m) {
    ch <- extra_term_choices(m)
    lab <- names(ch)[match(x, ch)]
    ifelse(is.na(lab), x, lab)
  }
  observeEvent(safe_get(meta()$covars), disappr_guard("safe_get(meta()$covars)", {
    ch <- extra_term_choices(safe_get(meta()))
    for (mid in MODEL_IDS) {
      cur <- isolate(input[[paste0("extra_", mid)]]) %||% character(0)
      updateSelectizeInput(session, paste0("extra_", mid), choices = ch, selected = intersect(cur, ch))
    }
  }), ignoreNULL = FALSE)
  lapply(MODEL_IDS, function(mid) {
    observeEvent(input[[paste0("info_model_", mid)]], disappr_guard("input[[paste0('info_model_', mid)]]", {
      s <- model_settings()
      m <- safe_get(meta())
      d <- safe_get(dat())
      z <- tryCatch(model_info_content(mid, s, m), error = function(e) NULL)
      if (is.null(z) || !is.function(session$sendModal)) return(invisible(NULL))
      lv <- list()
      if (!is.null(m) && !is.null(d) && length(m$covars)) {
        cat_cov <- m$covars[vapply(m$covars, function(cv) !is.null(d[[cv]]) && !is.numeric(d[[cv]]), logical(1))]
        for (cv in cat_cov) lv[[cv]] <- sort(unique(as.character(stats::na.omit(d[[cv]]))))
      }
      eqn <- tryCatch(model_equation(mid, s, m, lv, d), error = function(e) list(available = FALSE, message = paste("The equation could not be built:", conditionMessage(e))))
      body <- tagList(
        p(z$meaning),
        tags$table(class = "table table-condensed",
                   tags$tbody(lapply(names(z$attributes), function(k) tags$tr(tags$th(k), tags$td(z$attributes[[k]]))))),
        p(strong("Error family: "), z$family, br(), strong("Ageing function: "), z$ageing,
          if (length(z$extra)) tagList(br(), strong("Terms added to this model: "), paste(extra_term_label(z$extra, m), collapse = ", ")) else NULL),
        h5(strong("Model equation")),
        if (isTRUE(eqn$available)) tagList(
          div(class = "model-eq", lapply(eqn$equation, function(x) div(HTML(x)))),
          tags$ul(class = "small-note", lapply(eqn$distributions, function(x) tags$li(HTML(x)))),
          h5(strong("Coefficients (each term has its own coefficient)")),
          tags$table(class = "table table-condensed coef-eq",
                     tags$thead(tags$tr(tags$th("Coefficient"), tags$th("Multiplies"), tags$th("Name in the fitted R model"))),
                     tags$tbody(lapply(seq_len(nrow(eqn$coefficients)), function(k) {
                       tags$tr(tags$td(HTML(eqn$coefficients$Coefficient[[k]])), tags$td(HTML(eqn$coefficients$Multiplies[[k]])),
                               tags$td(tags$code(eqn$coefficients$R_term[[k]])))
                     }))),
          h5(strong("R syntax")),
          tags$pre(class = "code-out", paste(c(if (!is.null(eqn$basis)) paste("#", eqn$basis),
                                               "# as written in the app", eqn$r_compact,
                                               if (!is.null(eqn$r_expanded)) c("", "# fully expanded (A * B = A + B + A:B)", eqn$r_expanded),
                                               "", "# fitting call", eqn$r_call), collapse = "\n"))
        ) else p(eqn$message %||% "The equation could not be built."),
        p(class = "small-note", "i = individual, j = sampling occasion, k = group. age, ALR (age at last record), AFR (age at first observation) and LS (known lifespan) are centred and scaled when 'Standardise' is ticked. mean(age) is an individual's mean age over its records and \u0394age = age \u2212 mean(age). [x = level] is 1 when a categorical covariate takes that level and 0 otherwise; the first level is the reference. In R, f1 is the (transformed) age term, f2 and f3 its square and cube, and cv_ marks covariates."),
        p(class = "small-note", "Covariates and interactions are set on the Data tab; the other settings on this tab.")
      )
      showModal(modalDialog(title = z$title, body, easyClose = TRUE, footer = modalButton("Close"), size = "l"))
    }), ignoreInit = TRUE)
  })


  output$family_hint <- renderUI({
    fam <- safe_get(integrity()$family)
    s <- model_settings()
    msgs <- list()
    if (!is.null(fam)) {
      if (!identical(fam$family, "gaussian") && identical(s$family, "gaussian")) {
        msgs <- c(msgs, list(div(class = "diagnosis-card", div(class = "diagnosis-detail",
          strong("Integer data treated as Gaussian. "), "Check that the trait is continuous rather than a count, to avoid a misfit family."))))
      }
      if (identical(fam$family, "gaussian") && !identical(s$family, "gaussian")) {
        msgs <- c(msgs, list(div(class = "diagnosis-card", div(class = "diagnosis-detail", "The trait is not a non-negative integer count: use the Gaussian family."))))
      }
    }
    if (!identical(s$family, "gaussian") && !HAS_GLMMTMB) {
      msgs <- c(msgs, list(div(class = "diagnosis-card", div(class = "diagnosis-detail", "Install glmmTMB to fit count families."))))
    }
    tagList(msgs)
  })

  output$structure_note <- renderUI({
    m <- safe_get(meta())
    if (is.null(m)) return(NULL)
    s <- model_settings()
    rstr <- random_display(m, s$random_slope)
    cov <- if (length(m$covars)) {
      paste(vapply(m$covars, function(cv) paste0(m$cov_labels[[cv]], if (cv %in% (m$cov_age %||% character(0))) " (\u00d7 age)" else ""), character(1)), collapse = ", ")
    } else "none"
    div(class = "truth-card small-note",
        strong("Random effects: "), rstr, if (isTRUE(m$has_group) && isTRUE(m$nested)) " (IDs nested in group)" else "", br(),
        strong("Covariates: "), cov, br(),
        strong("Among-individual terms: "), if (identical(s$among, "same")) "same polynomial order as the ageing function" else "linear", br(),
        if (length(s$extra)) tagList(strong("Extra terms: "), paste(paste0(model_label(names(s$extra)), ": ",
                                     vapply(s$extra, function(x) paste(extra_term_label(x, m), collapse = ", "), character(1))), collapse = "; "), br()) else NULL,
        "Change covariates and grouping on the Data tab.")
  })

  output$b2_settings <- renderUI({
    s <- model_settings()
    p(class = "small-note", strong("Current settings: "), family_label(s$family),
      if (s$family %in% c("zip", "zinb")) paste0(", zero-inflation ", s$zi) else "",
      ", ", slope_text(s$random_slope), ".")
  })

  # ---- B2: population-level function comparison ----
  b2_res <- reactiveVal(NULL)
  observeEvent(input$run_functions, disappr_guard("input$run_functions", {
    d <- safe_get(dat())
    if (is.null(d)) {
      notify("The data are not ready: check the Data tab.", type = "error")
      return(invisible(NULL))
    }
    s <- model_settings()
    fns <- intersect(input$b2_functions %||% character(0), MODEL_FUNCTIONS)
    if (A3_NONLINEAR %in% fns && !identical(s$family, "gaussian")) {
      fns <- setdiff(fns, A3_NONLINEAR)
      notify("With a count family the non-linear exponential equals the Linear function on the log scale, so it was left out.", type = "warning", duration = 6)
    }
    if (!length(fns)) {
      notify("Tick at least one ageing function.", type = "warning")
      return(invisible(NULL))
    }
    res <- run_with_progress("Comparing ageing functions", function(pr) {
      compare_population_functions(d, meta(), s$family, s$random_slope, s$standardise, s$zi, progress = pr,
                                   functions = fns, model = "M1", among = s$among, include_invalid = s$include_invalid)
    })
    res$signature <- paste(data_sig(), settings_sig(s))
    b2_res(res)
  }))
  b2_current <- reactive({
    r <- b2_res()
    s <- model_settings()
    validate(need(!is.null(r), "Press 'Compare ageing functions'."))
    validate(need(identical(r$signature, paste(data_sig(), settings_sig(s))),
                  "Data or settings changed: press 'Compare ageing functions' again."))
    r
  })
  # Results for the currently ticked functions (delta AIC recomputed among them)
  b2_shown <- reactive({
    r <- b2_current()
    keep <- intersect(input$b2_functions %||% character(0), r$table$Function)
    validate(need(length(keep) > 0, "Press 'Compare ageing functions' to fit the ticked functions."))
    tab <- r$table[r$table$Function %in% keep, , drop = FALSE]
    if (any(is.finite(tab$AIC))) tab$Delta_AIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
    cv <- if (nrow(r$curves)) r$curves[r$curves$Function %in% keep, , drop = FALSE] else r$curves
    list(table = tab, curves = cv)
  })
  output$b2_table <- renderTable({
    tab <- b2_shown()$table
    tab$AIC <- round(tab$AIC, 1)
    tab$Delta_AIC <- round(tab$Delta_AIC, 1)
    tab
  }, striped = TRUE, spacing = "xs")
  output$b2_plot <- renderPlot(b2_plot_obj())
  b2_plot_obj <- reactive({
    r <- b2_shown()
    validate(need(nrow(r$curves) > 0, "No function could be fitted."))
    obs <- observed_trajectory(dat())
    ggplot() +
      geom_point(data = obs, aes(age, fitted, size = n), shape = 21, fill = "grey55", colour = "black", stroke = 0.6, alpha = 0.85) +
      geom_line(data = r$curves, aes(age, fitted, colour = Function), linewidth = 1.2) +
      scale_colour_manual(values = FUNCTION_COLOURS) +
      scale_size_area(max_size = 4, guide = "none") +
      labs(x = "Age", y = current_map()$trait, subtitle = "Lines: Model 1 predictions for each function; outlined points: observed means") +
      theme_disappR(13)
  })

  # ---- B4: comparative models ----
  model_res <- reactiveVal(NULL)
  observeEvent(input$fit_models, disappr_guard("input$fit_models", {
    d <- safe_get(dat())
    if (is.null(d)) {
      notify("The data are not ready: check the Data tab.", type = "error")
      return(invisible(NULL))
    }
    s <- model_settings()
    if (!length(s$models)) {
      notify("Select at least one model.", type = "warning")
      return(invisible(NULL))
    }
    res <- run_with_progress("Fitting models", function(pr) {
      fit_model_suite(d, meta(), s$models, s$age_function, s$family, s$random_slope, s$standardise, s$zi, progress = pr,
                      among = s$among, include_invalid = s$include_invalid, extra = s$extra)
    })
    res$signature <- model_sig()
    model_res(res)
    dharma_res(NULL)
    if (!isTRUE(res$ok)) notify(res$message, type = "error", duration = 8)
  }))

  current_models <- reactive({
    r <- model_res()
    validate(need(!is.null(r), "Choose settings and press 'Fit models'."))
    validate(need(identical(r$signature, model_sig()), "Data or model settings changed since the last fit: press 'Fit models' to update."))
    r
  })
  fitted_models <- reactive({
    r <- current_models()
    validate(need(isTRUE(r$ok), r$message))
    r
  })

  output$model_fit_note <- renderUI({
    r <- current_models()
    if (!isTRUE(r$ok)) {
      return(div(class = "diagnosis-card", div(class = "diagnosis-title", "Models not fitted"), div(class = "diagnosis-detail", r$message)))
    }
    aic <- r$aic
    best <- aic$Model[[1]]
    best_id <- sub("^Model ", "M", best)
    elig <- aic[aic$Eligible, , drop = FALSE]
    d45 <- if (all(c("Model 4", "Model 5") %in% elig$Model)) {
      elig$AIC[elig$Model == "Model 5"] - elig$AIC[elig$Model == "Model 4"]
    } else NA_real_
    warn <- sum(unlist(r$validity) != "Valid")
    div(class = "diagnosis-card",
        div(class = "diagnosis-title", paste0("Lowest AIC among eligible fits: ", best, " (", aic$Fit[[1]], ")",
                                              if (nrow(elig) > 1) sprintf("; next model \u0394AIC = %.1f", elig$Delta_AIC[[2]]) else "")),
        div(class = "diagnosis-detail", paste0(family_label(r$family), "; ", r$age_function, " ageing; random effects ",
                                               if (isTRUE(r$nonlinear)) r$random else random_display(meta(), r$random_slope),
                                               "; N = ", aic$N[[1]], " observations from ", length(unique(r$data$id)), " individuals.")),
        div(class = "diagnosis-detail small-note", if (isTRUE(r$nonlinear)) paste0("Fitted model: ", r$formulas[[best_id]]) else
          paste0("Fitted formula of this model: trait ~ ", display_term(r$formulas[[best_id]]), " + ", display_term(r$random))),
        if (nzchar(r$random_note %||% "")) div(class = "diagnosis-detail small-note", r$random_note) else NULL,
        if (isTRUE(r$n_excluded > 0)) div(class = "diagnosis-detail", strong(sprintf("%d fit(s) with an invalid Hessian are excluded from the ranking and from automated interpretation.", r$n_excluded))) else NULL,
        if (is.finite(d45)) div(class = "diagnosis-detail", sprintf("AIC(Model 5) \u2212 AIC(Model 4) = %.1f: %s", d45,
          if (d45 > 2) "ALR interaction better supported than mean-age centring." else if (d45 < -2) "mean-age centring better supported than the ALR interaction." else "Models 4 and 5 are similarly supported.")) else NULL,
        if (r$n_dropped > 0) div(class = "diagnosis-detail", paste0(r$n_dropped, " rows with missing values in model variables were dropped from all models (", paste(r$drop_by, collapse = "; "), ").")) else NULL,
        if (warn > 0) div(class = "diagnosis-detail", strong(paste0(warn, " model(s) are classified Caution or Failed: see the fitting status below."))) else NULL)
  })

  output$aic_plot <- renderPlot({
    all_aic <- fitted_models()$aic
    aic <- all_aic[all_aic$Eligible, , drop = FALSE]
    validate(need(nrow(aic) > 0, "No eligible AIC values."))
    aic$Model <- factor(aic$Model, levels = rev(aic$Model))
    aic$lab <- paste0(sprintf("%.1f", aic$Delta_AIC), ifelse(aic$Fit == "Valid", "", paste0(" (", tolower(aic$Fit), ")")))
    excl <- all_aic$Model[!all_aic$Eligible]
    ggplot(aic, aes(Model, Delta_AIC, fill = Model, alpha = Fit)) + geom_col(width = 0.7) + coord_flip() +
      geom_hline(yintercept = 2, linetype = 2, colour = "grey45") +
      geom_text(aes(label = lab), hjust = -0.1, size = 3.5, colour = "#453A32", alpha = 1) +
      scale_fill_manual(values = METHOD_COLOURS, guide = "none") +
      scale_alpha_manual(values = c(Valid = 1, Caution = 0.4, Failed = 0.2), guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(x = NULL, y = "\u0394AIC", subtitle = paste0("Dashed line: \u0394AIC = 2; faded bars: Caution fits",
                                                       if (length(excl)) paste0("; excluded (invalid Hessian): ", paste(excl, collapse = ", ")) else "")) +
      theme_disappR(12)
  })

  output$aic_table <- renderTable({
    aic <- fitted_models()$aic
    data.frame(Model = aic$Model, Fit = ifelse(aic$Eligible, aic$Fit, paste(aic$Fit, "(excluded)")), AIC = round(aic$AIC, 1),
               dAIC = round(aic$Delta_AIC, 1), Weight = round(aic$Weight, 3), df = aic$df, logLik = round(aic$logLik, 1), N = aic$N)
  }, striped = TRUE, spacing = "xs")

  output$lrt_table <- renderTable({
    l <- fitted_models()$lrt
    if (!nrow(l)) return(data.frame(Note = "Fit nested pairs (e.g. Models 1, 2 and 4) to obtain likelihood-ratio tests."))
    l$Chisq <- round(l$Chisq, 2)
    l$P_value <- vapply(l$P_value, format_p, character(1))
    l
  }, striped = TRUE, spacing = "xs")

  output$status_table <- renderTable({
    st <- current_models()$status
    if (!length(st)) return(NULL)
    data.frame(Model = model_label(names(st)), Status = unlist(st, use.names = FALSE))
  }, striped = TRUE, spacing = "xs")

  varcomp_display <- reactive({
    r <- fitted_models()
    vc <- r$varcomp
    validate(need(is.data.frame(vc) && nrow(vc) > 0, "No random-effect variances available."))
    m <- meta()
    grp_lab <- function(g) {
      out <- g
      out[g == "id"] <- if (isTRUE(m$has_group) && isTRUE(m$nested)) paste0(m$map$group, ":", m$map$id) else m$map$id
      out[g == "group"] <- m$map$group %||% "group"
      rt <- g %in% names(m$random_labels)
      out[rt] <- unname(m$random_labels[g[rt]])
      out
    }
    term_lab <- function(t) gsub("f1", "age (first ageing term)", t, fixed = TRUE)
    data.frame(Model = vc$Model, Group = grp_lab(vc$Group), Term = term_lab(vc$Term), Type = vc$Type,
               Variance = signif(vc$Variance %||% rep(NA_real_, nrow(vc)), 4),
               `SD / correlation / dispersion` = signif(vc$Estimate, 4), check.names = FALSE, stringsAsFactors = FALSE)
  })
  output$varcomp_table <- renderTable(varcomp_display(), striped = TRUE, spacing = "xs")

  output$drop_note <- renderUI({
    r <- current_models()
    if (is.null(r$n_dropped) || r$n_dropped == 0) return(NULL)
    p(class = "small-note", paste0("Rows dropped for the common-data comparison: ", paste(r$drop_by, collapse = "; "), "."))
  })

  eligible_models <- reactive({
    r <- fitted_models()
    intersect(names(r$fits), sub("^Model ", "M", r$aic$Model[r$aic$Eligible]))
  })
  pred_curves_for <- function(r, models, by = NULL) {
    ages <- prediction_ages(r$data$age)
    lst <- lapply(models, function(m) {
      cv <- predict_population_curve(r$fits[[m]], r, ages, by = by)
      if (!nrow(cv)) return(NULL)
      cv$Method <- model_label(m)
      cv$Fit <- r$validity[[m]] %||% "Valid"
      cv
    })
    lst <- Filter(Negate(is.null), lst)
    if (length(lst)) do.call(rbind, lst) else data.frame(age = numeric(0), fitted = numeric(0), Method = character(0), Fit = character(0))
  }
  pred_curves <- reactive(pred_curves_for(fitted_models(), eligible_models()))
  # NULL = the tick boxes have not been used yet (draw every eligible model); character(0) = the user unticked all
  pred_sel <- reactiveVal(NULL)
  observeEvent(input$pred_models, disappr_guard("pred_models", pred_sel(as.character(input$pred_models %||% character(0)))),
               ignoreNULL = FALSE, ignoreInit = TRUE)
  pred_curves_shown <- reactive({
    r <- fitted_models()
    el <- eligible_models()
    chosen <- pred_sel()
    sel <- if (is.null(chosen)) el else intersect(chosen, el)
    by <- input$pred_by %||% ""
    pred_curves_for(r, sel, if (nzchar(by)) by else NULL)
  })
  output$pred_models_ui <- renderUI({
    r <- safe_get(fitted_models())
    if (is.null(r)) return(NULL)
    ch <- safe_get(eligible_models()) %||% character(0)
    if (!length(ch)) return(NULL)
    chosen <- isolate(pred_sel())
    sel <- if (is.null(chosen)) ch else intersect(chosen, ch)
    if (!is.null(chosen) && length(chosen) && !length(sel)) sel <- ch   # after refitting with other models
    checkboxGroupInput("pred_models", "Models to draw", choices = stats::setNames(ch, model_label(ch)), selected = sel, inline = TRUE)
  })
  output$pred_by_ui <- renderUI({
    r <- safe_get(fitted_models())
    if (is.null(r)) return(NULL)
    cv <- intersect(r$covars %||% character(0), names(r$data))
    fac <- cv[!vapply(r$data[cv], is.numeric, logical(1))]
    if (!length(fac)) return(NULL)
    m <- meta()
    cur <- isolate(input$pred_by) %||% ""
    selectInput("pred_by", "Show predictions by", choices = c("(pooled)" = "", stats::setNames(fac, unname(m$cov_labels[fac]))),
                selected = if (cur %in% fac) cur else "")
  })
  obs_traj <- reactive(observed_trajectory(dat()))
  decomp_traj <- reactive(decomposition_trajectory(dat()))

  output$pred_plot <- renderPlot(pred_plot_obj())
  pred_plot_obj <- reactive({
    pc <- pred_curves_shown()
    validate(need(nrow(pc) > 0, "Choose at least one model to draw (or predictions could not be computed for the chosen models)."))
    by_on <- "level" %in% names(pc)
    p <- ggplot()
    if (isTRUE(input$show_observed)) {
      ob <- if (by_on) {
        byv <- input$pred_by
        dd <- dat()
        parts <- lapply(unique(pc$level), function(lv) {
          o <- observed_trajectory(dd[!is.na(dd[[byv]]) & as.character(dd[[byv]]) == lv, , drop = FALSE])
          if (!nrow(o)) return(NULL)
          o$level <- lv
          o
        })
        parts <- Filter(Negate(is.null), parts)
        if (length(parts)) do.call(rbind, parts) else NULL
      } else obs_traj()
      if (!is.null(ob) && nrow(ob)) {
        ob$Method <- "Observed"
        p <- p + geom_point(data = ob, aes(age, fitted, size = n), shape = 21, fill = "grey55", colour = "black", stroke = 0.6, alpha = 0.85)
      }
    }
    p <- p + geom_line(data = pc, aes(age, fitted, colour = Method, linetype = Fit), linewidth = 1.1)
    if (isTRUE(input$show_a3) && !by_on) {
      z <- safe_get(a3_fit())
      if (!is.null(z) && nrow(z$mean_curve)) {
        mc <- z$mean_curve
        a3 <- data.frame(age = mc$age, fitted = mc$fitted, Method = "Individual fits: mean of coefficients")
        if (isTRUE(z$nonlinear)) a3 <- rbind(a3, data.frame(age = mc$age, fitted = mc$fitted_fun, Method = "Individual fits: mean of functions"))
        a3 <- a3[is.finite(a3$fitted), , drop = FALSE]
        if (nrow(a3)) p <- p + geom_line(data = a3, aes(age, fitted, colour = Method), linewidth = 1.2, linetype = "dotdash")
      }
    }
    if (isTRUE(input$show_decomp) && !by_on) {
      de <- decomp_traj()
      if (nrow(de)) {
        de$Method <- "Decomposition"
        p <- p + geom_line(data = de, aes(age, fitted, colour = Method), linewidth = 1.1, linetype = 2)
      }
    }
    tc <- truth()
    if (!is.null(tc)) {
      ag <- sort(unique(c(pc$age)))
      tr <- data.frame(age = ag, fitted = toy_true_curve(tc, ag), Method = "True (simulated)")
      p <- p + geom_line(data = tr, aes(age, fitted, colour = Method), linewidth = 1.6)
    }
    p + scale_colour_manual(values = METHOD_COLOURS) + scale_size_area(max_size = 4, guide = "none") +
      scale_linetype_manual(values = c(Valid = "solid", Caution = "dashed", Failed = "dotted"), guide = "none") +
      (if (by_on) facet_wrap(~ level) else NULL) +
      labs(x = "Age", y = current_map()$trait, colour = NULL,
           subtitle = paste0(if (!is.null(tc)) "Black: simulated typical-individual trajectory; " else "",
                             if (!by_on) "dashed purple: decomposition; dot-dash: reconstruction from individual fits (if shown); " else "panels: levels of the chosen covariate; ",
                             "points: observed means; dashed model lines: Caution fits")) +
      theme_disappR(13) + guides(colour = guide_legend(nrow = 2))
  })

  output$deviation_note <- renderUI({
    if (!is_toy()) {
      return(p(class = "small-note", "The true trajectory is unknown for empirical data. Use a simulated dataset to see how each model and the decomposition recover it."))
    }
    p(class = "small-note", "Relativised deviation D = 100 \u00d7 (estimate \u2212 truth) / truth (Methods Eq. 11), over ages with \u2265 10 observed individuals.",
      if (grepl("^count", toy_cfg()$trait)) " For counts the truth is the typical-individual (fixed-effect) curve, so observed means differ from it even without selection." else "")
  })

  output$deviation_table <- renderTable(deviation_tab(), striped = TRUE, spacing = "xs")
  deviation_tab <- reactive({
    dv <- deviation_raw()
    if (is.null(dv) || !nrow(dv)) return(NULL)
    dv$Mean_abs_D <- round(dv$Mean_abs_D, 1)
    dv$Max_abs_D <- round(dv$Max_abs_D, 1)
    dv$D_at_oldest_age <- round(dv$D_at_oldest_age, 1)
    names(dv) <- c("Method", "Mean |D| %", "Max |D| %", "D at oldest %", "Oldest age")
    dv
  })
  deviation_raw <- reactive({
    validate(need(is_toy(), ""))
    pc <- pred_curves()
    tc <- truth()
    ob <- obs_traj()
    parts <- list(pc[, c("age", "fitted", "Method"), drop = FALSE])
    ob2 <- ob[, c("age", "fitted")]
    ob2$Method <- "Observed"
    parts <- c(parts, list(ob2))
    de <- decomp_traj()
    if (nrow(de)) {
      de2 <- de[, c("age", "fitted")]
      de2$Method <- "Decomposition"
      parts <- c(parts, list(de2))
    }
    allc <- do.call(rbind, parts)
    deviation_from_truth(allc[, c("age", "fitted", "Method")], tc, data.frame(age = ob$age, n = ob$n))
  })
  output$consistency_ui <- renderUI({
    r <- safe_get(fitted_models())
    if (is.null(r)) return(NULL)
    dev <- if (isTRUE(safe_get(is_toy()))) safe_get(deviation_raw()) else NULL
    tr <- safe_get(a2_slope_trend(a2_slopes()))
    notes <- tryCatch(consistency_notes(r, dev, if (is.data.frame(tr) && nrow(tr) == 1) tr else NULL), error = function(e) character(0))
    tagList(
      h5(info_title(strong("Internal consistency"), "consistency")),
      if (length(notes)) {
        div(class = "diagnosis-card", lapply(notes, function(x) div(class = "diagnosis-detail", icon("exclamation-triangle"), " ", x)))
      } else {
        p(class = "small-note", "No inconsistencies flagged between AIC, likelihood-ratio tests, fit validity, random-slope support",
          if (isTRUE(safe_get(is_toy()))) ", recovery of the simulated truth" else "", " and the trend of the lifespan\u2013trait coefficient across age bins. This does not prove that the evidence is coherent.")
      })
  })

  output$coef_model_ui <- renderUI({
    r <- safe_get(fitted_models())
    if (is.null(r)) return(NULL)
    ch <- names(r$fits)
    cur <- isolate(input$coef_model)
    selectInput("coef_model", NULL, choices = stats::setNames(ch, model_label(ch)),
                selected = if (isTRUE(cur %in% ch)) cur else ch[[length(ch)]])
  })

  output$coef_table <- renderTable({
    r <- fitted_models()
    req(input$coef_model)
    x <- r$coefficients[r$coefficients$Model == model_label(input$coef_model), , drop = FALSE]
    if (!nrow(x)) return(data.frame(Note = "No coefficient table."))
    coef_display(x, r)
  }, striped = TRUE, spacing = "xs")
  coef_display <- function(x, r) {
    sc <- tryCatch(term_scaling(x$Raw_term, r), error = function(e) data.frame(Scale_factor = rep(NA_real_, nrow(x)), Per = ""))
    data.frame(Term = x$Term, Estimate = signif(x$Estimate, 4), SE = signif(x$SE, 3),
               Statistic = round(x$Statistic, 2), P = vapply(x$P_value, format_p, character(1)),
               `Scale factor` = signif(sc$Scale_factor, 4), `Estimate per original unit` = signif(x$Estimate / sc$Scale_factor, 4),
               `SE per original unit` = signif(x$SE / sc$Scale_factor, 3), Per = sc$Per, check.names = FALSE, stringsAsFactors = FALSE)
  }
  output$scaling_table <- renderTable({
    r <- fitted_models()
    sc <- scaling_constants(r)
    data.frame(Variable = sc$Variable, Centre = signif(sc$Centre, 5), `Scale (SD)` = signif(sc$Scale, 5), check.names = FALSE)
  }, striped = TRUE, spacing = "xs")
  coef_interpretation_text <- reactive({
    r <- fitted_models()
    req(input$coef_model)
    tryCatch(interpret_model_terms(r, input$coef_model), error = function(e) paste("Interpretation not available:", conditionMessage(e)))
  })
  output$coef_interpretation <- renderUI({
    txt <- safe_get(coef_interpretation_text())
    if (is.null(txt) || !length(txt)) return(NULL)
    div(class = "small-note", lapply(txt, p))
  })

  coef_re_tab <- reactive({
    r <- fitted_models()
    req(input$coef_model)
    vd <- safe_get(varcomp_display())
    if (is.null(vd)) return(data.frame(Note = "No random-effect estimates."))
    x <- vd[vd$Model == model_label(input$coef_model), setdiff(names(vd), "Model"), drop = FALSE]
    if (!nrow(x)) return(data.frame(Note = "No random-effect estimates."))
    rownames(x) <- NULL
    x
  })
  output$coef_re_table <- renderTable(coef_re_tab(), striped = TRUE, spacing = "xs")

  output$definition_table <- renderTable({
    s <- model_settings()
    m <- safe_get(meta())
    covs <- if (is.null(m)) character(0) else m$covars
    labs <- if (is.null(m)) character(0) else m$cov_labels
    rstr <- if (is.null(m)) "(1 | ID)" else random_display(m, s$random_slope)
    defs <- model_definition_table(s$age_function, covs, rstr, labs, among = s$among,
                                   cov_age = if (is.null(m)) character(0) else m$cov_age, cov_pairs = if (is.null(m)) character(0) else m$cov_pairs)
    defs[defs$Model %in% model_label(s$models), c("Model", "Name", "Fixed_effects", "Question"), drop = FALSE]
  }, striped = TRUE, spacing = "xs")

  # ---- B1: family check ----
  family_res <- reactiveVal(NULL)
  observeEvent(input$run_family_check, disappr_guard("input$run_family_check", {
    d <- safe_get(dat())
    if (is.null(d)) return(invisible(NULL))
    s <- model_settings()
    res <- run_with_progress("Comparing count families", function(pr) {
      fams <- if ((s$family %||% "gaussian") %in% BINOMIAL_FAMILIES) BINOMIAL_FAMILIES else COUNT_FAMILIES
      compare_families(d, meta(), model = input$family_check_model %||% "M4", age_function = s$age_function,
                       random_slope = s$random_slope, standardise = s$standardise, zi_str = s$zi, progress = pr,
                       among = s$among, include_invalid = s$include_invalid, families = fams)
    })
    res$signature <- paste(data_sig(), input$family_check_model, s$age_function, settings_sig(s))
    family_res(res)
  }))
  output$family_table <- renderTable({
    r <- family_res()
    validate(need(!is.null(r), "Press 'Compare count families'."))
    validate(need(isTRUE(r$ok), r$message))
    tab <- r$table
    tab$AIC <- round(tab$AIC, 1)
    tab$Delta_AIC <- round(tab$Delta_AIC, 1)
    tab
  }, striped = TRUE, spacing = "xs")


  # ======================================================================
  # 2b. A4-A7 disappearance diagnostics
  # ======================================================================
  a4_ia <- reactive(disappearance_data(dat(), meta(), FALSE))
  fmt_ci <- function(lo, hi) {
    ifelse(is.finite(lo) & is.finite(hi), paste0(vapply(lo, format_num, character(1)), " to ", vapply(hi, format_num, character(1))), "NA")
  }
  output$a4_status <- renderUI({
    ia <- safe_get(a4_ia())
    if (is.null(ia)) return(p(class = "small-note", "Too few records for the disappearance diagnostics."))
    known <- isTRUE(attr(ia, "lifespan_known"))
    m <- meta()
    div(class = "small-note",
        sprintf("%s records with a known outcome and %s disappearances. Disappearance = %s. %s",
                format(sum(is.finite(ia$event)), big.mark = ","), format(sum(ia$event == 1, na.rm = TRUE), big.mark = ","),
                if (known) "death before the next occasion (known lifespan)" else "last record (lifespan unknown: includes emigration and missed detections)",
                if (isTRUE(m$has_censor)) sprintf("%d individuals are censored and their final interval is excluded.", m$n_censored)
                else if (!known) "No censoring column is mapped, so individuals alive at the end of the study count as disappearing (except at the oldest sampled age)." else ""))
  })
  a5_data <- reactive(terminal_data(dat(), meta(), 3))
  a5_plot_obj <- reactive({
    z <- a5_data()
    validate(need(is.data.frame(z) && nrow(z) > 0, "Too few individuals with a known end (death or last record) at each occasion before death."))
    known <- isTRUE(attr(z, "lifespan_known"))
    u <- sort(unique(z$before))
    ggplot(z, aes(before, trait, colour = end_bin, group = end_bin)) +
      geom_line(linewidth = 1.1) + geom_point(aes(size = n)) +
      scale_x_reverse(breaks = if (length(u) <= 15) u else waiver()) +
      scale_colour_manual(values = ordered_colours(nlevels(factor(z$end_bin)))) + scale_size_area(max_size = 4.5, guide = "none") +
      labs(x = if (known) "Occasions before death (0 = last occasion before death)" else "Occasions before the last record (0 = last record)",
           y = paste("Mean", trait_label()),
           colour = if (known) "Lifespan" else "ALR",
           title = "Trait before death", subtitle = "Individuals aligned on time before death, by lifespan tercile; point size = individuals") +
      theme_disappR(13)
  })
  output$a5_plot <- renderPlot(a5_plot_obj())
  a5_lines <- reactive({
    z <- a5_data()
    if (!is.data.frame(z) || !nrow(z)) return(character(0))
    parts <- lapply(split(z, z$end_bin), function(x) {
      last <- x$resid[x$before == 0]
      earlier <- x$resid[x$before >= 2]
      if (!length(last) || !length(earlier)) return(NULL)
      dd <- last[[1]] - mean(earlier)
      sprintf("%s: %s%s", as.character(x$end_bin[[1]]), if (dd >= 0) "+" else "\u2212", format_num(abs(dd)))
    })
    parts <- unlist(Filter(Negate(is.null), parts))
    if (!length(parts)) return("Too few occasions before death to compare the final occasion with earlier ones.")
    paste0("Change at the last occasion relative to occasions two or more before death, using values relative to the mean at the same age (so the normal age trend is removed): ", paste(parts, collapse = "; "),
           ". Similar changes in every lifespan group suggest terminal effects of similar strength at all ages; clearly different changes suggest age-dependent selective disappearance.",
           if (!isTRUE(attr(z, "lifespan_known"))) " Lifespan is unknown, so the last record may not be close to death." else "")
  })
  output$a5_note <- renderUI({
    txt <- safe_get(a5_lines())
    if (is.null(txt) || !length(txt)) return(NULL)
    p(class = "small-note", txt)
  })

  a6_data <- reactive(selection_differentials(a4_ia()))
  a6_values <- reactive({
    s <- a6_data()
    contrast <- identical(input$a6_compare, "contrast")
    if (nrow(s)) {
      s$y <- if (contrast) s$Contrast else s$Differential
      s$se <- if (contrast) s$Contrast_SE else s$Differential_SE
    }
    s
  })
  a6_trend <- reactive({
    s <- a6_values()
    if (nrow(s) < 3) return(NULL)
    ok <- is.finite(s$y) & is.finite(s$se) & s$se > 0
    if (sum(ok) < 3) return(NULL)
    x <- s[ok, , drop = FALSE]
    fit <- stats::lm(y ~ Age, data = x, weights = 1 / x$se^2)
    cf <- suppressWarnings(summary(fit)$coefficients)
    if (!"Age" %in% rownames(cf)) return(NULL)
    list(slope = cf["Age", 1], p = cf["Age", 4], mean = stats::weighted.mean(x$y, 1 / x$se^2), n = nrow(x))
  })
  a6_plot_obj <- reactive({
    s <- a6_values()
    validate(need(nrow(s) > 0, "No age with at least three survivors and three disappearing individuals."))
    contrast <- identical(input$a6_compare, "contrast")
    s$w <- 1 / pmax(s$se, 1e-9)^2
    ggplot(s, aes(Age, y)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_errorbar(aes(ymin = y - 1.96 * se, ymax = y + 1.96 * se), width = 0, colour = "#0072B2", alpha = 0.7) +
      geom_point(aes(size = N), colour = "#0072B2") +
      geom_smooth(aes(weight = w), method = "lm", formula = y ~ x, se = FALSE, colour = "#D55E00", linewidth = 1, na.rm = TRUE) +
      scale_size_area(max_size = 5, guide = "none") +
      labs(x = "Age", y = if (contrast) "Survivors \u2212 disappearing (within-age SD)" else "Selection differential (within-age SD)",
           title = "Selection differentials by age", subtitle = "Points \u00b1 95% intervals; orange line: inverse-variance weighted trend across age") +
      theme_disappR(13)
  })
  output$a6_plot <- renderPlot(a6_plot_obj())
  a6_lines <- reactive({
    tr <- a6_trend()
    if (is.null(tr)) return("Fewer than three ages with enough survivors and disappearing individuals to describe a trend.")
    paste0(sprintf("Weighted mean %s: %s within-age SD; trend across age %s per unit age (p = %s). ",
                   if (identical(input$a6_compare, "contrast")) "survivor\u2013disappearing difference" else "selection differential",
                   format_num(tr$mean), format_num(tr$slope), format_p(tr$p)),
           if (is.finite(tr$p) && tr$p < 0.05) "The differential changes with age, consistent with age-dependent selective disappearance."
           else "No clear trend across age, consistent with age-independent selective disappearance (or too little information).")
  })
  output$a6_note <- renderUI({
    txt <- safe_get(a6_lines())
    if (is.null(txt)) return(NULL)
    p(class = "small-note", txt)
  })

  a7_data <- reactive(life_table(dat(), meta()))
  a7_plot_obj <- reactive({
    z <- a7_data()
    validate(need(isTRUE(z$ok), z$message %||% "No hazard could be estimated."))
    lt <- z$table
    ggplot(lt, aes(Age, Hazard)) +
      geom_hline(yintercept = z$overall, linetype = 2, colour = "#D55E00", linewidth = 1) +
      geom_errorbar(aes(ymin = Lower_95, ymax = Upper_95), width = 0, colour = "#0072B2", alpha = 0.7) +
      geom_line(colour = "#0072B2", alpha = 0.6) + geom_point(aes(size = At_risk), colour = "#0072B2") +
      scale_size_area(max_size = 5, guide = "none") +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
      labs(x = "Age", y = "Probability of disappearing before the next occasion",
           title = "Disappearance hazard", subtitle = "Points: age-specific hazard \u00b1 95% Wilson interval (size = at risk); dashed: overall, age-independent hazard") +
      theme_disappR(13)
  })
  output$a7_plot <- renderPlot(a7_plot_obj())
  a7_lines <- reactive({
    z <- a7_data()
    if (!isTRUE(z$ok)) return(z$message %||% "")
    c(sprintf("Overall (age-independent) hazard: %.1f%% per occasion (%s disappearances in %s individual-occasions at risk).",
              100 * z$overall, format(sum(z$table$Disappearances), big.mark = ","), format(sum(z$table$At_risk), big.mark = ",")),
      sprintf("Test of a constant hazard against age-specific hazards: p = %s%s", format_p(z$p_constant),
              if (is.finite(z$p_constant) && z$p_constant < 0.05) ": the hazard changes with age." else ": no clear change with age."),
      if (!isTRUE(z$lifespan_known)) "Lifespan is unknown: disappearance is the last record, which includes emigration and missed detections." else NULL)
  })
  output$a7_note <- renderUI({
    txt <- safe_get(a7_lines())
    if (is.null(txt)) return(NULL)
    div(class = "small-note", lapply(txt, p))
  })
  # ---- sampling caveat and A3 data support
  output$sampling_caveat <- renderUI({
    sch <- safe_get(integrity()$schedule)
    if (is.null(sch)) return(NULL)
    if (isTRUE(sch$irregular)) {
      div(class = "diagnosis-card", style = "border-left: 5px solid #A50026;",
          div(class = "diagnosis-title", icon("exclamation-triangle"), " Sampling is irregular: read the missingness grid with caution"),
          div(class = "diagnosis-detail", sprintf("%.0f%% of records are not on a common sampling schedule (step %s). The missingness grid assumes regular, scheduled occasions and is most reliable under regular sampling; here expected occasions are approximate, so missingness percentages, patterns and drivers are indicative only. Round ages to a sampling resolution on the Data tab for a more reliable grid.",
                                                  100 * sch$off_share, format_num(sch$step))))
    } else if (isTRUE(sch$n_schedules > 1)) {
      p(class = "small-note", sprintf("Individuals follow %d offset sampling schedules; expected occasions are anchored on each individual's first record. The grid is most reliable under regular sampling.", sch$n_schedules))
    } else {
      p(class = "small-note", "The missingness grid assumes regular, scheduled sampling occasions, which these data appear to follow.")
    }
  })
  output$a3_support <- renderUI({
    sup <- safe_get(individual_data_support(dat()))
    if (is.null(sup)) return(NULL)
    div(class = "truth-card small-note",
        strong("Data support for individual slopes"), br(),
        sprintf("Individuals: %s; trait observations: %s", format(sup$individuals, big.mark = ","), format(sup$observations, big.mark = ",")), br(),
        sprintf("Observations per individual: median %s (IQR %s\u2013%s)", format_num(sup$median_obs), format_num(sup$iqr_low), format_num(sup$iqr_high)), br(),
        sprintf("Time points: \u2265 2 in %.0f%%, \u2265 3 in %.0f%%, \u2265 4 in %.0f%% of individuals", sup$pct_2, sup$pct_3, sup$pct_4), br(),
        sprintf("Distinct ages: %d; median age span per individual: %s", sup$distinct_ages, format_num(sup$median_span)))
  })

  # ---- performance diagnostics
  perf_res <- reactiveVal(NULL)
  output$performance_model_ui <- renderUI({
    r <- safe_get(fitted_models())
    if (is.null(r)) return(p(class = "small-note", "Fit models first."))
    ch <- names(r$fits)
    selectInput("performance_model", "Model", choices = stats::setNames(ch, model_label(ch)), selected = sub("^Model ", "M", r$aic$Model[[1]]))
  })
  observeEvent(input$run_performance, disappr_guard("input$run_performance", {
    r <- safe_get(fitted_models())
    m <- input$performance_model
    if (is.null(r) || is.null(m) || is.null(r$fits[[m]])) {
      notify("Fit models first.", type = "warning")
      return(invisible(NULL))
    }
    res <- tryCatch(run_with_progress("Running performance checks", function(pr) run_performance(r$fits[[m]], r$family)),
                    error = function(e) list(ok = FALSE, message = paste("The performance checks failed for this model:", conditionMessage(e))))
    res$model <- m
    res$signature <- model_sig()
    perf_res(res)
    if (!isTRUE(res$ok)) notify(res$message, type = "error", duration = 8)
  }))
  output$performance_table <- renderTable({
    z <- perf_res()
    validate(need(!is.null(z), "Press 'Run performance checks'."))
    validate(need(isTRUE(z$ok), z$message %||% "The performance checks failed."))
    validate(need(identical(z$signature, model_sig()), "Models were refitted: run the checks again."))
    cbind(Model = model_label(z$model), z$table, stringsAsFactors = FALSE)
  }, striped = TRUE, spacing = "xs")
  output$performance_note <- renderUI({
    z <- perf_res()
    if (is.null(z) || !isTRUE(z$ok) || !nzchar(z$plot_message %||% "")) return(NULL)
    p(class = "small-note", z$plot_message)
  })
  draw_check_model <- function(cm) {
    function() {
      tryCatch(with_time_limit(print(plot(cm)), 90), error = function(e) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste("check_model() panels could not be drawn:", conditionMessage(e)), cex = 0.9)
      })
      invisible(NULL)
    }
  }
  output$performance_plot <- renderPlot({
    z <- perf_res()
    validate(need(!is.null(z) && isTRUE(z$ok) && !is.null(z$check_model), ""))
    validate(need(identical(z$signature, model_sig()), ""))
    draw_check_model(z$check_model)()
  })


  # ======================================================================
  # R code for each section: current settings + the app's own functions
  # ======================================================================
  dput_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "\n")
  script_header <- function(entries, packages = "ggplot2") {
    defs <- tryCatch(app_code_definitions(app_code_closure(entries)), error = function(e) paste("# helper functions could not be collected:", conditionMessage(e)))
    c("# disappR: R code for this section, generated with the current settings",
      paste0("# Created ", format(Sys.time(), "%Y-%m-%d %H:%M")),
      "# The helper functions are copied from disappR, so this script runs on its own.",
      paste0("library(", packages, ")"),
      "",
      "# ==== Helper functions from disappR (no need to edit) ====",
      defs,
      "# ==== End of helper functions ====",
      "")
  }
  script_data_lines <- function() {
    src <- input$data_source %||% "toy"
    lines <- "# ---- Data ----"
    if (identical(src, "upload")) {
      lines <- c(lines, paste0("raw <- read_user_csv(", dput_text(input$data_file$name %||% "your_data.csv"), ")$data  # the uploaded file, placed in the working directory"))
    } else if (identical(src, "example")) {
      lines <- c(lines, sprintf("raw <- utils::read.csv(system.file(\"app\", \"data\", %s, package = \"disappR\"), stringsAsFactors = FALSE, check.names = FALSE)",
                                dput_text(current_example()$file)))
    } else {
      lines <- c(lines, paste0("raw <- simulate_toy_data(", dput_text(toy_cfg()), ")"))
    }
    sv <- input$subset_var %||% ""
    if (nzchar(sv) && length(input$subset_levels)) {
      lines <- c(lines, sprintf("raw <- raw[!is.na(raw[[%s]]) & as.character(raw[[%s]]) %%in%% %s, , drop = FALSE]  # subset used in the app",
                                dput_text(sv), dput_text(sv), dput_text(as.character(input$subset_levels))))
    }
    c(lines,
      paste0("map <- ", dput_text(current_map())),
      paste0("prepared <- standardise_data(raw, map, ", dput_text(input$dup_action %||% "keep"), ")"),
      "dat <- prepared$data   # one row per individual x age",
      "meta <- prepared$meta",
      "im <- individual_metrics(dat)   # one row per individual: ALR, AFR, mean age, lifespan",
      "")
  }
  data_entries <- function(...) c("read_user_csv", "standardise_data", "individual_metrics", if (identical(input$data_source %||% "toy", "toy")) "simulate_toy_data", ...)

  section_code <- function(section) {
    lines <- switch(section,
      data = c(
        script_header(data_entries("data_integrity")),
        script_data_lines(),
        "# ---- Data integrity checks ----",
        "integrity <- data_integrity(dat, meta)",
        "print(integrity$table[, c(\"Check\", \"Result\", \"Status\")])",
        "",
        "# ---- Distribution of a variable ----",
        paste0("variable <- ", dput_text(input$dist_var %||% current_map()$trait)),
        "values <- raw[[variable]]",
        "num <- suppressWarnings(as.numeric(values))",
        "if (mean(is.finite(num[!is.na(values)])) > 0.9) {",
        "  print(ggplot(data.frame(x = num[is.finite(num)]), aes(x)) + geom_histogram(bins = 30, fill = \"#0072B2\", colour = \"white\") +",
        "          labs(x = variable, y = \"Rows\") + theme_minimal())",
        "} else {",
        "  tab <- sort(table(values), decreasing = TRUE)",
        "  print(ggplot(data.frame(level = names(tab), n = as.numeric(tab)), aes(stats::reorder(level, n), n)) + geom_col(fill = \"#0072B2\") +",
        "          coord_flip() + labs(x = variable, y = \"Rows\") + theme_minimal())",
        "}",
        "print(utils::head(raw, 10))"),
      visual = c(
        script_header(data_entries("mode_or_na", "proxy_values", "binned_trajectory", "bin_differences", "bin_difference_trends",
                                   "trait_by_age_bins", "a2_bin_slopes", "a2_slope_trend", "disappearance_data", "terminal_data",
                                   "selection_differentials", "life_table")),
        script_data_lines(),
        "# ---- Settings (copied from the app) ----",
        paste0("proxy <- ", dput_text(safe_get(vis_px()) %||% "ALR"), "   # \"ALR\", \"Mean age\", \"LS\" or \"AFR\""),
        paste0("n_bins <- ", dput_text(as.numeric(input$n_bins %||% 4))),
        paste0("bin_method <- ", dput_text(input$bin_method %||% "equal"), "   # \"equal\" or \"quantile\""),
        paste0("trait_scale <- ", dput_text(input$trait_scale %||% "raw"), "   # \"raw\" or \"log1p\""),
        paste0("facet_column <- ", dput_text(input$facet_var %||% ""), "   # \"\" for no panels"),
        "",
        "vis <- dat",
        "if (identical(trait_scale, \"log1p\")) vis$trait <- log1p(vis$trait)",
        "facet <- NULL",
        "if (nzchar(facet_column)) {",
        "  ids <- as.character(raw[[map$id]])",
        "  if (isTRUE(meta$has_group) && isTRUE(meta$nested)) {",
        "    g <- as.character(raw[[map$group]])",
        "    ids <- ifelse(is.na(g), ids, paste(g, ids, sep = \"/\"))",
        "  }",
        "  lv <- as.character(raw[[facet_column]])",
        "  ok <- !is.na(ids) & nzchar(ids) & !is.na(lv) & nzchar(lv)",
        "  per <- tapply(lv[ok], ids[ok], mode_or_na)",
        "  facet <- stats::setNames(paste0(facet_column, \": \", as.character(per)), names(per))",
        "}",
        "panels <- if (is.null(facet)) NULL else facet_wrap(~ facet)",
        "pvals <- proxy_values(im, proxy)",
        "",
        "# ---- Trait trajectory within bins (at least 3 individuals per point) ----",
        "traj <- binned_trajectory(vis, pvals, n_bins, bin_method, 3, facet = facet)",
        "print(ggplot(traj, aes(age, mean, colour = bin, group = bin)) + geom_line(linewidth = 1) + geom_point(aes(size = n)) + panels +",
        "        labs(x = \"Age\", y = \"Mean trait\", colour = paste(proxy, \"bin\"), size = \"Individuals\") + theme_minimal())",
        "",
        "# ---- Difference between consecutive bins at each age ----",
        "diffs <- bin_differences(traj, \"successive\")",
        "print(bin_difference_trends(diffs))",
        "print(ggplot(diffs, aes(age, difference, colour = pair)) + geom_hline(yintercept = 0, linetype = 2) + geom_point() +",
        "        geom_smooth(method = \"lm\", formula = y ~ x, se = FALSE) + panels + labs(x = \"Age\", y = \"Difference in mean trait\", colour = NULL) + theme_minimal())",
        "",
        "# ---- Trait against the grouping variable within age bins (at least 3 individuals per bin) ----",
        "a2 <- trait_by_age_bins(vis, pvals, n_bins, facet = facet)",
        "a2 <- a2[stats::ave(rep(1, nrow(a2)), a2$facet, a2$age_bin, FUN = length) >= 3, , drop = FALSE]",
        "slopes <- a2_bin_slopes(a2)",
        "print(slopes)",
        "print(a2_slope_trend(slopes))",
        "print(ggplot(a2, aes(proxy, trait, colour = age_bin)) + geom_point(alpha = 0.35) +",
        "        geom_smooth(aes(group = age_bin), method = \"lm\", formula = y ~ x, se = FALSE) + panels +",
        "        labs(x = proxy, y = \"Mean trait within the age bin\", colour = \"Age bin\") + theme_minimal())",
        "",
        "# ---- Disappearance diagnostics ----",
        "terminal <- terminal_data(dat, meta, 3)",
        "if (nrow(terminal)) print(ggplot(terminal, aes(before, trait, colour = end_bin, group = end_bin)) + geom_line() + geom_point(aes(size = n)) +",
        "                         scale_x_reverse() + labs(x = \"Occasions before death\", y = \"Mean trait\", colour = \"Lifespan\") + theme_minimal())",
        "sel <- selection_differentials(disappearance_data(dat, meta, FALSE))",
        "print(sel)",
        "if (nrow(sel)) print(ggplot(sel, aes(Age, Differential)) + geom_hline(yintercept = 0, linetype = 2) +",
        "                    geom_errorbar(aes(ymin = Differential - 1.96 * Differential_SE, ymax = Differential + 1.96 * Differential_SE), width = 0) +",
        "                    geom_point(aes(size = N)) + labs(y = \"Selection differential (within-age SD)\") + theme_minimal())",
        "hazard <- life_table(dat, meta)",
        "if (isTRUE(hazard$ok)) {",
        "  print(hazard$table)",
        "  cat(\"Overall hazard:\", round(hazard$overall, 3), \"; test of a constant hazard: p =\", signif(hazard$p_constant, 3), \"\\n\")",
        "  print(ggplot(hazard$table, aes(Age, Hazard)) + geom_hline(yintercept = hazard$overall, linetype = 2) +",
        "          geom_errorbar(aes(ymin = Lower_95, ymax = Upper_95), width = 0) + geom_point(aes(size = At_risk)) +",
        "          labs(y = \"Probability of disappearing before the next occasion\") + theme_minimal())",
        "}"),
      sampling = c(
        script_header(data_entries("build_missing_grid", "missingness_summary", "missingness_by_variable", "grid_display")),
        script_data_lines(),
        "# ---- Settings (copied from the app) ----",
        paste0("start_age <- ", dput_text(miss_start_age()), "  # age at first trait expression (NA unless missingness is counted from AFE)"),
        "start_mode <- if (is.finite(start_age)) \"same\" else \"afr\"",
        paste0("n_show <- ", dput_text(as.numeric(safe_get(heat_n()) %||% 150))),
        "life_known <- isTRUE(meta$has_life) && !isTRUE(meta$life_auto)",
        "",
        "# ---- Expected occasions and missingness ----",
        "grid <- build_missing_grid(dat, margin = 0, start_mode = start_mode, start_age = start_age)",
        "ms <- missingness_summary(dat, grid, 0.05, life_known)",
        "cat(sprintf(\"Missing expected occasions: %.1f%%\\n\", ms$percent))",
        "print(ms$drivers)",
        "cat(ms$guidance, \"\\n\")",
        "",
        "# ---- Sampling grid (individuals ordered by ALR) ----",
        "ids <- unique(grid$id)",
        "ids <- ids[order(im$alr[match(ids, im$id)])]",
        "if (length(ids) > n_show) ids <- ids[unique(round(seq(1, length(ids), length.out = n_show)))]",
        "cells <- grid_display(grid, ids)",
        "cells$id <- factor(cells$id, levels = rev(ids))",
        "print(ggplot(cells, aes(age, id, fill = status)) + geom_tile() +",
        "        scale_fill_manual(values = c(\"Observed\" = \"#7C8060\", \"Missed\" = \"#C0392B\", \"Death or ALR\" = \"#E67E22\", \"Not expected\" = \"#FFFFFF\")) +",
        "        labs(x = \"Age\", y = NULL, fill = NULL) + theme_minimal() + theme(axis.text.y = element_blank()))",
        "",
        "# ---- Missingness by age ----",
        "by_age <- missingness_by_variable(dat, grid, \"Age\")",
        "print(ggplot(by_age, aes(value, 100 * missing_rate)) + geom_line() + geom_point(aes(size = n)) +",
        "        labs(x = \"Age\", y = \"% of expected occasions missed\") + theme_minimal())",
        "",
        "# ---- Agreement between lifespan proxies ----",
        "cat(\"r(mean age, ALR) =\", round(stats::cor(im$mean_age, im$alr, use = \"complete.obs\"), 3), \"\\n\")",
        "print(ggplot(im, aes(alr, mean_age)) + geom_point(alpha = 0.5) + geom_smooth(method = \"lm\", formula = y ~ x) +",
        "        labs(x = \"ALR (age at last record)\", y = \"Mean age of the records\") + theme_minimal())",
        "if (life_known) {",
        "  cat(\"r(ALR, LS) =\", round(stats::cor(im$alr, im$lifespan, use = \"complete.obs\"), 3), \"\\n\")",
        "  print(ggplot(im, aes(lifespan, alr)) + geom_point(alpha = 0.5) + geom_abline(linetype = 2) + labs(x = \"Lifespan (LS)\", y = \"ALR\") + theme_minimal())",
        "}"),
      individual = {
        s <- model_settings()
        fns <- intersect(input$b2_functions %||% AGE_FUNCTIONS, MODEL_FUNCTIONS)
        c(
          script_header(data_entries("fit_individual_function", "compare_individual_functions", "compare_population_functions", "observed_trajectory"),
                        packages = c("ggplot2", "lme4")),
          script_data_lines(),
          "# ---- Settings (copied from the app) ----",
          paste0("fun <- ", dput_text(input$a3_function %||% "Quadratic"), "   # ageing function fitted to each individual"),
          paste0("functions <- ", dput_text(fns), "   # functions compared at the population level"),
          paste0("family <- ", dput_text(s$family)),
          paste0("random_slope <- ", dput_text(s$random_slope)),
          paste0("zero_inflation <- ", dput_text(s$zi)),
          "",
          "# ---- Individual fits ----",
          "fits <- fit_individual_function(dat, fun, 1, draw_ids = unique(dat$id))",
          "cat(\"Individuals fitted:\", nrow(fits$coefs), \"of\", fits$n_total, \"\\n\")",
          "print(fits$mean_curve)",
          "obs <- observed_trajectory(dat)",
          "print(ggplot() + geom_point(data = obs, aes(age, fitted, size = n), shape = 21, fill = \"grey55\") +",
          "        geom_ribbon(data = fits$mean_curve, aes(age, ymin = lo, ymax = hi), alpha = 0.2) +",
          "        geom_line(data = fits$mean_curve, aes(age, fitted), colour = \"#D55E00\", linewidth = 1.2) +",
          "        labs(x = \"Age\", y = \"Trait\", subtitle = \"Line: mean of the individual coefficients (95% band); points: observed means\") + theme_minimal())",
          "show <- utils::head(unique(fits$curves$id), 20)",
          "print(ggplot(fits$curves[fits$curves$id %in% show, ], aes(age, fitted)) + geom_line(colour = \"#0072B2\") +",
          "        geom_point(data = dat[dat$id %in% show, ], aes(age, trait), size = 1) + facet_wrap(~ id, scales = \"free\") + theme_minimal())",
          "",
          "# ---- Ageing functions compared across individuals (mean delta AICc) ----",
          "print(compare_individual_functions(dat, 1))",
          "",
          "# ---- Ageing functions compared at the population level (Model 1) ----",
          paste0("population <- compare_population_functions(dat, meta, family = family, random_slope = random_slope, standardise = ",
                 dput_text(isTRUE(s$standardise)), ", zi_str = zero_inflation, functions = functions, model = \"M1\", among = ", dput_text(s$among), ")"),
          "print(population$table)",
          "if (nrow(population$curves)) print(ggplot() + geom_point(data = obs, aes(age, fitted, size = n), shape = 21, fill = \"grey55\") +",
          "        geom_line(data = population$curves, aes(age, fitted, colour = Function), linewidth = 1.1) + labs(x = \"Age\", y = \"Trait\") + theme_minimal())")
      },
      character(0))
    paste(lines, collapse = "\n")
  }
  for (sec in c("data", "visual", "sampling", "individual")) {
    local({
      section <- sec
      output[[paste0("code_", section, "_text")]] <- renderText({
        tryCatch(section_code(section), shiny.silent.error = function(e) "",
                 error = function(e) paste("# The R code could not be generated:", conditionMessage(e)))
      })
      output[[paste0("download_code_", section)]] <- downloadHandler(
        filename = function() paste0("disappR_", section, "_section.R"),
        content = function(file) {
          txt <- tryCatch(section_code(section), error = function(e) paste("# The R code could not be generated:", conditionMessage(e)))
          con <- file(file, open = "w", encoding = "UTF-8")
          on.exit(close(con))
          writeLines(txt, con)
        })
    })
  }

  # ---- DHARMa ----
  dharma_res <- reactiveVal(NULL)
  output$dharma_model_ui <- renderUI({
    r <- safe_get(fitted_models())
    if (is.null(r)) return(p(class = "small-note", "Fit models first."))
    ch <- names(r$fits)
    best <- sub("Model ", "M", r$aic$Model[[1]])
    selectInput("dharma_model", "Model", choices = stats::setNames(ch, model_label(ch)), selected = best)
  })
  observeEvent(input$run_dharma, disappr_guard("input$run_dharma", {
    r <- safe_get(fitted_models())
    if (is.null(r) || is.null(input$dharma_model) || is.null(r$fits[[input$dharma_model]])) {
      notify("Fit models first.", type = "warning")
      return(invisible(NULL))
    }
    res <- run_with_progress("Simulating residuals", function(pr) run_dharma(r$fits[[input$dharma_model]], r$family))
    res$model <- input$dharma_model
    res$signature <- model_sig()
    dharma_res(res)
  }))
  output$dharma_table <- renderTable({
    z <- dharma_res()
    validate(need(!is.null(z), "Press 'Simulate residuals'."))
    validate(need(isTRUE(z$ok), z$message %||% "DHARMa failed."))
    validate(need(identical(z$signature, model_sig()), "Models were refitted: simulate residuals again."))
    data.frame(Model = model_label(z$model), Test = z$table$Test, P = vapply(z$table$P_value, format_p, character(1)))
  }, striped = TRUE, spacing = "xs")
  output$dharma_plot <- renderPlot({
    z <- dharma_res()
    validate(need(!is.null(z) && isTRUE(z$ok), ""))
    validate(need(identical(z$signature, model_sig()), ""))
    graphics::plot(z$sim)
  })

  # ---- R code ----
  code_text <- reactive({
    r <- fitted_models()
    src <- switch(input$data_source %||% "toy", upload = input$data_file$name %||% "your_data.csv",
                  example = current_example()$file, "disappR_simulated.csv")
    model_r_code(r, meta(), src)
  })
  output$code_ui <- renderUI(tags$pre(class = "code-out", code_text()))
  make_code_download <- function() {
    downloadHandler(
      filename = function() paste0("disappR_models_", Sys.Date(), ".R"),
      content = function(file) writeLines(code_text(), file)
    )
  }
  output$download_code <- make_code_download()
  output$download_code2 <- make_code_download()

  # ======================================================================
  # 7. Summary and report
  # ======================================================================
  summary_parts <- reactive({
    d <- dat()
    m <- meta()
    integ <- integrity()
    ms <- msum()
    im <- imet()
    warn_rows <- integ$table[integ$table$Status == "Warning", , drop = FALSE]
    a2_alr <- proxy_slopes_by_age(d, proxy_values(im, if (life_known()) "LS" else "ALR"), "r")
    a2_afr <- if (length(unique(im$entry[is.finite(im$entry)])) >= 2) proxy_slopes_by_age(d, proxy_values(im, "AFR"), "r") else NULL
    r <- model_res()
    models_ok <- !is.null(r) && isTRUE(r$ok) && identical(r$signature, model_sig())
    rec <- character(0)
    if (models_ok) {
      aic <- r$aic
      lrt <- r$lrt
      p_of <- function(cmp) if (nrow(lrt) && cmp %in% lrt$Comparison) lrt$P_value[lrt$Comparison == cmp] else NA_real_
      p24 <- p_of("Model 2 vs Model 4")
      p12 <- p_of("Model 1 vs Model 2")
      dA <- function(mod) if (mod %in% aic$Model) aic$Delta_AIC[aic$Model == mod] else NA_real_
      if (is.finite(p24) && p24 < 0.05 && isTRUE(dA("Model 4") + 2 < dA("Model 2"))) {
        rec <- c(rec, sprintf("The data are consistent with age-dependent selective disappearance: Model 4 fits better than Model 2 (LRT p = %s; \u0394AIC %.1f vs %.1f). Consider reporting Model 4 alongside Model 2, and check that the pattern also appears in the visual diagnosis plots.",
                              format_p(p24), dA("Model 4"), dA("Model 2")))
      } else if (is.finite(p24)) {
        rec <- c(rec, sprintf("The ALR \u00d7 age interaction is not clearly supported (Model 2 vs 4 LRT p = %s); this may mean no age-dependent selection, or too little information at old ages.", format_p(p24)))
        if (is.finite(p12) && p12 < 0.05) rec <- c(rec, sprintf("The additive ALR term is supported (Model 1 vs 2 p = %s), which is consistent with age-independent selective disappearance.", format_p(p12)))
      }
      if (all(c("Model 4", "Model 5") %in% aic$Model[aic$Eligible])) {
        d45 <- dA("Model 5") - dA("Model 4")
        rec <- c(rec, sprintf("Model 4 vs Model 5: \u0394AIC = %.1f in favour of %s%s. %s", abs(d45), if (d45 >= 0) "Model 4 (ALR)" else "Model 5 (mean age)",
                              if (abs(d45) < 2) " (similar support)" else "", ms$guidance))
      }
      if (all(c("Model 4", "Model 8") %in% aic$Model)) {
        p48 <- p_of("Model 4 vs Model 8")
        if (is.finite(p48)) rec <- c(rec, sprintf("Selective appearance (Model 4 vs 8): LRT p = %s.", format_p(p48)))
      }
      if (any(unlist(r$validity) != "Valid")) rec <- c(rec, "Some fits are classified Caution or Failed; inspect the fitting status and check whether a simpler random-effect structure gives the same conclusion.")
      if (isTRUE(r$n_excluded > 0)) rec <- c(rec, sprintf("%d fit(s) with an invalid Hessian were excluded from these statements.", r$n_excluded))
      rec <- c(rec, tryCatch(consistency_notes(r, if (is_toy()) safe_get(deviation_raw()) else NULL, NULL), error = function(e) character(0)),
               "These statements come from automated rules. They depend on the error family, random-effect structure, ageing function, proxy quality, sampling and fit validity, and are not proof of selection.")
    }
    fam <- integ$family
    fam_note <- if (!identical(fam$family, "gaussian") && models_ok && identical(r$family, "gaussian")) "Warning: a count trait was modelled as Gaussian; refit with a count family." else NULL
    list(d = d, m = m, integ = integ, warn_rows = warn_rows, ms = ms, a2 = a2_alr, a2_afr = a2_afr, r = r,
         models_ok = models_ok, rec = rec, fam_note = fam_note, b2 = b2_res(), a3 = a3_compare())
  })

  output$summary_ui <- renderUI({
    z <- summary_parts()
    tc <- truth()
    ms <- z$ms
    blocks <- list()
    if (!is.null(tc)) {
      txt <- toy_truth_text(tc$cfg)
      blocks <- c(blocks, list(div(class = "summary-block", h4("Simulated truth"), p(txt$sim), tags$ul(lapply(txt$expect, tags$li)))))
    }
    blocks <- c(blocks, list(
      div(class = "summary-block", h4("1 \u00b7 Data"),
          p(sprintf("%s rows from %s individuals; %d distinct ages.", format(nrow(z$d), big.mark = ","), format(length(unique(z$d$id)), big.mark = ","), length(unique(z$d$age)))),
          if (nrow(z$warn_rows)) tags$ul(lapply(paste0(z$warn_rows$Check, ": ", z$warn_rows$Result), tags$li)) else p("No integrity warnings."),
          p(z$integ$family$text)),
      div(class = "summary-block", h4("2 \u00b7 Visual diagnosis"),
          p(strong(if (life_known()) "LS: " else "ALR: "), a2_interpretation(z$a2, if (life_known()) "LS" else "ALR")),
          if (!is.null(z$a2_afr)) p(strong("AFR: "), a2_interpretation(z$a2_afr, "AFR")) else p("AFR does not vary: selective appearance not assessable.")),
      div(class = "summary-block", h4("3 \u00b7 Missingness and proxies"),
          p(sprintf("%.1f%% of expected occasions missing; p = %.2f.", ms$percent, ms$p), ms$type),
          p(ms$guidance)),
      div(class = "summary-block", h4("4 \u00b7 Individual and population trajectories"),
          p(if (nrow(z$a3)) paste0("Individual fits: best mean \u0394AICc = ", z$a3$Function[[1]], " (", z$a3$N_common[[1]], " individuals with all functions estimable).") else "Individual fits not estimable."),
          p(if (!is.null(z$b2) && nrow(z$b2$table)) paste0("Population level (last run): best AIC = ", z$b2$table$Function[[1]], ".") else "Population-level function comparison not yet run.")),
      div(class = "summary-block", h4("5 \u00b7 Modelling"),
          if (z$models_ok) tagList(
            p(sprintf("%s; lowest AIC among eligible fits: %s (%s). A lower AIC is relative support, not proof.", family_label(z$r$family), z$r$aic$Model[[1]], z$r$aic$Fit[[1]])),
            tags$ul(lapply(z$rec, tags$li)),
            if (!is.null(z$fam_note)) p(strong(z$fam_note)) else NULL,
            p(class = "small-note", "Decomposition: the manuscript recommends against using it to recover the latent trajectory unless there is no age-dependent selective disappearance, little missingness and many consecutively sampled individuals.")
          ) else p("Models not fitted for the current data and settings (Models tab \u2192 'Fit models')."))
    ))
    tagList(blocks)
  })

  # ======================================================================
  # Saved results: every "Save to summary" button adds a snapshot here
  # ======================================================================
  saved_results <- reactiveVal(list())
  save_counter <- reactiveVal(0L)
  nothing_to_save <- function(msg = "Nothing to save yet: this result is not available for the current data and settings.") {
    notify(msg, type = "warning")
    invisible(NULL)
  }
  data_description <- function() {
    src <- input$data_source %||% "toy"
    d <- safe_get(dat())
    n_txt <- if (is.null(d)) "" else sprintf(" (%s rows, %s individuals)", format(nrow(d), big.mark = ","), format(length(unique(d$id)), big.mark = ","))
    base <- switch(src, example = current_example()$label, upload = paste("Uploaded file", input$data_file$name %||% ""),
                   paste("Simulated:", toy_truth_text(toy_cfg())$sim))
    paste0(base, n_txt)
  }
  add_saved <- function(section, title, text = character(0), tables = list(), plots = list(), code = NULL) {
    k <- save_counter() + 1L
    save_counter(k)
    id <- paste0("saved", k)
    tables <- Filter(function(t) is.data.frame(t) && nrow(t) > 0, tables)
    plots <- Filter(Negate(is.null), plots)
    entry <- list(id = id, section = section, title = title, time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  data = data_description(), text = as.character(text), tables = tables, plots = plots, code = code)
    lst <- saved_results()
    lst[[id]] <- entry
    saved_results(lst)
    for (j in seq_along(plots)) {
      local({
        pl <- plots[[j]]
        output[[paste0(id, "_plot", j)]] <- renderPlot(draw_saved_plot(pl))
      })
    }
    for (j in seq_along(tables)) {
      local({
        tb <- tables[[j]]
        output[[paste0(id, "_table", j)]] <- renderTable(tb, striped = TRUE, spacing = "xs")
      })
    }
    notify(paste0("Saved to Summary & export: ", title), duration = 3)
  }
  fmt_r <- function(r) if (length(r) == 1 && is.finite(r)) sprintf("%.2f", r) else "NA"
  aic_display <- function(aic) {
    data.frame(Model = aic$Model, AIC = round(aic$AIC, 1), dAIC = round(aic$Delta_AIC, 1), Weight = round(aic$Weight, 3),
               df = aic$df, logLik = round(aic$logLik, 1), N = aic$N)
  }
  lrt_display <- function(l) {
    if (is.null(l) || !nrow(l)) return(data.frame())
    l$Chisq <- round(l$Chisq, 2)
    l$P_value <- vapply(l$P_value, format_p, character(1))
    l
  }
  start_text <- function() {
    if (isTRUE(safe_get(is_toy()))) return("set by the simulation")
    {
      afe <- miss_start_age()
      paste0(if (isTRUE(is.finite(afe))) paste0("missingness counted from age ", format_num(afe), " (first trait expression); ") else "",
             if (identical(input$afr_mode %||% "auto", "column")) "AFR from a mapped column" else "AFR from each individual's first record")
    }
  }

  observeEvent(input$save_integrity, disappr_guard("input$save_integrity", {
    ig <- safe_get(integrity())
    if (is.null(ig)) return(nothing_to_save())
    add_saved("1 Data", "Data integrity checks",
              text = c(ig$family$text, if (length(ig$bad_ls_ids)) paste("IDs with LS < last record:", paste(utils::head(ig$bad_ls_ids, 30), collapse = ", "))),
              tables = list(`Integrity checks` = ig$table))
  }))
  observeEvent(input$save_distribution, disappr_guard("input$save_distribution", {
    pl <- safe_get(dist_plot_obj())
    if (is.null(pl)) return(nothing_to_save())
    add_saved("1 Data", paste("Distribution of", input$dist_var),
              text = paste("One value per", if (identical(input$dist_unit, "individual")) "individual" else "row"), plots = list(pl))
  }))
  observeEvent(input$save_a1, disappr_guard("input$save_a1", {
    pl <- safe_get(a1_plot_obj())
    if (is.null(pl)) return(nothing_to_save())
    px <- safe_get(a1_px()) %||% ""
    add_saved("2 Visual diagnosis", paste("Trait trajectory by", px, "bin"),
              text = paste0(sprintf("%s %s bins (%s boundaries); at least 3 individuals per point; trait scale: %s",
                                    input$n_bins %||% 4, px, input$bin_method %||% "equal", input$trait_scale %||% "raw"),
                            facet_text()),
              plots = list(pl))
  }))
  observeEvent(input$save_a1_diff, disappr_guard("input$save_a1_diff", {
    pl <- safe_get(bin_diff_obj())
    dz <- safe_get(bin_diff_data())
    if (is.null(pl) || is.null(dz)) return(nothing_to_save())
    faceted <- !is.null(facet_map())
    tab <- data.frame(Pair = as.character(dz$pair), Age = dz$age, Difference = signif(dz$difference, 4))
    if (faceted) tab <- cbind(Panel = dz$facet, tab, stringsAsFactors = FALSE)
    tr <- safe_get(bin_diff_trends())
    if (!is.null(tr) && nrow(tr)) {
      tr <- data.frame(Panel = tr$Facet, Pair = tr$Pair, Ages = tr$N_ages, From = tr$Age_from, To = tr$Age_to,
                       `Slope per unit age` = signif(tr$Slope_per_age, 3), P = vapply(tr$P, format_p, character(1)),
                       check.names = FALSE, stringsAsFactors = FALSE)
      if (!faceted) tr$Panel <- NULL
    }
    pooled <- NULL
    if (identical(input$diff_lines, "pooled")) {
      pz <- safe_get(bin_diff_pooled())
      if (!is.null(pz) && nrow(pz)) {
        pooled <- data.frame(Panel = pz$facet, `Slope per unit age (all pairs)` = signif(pz$Slope_per_age, 3), P = vapply(pz$P, format_p, character(1)),
                             Points = pz$N_points, check.names = FALSE, stringsAsFactors = FALSE)
        if (!faceted) pooled$Panel <- NULL
      }
    }
    add_saved("2 Visual diagnosis", paste("Difference between", safe_get(a1_px()) %||% "", "bins at each age"),
              text = paste0("Consecutive bin pairs; trend lines: ", if (identical(input$diff_lines, "pooled")) "one across all pairs" else "one per pair", facet_text()),
              tables = list(`Trend across all pairs` = pooled, `Trend of each difference against age` = tr, `Differences (higher minus lower bin)` = tab), plots = list(pl))
  }))
  observeEvent(input$save_a2, disappr_guard("input$save_a2", {
    pl <- safe_get(a2_plot_obj())
    if (is.null(pl)) return(nothing_to_save())
    px <- safe_get(a2_px()) %||% "ALR"
    stat <- safe_get(proxy_slopes_by_age(visual_dat(), proxy_values(imet(), px), "r"))
    sl <- safe_get(a2_slopes())
    add_saved("2 Visual diagnosis", paste("Trait against", px, "within age bins"),
              text = c(paste0(sprintf("%s on the x-axis; %s age bins; trait scale: %s", px, input$n_bins %||% 4, input$trait_scale %||% "raw"), facet_text()),
                       safe_get(a2_trend_lines()),
                       if (!is.null(stat)) paste("Correlation at each age:", a2_interpretation(stat, px))),
              tables = list(`Regression coefficient in each age bin` = if (!is.null(sl) && nrow(sl)) a2_slope_display(sl) else NULL),
              plots = list(pl))
  }))
  observeEvent(input$save_a2_table, disappr_guard("input$save_a2_table", {
    sl <- safe_get(a2_slopes())
    if (is.null(sl) || !nrow(sl)) return(nothing_to_save())
    px <- safe_get(a2_px()) %||% "ALR"
    add_saved("2 Visual diagnosis", paste("Regression coefficient of the trait on", px, "across age bins"),
              text = c(paste0(input$n_bins %||% 4, " age bins; trait scale: ", input$trait_scale %||% "raw", facet_text()),
                       safe_get(a2_trend_lines())),
              tables = list(`Regression coefficient in each age bin` = a2_slope_display(sl)))
  }))
  observeEvent(input$save_heatmap, disappr_guard("input$save_heatmap", {
    pl <- safe_get(heatmap_obj())
    if (is.null(pl)) return(nothing_to_save())
    ms <- safe_get(msum())
    add_saved("3 Missingness and proxies", "Sampling grid",
              text = c(if (!is.null(ms)) sprintf("%.1f%% of expected occasions missed; %s", ms$percent, ms$pattern),
                       paste("Start of sampling:", start_text()),
                       paste("Rows ordered by", switch(input$heat_order %||% "alr", alr = "ALR", afr = "AFR", trait = "mean trait value", id = "ID", "random order"))),
              plots = list(pl))
  }))
  observeEvent(input$save_sampling, disappr_guard("input$save_sampling", {
    ms <- safe_get(msum())
    if (is.null(ms)) return(nothing_to_save())
    add_saved("3 Missingness and proxies", "Sampling summary and proxy guidance",
              text = c(sprintf("Missed expected occasions: %.1f%%; detection p = %.2f", ms$percent, ms$p), ms$type, ms$guidance,
                       paste("Missed occasions:", ms$pattern),
                       sprintf("r(mean age, ALR) = %s; r(ALR, LS) = %s; r(mean age, LS) = %s", fmt_r(ms$r_mean_alr), fmt_r(ms$r_alr_ls), fmt_r(ms$r_mean_ls)),
                       paste("Start of sampling:", start_text())),
              tables = list(`Associations with missingness` = ms$drivers))
  }))
  observeEvent(input$save_missing_age, disappr_guard("input$save_missing_age", {
    pl <- safe_get(missing_by_age_obj())
    if (is.null(pl)) return(nothing_to_save())
    add_saved("3 Missingness and proxies", "Missingness and sample size by age",
              tables = list(`Coverage by age` = safe_get(coverage_by_age(grid_r()))), plots = list(pl))
  }))
  observeEvent(input$save_missing_var, disappr_guard("input$save_missing_var", {
    pl <- safe_get(missing_vs_var_obj())
    ms <- safe_get(msum())
    if (is.null(pl) && is.null(ms)) return(nothing_to_save())
    add_saved("3 Missingness and proxies", paste("Missingness against", input$miss_var %||% "ALR"),
              tables = list(`Associations with missingness` = if (!is.null(ms)) ms$drivers else NULL), plots = list(pl))
  }))
  observeEvent(input$save_proxies, disappr_guard("input$save_proxies", {
    lk <- isTRUE(safe_get(life_known()))
    pls <- list(safe_get(proxy_alr_mean_obj()), if (lk) safe_get(proxy_alr_ls_obj()), if (lk) safe_get(proxy_mean_ls_obj()))
    pls <- Filter(Negate(is.null), pls)
    if (!length(pls)) return(nothing_to_save())
    ms <- safe_get(msum())
    add_saved("3 Missingness and proxies", "Agreement between ALR, mean age and lifespan",
              text = if (!is.null(ms)) sprintf("r(mean age, ALR) = %s; r(ALR, LS) = %s; r(mean age, LS) = %s",
                                              fmt_r(ms$r_mean_alr), fmt_r(ms$r_alr_ls), fmt_r(ms$r_mean_ls)),
              plots = pls)
  }))
  observeEvent(input$save_a3, disappr_guard("input$save_a3", {
    pl <- safe_get(a3_mean_obj())
    z <- safe_get(a3_fit())
    if (is.null(pl) || is.null(z)) return(nothing_to_save())
    fn <- input$a3_function %||% "Quadratic"
    im <- safe_get(imet())
    txt <- c(sprintf("%s function fitted to %d of %d individuals", fn, length(z$fitted_ids), z$n_total),
             paste("Population curve:", switch(input$a3_recon %||% "both", coef = "mean of coefficients",
                                               fun = "mean of individual functions", "mean of coefficients and mean of individual functions"),
                   if (isTRUE(z$nonlinear)) "(non-linear function: the two differ)" else "(linear in parameters: the two are identical)"))
    if (!is.null(im) && length(z$fitted_ids)) {
      txt <- c(txt, sprintf("Mean ALR of fitted individuals %s vs %s for all individuals",
                            format_num(mean(im$alr[im$id %in% z$fitted_ids], na.rm = TRUE)), format_num(mean(im$alr, na.rm = TRUE))))
    }
    cf <- z$coefs
    nms <- setdiff(names(cf), "id")
    tab <- data.frame(Coefficient = nms,
                      Mean = signif(vapply(cf[nms], function(v) mean(v, na.rm = TRUE), numeric(1)), 5),
                      SD = signif(vapply(cf[nms], function(v) stats::sd(v, na.rm = TRUE), numeric(1)), 4))
    add_saved("4 Individual and population trajectories", "Mean-coefficient trajectory", text = txt, tables = list(`Individual coefficients` = tab), plots = list(pl))
  }))
  observeEvent(input$save_a3_compare, disappr_guard("input$save_a3_compare", {
    tab <- safe_get(a3_compare())
    if (is.null(tab) || !nrow(tab)) return(nothing_to_save())
    tab[] <- lapply(tab, function(v) if (is.numeric(v)) signif(v, 4) else v)
    add_saved("4 Individual and population trajectories", "Individual-level function comparison", tables = list(`Function comparison` = tab))
  }))
  observeEvent(input$save_b2, disappr_guard("input$save_b2", {
    r <- safe_get(b2_shown())
    if (is.null(r)) return(nothing_to_save("Run 'Compare ageing functions' first."))
    s <- model_settings()
    add_saved("4 Individual and population trajectories", "Population-level ageing functions",
              text = paste0("Model 1 structure; ", family_label(s$family), "; ", slope_text(s$random_slope)),
              tables = list(`Function comparison` = r$table), plots = list(safe_get(b2_plot_obj())))
  }))
  observeEvent(input$save_models, disappr_guard("input$save_models", {
    r <- safe_get(fitted_models())
    if (is.null(r)) return(nothing_to_save("Fit models first."))
    aic <- r$aic
    best_id <- sub("^Model ", "M", aic$Model[[1]])
    elig <- aic[aic$Eligible, , drop = FALSE]
    d45 <- if (all(c("Model 4", "Model 5") %in% elig$Model)) elig$AIC[elig$Model == "Model 5"] - elig$AIC[elig$Model == "Model 4"] else NA_real_
    txt <- c(sprintf("%s; %s ageing function; random effects %s", family_label(r$family), r$age_function,
                     if (isTRUE(r$nonlinear)) r$random else random_display(meta(), r$random_slope)),
             sprintf("Lowest AIC among eligible fits: %s (%s)%s", aic$Model[[1]], aic$Fit[[1]], if (nrow(elig) > 1) sprintf("; next model \u0394AIC = %.1f", elig$Delta_AIC[[2]]) else ""),
             if (isTRUE(r$nonlinear)) paste("Model:", r$formulas[[best_id]]) else paste0("Formula: trait ~ ", display_term(r$formulas[[best_id]]), " + ", display_term(r$random)),
             if (nzchar(r$random_note %||% "")) r$random_note,
             if (isTRUE(r$n_excluded > 0)) sprintf("%d fit(s) with an invalid Hessian excluded from the ranking", r$n_excluded),
             if (is.finite(d45)) sprintf("AIC(Model 5) \u2212 AIC(Model 4) = %.1f", d45),
             if (r$n_dropped > 0) sprintf("%d rows dropped from all models (%s)", r$n_dropped, paste(r$drop_by, collapse = "; ")),
             safe_get(summary_parts()$rec))
    st <- data.frame(Model = model_label(names(r$status)), Status = unlist(r$status, use.names = FALSE))
    add_saved("5 Modelling", "Model comparison", text = txt,
              tables = list(AIC = aic_display(aic), `Likelihood-ratio tests` = lrt_display(r$lrt), `Fitting status` = st,
                            `Random-effect variances` = safe_get(varcomp_display())))
  }))
  observeEvent(input$save_predictions, disappr_guard("input$save_predictions", {
    pl <- safe_get(pred_plot_obj())
    if (is.null(pl)) return(nothing_to_save("Fit models first."))
    dv <- if (isTRUE(safe_get(is_toy()))) safe_get(deviation_tab()) else NULL
    add_saved("5 Modelling", "Population-level ageing trajectories",
              text = paste0("Models drawn: ", paste(model_label(intersect(input$pred_models %||% safe_get(eligible_models()), safe_get(eligible_models()) %||% character(0))), collapse = ", "),
                            if (nzchar(input$pred_by %||% "")) paste0("; by ", input$pred_by) else "",
                            if (isTRUE(input$show_a3)) paste0("; reconstruction from individual fits with the ", input$a3_function %||% "Quadratic", " function") else ""),
              tables = list(`Deviation from the simulated truth` = dv), plots = list(pl))
  }))
  observeEvent(input$save_coefs, disappr_guard("input$save_coefs", {
    r <- safe_get(fitted_models())
    if (is.null(r)) return(nothing_to_save("Fit models first."))
    m <- input$coef_model %||% sub("Model ", "M", r$aic$Model[[1]])
    x <- r$coefficients[r$coefficients$Model == model_label(m), , drop = FALSE]
    if (!nrow(x)) return(nothing_to_save())
    tab <- coef_display(x, r)
    sc <- safe_get(scaling_constants(r))
    add_saved("5 Modelling", paste("Coefficients of", model_label(m)),
              text = c(if (isTRUE(r$standardise)) "Estimates on the standardised scale; 'per original unit' columns divide by the scale factor." else "Coefficients on the original scale.",
                       safe_get(interpret_model_terms(r, m))),
              tables = list(`Fixed effects` = tab, `Scaling constants` = if (!is.null(sc)) data.frame(Variable = sc$Variable, Centre = signif(sc$Centre, 5), Scale = signif(sc$Scale, 5)) else NULL,
                            `Random effects` = safe_get(coef_re_tab())))
  }))
  observeEvent(input$save_family, disappr_guard("input$save_family", {
    fr <- family_res()
    if (is.null(fr) || !isTRUE(fr$ok)) return(nothing_to_save("Run the error-family check first."))
    add_saved("5 Modelling", paste("Error-family check with", model_label(fr$model)), tables = list(`Count families` = fr$table))
  }))
  observeEvent(input$save_dharma, disappr_guard("input$save_dharma", {
    z <- dharma_res()
    if (is.null(z) || !isTRUE(z$ok)) return(nothing_to_save("Simulate residuals first."))
    sim <- z$sim
    add_saved("5 Modelling", paste("Residual checks for", model_label(z$model)),
              tables = list(`DHARMa tests` = data.frame(Test = z$table$Test, P = vapply(z$table$P_value, format_p, character(1)))),
              plots = list(function() graphics::plot(sim)))
  }))
  observeEvent(input$save_a5, disappr_guard("input$save_a5", {
    pl <- safe_get(a5_plot_obj())
    if (is.null(pl)) return(nothing_to_save())
    add_saved("2 Visual diagnosis", "Trait before death", text = safe_get(a5_lines()), plots = list(pl))
  }))
  observeEvent(input$save_a6, disappr_guard("input$save_a6", {
    pl <- safe_get(a6_plot_obj())
    s <- safe_get(a6_data())
    if (is.null(pl) || is.null(s) || !nrow(s)) return(nothing_to_save())
    add_saved("2 Visual diagnosis", "Selection differentials by age", text = safe_get(a6_lines()),
              tables = list(`Selection differentials` = data.frame(Age = s$Age, N = s$N, Disappearing = s$Disappearing,
                                                                   Differential = signif(s$Differential, 3), SE = signif(s$Differential_SE, 3),
                                                                   `Survivors - disappearing` = signif(s$Contrast, 3), check.names = FALSE)),
              plots = list(pl))
  }))
  observeEvent(input$save_a7, disappr_guard("input$save_a7", {
    pl <- safe_get(a7_plot_obj())
    if (is.null(pl)) return(nothing_to_save())
    add_saved("2 Visual diagnosis", "Disappearance hazard", text = safe_get(a7_lines()), plots = list(pl))
  }))
  observeEvent(input$save_performance, disappr_guard("input$save_performance", {
    z <- perf_res()
    if (is.null(z) || !isTRUE(z$ok)) return(nothing_to_save("Run the performance checks first."))
    cm <- z$check_model
    add_saved("5 Modelling", paste("performance checks for", model_label(z$model)), text = if (nzchar(z$plot_message %||% "")) z$plot_message else NULL,
              tables = list(`performance checks` = z$table), plots = if (!is.null(cm)) list(draw_check_model(cm)) else list())
  }))
  observeEvent(input$save_code, disappr_guard("input$save_code", {
    code <- safe_get(code_text())
    if (is.null(code)) return(nothing_to_save("Fit models first."))
    add_saved("5 Modelling", "Reproducible R code", code = code)
  }))

  overview_text <- reactive({
    z <- summary_parts()
    ms <- z$ms
    lk <- life_known()
    c(sprintf("Data: %s rows from %s individuals; %d distinct ages.", format(nrow(z$d), big.mark = ","),
              format(length(unique(z$d$id)), big.mark = ","), length(unique(z$d$age))),
      if (nrow(z$warn_rows)) paste0("Integrity warning \u2014 ", z$warn_rows$Check, ": ", z$warn_rows$Result),
      z$integ$family$text,
      paste(if (lk) "LS:" else "ALR:", a2_interpretation(z$a2, if (lk) "LS" else "ALR")),
      if (!is.null(z$a2_afr)) paste("AFR:", a2_interpretation(z$a2_afr, "AFR")),
      sprintf("Sampling: %.1f%% of expected occasions missed; p = %.2f. %s", ms$percent, ms$p, ms$type),
      ms$guidance,
      if (nrow(z$a3)) paste0("Individual fits: best mean \u0394AICc = ", z$a3$Function[[1]], "."),
      if (z$models_ok) c(sprintf("Models: %s; lowest AIC among eligible fits: %s (%s).", family_label(z$r$family), z$r$aic$Model[[1]], z$r$aic$Fit[[1]]), z$rec)
      else "Models not fitted for the current data and settings.")
  })
  observeEvent(input$save_overview, disappr_guard("input$save_overview", {
    txt <- safe_get(overview_text())
    if (is.null(txt)) return(nothing_to_save())
    add_saved("6 Summary and report", "Automatic overview", text = txt)
  }))

  output$saved_controls <- renderUI({
    lst <- saved_results()
    if (!length(lst)) return(p(class = "small-note", "No saved results yet."))
    ch <- stats::setNames(names(lst), vapply(lst, function(e) paste0(e$section, " \u00b7 ", e$title, " (", e$time, ")"), character(1)))
    tagList(
      selectInput("saved_remove_id", "Remove a saved result", choices = ch),
      actionButton("saved_remove", "Remove", icon = icon("trash"), class = "btn-default btn-sm"),
      actionButton("saved_clear", "Clear all", icon = icon("broom"), class = "btn-default btn-sm")
    )
  })
  observeEvent(input$saved_remove, disappr_guard("input$saved_remove", {
    lst <- saved_results()
    id <- input$saved_remove_id
    if (!is.null(id) && id %in% names(lst)) {
      lst[[id]] <- NULL
      saved_results(lst)
    }
  }))
  observeEvent(input$saved_clear, disappr_guard("input$saved_clear", saved_results(list())))

  output$saved_ui <- renderUI({
    lst <- saved_results()
    if (!length(lst)) {
      return(div(class = "truth-card", "Nothing saved yet. Use the 'Save to summary' buttons on the other tabs; saved results appear here and make up the downloaded reports."))
    }
    tagList(lapply(lst, function(e) {
      div(class = "summary-block",
          h4(paste0(e$section, " \u00b7 ", e$title)),
          p(class = "small-note", paste0("Saved ", e$time, " \u00b7 ", e$data)),
          if (length(e$text)) tags$ul(lapply(e$text, tags$li)) else NULL,
          tagList(lapply(seq_along(e$tables), function(j) {
            cap <- names(e$tables)[j]
            tagList(if (!is.null(cap) && !is.na(cap) && nzchar(cap)) h5(strong(cap)) else NULL,
                    tableOutput(paste0(e$id, "_table", j)))
          })),
          tagList(lapply(seq_along(e$plots), function(j) plotOutput(paste0(e$id, "_plot", j), height = 380))),
          if (!is.null(e$code)) tags$pre(class = "code-out", e$code) else NULL)
    }))
  })

  output$download_report_html <- downloadHandler(
    filename = function() paste0("disappR_report_", Sys.Date(), ".html"),
    content = function(file) {
      con <- file(file, open = "w", encoding = "UTF-8")
      on.exit(close(con))
      writeLines(html_report(saved_results(), "disappR report: saved results"), con)
    }
  )
  output$download_report <- downloadHandler(
    filename = function() paste0("disappR_report_", Sys.Date(), ".txt"),
    content = function(file) {
      con <- file(file, open = "w", encoding = "UTF-8")
      on.exit(close(con))
      writeLines(text_report(saved_results(), "disappR report: saved results"), con)
    }
  )
}
