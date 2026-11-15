# Step 7 - GAMMs (optional). The tab's interface, built by gamm_tab_ui() from ui.R. Every input and output id starts
# gamm_ (save buttons save_gamm_, help keys info_gamm_), so the tab shares nothing with steps 1-6 except the data.

gamm_tab_ui <- function() {
  tabItem(tabName = "gamm",
    page_title("Does a smooth model of age (GAMM) tell the same ageing story?"),
    guide_box("gamm",
      goal = "Estimate the ageing curve without choosing its shape (a smooth of age), and check whether step 5's conclusions hold. Uses step 1's data; changes nothing in steps 1\u20136.",
      start_with = "press Fit GAMM with the default settings and read the curve and summary.",
      then = "follow the GAMM guide on Model output: is a curve needed, do ALR or AFR matter (Selective disappearance), do groups differ (Groups), is the model adequate (Checks), and where does the trait change (Ageing rate).",
      tip = "Save results as you go; they appear in step 6's report and methods draft.",
      further = "fit a random slope or random smooth per individual before trusting an age-dependent proxy effect, or add a terminal-record effect."),
    fluidRow(box(width = 12, status = "info", solidHeader = FALSE, uiOutput("gamm_data_note"))),
    tabBox(width = 12, id = "gamm_tabs",
      # ------------------------------------------------------------------ Fit
      tabPanel("Fit", icon = icon("play"),
        fluidRow(
          box(width = 4, title = info_title("GAMM settings", "gamm_settings"), status = "primary", solidHeader = TRUE,
            selectInput("gamm_family", "Error family", GAMM_FAMILIES, "gaussian"),
            uiOutput("gamm_family_hint"),
            conditionalPanel("input.gamm_family == 'poisson' || input.gamm_family == 'binomial'",
              checkboxInput("gamm_olre", "Observation-level random effect (overdispersed counts or proportions)", FALSE)),
            numericInput("gamm_k", "Basis size k (maximum wiggliness)", 10, min = 3, max = 40, step = 1),
            selectInput("gamm_individual", "Differences between individuals", GAMM_INDIVIDUAL, "ri"),
            uiOutput("gamm_individual_note"),
            h5(strong(info_title("Selective disappearance and appearance (how do I choose?)", "gamm_choose"))),
            selectInput("gamm_proxy", "Lifespan proxy", c("Age at last reproduction or record (ALR)" = "ALR", "Known lifespan (LS)" = "LS", "None" = "none"), "ALR"),
            conditionalPanel("input.gamm_proxy != 'none'",
              selectInput("gamm_proxy_form", "\u2026entered as", GAMM_PROXY_FORMS[-1], "linear")),
            uiOutput("gamm_proxy_hint"),
            selectInput("gamm_afr_form", "Age at first reproduction or record (AFR), entered as", GAMM_PROXY_FORMS, "none"),
            checkboxInput("gamm_terminal", "Terminal-record effect (last record before death)", FALSE),
            uiOutput("gamm_group_ui"),
            tags$details(class = "advanced-block",
              tags$summary("Advanced options (usually leave as they are)"),
              selectInput("gamm_engine", "Fitting engine", GAMM_ENGINES, "auto"),
              checkboxInput("gamm_covars", "Include the covariates mapped on the Data tab", TRUE),
              checkboxInput("gamm_cov_smooth", "Numeric covariates as smooths", FALSE),
              checkboxInput("gamm_random", "Include the random intercepts mapped on the Data tab", TRUE),
              checkboxInput("gamm_ar1", "AR1 errors within individuals (Gaussian traits; uses bam)", FALSE)),
            actionButton("gamm_fit", "Fit GAMM", class = "btn-primary", icon = icon("play")),
            uiOutput("gamm_formula_ui")),
          box(width = 8, title = info_title("What does the fitted ageing trajectory look like?", "gamm_curve"), status = "info", solidHeader = TRUE,
            uiOutput("gamm_reading_ui"),
            uiOutput("gamm_cautions_ui"),
            uiOutput("gamm_next_ui"),
            fluidRow(
              column(4, selectInput("gamm_show_data", "Show the data as", c("Observed means at each age" = "means", "Observed means in age bins" = "bins",
                                                                          "Individual records" = "points", "Nothing" = "none"), "means")),
              column(3, conditionalPanel("input.gamm_show_data == 'bins'", numericInput("gamm_age_bins", "Number of age bins", 8, min = 3, max = 25, step = 1))),
              column(3, uiOutput("gamm_n_proxy_ui")),
              column(2, uiOutput("gamm_truth_ui"))),
            uiOutput("gamm_scale_ui"),
            plotOutput("gamm_curve_plot", height = 460),
            uiOutput("gamm_curve_note"),
            plotOutput("gamm_support_plot", height = 120),
            look_for("the overall shape and where the interval is narrow; ignore wiggles where few individuals remain (the bars)."),
            save_button("save_gamm_curve"))
        )
      ),
      # ------------------------------------------------------------------ Model output
      tabPanel("Model output", icon = icon("table-list"),
        fluidRow(
          box(width = 12, title = info_title("Reading this GAMM, step by step", "gamm_guide"), status = "primary", solidHeader = TRUE,
            uiOutput("gamm_guide_ui"))
        ),
        fluidRow(
          box(width = 12, title = info_title("What does each term of the GAMM estimate?", "gamm_output"), status = "primary", solidHeader = TRUE,
            uiOutput("gamm_output_reading"),
            fluidRow(
              column(4, h5(strong("Model statistics")), tableOutput("gamm_out_stats")),
              column(8,
                h5(strong("Smooth terms (fixed)")), tableOutput("gamm_out_smooth"),
                h5(strong("Parametric terms (fixed)")), tableOutput("gamm_out_param"))),
            p(class = "small-note", uiOutput("gamm_out_formula", inline = TRUE)),
            save_button("save_gamm_output"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Is a smooth needed? No age, straight line, polynomial or smooth", "gamm_poly"), status = "primary", solidHeader = TRUE,
            p(class = "small-note", "Does the smooth beat no age, linear, quadratic or cubic age (step 5's functions)? Compared by AIC on the same records."),
            actionButton("gamm_run_poly", "Compare: no age, linear, quadratic, cubic and smooth", class = "btn-default", icon = icon("scale-balanced")),
            uiOutput("gamm_engine_note_poly"),
            uiOutput("gamm_poly_reading"),
            tableOutput("gamm_poly_table"),
            save_button("save_gamm_poly"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Random effects", "gamm_random"), status = "primary", solidHeader = TRUE,
            tableOutput("gamm_out_random"),
            actionButton("gamm_run_re_test", "Test the random effects (mgcv; can take a while)", class = "btn-default", icon = icon("vial")),
            uiOutput("gamm_re_test_note"),
            tableOutput("gamm_re_test_table"))
        )
      ),
      # ------------------------------------------------------------------ Ageing rate
      tabPanel("Ageing rate", icon = icon("arrow-trend-down"),
        fluidRow(
          box(width = 12, title = info_title("Where does the trait change with age?", "gamm_deriv"), status = "info", solidHeader = TRUE,
            p(class = "small-note", "Where the 95% interval excludes zero, the trait changes significantly with age."),
            fluidRow(
              column(7, plotOutput("gamm_deriv_plot", height = 380)),
              column(5, h5(strong("Periods of significant change")), tableOutput("gamm_periods_table"),
                     h5(strong(info_title("Peak age and onset of the final decline", "peak_onset"))), tableOutput("gamm_peak_table"))),
            look_for("runs of red (significant decline) and blue (significant increase) points; the onset of senescence is where the final decline begins."),
            save_button("save_gamm_deriv"))
        )
      ),
      # ------------------------------------------------------------------ Selective disappearance
      tabPanel("Selective disappearance", icon = icon("user-minus"),
        fluidRow(
          box(width = 5, title = info_title("Does selective disappearance or appearance shape the curve?", "gamm_proxy"), status = "info", solidHeader = TRUE,
            h5(strong("Proxy terms in the fitted model")), tableOutput("gamm_proxy_terms_table"),
            actionButton("gamm_run_proxy", "Compare proxy forms (ML, same records)", class = "btn-default", icon = icon("scale-balanced")),
            uiOutput("gamm_engine_note_proxy"),
            uiOutput("gamm_proxy_verdict"),
            uiOutput("gamm_proxy_tables"),
            save_button("save_gamm_proxy")),
          box(width = 7, title = "How far do the proxy terms move the curve?", status = "info", solidHeader = TRUE,
            p(class = "small-note", "The curve without ALR, lifespan or AFR terms (a naive analysis) against the fitted model's curve. A gap means selective disappearance or appearance distorts the naive curve. A linear proxy often changes it little, because individual random intercepts already absorb level differences."),
            plotOutput("gamm_naive_plot", height = 380),
            look_for("whether the two curves differ in shape (age-dependent selective disappearance) or only in level (age-independent)."))
        )
      ),
      # ------------------------------------------------------------------ Groups
      tabPanel("Groups", icon = icon("people-group"),
        fluidRow(
          box(width = 12, title = info_title("Do groups age differently?", "gamm_groups"), status = "info", solidHeader = TRUE,
            uiOutput("gamm_groups_note"),
            tags$details(class = "info-more", open = NA, tags$summary("Which group structure? G, GS, GI, S and I explained"),
              div(class = "info-more-body", tags$table(class = "table table-condensed small-note",
                tags$thead(tags$tr(tags$th("Structure"), tags$th("What it assumes"), tags$th("In mgcv"), tags$th("Choose it when"))),
                tags$tbody(
                  tags$tr(tags$td(strong("G")), tags$td("One ageing curve for all groups; groups differ only in overall level."), tags$td(code("group + s(age)")), tags$td("Groups age the same way (as group + age in step 5).")),
                  tags$tr(tags$td(strong("GS")), tags$td("A common curve plus a deviation for each group; all deviations share one wiggliness."), tags$td(code('s(age) + s(age, group, bs = "sz")')), tags$td("Groups share a pattern and differ modestly; good for small groups.")),
                  tags$tr(tags$td(strong("GI")), tags$td("A common curve plus a deviation for each group, each deviation with its own wiggliness."), tags$td(code("s(age) + s(age, by = group, m = 1)")), tags$td("Groups share a pattern, but one departs from it more than another.")),
                  tags$tr(tags$td(strong("S")), tags$td("A separate curve per group, no common curve, all with the same wiggliness."), tags$td(code("s(age, by = group, id = 1)")), tags$td("No shared pattern, but similarly smooth curves.")),
                  tags$tr(tags$td(strong("I")), tags$td("Fully separate curves, each with its own wiggliness."), tags$td(code("s(age, by = group)")), tags$td("Groups may age in completely different ways (as age \u00d7 group in step 5)."))
                )),
                p(class = "small-note", "Compare all five (button below). Within 2 AIC, prefer the simpler structure (G, GS, S, GI, I). If G is within 2 AIC of the best, report one curve.")
              )),
            actionButton("gamm_run_groups", "Compare group structures (ML, same records)", class = "btn-default", icon = icon("layer-group")),
            uiOutput("gamm_engine_note_groups"),
            uiOutput("gamm_groups_verdict"),
            fluidRow(
              column(6, tableOutput("gamm_groups_table"),
                     h5(strong("Do the shapes differ from the reference group?")), tableOutput("gamm_diff_test_table"),
                     tableOutput("gamm_lrt_table")),
              column(6, plotOutput("gamm_diff_plot", height = 360))),
            save_button("save_gamm_groups"))
        )
      ),
      # ------------------------------------------------------------------ Checks
      tabPanel("Checks", icon = icon("stethoscope"),
        fluidRow(
          box(width = 12, title = info_title("Is the model adequate?", "gamm_check"), status = "info", solidHeader = TRUE,
            uiOutput("gamm_check_overall"),
            tableOutput("gamm_check_summary"),
            tags$details(class = "info-more", tags$summary("Details: basis dimension (k) check"),
              div(class = "info-more-body",
                p(class = "small-note", "k may be too low only when all three hold: edf close to k', k-index below 1 and a small p-value."),
                tableOutput("gamm_kcheck_table"))),
            tags$details(class = "info-more", tags$summary("Details: concurvity between fixed terms"),
              div(class = "info-more-body",
                p(class = "small-note", "How far each smooth can be mimicked by another term (0 = not at all, 1 = completely). Above 0.8, their separate effects are poorly identified."),
                tableOutput("gamm_concurvity_table"))),
            save_button("save_gamm_check"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Do the residuals look right?", "gamm_resid"), status = "info", solidHeader = TRUE,
            uiOutput("gamm_resid_note"),
            plotOutput("gamm_resid_plot", height = 520))
        )
      ),
      # ------------------------------------------------------------------ Code
      tabPanel("R code", icon = icon("code"),
        fluidRow(
          box(width = 12, title = info_title("Reproduce this GAMM in R", "gamm_code"), status = "info", solidHeader = TRUE,
            div(class = "code-toolbar",
                tags$button(type = "button", class = "btn btn-default btn-sm", onclick = "disapprCopy('gamm_code', this)", icon("copy"), " Copy"),
                downloadButton("download_gamm_code", "Download .R", class = "btn-default btn-sm"),
                downloadButton("download_gamm_data", "Download the model's data (.csv)", class = "btn-default btn-sm")),
            div(class = "section-code", verbatimTextOutput("gamm_code")))
        )
      )
    )
  )
}
