
ui <- dashboardPage(
  title = "disappR \u00b7 selective disappearance and ageing",
  skin = "blue",
  dashboardHeader(
    title = tags$span(tags$img(src = "logo.png", height = "32px", style = "margin:-4px 8px 0 0;"), "disappR"),
    titleWidth = 250
  ),
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = "tabs",
      menuItem("Start here", tabName = "overview", icon = icon("compass")),
      menuItem("1 \u00b7 Data", tabName = "data", icon = icon("table")),
      menuItem("2 \u00b7 Visual diagnosis", tabName = "visual", icon = icon("eye")),
      menuItem("3 \u00b7 Missingness and proxies", tabName = "sampling", icon = icon("table-cells")),
      menuItem("4 \u00b7 Individual and population trajectories", tabName = "individual", icon = icon("chart-line")),
      menuItem("5 \u00b7 Modelling", tabName = "models", icon = icon("layer-group")),
      menuItem("6 \u00b7 Summary and report", tabName = "summary", icon = icon("square-check")),
      menuItem("Help and videos", tabName = "help", icon = icon("circle-question"))
    ),
    tags$hr(),
    div(class = "sidebar-note", strong("Diagnose visually, then model."), br(),
        "Visual and descriptive diagnostics first, then sampling and ageing functions, then the comparative models."),
    tags$hr(),
    div(class = "sidebar-cite", strong("Citation"),
        tags$ol(
          tags$li("Sanghvi, K., Ivimey-Cook, E. 2026. disappR: a shiny app to model ageing and selective [dis]appearance."),
          tags$li("Sanghvi, K., Ivimey-Cook, E.R., Bouwhuis, S., Sepil, I. and van de Pol, M., 2026. A comparison of methods to assess selective disappearance and quantify ageing. ",
                  tags$em("EcoEvoRxiv"))))
  ),
  dashboardBody(
    tags$head(tags$script(HTML("function disapprCopy(id){var el=document.getElementById(id);if(!el){return;}var txt=el.innerText||el.textContent;if(navigator.clipboard&&window.isSecureContext){navigator.clipboard.writeText(txt);}else{var ta=document.createElement('textarea');ta.value=txt;document.body.appendChild(ta);ta.select();try{document.execCommand('copy');}catch(e){}document.body.removeChild(ta);}}")), tags$style(HTML("
      .section-code pre { max-height: 440px; overflow: auto; font-size: 11.5px; }
      .code-toolbar { display: flex; gap: 6px; margin-bottom: 6px; }
      .model-eq { font-family: 'Cambria Math', 'Times New Roman', serif; font-size: 15px; line-height: 1.8; margin: 4px 0 8px 0; overflow-x: auto; }
      .coef-eq td { vertical-align: top; }
      body, .content-wrapper, .right-side { font-family: Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; }
      .content-wrapper, .right-side { background:#F2EADF !important; }
      .main-header .logo { background:#685243 !important; color:#FFF9F2 !important; font-weight:700; }
      .model-row { display:flex; align-items:center; gap:8px; border-bottom:1px solid #EFE6DA; padding:2px 0; }
      .model-row-check { flex:0 0 42%; } .model-row-check .checkbox { margin:4px 0; }
      .model-row-extra { flex:1 1 auto; } .model-row-extra .form-group { margin:3px 0; } .model-row-extra .selectize-control { margin:0; }
      .main-header .navbar { background:#91634E !important; }
      .skin-blue .main-sidebar { background:#493D34 !important; }
      .skin-blue .sidebar-menu>li>a { color:#EFE4D6 !important; border-left:3px solid transparent; }
      .skin-blue .sidebar-menu>li:hover>a, .skin-blue .sidebar-menu>li.active>a { background:#5B4B40 !important; color:#FFF9F2 !important; border-left-color:#D6A15E !important; }
      .sidebar-note { padding:2px 18px 18px; color:#DCCFC0; font-size:12px; line-height:1.55; }
      .sidebar-cite { padding:2px 18px 18px; color:#DCCFC0; font-size:11.5px; line-height:1.45; }
      .purpose .purpose-lead { color:#4D4036; font-size:14.5px; line-height:1.6; margin:0 0 8px; }
      .purpose ul { padding-left:18px; margin:4px 0 10px; }
      .purpose li { color:#5E4E42; font-size:13.5px; line-height:1.55; margin-bottom:6px; }
      .sidebar-cite ol { padding-left:16px; margin:6px 0 0; }
      .sidebar-cite li { margin-bottom:7px; }
      .main-sidebar hr { border-color:rgba(255,255,255,.10); margin:14px 18px; }
      .box { border-radius:12px; border-top:0 !important; box-shadow:0 3px 14px rgba(79,60,44,.08); background:#FFFDF9; }
      .box.box-solid.box-primary>.box-header, .box.box-solid.box-info>.box-header, .box.box-solid.box-warning>.box-header { background:#E5D2BE !important; color:#463B33 !important; }
      .box-title { font-weight:700 !important; }
      .box-body { overflow-x:auto; }
      .metric-card { background:#FFFDF9; border:1px solid #E6D8C8; border-radius:12px; padding:13px 15px; min-height:96px; margin-bottom:12px; }
      .metric-label { color:#8A7564; font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:.06em; }
      .metric-value { color:#483D35; font-size:24px; font-weight:700; margin-top:4px; line-height:1.1; }
      .metric-sub { color:#78685B; font-size:12px; margin-top:5px; line-height:1.4; }
      .diagnosis-card { border-left:5px solid #B56C54; background:#F8EEE4; padding:13px 16px; border-radius:10px; margin-bottom:12px; }
      .diagnosis-title { font-size:16px; font-weight:700; color:#5A4033; }
      .diagnosis-detail { margin-top:6px; color:#66564A; line-height:1.5; }
      .truth-card { background:#EFE8D8; border:1px solid #D8C8AC; border-radius:10px; padding:12px 14px; color:#5B4C3F; line-height:1.5; margin:10px 0; }
      .truth-card ul { padding-left:18px; margin:6px 0 0; }
      .small-note { color:#7B6B5D; font-size:12.5px; line-height:1.5; }
      .section-lead { color:#66574B; font-size:14px; line-height:1.6; margin-bottom:12px; }
      .step-card { background:#FAF3EA; border:1px solid #E8DACB; border-radius:12px; padding:12px 14px; margin-bottom:10px; }
      .step-kicker { font-size:11px; text-transform:uppercase; letter-spacing:.08em; font-weight:700; color:#A26750; }
      .step-title { font-weight:700; color:#4D4036; margin:2px 0; }
      .step-text { color:#756458; font-size:12.5px; line-height:1.45; }
      .summary-block { background:#FFFDF9; border:1px solid #E5D6C6; border-radius:12px; padding:14px 16px; margin-bottom:12px; }
      .summary-block h4 { margin-top:0; color:#58463A; font-weight:700; }
      .btn-primary { background:#A8644E !important; border-color:#A8644E !important; color:white !important; }
      .btn-primary:hover { background:#8E5442 !important; }
      .btn-default { background:#F4E8DA !important; border-color:#D8C3AF !important; color:#5E4A3B !important; }
      .btn-block-space { margin:6px 0; width:100%; }
      table { color:#584B41; }
      th, td { white-space:normal !important; }
      pre.code-out { background:#2E2620; color:#F3E9DC; border-radius:10px; font-size:12px; max-height:420px; overflow:auto; }
      details summary { cursor:pointer; color:#7B5A45; font-weight:600; margin:8px 0; }
      .info-link { margin-left:7px; color:#A8644E !important; font-size:15px; }
      .info-link:hover { color:#6D3F31 !important; }
      .save-btn { margin-top:8px; }
      .saved-controls { background:#FAF3EA; border:1px solid #E8DACB; border-radius:10px; padding:10px 12px; margin-bottom:12px; }
    "))),

    tabItems(
      # ---------------------------------------------------------------- overview
      tabItem(tabName = "overview",
        fluidRow(
          box(width = 12, title = "Purpose", status = "primary", solidHeader = TRUE,
            div(class = "purpose",
              p(class = "purpose-lead",
                "Longitudinal studies of ageing do not follow a fixed set of individuals: some die or leave early (selective disappearance) and some are first recorded late (selective appearance). When these individuals differ in their trait, the ageing pattern seen in the data can differ markedly from how individuals actually age."),
              tags$ul(
                tags$li(strong("Pervasive, but the remedies are recent. "),
                        "Selective disappearance is found across traits and taxa, yet the methods that account for it date from the last two decades and were built for simple cases: linear ageing, complete sampling, and lifespan linked only to an individual's average trait level."),
                tags$li(strong("Missing data bias some methods. "),
                        "With incomplete sampling, mean age becomes a poor proxy of lifespan, and the decomposition runs short of individuals sampled at successive occasions."),
                tags$li(strong("Age-dependent selection is widely overlooked. "),
                        "Lifespan, or the age at first record, can be linked to how fast or in what shape individuals age, not only to their average level: longer-lived individuals may senesce more slowly, and late starters may age differently from early ones. The bias then grows with age, and additive lifespan terms, mean-age centring, random slopes and the decomposition do not remove it, yet most analyses assume without testing that selection acts on the level alone."),
                tags$li(strong("Many options, little guidance. "),
                        "It is still unclear which of the many available models to fit, and how to interpret them, when selection is present.")),
              p(class = "purpose-lead",
                "disappR diagnoses these problems in your data first, visually and descriptively, checks sampling and missingness, and then fits and compares the models side by side, with simulations whose true answer is known.")))
        ),
        fluidRow(
          box(width = 8, title = info_title("A visual-first workflow for selective disappearance and appearance", "overview"),
              status = "primary", solidHeader = TRUE,
            div(class = "section-lead",
                "Longitudinal ageing trajectories are biased when individuals that die (or enter) early differ from the others in their trait values or in how the trait changes with age. The tabs follow a diagnose-then-model workflow: first diagnose selective disappearance and appearance visually and descriptively, then check sampling and ageing shapes, then fit and compare models. Click the ", icon("circle-info"), " symbols for help on any section."),
            fluidRow(
              column(6,
                div(class = "step-card", div(class = "step-kicker", "1 \u00b7 Data"),
                    div(class = "step-title", "Load, map, subset and check your data"),
                    div(class = "step-text", "Use simulated data (known answer), one of the published empirical datasets, or your own CSV. Map continuous and categorical covariates and their interactions, and check that ages, lifespans, AFR (age at first observation) and censoring are read correctly.")),
                div(class = "step-card", div(class = "step-kicker", "2 \u00b7 Visual diagnosis"),
                    div(class = "step-title", "Diagnose selective disappearance and appearance visually"),
                    div(class = "step-text", "Trait trajectories of lifespan (ALR, LS, mean age) or AFR groups and their differences across age, the trait against lifespan within age bins, the trait before death, selection differentials by age and the disappearance hazard. Constant differences suggest age-independent selection; differences that change with age suggest age-dependent selection.")),
                div(class = "step-card", div(class = "step-kicker", "3 \u00b7 Missingness and proxies"),
                    div(class = "step-title", "How complete is sampling, and which lifespan proxy to use?"),
                    div(class = "step-text", "Sampling grid, missingness by age and against other variables, and agreement between ALR, mean age and lifespan. Guides the choice between ALR and mean age in the models."))
              ),
              column(6,
                div(class = "step-card", div(class = "step-kicker", "4 \u00b7 Individual and population trajectories"),
                    div(class = "step-title", "Diagnose ageing shapes visually, then compare functions"),
                    div(class = "step-text", "Fits a function (linear to exponential) to each individual to see how individuals age, reconstructs the population curve from individual fits, and compares ageing functions at the population level. Choose the ageing function for the mixed models here.")),
                div(class = "step-card", div(class = "step-kicker", "5 \u00b7 Modelling"),
                    div(class = "step-title", "Quantify selective disappearance and appearance"),
                    div(class = "step-text", "Mixed models with and without lifespan and AFR terms, random slopes, count families, fit validity, population-level predictions, coefficients in original units, residual diagnostics and downloadable R code for the models and figures.")),
                div(class = "step-card", div(class = "step-kicker", "6 \u00b7 Summary and report"),
                    div(class = "step-title", "Collect and export results"),
                    div(class = "step-text", "Every 'Save to summary' result is gathered here with an automatic, cautious overview, and exported as an HTML or text report."))
              )
            )
          ),
          box(width = 4, title = info_title("Quick start", "quickstart"), status = "info", solidHeader = TRUE,
            p("Start with a simulated dataset where the answer is known and the patterns are deliberately obvious."),
            actionButton("go_toy", "Simulated teaching data", class = "btn-primary btn-block-space", icon = icon("flask")),
            p(style = "margin-top:10px;", "Or explore published data: ten empirical examples from laboratory systems (fruit flies, seed beetles, leafcutting bees) and wild populations (common terns, painted turtles, great tits, eastern chipmunks, Soay sheep), each pre-set to the analysis reported in its paper."),
            actionButton("go_fly", "Empirical examples (published data)", class = "btn-default btn-block-space", icon = icon("bug")),
            tags$hr(),
            p(class = "small-note", "Use 'Save to summary' on any result you want in the exported report (tab 7)."),
            uiOutput("package_status")
          )
        )
      ),

      # ---------------------------------------------------------------- data
      tabItem(tabName = "data",
        fluidRow(
          box(width = 4, title = info_title("Data source", "data_source"), status = "primary", solidHeader = TRUE,
            radioButtons("data_source", NULL, choices = c("Simulated teaching data" = "toy", "Empirical examples (published data)" = "example", "Upload my CSV" = "upload"), selected = "toy"),
            conditionalPanel("input.data_source == 'toy'",
              selectInput("toy_trait", "Trait and distribution", choices = TOY_TRAITS, selected = "mass"),
              conditionalPanel("input.toy_trait != 'paper'",
                selectInput("toy_form", "Ageing form", choices = AGE_FUNCTIONS, selected = "Quadratic"),
                selectInput("toy_shape", "Shape within this form (biological scenario)", choices = toy_shape_choices("Quadratic"), selected = "default"),
                radioButtons("toy_strength", "Pattern strength", inline = TRUE,
                             choices = c("Dramatic (teaching)" = "dramatic", "Moderate" = "moderate"), selected = "dramatic"),
                sliderInput("toy_mean_ls", "Mean lifespan (number of sampling occasions)", min = 3, max = 30, value = 20, step = 1),
                radioButtons("toy_rate_var", "Individual differences in ageing rates", inline = TRUE,
                             choices = c("Small" = "low", "Large" = "high"), selected = "low")
              ),
              radioButtons("toy_sd_type", "Selective disappearance", choices = SELECTION_TYPES, selected = "both"),
              conditionalPanel("input.toy_sd_type != 'none'",
                radioButtons("toy_sd_dir", "Direction", inline = TRUE,
                             choices = c("Positive: longer-lived individuals have higher values" = "1",
                                         "Negative: lower values" = "-1"), selected = "1")),
              box_note("Age-independent: lifespan is linked to the level of the trajectory. Age-dependent: lifespan is linked to the rate of ageing."),
              selectInput("toy_missingness", "Sampling design", choices = c(
                "Complete sampling" = "complete", "Missing completely at random (MCAR, ~25% missing)" = "mcar",
                "Missing when old (MWO)" = "mwo", "Missing when young (MWY)" = "mwy",
                "Trait-dependent missingness" = "trait", "Condition-dependent missingness" = "condition"), selected = "complete"),
              radioButtons("toy_afr_mode", "Age at first observation (AFR)",
                           choices = c("Same for everyone (age 1)" = "same", "Individual-specific" = "individual"), selected = "same"),
              conditionalPanel("input.toy_afr_mode == 'individual'",
                radioButtons("toy_sa_type", "Selective appearance",
                             choices = c("None (AFR varies at random)" = "none", "Age-independent" = "independent",
                                         "Age-dependent" = "dependent", "Age-independent and age-dependent" = "both"),
                             selected = "none"),
                conditionalPanel("input.toy_sa_type != 'none'",
                  radioButtons("toy_sa_dir", "Direction", inline = TRUE,
                               choices = c("Positive: later starters have higher values" = "1",
                                           "Negative: lower values" = "-1"), selected = "1"))),
              checkboxInput("toy_diet", "Add a diet covariate (poor diet lowers the trait; no effect on lifespan)", FALSE),
              checkboxInput("toy_groups", "Add families (nested random effect: individuals within 25 families)", FALSE),
              fluidRow(
                column(6, numericInput("toy_n", "Individuals", value = 300, min = 30, max = 2000, step = 50)),
                column(6, numericInput("toy_seed", "Seed", value = 1, min = 1, step = 1))
              ),
              actionButton("toy_commit", "Simulate & use this dataset", class = "btn-primary btn-block-space", icon = icon("play")),
              uiOutput("toy_status"),
              downloadButton("download_toy", "Download simulated CSV", class = "btn-default btn-block-space")
            ),
            conditionalPanel("input.data_source == 'example'",
              selectInput("example_id", "Published dataset", choices = EXAMPLE_CHOICES, selected = "fly"),
              uiOutput("example_note"),
              box_note("Column mapping, covariates, nesting, random effects, the error family and the model selection are pre-filled to match the published analysis. Change any of them to explore. How each bundled file relates to its archived original is documented in data/PROVENANCE.md, and REPLICATION.md in the package root sets each example against its published analysis.")
            ),
            conditionalPanel("input.data_source == 'upload'",
              fileInput("data_file", "CSV file (one row per individual \u00d7 age)", accept = ".csv"),
              box_note("Required: individual ID, age (numeric: whole numbers or decimals, not categories), trait. Optional: known lifespan, ALR, AFR, condition, covariates, grouping and random-effect variables, a censoring indicator. Missing values can be blank cells or NA (both are read as missing, as are '.', '-', 'NaN', 'N/A', '#N/A' and 'NULL').")
            ),
            tags$hr(),
            radioButtons("dup_action", "Repeated ID \u00d7 age records", inline = TRUE,
                         choices = c("Keep all" = "keep", "Average them" = "mean"), selected = "keep")
          ),
          box(width = 8, title = info_title("Map columns", "mapping"), status = "primary", solidHeader = TRUE,
            uiOutput("mapping_ui"),
            box_note(strong("ALR"), " = age at last record (automatic by default). ",
                     strong("LS"), " = known lifespan; an automatic last-age LS is only a proxy and cannot serve as the Model 6 positive control. ",
                     strong("Nesting"), " identifies individuals by group + ID, i.e. (1 | group) + (1 | group:ID); with a next level up (e.g. individual within father within family) it adds (1 | family) and uses (1 | family:father) + (1 | family:father:ID). ",
                     strong("Additional random intercepts"), " enter every model as + (1 | X). ",
                     strong("Censored"), " individuals keep their records but their LS is treated as unknown.")
          )
        ),
        fluidRow(
          box(width = 12, title = info_title("Subset the data (optional)", "subset"), status = "info", solidHeader = TRUE, collapsible = TRUE,
              uiOutput("subset_ui"))
        ),
        fluidRow(uiOutput("data_metrics")),
        fluidRow(
          box(width = 12, title = info_title("Distribution of a variable", "distribution"), status = "info", solidHeader = TRUE, collapsible = TRUE,
            fluidRow(
              column(3, uiOutput("dist_var_ui"),
                     radioButtons("dist_unit", "One value per", choices = c("Row" = "row", "Individual (mean or most common value)" = "individual"),
                                  selected = "row"),
                     save_button("save_distribution")),
              column(9, plotOutput("dist_plot", height = 320))
            ))
        ),
        fluidRow(box(width = 12, title = "Preview (first 10 rows)", status = "info", solidHeader = TRUE, collapsible = TRUE,
                     tableOutput("data_preview"))),
        fluidRow(
          box(width = 12, title = info_title("Data integrity checks", "integrity"), status = "warning", solidHeader = TRUE,
              tableOutput("integrity_table"), uiOutput("integrity_extra"), save_button("save_integrity"))
        ),
        section_code_box("data")
      ),

      # ---------------------------------------------------------------- visual
      tabItem(tabName = "visual",
        fluidRow(
          box(width = 12, title = info_title("Settings for the lifespan-group figures", "visual_settings"), status = "primary", solidHeader = TRUE,
            fluidRow(
              column(3, uiOutput("visual_proxy_ui")),
              column(3, sliderInput("n_bins", "Number of bins (lifespan groups; age bins)", min = 3, max = 6, value = 4, step = 1)),
              column(3, uiOutput("facet_ui")),
              column(3, radioButtons("trait_scale", "Trait scale", inline = TRUE, choices = c("Raw" = "raw", "log(trait + 1)" = "log1p")))
            ),
            box_note("These settings apply to both lifespan-group figures below. Every bin, bin \u00d7 age point and age bin needs at least 3 individuals; with no more distinct values than bins, each value is its own bin. Choose AFR to diagnose selective appearance."),
            uiOutput("toy_card_visual")
          )
        ),
        fluidRow(
          box(width = 12, title = info_title("Trait trajectory within bins", "a1"), status = "primary", solidHeader = TRUE,
            fluidRow(
              column(4, radioButtons("bin_method", "Bin boundaries", inline = TRUE, choices = c("Equal width" = "equal", "Quantiles" = "quantile"))),
              column(3, checkboxInput("show_se", "Show \u00b11 SE", FALSE)),
              column(5, uiOutput("a1_facet_ui"))
            ),
            plotOutput("a1_plot", height = 430),
            box_note("Bins are formed among individuals (each individual counts once); points are means of individual \u00d7 age means, sized by the number of individuals. 'Panels by' splits this figure, and the bin differences below, by a categorical variable."),
            save_button("save_a1")
          )
        ),
        fluidRow(
          box(width = 6, title = info_title("Difference between consecutive bins at each age", "a1_diff"), status = "primary", solidHeader = TRUE,
            radioButtons("diff_lines", "Trend lines", inline = TRUE, choices = c("One line per bin pair" = "pairs", "One line across all pairs" = "pooled"), selected = "pairs"),
            plotOutput("bin_diff_plot", height = 400), save_button("save_a1_diff")),
          box(width = 6, title = info_title("Trait against lifespan within age bins", "a2"), status = "primary", solidHeader = TRUE,
            checkboxInput("a2_points", "Show individual points (jittered)", TRUE),
            plotOutput("a2_plot", height = 400), save_button("save_a2"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Regression coefficient across consecutive age bins", "a2_table"), status = "primary", solidHeader = TRUE,
            uiOutput("a2_trend_note"),
            tableOutput("a2_slope_table"),
            save_button("save_a2_table"))
        ),
        fluidRow(
          box(width = 12, status = "info", solidHeader = FALSE,
              h4("Disappearance diagnostics"),
              div(class = "section-lead",
                  "These figures follow individuals up to their disappearance: the trait before death, whether survivors differ from those that disappear at each age, and how the chance of disappearing changes with age. Disappearance means death when lifespan is known, otherwise the last record; map a censoring column so that individuals alive at the end of the study are not counted as disappearing."),
              uiOutput("a4_status"))
        ),
        fluidRow(
          box(width = 6, title = info_title("Trait before death", "a5_terminal"), status = "primary", solidHeader = TRUE,
              plotOutput("a5_plot", height = 360), uiOutput("a5_note"), save_button("save_a5")),
          box(width = 6, title = info_title("Selection differentials by age", "a6_selection"), status = "primary", solidHeader = TRUE,
              radioButtons("a6_compare", "Compare survivors with", inline = TRUE, choices = c("All individuals present" = "population", "Individuals that disappear" = "contrast")),
              plotOutput("a6_plot", height = 320), uiOutput("a6_note"), save_button("save_a6"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Disappearance hazard: age-independent and age-specific", "a7_hazard"), status = "primary", solidHeader = TRUE,
              fluidRow(column(8, plotOutput("a7_plot", height = 340)), column(4, uiOutput("a7_note"))),
              save_button("save_a7"))
        ),
        section_code_box("visual")
      ),

      # ---------------------------------------------------------------- sampling
      tabItem(tabName = "sampling",
        fluidRow(column(12, uiOutput("sampling_caveat"))),
        fluidRow(uiOutput("sampling_metrics")),
        fluidRow(
          box(width = 8, title = info_title("Sampling grid (each row is an individual)", "heatmap"), status = "primary", solidHeader = TRUE,
              plotOutput("heatmap", height = "auto"), save_button("save_heatmap")),
          box(width = 4, title = info_title("Settings and proxy guidance", "sampling_guidance"), status = "info", solidHeader = TRUE,
              radioButtons("miss_start", info_title("Count missed occasions from", "afr_expression"),
                           choices = c("Age at first record (AFR)" = "afr",
                                       "Age at first trait expression (AFE)" = "afe"), selected = "afr"),
              conditionalPanel("input.miss_start == 'afe'",
                numericInput("afe_age", "Age at first trait expression (AFE)", value = NA)),
              box_note("Expected occasions end at the known lifespan (or the last record). This setting changes the missingness figures on this tab only: every model and statistic in the app uses AFR, as mapped on the Data tab."),
              sliderInput("heat_n", "Individuals to show in the grid", min = 20, max = 500, value = 150, step = 10),
              selectInput("heat_order", "Order individuals by", choices = c("ALR" = "alr", "AFR" = "afr", "Mean trait value" = "trait", "ID" = "id", "Random" = "random")),
              uiOutput("sampling_guidance"),
              save_button("save_sampling"))
        ),
        fluidRow(
          box(width = 6, title = info_title("Missingness and sample size by age", "missing_by_age"), status = "primary", solidHeader = TRUE,
              plotOutput("missing_by_age", height = 320), tableOutput("coverage_table"), save_button("save_missing_age")),
          box(width = 6, title = info_title("Missingness against other variables", "missing_vs_var"), status = "primary", solidHeader = TRUE,
              selectInput("miss_var", NULL, choices = c("ALR", "LS", "Condition", "Trait")),
              plotOutput("missing_vs_var", height = 260),
              h5(strong("Observed-data associations with missingness")), tableOutput("drivers_table"),
              uiOutput("missing_trait_note"),
              save_button("save_missing_var"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Agreement between ALR, mean age and lifespan", "proxy_agreement"), status = "warning", solidHeader = TRUE,
              uiOutput("proxy_plots_ui"),
              h5(strong(info_title("Agreement between AFR and ALR", "afr_alr"))),
              fluidRow(column(7, plotOutput("afr_alr_plot", height = 300)),
                       column(5, tableOutput("afr_alr_table"))),
              save_button("save_proxies"))
        ),
        section_code_box("sampling")
      ),

      # ---------------------------------------------------------------- A3 individual fits
      tabItem(tabName = "individual",
        fluidRow(
          box(width = 3, title = info_title("Individual fit settings", "a3_settings"), status = "primary", solidHeader = TRUE,
            selectInput("a3_function", "Individual ageing function", choices = A3_FUNCTIONS, selected = "Quadratic"),
            radioButtons("a3_recon", "Population curve from", choices = c("Mean of coefficients" = "coef", "Mean of individual functions" = "fun", "Both" = "both"),
                         selected = "both"),
            sliderInput("a3_n", "Individuals to draw", min = 5, max = 60, value = 20, step = 5),
            selectInput("a3_order", "Which individuals", choices = c("Random sample" = "random", "Longest ALR" = "longest", "Shortest ALR" = "shortest")),
            uiOutput("use_a3_function_ui"),
            uiOutput("a3_metrics"),
            uiOutput("a3_support")
          ),
          box(width = 9, title = info_title("Individual trajectories and fitted functions", "a3_fits"), status = "primary", solidHeader = TRUE,
              plotOutput("a3_plot", height = "auto"))
        ),
        fluidRow(
          box(width = 7, title = info_title("Mean-coefficient and mean-function trajectories", "a3_mean"), status = "info", solidHeader = TRUE,
              plotOutput("a3_mean_plot", height = 380), save_button("save_a3")),
          box(width = 5, title = info_title("Which function fits individuals best?", "a3_compare"), status = "info", solidHeader = TRUE,
              tableOutput("a3_compare_table"),
              box_note("\u0394AICc uses only individuals for which every function is estimable (N_common)."),
              save_button("save_a3_compare"))
        ),
        fluidRow(box(width = 12, title = info_title("Individual coefficient summary", "a3_coefs"), status = "info", solidHeader = TRUE, collapsible = TRUE,
                     tableOutput("a3_coef_table"))),
        fluidRow(
          box(width = 4, title = info_title("Population-level ageing function", "b2"), status = "primary", solidHeader = TRUE,
              div(class = "section-lead", "Evaluated at the population level, not per individual: each function is fitted as a mixed model with the Model 1 structure, using the error family, covariates and random effects chosen on the Modelling tab."),
              uiOutput("b2_settings"),
              checkboxGroupInput("b2_functions", "Functions to compare", choices = MODEL_FUNCTIONS, selected = AGE_FUNCTIONS, inline = TRUE),
              actionButton("run_functions", "Compare ageing functions", class = "btn-primary btn-block-space", icon = icon("play")),
              tableOutput("b2_table"),
              save_button("save_b2")),
          box(width = 8, title = "Fitted population trajectories", status = "primary", solidHeader = TRUE,
              plotOutput("b2_plot", height = 460))
        ),
        section_code_box("individual")
      ),

      # ---------------------------------------------------------------- models
      tabItem(tabName = "models",
        fluidRow(
          box(width = 12, title = info_title("Model settings", "model_settings"), status = "primary", solidHeader = TRUE,
            fluidRow(
              column(4,
                selectInput("model_family", "Error family", choices = MODEL_FAMILIES, selected = "gaussian"),
                uiOutput("family_hint"),
                conditionalPanel("input.model_family == 'zip' || input.model_family == 'zinb' || input.model_family == 'zinb1'", uiOutput("zi_ui")),
                conditionalPanel("input.model_family == 'binomial' || input.model_family == 'betabinomial'", uiOutput("trials_ui")),
                selectInput("model_age_function", "Ageing function", choices = c(AGE_FUNCTIONS, A3_NONLINEAR), selected = "Quadratic"),
                uiOutput("function_default_note"),
                selectInput("random_structure", "Random effects for individuals", choices = RANDOM_STRUCTURES, selected = "none"),
                uiOutput("random_support_note"),
                radioButtons("among_order", info_title("Among-individual terms (ALR, LS, AFR, mean age)", "among_order"),
                             choices = c("Linear (default)" = "linear", "Same polynomial order as the ageing function" = "same"), selected = "linear"),
                checkboxInput("standardise", "Standardise age and proxies (recommended; AIC unchanged)", TRUE),
                checkboxInput("include_invalid", "Include fits with invalid Hessians in the ranking (inspection only)", FALSE)
              ),
              column(8,
                h5(strong("Models to compare"), span(class = "small-note", " \u2014 tick models; click ", icon("circle-info"), " for what a model tests and its exact specification; add terms to a model in its menu, or click ", icon("wrench"), " to build a term from chosen terms joined by + or \u00d7")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M1", paste0("Model 1 \u00b7 ", MODEL_MEANING[["M1"]]$name), TRUE)),
                    actionLink("info_model_M1", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M1", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M1", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M2", paste0("Model 2 \u00b7 ", MODEL_MEANING[["M2"]]$name), TRUE)),
                    actionLink("info_model_M2", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M2", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M2", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M3", paste0("Model 3 \u00b7 ", MODEL_MEANING[["M3"]]$name), TRUE)),
                    actionLink("info_model_M3", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M3", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M3", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M4", paste0("Model 4 \u00b7 ", MODEL_MEANING[["M4"]]$name), TRUE)),
                    actionLink("info_model_M4", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M4", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M4", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M5", paste0("Model 5 \u00b7 ", MODEL_MEANING[["M5"]]$name), TRUE)),
                    actionLink("info_model_M5", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M5", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M5", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M6", paste0("Model 6 \u00b7 ", MODEL_MEANING[["M6"]]$name), FALSE)),
                    actionLink("info_model_M6", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M6", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M6", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M7", paste0("Model 7 \u00b7 ", MODEL_MEANING[["M7"]]$name), FALSE)),
                    actionLink("info_model_M7", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M7", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M7", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M8", paste0("Model 8 \u00b7 ", MODEL_MEANING[["M8"]]$name), FALSE)),
                    actionLink("info_model_M8", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M8", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M8", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M9", paste0("Model 9 \u00b7 ", MODEL_MEANING[["M9"]]$name), FALSE)),
                    actionLink("info_model_M9", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M9", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M9", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                div(class = "model-row",
                    div(class = "model-row-check", checkboxInput("use_M10", paste0("Model 10 \u00b7 ", MODEL_MEANING[["M10"]]$name), FALSE)),
                    actionLink("info_model_M10", icon("circle-info"), class = "info-link"),
                    div(class = "model-row-extra", selectizeInput("extra_M10", NULL, choices = EXTRA_TERM_KEYS, multiple = TRUE,
                                                                 options = list(placeholder = "add terms to this model (optional)"))),
                    actionLink("build_M10", icon("wrench"), class = "info-link", title = "Build a term: chosen terms joined by + or \u00d7 (up to three-way interactions)")),
                uiOutput("structure_note"),
                actionButton("fit_models", "Fit models", class = "btn-primary btn-block-space", icon = icon("play"))
              )
            )
          )
        ),
        fluidRow(
          box(width = 12, title = info_title("Model support", "model_support"), status = "info", solidHeader = TRUE,
            uiOutput("model_fit_note"),
            fluidRow(
              column(6, plotOutput("aic_plot", height = 300)),
              column(6, tableOutput("aic_table"))
            ),
            h5(strong("Nested likelihood-ratio tests")), tableOutput("lrt_table"),
            uiOutput("consistency_ui"),
            tags$details(tags$summary("Fitting status, convergence and dropped rows"), tableOutput("status_table"), uiOutput("drop_note")),
            save_button("save_models")
          )
        ),
        fluidRow(
          box(width = 12, title = info_title("Random-effect variances (check random slopes here)", "varcomp"), status = "info", solidHeader = TRUE, collapsible = TRUE,
              tableOutput("varcomp_table"))
        ),
        fluidRow(
          box(width = 7, title = info_title("Population-level ageing trajectories", "predictions"), status = "primary", solidHeader = TRUE,
            fluidRow(
              column(3, checkboxInput("show_observed", "Observed means", TRUE)),
              column(3, checkboxInput("show_raw_points", "Individual records (jittered)", FALSE)),
              column(3, checkboxInput("show_decomp", "Decomposition", TRUE)),
              column(3, checkboxInput("show_a3", "Reconstruction from individual fits", FALSE))
            ),
            fluidRow(column(8, uiOutput("pred_models_ui")), column(4, uiOutput("pred_by_ui"))),
            plotOutput("pred_plot", height = 430),
            uiOutput("decomp_note"),
            box_note("Random effects excluded; numeric covariates at their mean; factor covariates marginalised over observed level combinations (weighted by individuals); ALR, LS, AFR and mean age at their individual-level means; response scale. Failed fits are not drawn; Caution fits are dashed."),
            save_button("save_predictions")),
          box(width = 5, title = info_title("Accuracy against the simulated truth", "accuracy"), status = "info", solidHeader = TRUE,
            uiOutput("deviation_note"), tableOutput("deviation_table"))
        ),
        fluidRow(
          box(width = 6, title = info_title("Coefficients", "coefficients"), status = "primary", solidHeader = TRUE,
              uiOutput("coef_model_ui"), tableOutput("coef_table"),
              h5(strong("Scaling constants")), tableOutput("scaling_table"),
              h5(strong("Random effects")), tableOutput("coef_re_table"),
              h5(strong("What the selection terms suggest")), uiOutput("coef_interpretation"),
              box_note("Estimate per original unit = estimate / scale factor (per unit of the centred original variables). Covariate names are shown without the internal prefix."),
              save_button("save_coefs")),
          box(width = 6, title = info_title("Model definitions", "definitions"), status = "primary", solidHeader = TRUE, tableOutput("definition_table"))
        ),
        fluidRow(
          box(width = 6, title = info_title("Error-family check (count traits)", "family_check"), status = "warning", solidHeader = TRUE,
              selectInput("family_check_model", "Fit this model with each count family", choices = stats::setNames(MODEL_IDS, model_label(MODEL_IDS)), selected = "M4"),
              actionButton("run_family_check", "Compare count families", class = "btn-default btn-block-space", icon = icon("scale-balanced")),
              tableOutput("family_table"),
              save_button("save_family")),
          box(width = 6, title = info_title("Residual checks (DHARMa)", "dharma"), status = "warning", solidHeader = TRUE,
              uiOutput("dharma_model_ui"),
              actionButton("run_dharma", "Simulate residuals", class = "btn-default btn-block-space", icon = icon("stethoscope")),
              tableOutput("dharma_table"), plotOutput("dharma_plot", height = 300),
              save_button("save_dharma"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Model diagnostics (performance)", "performance"), status = "warning", solidHeader = TRUE, collapsible = TRUE,
            fluidRow(
              column(3, uiOutput("performance_model_ui"),
                     actionButton("run_performance", "Run performance checks", class = "btn-default btn-block-space", icon = icon("clipboard-check")),
                     save_button("save_performance")),
              column(9, tableOutput("performance_table"), uiOutput("performance_note"), plotOutput("performance_plot", height = 760))
            ))
        ),
        fluidRow(
          box(width = 12, title = info_title("Reproducible R code for this comparison", "rcode"), status = "info", solidHeader = TRUE, collapsible = TRUE,
              downloadButton("download_code", "Download R script", class = "btn-default"),
              save_button("save_code", "Save code to summary"),
              uiOutput("code_ui"))
        )
      ),

      # ---------------------------------------------------------------- summary
      tabItem(tabName = "summary",
        fluidRow(
          box(width = 12, title = info_title("Saved results (these are exported)", "saved"), status = "primary", solidHeader = TRUE,
            div(class = "saved-controls",
              fluidRow(
                column(5, uiOutput("saved_controls")),
                column(7,
                  downloadButton("download_report_html", "Download HTML report (with plots)", class = "btn-primary"),
                  downloadButton("download_report", "Download text report", class = "btn-default"),
                  downloadButton("download_code2", "Download R script", class = "btn-default"))
              )
            ),
            uiOutput("saved_ui"))
        ),
        fluidRow(
          box(width = 12, title = info_title("Automatic overview of the current data and settings (not exported unless saved)", "overview_auto"),
              status = "info", solidHeader = TRUE, collapsible = TRUE,
            uiOutput("summary_ui"),
            save_button("save_overview"))
        )
      ),

      # ---------------------------------------------------------------- help
      tabItem(tabName = "help",
        fluidRow(
          box(width = 8, title = "Video walkthroughs", status = "primary", solidHeader = TRUE,
            div(class = "section-lead", "Two short screen recordings walk through the whole workflow. Each opens in a popup; press Close or click outside it to stop."),
            fluidRow(lapply(names(HELP_VIDEOS), function(k) {
              v <- HELP_VIDEOS[[k]]
              column(6, div(class = "step-card video-card",
                div(class = "step-kicker", paste0("Video · ", v$minutes)),
                div(class = "step-title", v$title),
                div(class = "step-text", v$text),
                actionButton(paste0("play_", k), "Play video", class = "btn-primary btn-block-space", icon = icon("circle-play"))))
            }))
          ),
          box(width = 4, title = "More help", status = "info", solidHeader = TRUE,
            p("Every section has an ", icon("circle-info"), " symbol that explains what it shows and how to read it."),
            p("Data, code and the manuscript's supplementary material are on OSF: ",
              tags$a(href = "https://osf.io/kevnm/", target = "_blank", "osf.io/kevnm"), "."),
            tags$hr(),
            p(strong("Citation")),
            tags$ol(class = "small-note", style = "padding-left:18px; overflow-wrap:anywhere;",
              tags$li("Sanghvi, K., Ivimey-Cook, E. 2026. disappR: a shiny app to model ageing and selective [dis]appearance."),
              tags$li("Sanghvi, K., Ivimey-Cook, E.R., Bouwhuis, S., Sepil, I. and van de Pol, M., 2026. A comparison of methods to assess selective disappearance and quantify ageing. ",
                      tags$em("EcoEvoRxiv"))))
        )
      )
    )
  )
)
