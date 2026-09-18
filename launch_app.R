# Run the app without installing: open R in the package root (the folder with DESCRIPTION)
needed <- c("shiny", "shinydashboard", "ggplot2", "lme4")
optional <- c("glmmTMB", "lmerTest", "DHARMa")
missing <- c(needed, optional)[!vapply(c(needed, optional), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
shiny::runApp(file.path("inst", "app"), launch.browser = TRUE)
