#' Launch the disappR workflow
#'
#' Starts the Shiny dashboard for diagnosing and modelling selective
#' disappearance and selective appearance in longitudinal ageing data.
#'
#' @param ... Arguments passed to [shiny::runApp()], such as `port`,
#'   `launch.browser` or `host`.
#'
#' @return Runs the Shiny application; does not return a value.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   run_app()
#' }
run_app <- function(...) {
  needed <- c("shiny", "shinydashboard", "ggplot2", "lme4")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "),
         ". Install them with install.packages() and try again.", call. = FALSE)
  }
  optional <- c(glmmTMB = "count families", lmerTest = "Gaussian p-values", DHARMa = "residual checks")
  absent <- names(optional)[!vapply(names(optional), requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    message("Optional packages not installed (", paste(paste0(absent, ": ", optional[absent]), collapse = "; "),
            "). install.packages(c(", paste0('"', absent, '"', collapse = ", "), "))")
  }
  app_dir <- system.file("app", package = "disappR")
  if (!nzchar(app_dir)) stop("Could not find the app directory. Re-install disappR and try again.", call. = FALSE)
  shiny::runApp(app_dir, ...)
}
