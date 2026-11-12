#' Launch the interactive Shiny dashboard for agentcrabr
#'
#' Starts the Shiny web application providing a 4-tab interface for simulation
#' execution, demographic analysis, stock-recruitment calibration, and spatial
#' area-of-interest setup.
#'
#' @param port Integer port number for the local Shiny server (default 3840).
#' @param launch.browser Logical indicating whether to open the system browser.
#' @return Invisible Shiny app object.
#' @export
run_crab_app <- function(port = 3840, launch.browser = TRUE) {
  app_dir <- system.file("app", package = "agentcrabr")
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R"))) {
    app_dir <- getwd()
  }
  shiny::runApp(app_dir, port = port, launch.browser = launch.browser)
}

