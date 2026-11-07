script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
app_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(app_dir)
shiny::runApp(app_dir, port = 3840, launch.browser = TRUE)
