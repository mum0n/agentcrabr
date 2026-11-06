# Launcher script for agentcrabr Shiny app
app_dir <- "C:/Users/Nathalie/projects/agentcrabr"
setwd(app_dir)
source(file.path(app_dir, "R", "model.R"))
shiny::runApp(file.path(app_dir, "app.R"), port = 3840, launch.browser = TRUE)
