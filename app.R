library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(maps)
library(mapdata)

options(shiny.sanitize.errors = FALSE)

source("C:/Users/Nathalie/projects/agentcrabr/R/model.R")

ui <- fluidPage(
  titlePanel("agentcrabr - Crab Agent-Based Model (TethysBio-faithful)"),
  sidebarLayout(
    sidebarPanel(width = 4,
      div(class = "panel",
        h4("Dynamic Update"),
        fluidRow(
          column(6, numericInput("update_interval", "Interval (sec)", 2, min = 0.5, max = 60, step = 0.5)),
          column(6, checkboxInput("auto_update", "Auto-update", value = FALSE))
        ),
        div(style = "font-size: 11px; color: #64748b;", textOutput("update_status"))
      ),
      div(class = "panel",
        h4("Location Selection (click map)"),
        plotOutput("loc_map", height = "250px", click = "map_click"),
        fluidRow(
          column(6, numericInput("inp_lon", "Longitude", -63.5, min = -66.5, max = -57.5, step = 0.01)),
          column(6, numericInput("inp_lat", "Latitude", 44.5, min = 42.0, max = 48.5, step = 0.01))
        ),
        div(style = "font-size: 11px; color: #64748b;", textOutput("coord_display"))
      ),
      div(class = "panel",
        h4("Environmental Data"),
        tableOutput("env_table")
      ),
      div(class = "panel",
        h4("Model Parameters"),
        fluidRow(
          column(6, numericInput("inp_n_crabs", "N Crabs", 50, min = 1, max = 500, step = 1)),
          column(6, numericInput("inp_sim_time", "Sim Time (days)", 30, min = 1, max = 365, step = 1))
        ),
        fluidRow(
          column(6, numericInput("inp_energy_init", "Energy Init", 10, min = 0, max = 100, step = 0.1)),
          column(6, numericInput("inp_energy_max", "Energy Max", 20, min = 1, max = 100, step = 0.1))
        ),
        fluidRow(
          column(6, numericInput("inp_foraging_gain", "Foraging Gain", 0.3, min = 0, max = 10, step = 0.01)),
          column(6, numericInput("inp_resting_cost", "Resting Cost", 0.05, min = 0, max = 10, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_molting_cost", "Molting Cost", 3, min = 0, max = 20, step = 0.1)),
          column(6, numericInput("inp_molt_growth_male", "Molt Growth Male", 0.12, min = 0, max = 1, step = 0.01)),
          column(6, numericInput("inp_molt_growth_female", "Molt Growth Female", 0.08, min = 0, max = 1, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_foraging_rate", "Foraging Rate", 0.5, min = 0.01, max = 10, step = 0.01)),
          column(6, numericInput("inp_resting_rate", "Resting Rate", 0.3, min = 0.01, max = 10, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_molting_rate", "Molting Rate", 0.05, min = 0.001, max = 10, step = 0.001)),
          column(6, numericInput("inp_mortality_rate", "Mortality Rate", 0.02, min = 0.001, max = 1, step = 0.001))
        ),
        fluidRow(
          column(6, numericInput("inp_size_min", "Size Init Min (mm)", 20, min = 5, max = 100, step = 1)),
          column(6, numericInput("inp_size_max", "Size Init Max (mm)", 40, min = 5, max = 150, step = 1))
        ),
        fluidRow(
          column(6, numericInput("inp_legal_size", "Legal Size Male (mm)", 95, min = 50, max = 150, step = 1)),
          column(6, numericInput("inp_maturity_l50_m", "Maturity L50 Male (mm)", 91.8, min = 50, max = 150, step = 0.1)),
          column(6, numericInput("inp_maturity_slope_m", "Maturity Slope Male", 6.9957, min = 0, max = 50, step = 0.001))
        ),
        fluidRow(
          column(6, numericInput("inp_maturity_l50_f", "Maturity L50 Female (mm)", 51.6, min = 20, max = 100, step = 0.1)),
          column(6, numericInput("inp_maturity_slope_f", "Maturity Slope Female", 12.9812, min = 0, max = 50, step = 0.001)),
          column(6, numericInput("inp_max_size_m", "Max Size Male (mm)", 160, min = 100, max = 200, step = 1))
        ),
        fluidRow(
          column(6, numericInput("inp_max_size_f", "Max Size Female (mm)", 95, min = 50, max = 150, step = 1)),
          column(6, numericInput("inp_fishing_prob", "Fishing Probability", 0.05, min = 0, max = 1, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_growth_beta0_m", "Growth Beta0 Male", 0.30153, min = 0, max = 5, step = 0.001)),
          column(6, numericInput("inp_growth_beta0_f", "Growth Beta0 Female", 0.1931, min = 0, max = 5, step = 0.001))
        ),
        fluidRow(
          column(6, numericInput("inp_growth_beta1_m", "Growth Beta1 Male", 1.00027, min = 0, max = 5, step = 0.00001)),
          column(6, numericInput("inp_growth_beta1_f", "Growth Beta1 Female", 1.0372, min = 0, max = 5, step = 0.00001))
        ),
        fluidRow(
          column(6, numericInput("inp_growth_slp_m", "Growth Slope Male", 1.00027, min = 0, max = 5, step = 0.00001)),
          column(6, numericInput("inp_growth_slp_f", "Growth Slope Female", 1.0372, min = 0, max = 5, step = 0.00001))
        ),
        fluidRow(
          column(6, numericInput("inp_growth_beta2_m", "Growth Beta2 Male", -0.00643655, min = -1, max = 1, step = 0.000001)),
          column(6, numericInput("inp_growth_beta2_f", "Growth Beta2 Female", -0.017376, min = -1, max = 1, step = 0.000001))
        ),
        fluidRow(
          column(6, numericInput("inp_growth_beta3_m", "Growth Beta3 Male", 0.0017492, min = -1, max = 1, step = 0.000001)),
          column(6, numericInput("inp_growth_beta3_f", "Growth Beta3 Female", 0.005068, min = -1, max = 1, step = 0.000001))
        ),
        fluidRow(
          column(6, numericInput("inp_growth_temp_coeff", "Growth Temp Coeff", -0.004, min = -1, max = 1, step = 0.001)),
          column(6, numericInput("inp_growth_temp_ref", "Growth Temp Ref", 4.0, min = 0, max = 20, step = 0.1))
        ),
        fluidRow(
          column(6, numericInput("inp_fec_int_prim", "Fecund Intercept Primiparous", -0.7229, min = -5, max = 5, step = 0.001)),
          column(6, numericInput("inp_fec_slope_prim", "Fecund Slope Primiparous", 2.94, min = 0, max = 10, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_fec_int_multi", "Fecund Intercept Multiparous", 0.0158, min = -5, max = 5, step = 0.001)),
          column(6, numericInput("inp_fec_slope_multi", "Fecund Slope Multiparous", 2.5811, min = 0, max = 10, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_fec_sen_thresh", "Fecund Senescence Thresh (yr)", 3.0, min = 0, max = 10, step = 0.1)),
          column(6, numericInput("inp_fec_sen_floor", "Fecund Senescence Floor", 0.6, min = 0, max = 1, step = 0.01))
        ),
        fluidRow(
          column(6, numericInput("inp_fec_sen_rate", "Fecund Senescence Rate", 0.08, min = 0, max = 1, step = 0.01)),
          column(6, numericInput("inp_fec_min_clutch", "Min Clutch Size", 1000, min = 1, max = 10000, step = 1))
        ),
        actionButton("run_sim", "Run Simulation", class = "btn btn-primary btn-block", style = "margin-top:10px;")
      )
    ),
    mainPanel(width = 8,
      fluidRow(
        column(3, div(class = "well", h4(textOutput("kpi_total")), p("Total"))),
        column(3, div(class = "well", h4(textOutput("kpi_foraging")), p("Foraging"))),
        column(3, div(class = "well", h4(textOutput("kpi_mature_m")), p("Mature M"))),
        column(3, div(class = "well", h4(textOutput("kpi_legal")), p("Legal M")))
      ),
      fluidRow(
        column(3, div(class = "well", h4(textOutput("kpi_dead")), p("Dead"))),
        column(3, div(class = "well", h4(textOutput("kpi_harvested")), p("Harvested"))),
        column(3, div(class = "well", h4(textOutput("kpi_reproducing")), p("Reproducing"))),
        column(3, div(class = "well", h4(textOutput("kpi_mature_f")), p("Mature F")))
      ),
      tabsetPanel(
        tabPanel("Population Dynamics", plotlyOutput("pop_plot", height = "400px")),
        tabPanel("Abundance (log10)", plotlyOutput("pop_log_plot", height = "400px")),
        tabPanel("Spatial Map", plotOutput("spatial_map", height = "450px")),
        tabPanel("Stage Structure", plotlyOutput("stage_plot", height = "400px")),
        tabPanel("Size by Sex", plotlyOutput("size_plot", height = "400px")),
        tabPanel("Maturity & Selectivity", plotlyOutput("mat_sel_plot", height = "400px")),
        tabPanel("Harvest & Reproduction", plotlyOutput("harvest_plot", height = "400px")),
        tabPanel("Environmental Drivers", plotlyOutput("env_plot", height = "400px")),
        tabPanel("Stock-Recruitment", plotlyOutput("sr_plot", height = "400px")),
        tabPanel("Raw Data", DTOutput("raw_table"))
      ),
      div(style = "margin-top: 10px;", textOutput("status_text"))
    )
  )
)

server <- function(input, output, session) {
  sim_results <- reactiveVal(NULL)
  status_msg <- reactiveVal("Ready. Click 'Run Simulation'.")

  observeEvent(input$run_sim, {
    status_msg("Running...")
    tryCatch({
      sim <- create_crab_simulation(
        n_crabs = as.integer(input$inp_n_crabs),
        params = list(
          energy_init = as.numeric(input$inp_energy_init),
          energy_max = as.numeric(input$inp_energy_max),
          energy_gain_foraging = as.numeric(input$inp_foraging_gain),
          energy_cost_resting = as.numeric(input$inp_resting_cost),
          energy_cost_molting = as.numeric(input$inp_molting_cost),
          molt_growth_male = as.numeric(input$inp_molt_growth_male),
          molt_growth_female = as.numeric(input$inp_molt_growth_female),
          foraging_rate = as.numeric(input$inp_foraging_rate),
          resting_rate = as.numeric(input$inp_resting_rate),
          molting_rate = as.numeric(input$inp_molting_rate),
          mortality_rate = as.numeric(input$inp_mortality_rate),
          size_init_min = as.numeric(input$inp_size_min),
          size_init_max = as.numeric(input$inp_size_max),
          legal_size_male = as.numeric(input$inp_legal_size),
          maturity_l50_male = as.numeric(input$inp_maturity_l50_m),
          maturity_l50_female = as.numeric(input$inp_maturity_l50_f),
          maturity_slope_male = as.numeric(input$inp_maturity_slope_m),
          maturity_slope_female = as.numeric(input$inp_maturity_slope_f),
          maturity_slope_male = as.numeric(input$inp_maturity_slope_m),
          maturity_slope_female = as.numeric(input$inp_maturity_slope_f),
          max_size_male = as.numeric(input$inp_max_size_m),
          max_size_female = as.numeric(input$inp_max_size_f),
          fishing_probability = as.numeric(input$inp_fishing_prob),
          growth_beta0_male = as.numeric(input$inp_growth_beta0_m),
          growth_beta0_female = as.numeric(input$inp_growth_beta0_f),
          growth_beta1_male = as.numeric(input$inp_growth_beta1_m),
          growth_beta1_female = as.numeric(input$inp_growth_beta1_f),
          growth_slp_male = as.numeric(input$inp_growth_slp_m),
          growth_slp_female = as.numeric(input$inp_growth_slp_f),
          growth_beta2_male = as.numeric(input$inp_growth_beta2_m),
          growth_beta2_female = as.numeric(input$inp_growth_beta2_f),
          growth_beta3_male = as.numeric(input$inp_growth_beta3_m),
          growth_beta3_female = as.numeric(input$inp_growth_beta3_f),
          growth_temperature_coefficient = as.numeric(input$inp_growth_temp_coeff),
          growth_temperature_reference = as.numeric(input$inp_growth_temp_ref),
          fecundity_primiparous_intercept = as.numeric(input$inp_fec_int_prim),
          fecundity_primiparous_slope = as.numeric(input$inp_fec_slope_prim),
          fecundity_multiparous_intercept = as.numeric(input$inp_fec_int_multi),
          fecundity_multiparous_slope = as.numeric(input$inp_fec_slope_multi),
          fecundity_senescence_threshold = as.numeric(input$inp_fec_sen_thresh),
          fecundity_senescence_floor = as.numeric(input$inp_fec_sen_floor),
          fecundity_senescence_rate = as.numeric(input$inp_fec_sen_rate),
          minimum_clutch_size = as.numeric(input$inp_fec_min_clutch)
        )
      )
      times <- seq(0, as.numeric(input$inp_sim_time), by = 1)
      res <- run_crab_simulation(sim, times)
      sim_results(res)
      status_msg(paste("Done! ", nrow(res), "time steps."))
    }, error = function(e) {
      calls <- sys.calls()
      frames <- vapply(calls, function(c) paste(deparse(c), collapse = " "), "")
      hits <- grep("create_crab_simulation|run_crab_simulation|agentcrabr/R/model\\.R|agentcrabr/app\\.R", frames)
      sel <- unique(sort(unlist(lapply(hits, function(h) seq(max(1, h - 2), min(length(frames), h + 8))))))
      msg <- paste0("ERROR: ", conditionMessage(e),
             "\nCALL: ", deparse(conditionCall(e)),
             "\n--- relevant frames ---\n",
             paste(frames[sel], collapse = "\n"))
      message(msg)
      cat(msg, "\n", file = "error_log.txt", append = TRUE)
      msg
    })
  })

  env_data <- reactive({
    lat <- as.numeric(input$inp_lat)
    lon <- as.numeric(input$inp_lon)
    temp <- 8 + (lat - 42) * 0.5 + sin(lon + 63) * 2
    depth <- 50 + abs(lat - 44) * 30 + abs(lon + 63) * 10
    food <- 0.3 + cos(lat * 0.5) * 0.2 + sin(lon) * 0.1
    food <- max(0, min(1, food))
    hsi <- aoi_habitat_suitability(temp, depth)
    list(temp = temp, depth = depth, food = food, hsi = hsi)
  })

  output$loc_map <- renderPlot({
    # High-resolution bathymetry grid (300x300)
    lon_seq <- seq(-66.5, -57.5, length.out = 300)
    lat_seq <- seq(42.0, 48.5, length.out = 300)
    bathy_grid <- expand.grid(lon = lon_seq, lat = lat_seq)
    # Realistic bathymetry: continental shelf, slope, ridges, seamounts
    bathy_grid$depth <- 50 + abs(bathy_grid$lat - 44) * 25 + abs(bathy_grid$lon + 63) * 8
    # Add seamounts and ridges using sinusoidal variations
    bathy_grid$depth <- bathy_grid$depth + 
      80 * sin(bathy_grid$lon * 4) * cos(bathy_grid$lat * 3) +
      40 * sin(bathy_grid$lon * 9) * sin(bathy_grid$lat * 7) +
      25 * cos(bathy_grid$lon * 6) * sin(bathy_grid$lat * 4)
    bathy_grid$depth <- pmax(20, pmin(bathy_grid$depth, 800))
    coast <- map_data("worldHires")
    coast_region <- coast[coast$long >= -67 & coast$long <= -57 & coast$lat >= 41.5 & coast$lat <= 49, ]
    ggplot() +
      geom_raster(data = bathy_grid, aes(x = lon, y = lat, fill = depth), alpha = 0.75) +
      scale_fill_gradientn(colors = c("#08306b", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7", "#f7fbff"),
                           name = "Depth (m)", trans = "reverse") +
      geom_polygon(data = coast_region, aes(x = long, y = lat, group = group),
                   fill = "grey70", color = "grey30", linewidth = 0.3) +
      geom_point(aes(x = input$inp_lon, y = input$inp_lat), color = "red", size = 3) +
      coord_cartesian(xlim = c(-67, -57), ylim = c(41.5, 49)) +
      labs(x = "Longitude", y = "Latitude", title = "Study Area (click to select)") +
      theme_minimal() + theme(legend.position = "right")
  })

  observeEvent(input$map_click, {
    click <- input$map_click
    if (!is.null(click$x) && !is.null(click$y)) {
      lon <- max(-66.5, min(-57.5, click$x))
      lat <- max(42.0, min(48.5, click$y))
      updateNumericInput(session, "inp_lon", value = round(lon, 2))
      updateNumericInput(session, "inp_lat", value = round(lat, 2))
    }
  })

  output$coord_display <- renderText({
    paste0("Selected: ", input$inp_lat, " N, ", abs(input$inp_lon), " W")
  })

  output$env_table <- renderTable({
    e <- env_data()
    data.frame(
      Parameter = c("Bottom Temperature", "Depth", "Food Availability", "Habitat Quality (HSI)"),
      Value = c(round(e$temp, 1), round(e$depth, 0), round(e$food, 2), round(e$hsi, 2)),
      Unit = c("C", "m", "index", "0-1")
    )
  })

  output$update_status <- renderText({
    if (isTRUE(input$auto_update)) paste0("Auto-updating every ", input$update_interval, "s")
    else "Auto-update off"
  })

  # Helper to safely render KPIs
  safe_kpi <- function(expr) {
    tryCatch({
      expr
    }, error = function(e) {
      paste("Error:", conditionMessage(e))
    })
  }

  output$kpi_total <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(res$foraging[1] + res$resting[1] + res$molting[1]))
  })
  output$kpi_foraging <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$foraging, 1)))
  })
  output$kpi_dead <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$dead, 1)))
  })
  output$kpi_mean_size <- renderText("--")
  output$kpi_mature_m <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$mature_male, 1)))
  })
  output$kpi_mature_f <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$mature_female, 1)))
  })
  output$kpi_legal <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$legal_male, 1)))
  })
  output$kpi_harvested <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$harvested, 1)))
  })
  output$kpi_reproducing <- renderText({
    res <- sim_results()
    if (is.null(res)) return("--")
    safe_kpi(as.character(tail(res$reproducing, 1)))
  })

  output$status_text <- renderText({ status_msg() })

  output$pop_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      plot_ly(res, x = ~times) %>%
        add_lines(y = ~foraging, name = "Foraging", line = list(color = "#2ecc71")) %>%
        add_lines(y = ~resting, name = "Resting", line = list(color = "#3498db")) %>%
        add_lines(y = ~molting, name = "Molting", line = list(color = "#f39c12")) %>%
        add_lines(y = ~dead, name = "Dead", line = list(color = "#e74c3c")) %>%
        layout(title = "Population Dynamics", xaxis = list(title = "Time"), yaxis = list(title = "Count"))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$pop_log_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      plot_ly(res, x = ~times) %>%
        add_lines(y = ~log10(foraging + 1), name = "Foraging (log10)", line = list(color = "#2ecc71")) %>%
        add_lines(y = ~log10(resting + 1), name = "Resting (log10)", line = list(color = "#3498db")) %>%
        add_lines(y = ~log10(molting + 1), name = "Molting (log10)", line = list(color = "#f39c12")) %>%
        add_lines(y = ~log10(dead + 1), name = "Dead (log10)", line = list(color = "#e74c3c")) %>%
        layout(title = "Abundance Trajectories (log10)", xaxis = list(title = "Time"), yaxis = list(title = "log10(Count + 1)"))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$stage_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      last <- tail(res, 1)
      plot_ly(x = c("Foraging", "Resting", "Molting", "Dead"),
              y = c(last$foraging, last$resting, last$molting, last$dead),
              type = "bar", marker = list(color = c("#2ecc71", "#3498db", "#f39c12", "#e74c3c"))) %>%
        layout(title = "Living Stage Structure (Final)", xaxis = list(title = "State"), yaxis = list(title = "Count"))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$size_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      pos <- get("last_positions", envir = .GlobalEnv, inherits = TRUE)
      if (is.null(pos) || nrow(pos) == 0) return(plotly_empty())
      plot_ly(pos, x = ~size, color = ~sex, type = "histogram",
              nbinsx = 20,
              marker = list(color = c("#3498db" = "#3498db", "#e74c3c" = "#e74c3c"))) %>%
        layout(title = "Size Frequency Spectra by Sex", xaxis = list(title = "Size (mm)"), yaxis = list(title = "Count"),
               barmode = "overlay", opacity = 0.7)
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$mat_sel_plot <- renderPlotly({
    tryCatch({
      sizes <- seq(20, 160, by = 1)
      mat_m <- compute_maturity_ogive(sizes, as.numeric(input$inp_maturity_l50_m), 6.9957)
      mat_f <- compute_maturity_ogive(sizes, as.numeric(input$inp_maturity_l50_f), 12.9812)
      sel <- compute_selectivity(sizes, as.numeric(input$inp_legal_size), 5.0)
      plot_ly(x = sizes) %>%
        add_lines(y = mat_m, name = "Maturity (Male)", line = list(color = "#2980b9")) %>%
        add_lines(y = mat_f, name = "Maturity (Female)", line = list(color = "#e91e63")) %>%
        add_lines(y = sel, name = "Selectivity", line = list(color = "#27ae60", dash = "dash")) %>%
        layout(title = "Maturity Ogives & Selectivity", xaxis = list(title = "Size (mm)"), yaxis = list(title = "Probability", range = c(0, 1)))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$harvest_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      plot_ly(res, x = ~times) %>%
        add_lines(y = ~harvested, name = "Harvested", line = list(color = "#e74c3c")) %>%
        add_lines(y = ~reproducing, name = "Reproducing", line = list(color = "#9b59b6")) %>%
        layout(title = "Harvest & Reproduction Dynamics", xaxis = list(title = "Time"), yaxis = list(title = "Count"))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$env_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      plot_ly(res, x = ~times) %>%
        add_lines(y = ~foraging, name = "Foraging", line = list(color = "#2ecc71")) %>%
        add_lines(y = ~(foraging + resting + molting), name = "Alive", line = list(color = "#3498db", dash = "dash")) %>%
        layout(title = "Environmental Driver Response", xaxis = list(title = "Time"), yaxis = list(title = "Count"))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$sr_plot <- renderPlotly({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      plot_ly(res, x = ~mature_female, y = ~molting,
              type = "scatter", mode = "markers",
              marker = list(color = "#9b59b6", size = 8)) %>%
        add_lines(y = ~fitted(loess(molting ~ mature_female)),
                  name = "Loess fit", line = list(color = "#e74c3c")) %>%
        layout(title = "Stock-Recruitment Relationship",
               xaxis = list(title = "Mature Females (Spawner)"),
               yaxis = list(title = "Molting (Recruits)"))
    }, error = function(e) {
      plotly_empty() %>% layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$spatial_map <- renderPlot({
    res <- sim_results()
    if (is.null(res)) {
      plot(0, 0, type = "n", main = "Run simulation to see spatial map")
      return()
    }
    tryCatch({
      pos <- get("last_positions", envir = .GlobalEnv, inherits = TRUE)
      e <- env_data()
      x_seq <- seq(0, 1, length.out = 200)
      y_seq <- seq(0, 1, length.out = 200)
      grid <- expand.grid(x = x_seq, y = y_seq)
      grid$hsi <- e$hsi * (0.5 + sin(grid$x * 5) * 0.15 + cos(grid$y * 5) * 0.15 +
                          0.1 * sin(grid$x * 12) * cos(grid$y * 8))
      grid$hsi <- pmax(0, pmin(1, grid$hsi))
      ggplot() +
        geom_raster(data = grid, aes(x = x, y = y, fill = hsi), alpha = 0.6) +
        scale_fill_gradientn(colors = c("#08306b", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7"),
                             name = "HSI") +
        geom_point(data = pos, aes(x = x, y = y, color = state, shape = sex), size = 3, alpha = 0.8) +
        scale_color_manual(values = c("foraging" = "#2ecc71", "resting" = "#3498db",
                                     "molting" = "#f39c12", "dead" = "#e74c3c",
                                     "harvested" = "#95a5a6", "reproducing" = "#9b59b6")) +
        labs(title = "Spatial Distribution of Crabs on HSI Background", x = "X", y = "Y") +
        theme_minimal() + coord_fixed()
    }, error = function(e) {
      plot(0, 0, type = "n", main = paste("Error:", conditionMessage(e)))
    })
  })

  output$raw_table <- renderDT({
    res <- sim_results(); if (is.null(res)) return(NULL)
    tryCatch({
      datatable(res, options = list(pageLength = 10, scrollX = TRUE))
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })
}

shinyApp(ui, server)
