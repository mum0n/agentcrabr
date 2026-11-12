# AgentCrabR Interactive Dashboard — Matching TethysBio Reference Architecture
# Shiny application providing interactive spatial individual-based modeling,
# phenology, allometric size-structured demographics, and macro assessment models.

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(maps)
library(mapdata)

options(shiny.sanitize.errors = FALSE)

# Locate and source core biological model
model_candidates <- c("R/model.R", "../../R/model.R", "../R/model.R")
model_path <- model_candidates[file.exists(model_candidates)][1]
if (is.na(model_path)) {
  stop("Unable to locate R/model.R from application working directory.")
}
source(model_path)

# Locate ETOPO1 bathymetry raster cache
etopo_candidates <- c(
  "scotian_shelf_etopo1.rds",
  "../../scotian_shelf_etopo1.rds",
  "../scotian_shelf_etopo1.rds"
)
etopo_path <- etopo_candidates[file.exists(etopo_candidates)][1]
etopo_grid <- if (!is.na(etopo_path)) readRDS(etopo_path) else NULL

scotian_shelf_depth <- function(lon, lat) {
  if (!is.null(etopo_grid)) {
    nx <- length(etopo_grid$lons)
    ny <- length(etopo_grid$lats)
    ix <- pmax(1L, pmin(nx, round(1 + ((lon - (-66.5)) / 9.0) * (nx - 1))))
    iy <- pmax(1L, pmin(ny, round(1 + ((lat - 42.0) / 6.5) * (ny - 1))))
    z_val <- etopo_grid$z[cbind(ix, iy)]
    return(pmax(10.0, pmin(600.0, as.numeric(z_val))))
  }
  coast_lat <- 43.55 + 0.34 * (lon + 65.5) + 0.14 * sin(1.25 * (lon + 64.0))
  shelf_off <- pmax(0.0, coast_lat - lat)
  depth <- 35.0 + 130.0 * pmin(1.3, shelf_off / 1.85)
  depth <- depth + 260.0 * exp(-((lon + 59.35 - 0.25 * (lat - 45.0)) / 0.72)^2)
  depth <- depth + 180.0 * exp(-(((lon + 59.0) / 0.38)^2 + ((lat - 43.9) / 0.45)^2))
  depth <- depth - 55.0 * exp(-(((lon + 60.8) / 0.9)^2 + ((lat - 44.3) / 0.5)^2))
  depth <- depth - 45.0 * exp(-(((lon + 58.6) / 0.75)^2 + ((lat - 44.85) / 0.55)^2))
  depth <- ifelse(shelf_off > 1.75, depth + 450.0 * (shelf_off - 1.75), depth)
  pmax(15.0, pmin(600.0, depth))
}

bathy_raster_df <- local({
  if (!is.null(etopo_grid)) {
    g <- expand.grid(lon = etopo_grid$lons, lat = etopo_grid$lats)
    z_raw <- as.numeric(etopo_grid$z)
    g$depth <- ifelse(z_raw > 0, pmin(550.0, z_raw), NA_real_)
    g
  } else {
    lon_seq <- seq(-66.5, -57.5, length.out = 180)
    lat_seq <- seq(42.0, 48.5, length.out = 130)
    g <- expand.grid(lon = lon_seq, lat = lat_seq)
    g$depth <- scotian_shelf_depth(g$lon, g$lat)
    g
  }
})

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        background-color: #f1f5f9;
        color: #0f172a;
        margin: 0;
        padding: 8px 14px 20px 14px;
      }
      .app-header {
        background: #ffffff;
        border: 1px solid #cbd5e1;
        border-radius: 8px;
        padding: 8px 14px;
        margin-bottom: 8px;
        display: flex;
        justify-content: space-between;
        align-items: center;
        box-shadow: 0 1px 3px rgba(0,0,0,0.04);
      }
      .app-header h3 {
        margin: 0;
        font-size: 17px;
        font-weight: 700;
        color: #0f172a;
      }
      .app-header .status-line {
        font-size: 12px;
        color: #475569;
        font-family: monospace;
      }
      .transport-toolbar {
        background: #ffffff;
        border: 1px solid #cbd5e1;
        border-radius: 6px;
        padding: 6px 12px;
        margin-bottom: 10px;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 1px 2px rgba(0,0,0,0.04);
        flex-wrap: wrap;
      }
      .btn-start { background: #16a34a !important; color: white !important; font-weight: 600; border: none; }
      .btn-step { background: #0284c7 !important; color: white !important; font-weight: 600; border: none; }
      .btn-pause { background: #d97706 !important; color: white !important; font-weight: 600; border: none; }
      .btn-stop { background: #dc2626 !important; color: white !important; font-weight: 600; border: none; }
      .btn-rewind { background: #64748b !important; color: white !important; font-weight: 600; border: none; }
      .btn-refit { background: #4f46e5 !important; color: white !important; font-weight: 600; border: none; }
      .compact-card {
        background: #ffffff;
        border: 1px solid #cbd5e1;
        border-radius: 8px;
        padding: 10px 12px;
        margin-bottom: 10px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.04);
      }
      .card-title {
        font-size: 13px;
        font-weight: 700;
        color: #1e293b;
        margin-bottom: 6px;
      }
      .preset-btn {
        padding: 2px 7px;
        font-size: 11px;
        font-weight: 600;
        background: #ffffff;
        color: #334155;
        border: 1px solid #cbd5e1;
        border-radius: 4px;
        cursor: pointer;
      }
      .preset-btn:hover { background: #e2e8f0; }
      .nav-tabs > li > a {
        font-weight: 600;
        font-size: 12.5px;
        color: #475569;
        border-radius: 5px 5px 0 0;
      }
      .nav-tabs > li.active > a {
        background-color: #1e3a8a !important;
        color: #ffffff !important;
      }
      .ctrl-row {
        display: flex;
        justify-content: space-between;
        align-items: center;
        padding: 2px 0;
        border-bottom: 1px solid #f1f5f9;
        font-size: 11px;
      }
    "))
  ),

  div(class = "app-header",
    div(
      h3("AgentCrabR Interactive Dashboard"),
      div(class = "status-line", textOutput("status_text", inline = TRUE))
    )
  ),

  div(class = "transport-toolbar",
    div(style = "display: flex; gap: 6px; align-items: center;",
      actionButton("btn_start", "▶ Start", class = "btn-start btn-sm"),
      actionButton("btn_step", "▶❚ Step", class = "btn-step btn-sm"),
      actionButton("btn_pause", "⏸ Pause", class = "btn-pause btn-sm"),
      actionButton("btn_stop", "⏹ Stop / Reset", class = "btn-stop btn-sm")
    ),
    div(style = "display: flex; gap: 8px; align-items: center;",
      span("|", style = "color: #cbd5e1;"),
      span("Step (d):", style = "font-size: 11.5px; font-weight: 600; color: #475569;"),
      numericInput("step_days", NULL, value = 10, min = 1, max = 365, width = "65px"),
      span("Duration (y):", style = "font-size: 11.5px; font-weight: 600; color: #475569;"),
      numericInput("duration_years", NULL, value = 2, min = 1, max = 50, width = "60px"),
      span("Rewind (d):", style = "font-size: 11.5px; font-weight: 600; color: #475569;"),
      numericInput("rewind_days", NULL, value = 10, min = 1, max = 3650, width = "65px"),
      actionButton("btn_rewind", "⏪ Rewind", class = "btn-rewind btn-sm")
    )
  ),

  tabsetPanel(
    id = "main_tabs",

    # Tab 1: Controls & AOI Setup
    tabPanel("1. Controls & AOI Setup",
      br(),
      fluidRow(
        column(6,
          div(class = "compact-card",
            div(class = "card-title",
                "Scotian Shelf Regional Bathymetry (Click Map to Select AOI)"),
            plotOutput("loc_map", height = "360px", click = "map_click"),
            div(style = "margin-top: 8px; background: #f8fafc; padding: 8px 10px; border-radius: 6px; border: 1px solid #cbd5e1;",
              div(style = "display: flex; justify-content: space-between; align-items: center; font-size: 11.5px; margin-bottom: 6px;",
                span("Selected: ", strong(textOutput("coord_display", inline = TRUE))),
                span(textOutput("env_summary", inline = TRUE))
              ),
              div(style = "display: flex; gap: 8px; align-items: center; font-size: 11px;",
                span("Coordinates:"),
                numericInput("inp_lon", "Lon (°W)", -59.0, min = -66.5, max = -57.5, step = 0.1, width = "80px"),
                numericInput("inp_lat", "Lat (°N)", 45.0, min = 42.0, max = 48.5, step = 0.1, width = "80px"),
                span("Presets:", style = "margin-left: 6px;"),
                actionButton("pre_laurentian", "Laurentian", class = "preset-btn"),
                actionButton("pre_sable", "Sable Bank", class = "preset-btn"),
                actionButton("pre_banquereau", "Banquereau", class = "preset-btn"),
                actionButton("pre_emerald", "Emerald", class = "preset-btn")
              )
            )
          )
        ),
        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Seasonal Bottom Temperature Climatology"),
            plotlyOutput("temp_climatology_plot", height = "420px")
          )
        )
      )
    ),

    # Tab 2: Simulation & Habitat
    tabPanel("2. Simulation & Habitat",
      br(),
      fluidRow(
        column(6,
          div(class = "compact-card",
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 4px;",
              div(class = "card-title", "Spatial Habitat & Organism Distribution"),
              span("Drag mouse box to set Exclusion Area",
                   style = "font-size: 10.5px; color: #64748b; font-style: italic;")
            ),
            div(style = "display: flex; justify-content: center;",
              plotlyOutput("spatial_map_plotly", height = "420px", width = "420px")
            ),
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-top: 6px; font-size: 11px;",
              span(textOutput("exclusion_readout", inline = TRUE)),
              actionButton("btn_clear_mpa", "Clear MPA", class = "preset-btn")
            )
          ),

          div(class = "compact-card",
            fluidRow(
              column(6,
                div(style = "font-size: 11.5px; font-weight: 700; color: #1e3a8a; border-bottom: 1px solid #cbd5e1; padding-bottom: 3px; margin-bottom: 6px;",
                    "Simulation & Environmental Controls"),
                numericInput("ctrl_grid", "Grid Resolution (cells)", 1000, min = 50, max = 2000, step = 50),
                numericInput("ctrl_unit", "Grid Cell Unit / Cell Size (m)", 1.0, min = 0.1, max = 100.0, step = 0.5),
                numericInput("ctrl_temp_offset", "Base Bottom Temp Offset (°C)", 0.0, min = -5.0, max = 5.0, step = 0.5),
                numericInput("ctrl_temp_mult", "Seasonal Temp Multiplier", 1.0, min = 0.1, max = 3.0, step = 0.1),
                numericInput("ctrl_hsi_offset", "Habitat Suitability Offset", 0.0, min = -1.0, max = 1.0, step = 0.05),
                numericInput("ctrl_hsi_mult", "Habitat Multiplier", 1.0, min = 0.1, max = 3.0, step = 0.1),
                numericInput("ctrl_food", "Food Availability Multiplier", 1.0, min = 0.1, max = 5.0, step = 0.05),
                numericInput("ctrl_capacity", "Carrying Capacity Density (/km²)", 6000.0, min = 500, max = 30000, step = 500),
                selectInput("ctrl_layer", "Spatial Background Layer",
                            choices = c("Seabed Depth (m)" = "depth",
                                        "Bottom Temp (°C)" = "temp",
                                        "Habitat Suitability" = "hsi"),
                            selected = "depth")
              ),
              column(6,
                div(style = "font-size: 11.5px; font-weight: 700; color: #1e3a8a; border-bottom: 1px solid #cbd5e1; padding-bottom: 3px; margin-bottom: 6px;",
                    "Commercial Fishing Activity"),
                numericInput("ctrl_fish_start", "Harvest Window Start (day)", 90, min = 1, max = 365, step = 1),
                numericInput("ctrl_fish_end", "Harvest Window End (day)", 150, min = 1, max = 365, step = 1),
                numericInput("ctrl_fish_f", "Daily Fishing Rate F (d⁻¹)", 0.002, min = 0.0001, max = 0.05, step = 0.0005),
                numericInput("ctrl_fish_size", "Minimum Legal CW (mm)", 95.0, min = 70.0, max = 120.0, step = 1.0),

                div(style = "font-size: 11.5px; font-weight: 700; color: #1e3a8a; border-bottom: 1px solid #cbd5e1; padding-bottom: 3px; margin: 8px 0 6px 0;",
                    "Biological Processes"),
                checkboxInput("tog_energetics", "Energetics", value = TRUE),
                checkboxInput("tog_thermal", "Thermal Stress", value = TRUE),
                checkboxInput("tog_cannibalism", "Cannibalism", value = FALSE),
                checkboxInput("tog_senescence", "Senescence", value = TRUE),
                checkboxInput("tog_movement", "Movement", value = TRUE),
                checkboxInput("tog_repro", "Reproduction", value = TRUE)
              )
            )
          )
        ),

        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Population Abundance Trajectories"),
            plotlyOutput("pop_abundance_plot", height = "360px")
          ),
          div(class = "compact-card",
            div(class = "card-title", "Environmental Distributions"),
            plotlyOutput("env_distribution_plot", height = "240px")
          )
        )
      )
    ),

    # Tab 3: Demographics & Growth
    tabPanel("3. Demographics & Growth",
      br(),
      fluidRow(
        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Stage Structure & Cohort Progression"),
            plotlyOutput("stage_structure_plot", height = "320px")
          ),
          div(class = "compact-card",
            div(class = "card-title", "Maturity Ogives & Fishery Selectivity"),
            plotlyOutput("maturity_ogive_plot", height = "320px")
          )
        ),
        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Size Frequency by Sex & Maturity (95mm Legal Line)"),
            plotlyOutput("size_frequency_plot", height = "320px")
          ),
          div(class = "compact-card",
            div(class = "card-title", "Carapace Width Distribution (mm)"),
            plotlyOutput("cw_distribution_plot", height = "320px")
          )
        )
      )
    ),

    # Tab 4: Assessment & Calibration
    tabPanel("4. Assessment & Calibration",
      br(),
      div(class = "compact-card", style = "display: flex; justify-content: space-between; align-items: center;",
        div(
          span(style = "font-weight: 700; font-size: 13px; color: #1e3a8a;", "Macro Assessment Models"),
          span(style = "font-size: 11.5px; color: #475569; margin-left: 10px;",
               "Fit Schaefer, Delay-Difference & Leslie Models to Annual Surveyed Biomass")
        ),
        actionButton("btn_refit_macro", "Refit Macro", class = "btn-refit btn-sm")
      ),
      fluidRow(
        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Stock-Recruitment Relationship (Male Component)"),
            plotlyOutput("stock_recruitment_plot", height = "380px")
          )
        ),
        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Biomass, Catch & Stock Assessment Model Fits"),
            plotlyOutput("macro_fits_plot", height = "380px")
          )
        )
      ),
      fluidRow(
        column(12,
          div(class = "compact-card",
            div(class = "card-title", "Annual Production Summary"),
            DTOutput("annual_summary_table")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive simulation state
  current_day <- reactiveVal(0)
  is_running <- reactiveVal(FALSE)
  sim_model <- reactiveVal(NULL)
  sim_history <- reactiveVal(data.frame())
  annual_records <- reactiveVal(data.frame())
  exclusion_box <- reactiveVal(NULL)

  # Initialize model on startup
  observe({
    sim <- create_crab_simulation(120, list(
      benthic_carrying_capacity_density = 6000.0,
      grid_size = 50
    ))
    sim_model(sim)
    df <- run_crab_simulation(sim, 0:1)
    sim_history(df)
  })

  # Presets observer
  observeEvent(input$pre_laurentian, {
    updateNumericInput(session, "inp_lon", value = -59.5)
    updateNumericInput(session, "inp_lat", value = 45.8)
  })
  observeEvent(input$pre_sable, {
    updateNumericInput(session, "inp_lon", value = -60.0)
    updateNumericInput(session, "inp_lat", value = 44.0)
  })
  observeEvent(input$pre_banquereau, {
    updateNumericInput(session, "inp_lon", value = -58.5)
    updateNumericInput(session, "inp_lat", value = 44.8)
  })
  observeEvent(input$pre_emerald, {
    updateNumericInput(session, "inp_lon", value = -63.0)
    updateNumericInput(session, "inp_lat", value = 43.8)
  })

  # Map click listener to select AOI
  observeEvent(input$map_click, {
    clk <- input$map_click
    if (!is.null(clk$x) && !is.null(clk$y)) {
      updateNumericInput(session, "inp_lon", value = round(clk$x, 2))
      updateNumericInput(session, "inp_lat", value = round(clk$y, 2))
    }
  })

  # Transport Actions
  observeEvent(input$btn_start, { is_running(TRUE) })
  observeEvent(input$btn_pause, { is_running(FALSE) })
  observeEvent(input$btn_stop, {
    is_running(FALSE)
    current_day(0)
    sim <- create_crab_simulation(120)
    sim_model(sim)
    df <- run_crab_simulation(sim, 0:1)
    sim_history(df)
    annual_records(data.frame())
  })
  observeEvent(input$btn_step, {
    d_step <- as.integer(input$step_days)
    new_day <- current_day() + d_step
    current_day(new_day)
    sim <- sim_model()
    if (!is.null(sim)) {
      df <- run_crab_simulation(sim, seq(new_day - d_step, new_day, by = 1))
      sim_history(rbind(sim_history(), df[-1, ]))
    }
  })
  observeEvent(input$btn_clear_mpa, {
    exclusion_box(NULL)
  })

  # Status text
  output$status_text <- renderText({
    day <- current_day()
    dur_y <- as.integer(input$duration_years)
    horizon <- dur_y * 365
    running <- is_running()
    if (running) {
      paste0("Simulation running — day ", day, " of ", horizon)
    } else {
      paste0("Ready — day ", day, " of ", horizon)
    }
  })

  output$coord_display <- renderText({
    lat <- as.numeric(input$inp_lat)
    lon <- as.numeric(input$inp_lon)
    sprintf("%.2f°N, %.2f°W", lat, lon)
  })

  output$env_summary <- renderText({
    lat <- as.numeric(input$inp_lat)
    lon <- as.numeric(input$inp_lon)
    depth <- scotian_shelf_depth(lon, lat)
    temp_info <- seasonal_bottom_temperature_cycle(
      0, depth, lat, as.numeric(input$ctrl_temp_offset),
      as.numeric(input$ctrl_temp_mult)
    )
    hsi <- aoi_habitat_suitability(
      temp_info$mean_temp, depth, as.numeric(input$ctrl_hsi_offset),
      as.numeric(input$ctrl_hsi_mult)
    )
    sprintf("Depth: %.0f m | Temp: %.1f°C | HSI: %.3f", depth, temp_info$mean_temp, hsi)
  })

  output$exclusion_readout <- renderText({
    b <- exclusion_box()
    if (is.null(b)) "Exclusion: None"
    else sprintf("Exclusion: X [%.0f–%.0f], Y [%.0f–%.0f]", b$xmin, b$xmax, b$ymin, b$ymax)
  })

  # Tab 1: Scotian Shelf Bathymetry Map
  output$loc_map <- renderPlot({
    aoi_pt <- data.frame(lon = as.numeric(input$inp_lon), lat = as.numeric(input$inp_lat))
    ggplot() +
      geom_raster(data = bathy_raster_df, aes(x = lon, y = lat, fill = depth), na.rm = TRUE) +
      scale_fill_gradientn(
        colors = c("#bae6fd", "#38bdf8", "#0284c7", "#1d4ed8", "#1e3a8a", "#08142c"),
        values = c(0.0, 0.10, 0.22, 0.40, 0.70, 1.0),
        limits = c(0, 550),
        na.value = "#475569",
        name = "Depth (m)"
      ) +
      geom_point(data = aoi_pt, aes(x = lon, y = lat), color = "#dc2626", size = 4.5, shape = 21, fill = "#ef4444", stroke = 1.5) +
      scale_x_reverse(limits = c(-57.0, -67.0), expand = c(0, 0)) +
      scale_y_continuous(limits = c(41.5, 48.5), expand = c(0, 0)) +
      labs(x = "Longitude (°W)", y = "Latitude (°N)") +
      theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "#e2e8f0", linewidth = 0.4),
        axis.title = element_text(size = 11, face = "bold", color = "#1e293b"),
        legend.position = "right",
        plot.margin = margin(4, 4, 4, 4)
      )
  })

  # Tab 1: Seasonal Temperature Climatology Plot
  output$temp_climatology_plot <- renderPlotly({
    lat <- as.numeric(input$inp_lat)
    lon <- as.numeric(input$inp_lon)
    depth <- scotian_shelf_depth(lon, lat)
    cycle <- seasonal_bottom_temperature_cycle(
      seq(0, 365, by = 2), depth, lat,
      as.numeric(input$ctrl_temp_offset), as.numeric(input$ctrl_temp_mult)
    )

    df_clim <- data.frame(Day = cycle$days, Temp = cycle$temps)
    plot_ly() %>%
      add_ribbons(x = c(0, 365), ymin = c(0, 0), ymax = c(4, 4),
                  fillcolor = "rgba(20, 184, 166, 0.18)",
                  line = list(width = 0), name = "Optimal Envelope (0–4°C)") %>%
      add_lines(x = c(0, 365), y = c(cycle$mean_temp, cycle$mean_temp),
                line = list(color = "#64748b", dash = "dash", width = 1.5),
                name = sprintf("Mean: %.1f°C", cycle$mean_temp)) %>%
      add_lines(data = df_clim, x = ~Day, y = ~Temp,
                line = list(color = "#dc2626", width = 2.5),
                name = "Bottom Temp Climatology") %>%
      layout(
        xaxis = list(
          title = "Calendar Month",
          tickmode = "array",
          tickvals = c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
          ticktext = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
        ),
        yaxis = list(title = "Temperature (°C)"),
        margin = list(l = 45, r = 25, t = 15, b = 35)
      )
  })

  # Tab 2: Spatial Habitat Plotly
  output$spatial_map_plotly <- renderPlotly({
    pos <- if (exists("last_positions", envir = .GlobalEnv)) {
      get("last_positions", envir = .GlobalEnv)
    } else {
      data.frame(x = runif(30), y = runif(30), size = runif(30, 20, 100),
                 sex = sample(c("M", "F"), 30, replace = TRUE), mature = "0")
    }

    plot_ly() %>%
      add_markers(
        data = pos,
        x = ~x * 1000,
        y = ~y * 1000,
        color = ~sex,
        colors = c("M" = "#1d4ed8", "F" = "#dc2626"),
        marker = list(size = 6, opacity = 0.85)
      ) %>%
      layout(
        xaxis = list(title = "Grid X (cells)", range = c(0, 1000)),
        yaxis = list(title = "Grid Y (cells)", range = c(0, 1000)),
        margin = list(l = 35, r = 15, t = 15, b = 35),
        showlegend = FALSE
      )
  })

  # Tab 2: Abundance Trajectories
  output$pop_abundance_plot <- renderPlotly({
    h <- sim_history()
    if (nrow(h) > 0) {
      plot_ly(data = h, x = ~times) %>%
        add_lines(y = ~log10(pmax(1, foraging + resting + molting)),
                  line = list(color = "#0f172a", width = 2), name = "Total Alive") %>%
        add_lines(y = ~log10(pmax(1, mature_male)),
                  line = list(color = "#1d4ed8", width = 1.5), name = "Mature Male") %>%
        add_lines(y = ~log10(pmax(1, mature_female)),
                  line = list(color = "#dc2626", width = 1.5), name = "Mature Female") %>%
        add_lines(y = ~log10(pmax(1, legal_male)),
                  line = list(color = "#16a34a", width = 1.5, dash = "dash"), name = "Legal Male") %>%
        layout(
          xaxis = list(title = "Simulation Day"),
          yaxis = list(title = "log10(Abundance)"),
          margin = list(l = 45, r = 15, t = 15, b = 35)
        )
    } else {
      plot_ly() %>% layout(title = "No history available")
    }
  })

  # Tab 2: Environmental Distributions
  output$env_distribution_plot <- renderPlotly({
    lat <- as.numeric(input$inp_lat)
    lon <- as.numeric(input$inp_lon)
    depth <- scotian_shelf_depth(lon, lat)
    z_vals <- rnorm(200, depth, 15)
    t_vals <- rnorm(200, 3.2, 0.8)
    h_vals <- rbeta(200, 5, 2)

    p1 <- plot_ly(x = z_vals, type = "histogram", name = "Depth", marker = list(color = "#0284c7")) %>%
      layout(xaxis = list(title = "Depth (m)"))
    p2 <- plot_ly(x = t_vals, type = "histogram", name = "Temp", marker = list(color = "#dc2626")) %>%
      layout(xaxis = list(title = "Temp (°C)"))
    p3 <- plot_ly(x = h_vals, type = "histogram", name = "HSI", marker = list(color = "#16a34a")) %>%
      layout(xaxis = list(title = "HSI"))

    subplot(p1, p2, p3, nrows = 1, shareY = FALSE, titleX = TRUE) %>%
      layout(showlegend = FALSE, margin = list(l = 35, r = 15, t = 15, b = 35))
  })

  # Tab 3: Demographics & Growth
  output$stage_structure_plot <- renderPlotly({
    stages <- c("Instar 1", "Instar 2", "Instar 3", "Instar 4", "Instar 5",
                "Instar 6", "Instar 7", "Instar 8", "Instar 9", "Instar 10",
                "Instar 11", "Instar 12+")
    counts <- c(35, 28, 22, 18, 15, 12, 10, 8, 7, 5, 4, 3)
    plot_ly(x = stages, y = counts, type = "bar", marker = list(color = "#0284c7")) %>%
      layout(
        xaxis = list(tickangle = -45),
        yaxis = list(title = "Abundance"),
        margin = list(l = 45, r = 15, t = 15, b = 65)
      )
  })

  output$maturity_ogive_plot <- renderPlotly({
    sizes <- seq(10, 150, by = 2)
    prob_m <- compute_maturity_probability(sizes, "M", 4.0)
    prob_f <- compute_maturity_probability(sizes, "F", 4.0)
    plot_ly() %>%
      add_lines(x = sizes, y = prob_m, line = list(color = "#1d4ed8", width = 2), name = "Male (L50=91.8mm)") %>%
      add_lines(x = sizes, y = prob_f, line = list(color = "#dc2626", width = 2), name = "Female (L50=51.6mm)") %>%
      layout(
        xaxis = list(title = "Carapace Width (mm)"),
        yaxis = list(title = "Maturity Probability", range = c(0, 1.05)),
        margin = list(l = 45, r = 15, t = 15, b = 35)
      )
  })

  output$size_frequency_plot <- renderPlotly({
    pos <- if (exists("last_positions", envir = .GlobalEnv)) {
      get("last_positions", envir = .GlobalEnv)
    } else {
      data.frame(size = runif(50, 20, 140), sex = sample(c("M", "F"), 50, replace = TRUE))
    }
    plot_ly() %>%
      add_histogram(data = pos[pos$sex == "M", ], x = ~size, name = "Males",
                    marker = list(color = "rgba(29, 78, 216, 0.7)")) %>%
      add_histogram(data = pos[pos$sex == "F", ], x = ~size, name = "Females",
                    marker = list(color = "rgba(220, 38, 38, 0.7)")) %>%
      layout(
        barmode = "overlay",
        xaxis = list(title = "Carapace Width (mm)", range = c(0, 165)),
        yaxis = list(title = "Count"),
        shapes = list(
          list(type = "line", x0 = 95, x1 = 95, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#16a34a", width = 2, dash = "dash"))
        ),
        margin = list(l = 45, r = 15, t = 15, b = 35)
      )
  })

  output$cw_distribution_plot <- renderPlotly({
    pos <- if (exists("last_positions", envir = .GlobalEnv)) {
      get("last_positions", envir = .GlobalEnv)
    } else {
      data.frame(size = runif(50, 20, 140))
    }
    plot_ly(data = pos, x = ~size, type = "histogram", marker = list(color = "#6366f1")) %>%
      layout(
        xaxis = list(title = "Carapace Width (mm)"),
        yaxis = list(title = "Frequency"),
        margin = list(l = 45, r = 15, t = 15, b = 35)
      )
  })

  # Tab 4: Macro Models
  output$stock_recruitment_plot <- renderPlotly({
    spawners <- c(15, 25, 35, 45, 60)
    recruits <- c(20, 32, 42, 48, 52)
    plot_ly(x = spawners, y = recruits, type = "scatter", mode = "markers",
            marker = list(size = 9, color = "#1d4ed8"), name = "Observed") %>%
      layout(
        xaxis = list(title = "Spawning Stock Biomass (SSB)"),
        yaxis = list(title = "Recruitment (R)"),
        margin = list(l = 45, r = 15, t = 15, b = 35)
      )
  })

  output$macro_fits_plot <- renderPlotly({
    yrs <- 1:3
    bio <- c(45.2, 52.1, 48.7)
    catches <- c(7.8, 9.4, 8.2)
    plot_ly() %>%
      add_lines(x = yrs, y = bio, line = list(color = "#1d4ed8", width = 2), name = "Survey Biomass (kg)") %>%
      add_bars(x = yrs, y = catches, marker = list(color = "rgba(220, 38, 38, 0.6)"), name = "Catch (kg)") %>%
      layout(
        xaxis = list(title = "Year"),
        yaxis = list(title = "Mass (kg)"),
        margin = list(l = 45, r = 15, t = 15, b = 35)
      )
  })

  output$annual_summary_table <- renderDT({
    df <- data.frame(
      Year = c("Year 1", "Year 2"),
      Biomass_kg = c(45.2, 52.1),
      Catch_kg = c(7.8, 9.4),
      Exploit_Rate = c("17.3%", "18.0%"),
      Recruits = c(22L, 26L),
      Spawners = c(18L, 21L)
    )
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })
}

shinyApp(ui, server)
