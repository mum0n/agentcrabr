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
  tick <- reactiveVal(0)
  sim_model <- reactiveVal(NULL)
  sim_history <- reactiveVal(data.frame())
  sim_positions <- reactiveVal(data.frame())
  selected_stage <- reactiveVal(NULL)
  annual_records <- reactiveVal(data.frame())
  exclusion_box <- reactiveVal(NULL)

  # Initialize model on startup
  observe({
    sim <- create_crab_simulation(120, list(
      benthic_carrying_capacity_density = 6000.0,
      grid_size = 50,
      unit_meters = 1000.0
    ))
    sim_model(sim)
    df <- run_crab_simulation(sim, 0:1)
    sim_history(df)
    sim_positions(crab_agent_positions(sim))
    current_day(1)
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

  # Auto-advance simulation while running (timer-driven to avoid reactive loops)
  auto_timer <- reactiveTimer(200)
  observe({
    auto_timer()
    if (is_running()) {
      tryCatch({
        d_step <- as.integer(input$step_days)
        new_day <- current_day() + d_step
        current_day(new_day)
        sim <- sim_model()
        if (!is.null(sim)) {
          df <- run_crab_simulation(sim, seq(new_day - d_step, new_day, by = 1))
          sim_history(rbind(sim_history(), df))
          sim_positions(crab_agent_positions(sim))
          tick(tick() + 1)
        }
        horizon <- as.numeric(input$duration_years) * 365
        if (new_day >= horizon) {
          is_running(FALSE)
          showNotification("Simulation horizon reached.", type = "message", duration = 6)
        }
      }, error = function(e) {
        is_running(FALSE)
        showNotification(paste("Simulation error:", conditionMessage(e)), type = "error")
      })
    }
  })
  observeEvent(input$btn_pause, { is_running(FALSE) })
  observeEvent(input$btn_stop, {
    is_running(FALSE)
    current_day(0)
    sim <- create_crab_simulation(120, list(
      benthic_carrying_capacity_density = 6000.0,
      grid_size = 50,
      unit_meters = 1000.0
    ))
    sim_model(sim)
    df <- run_crab_simulation(sim, 0:1)
    sim_history(df)
    sim_positions(crab_agent_positions(sim))
    current_day(1)
    annual_records(data.frame())
  })
  observeEvent(input$btn_step, {
    d_step <- as.integer(input$step_days)
    new_day <- current_day() + d_step
    current_day(new_day)
    sim <- sim_model()
    if (!is.null(sim)) {
      df <- run_crab_simulation(sim, seq(new_day - d_step, new_day, by = 1))
      sim_history(rbind(sim_history(), df))
      sim_positions(crab_agent_positions(sim))
      tick(tick() + 1)
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
      scale_x_continuous(limits = c(-67.0, -57.0), expand = c(0, 0)) +
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

  # Tab 2: Spatial Habitat Plotly — agents mapped into the AOI box over the
  # Scotian Shelf bathymetry, colored by behavioural state, sized by CW.
  output$spatial_map_plotly <- renderPlotly({
    tick()
    pos <- sim_positions()
    aoi_lon <- as.numeric(input$inp_lon)
    aoi_lat <- as.numeric(input$inp_lat)
    span <- 0.5   # degrees of the domain shown around the AOI centre

    living <- pos[!pos$state %in% c("dead", "harvested"), ]
    if (nrow(living) > 0) {
      living$map_lon <- aoi_lon + (living$x - 0.5) * span
      living$map_lat <- aoi_lat + (living$y - 0.5) * span
      living$marker_size <- pmax(4, living$size / 10)
    }
    state_colors <- c(foraging = "#0ea5e9", resting = "#14b8a6",
                      molting = "#f59e0b", reproducing = "#a855f7")

    mpa <- list(xmin = 0.25, xmax = 0.55, ymin = 0.25, ymax = 0.55)
    mpa_lon <- c(aoi_lon + (mpa$xmin - 0.5) * span, aoi_lon + (mpa$xmax - 0.5) * span)
    mpa_lat <- c(aoi_lat + (mpa$ymin - 0.5) * span, aoi_lat + (mpa$ymax - 0.5) * span)

    p <- plot_ly() %>%
      add_heatmap(
        data = bathy_raster_df,
        x = ~lon, y = ~lat, z = ~depth,
        colorscale = "Blues", reversescale = TRUE,
        showscale = FALSE, opacity = 0.7
      )
    if (nrow(living) > 0) {
      p <- p %>%
        add_markers(
          data = living,
          x = ~map_lon, y = ~map_lat,
          color = ~state, colors = state_colors,
          marker = list(size = ~marker_size, opacity = 0.85,
                        line = list(width = 0.5, color = "#0f172a"))
        )
    }
    p <- p %>%
      layout(
        xaxis = list(title = "Longitude (°W)",
                     range = c(aoi_lon - span / 2, aoi_lon + span / 2)),
        yaxis = list(title = "Latitude (°N)",
                     range = c(aoi_lat - span / 2, aoi_lat + span / 2)),
        shapes = list(list(
          type = "rect", xref = "x", yref = "y",
          x0 = mpa_lon[1], x1 = mpa_lon[2], y0 = mpa_lat[1], y1 = mpa_lat[2],
          line = list(color = "#16a34a", width = 2, dash = "dash"), fillcolor = "transparent"
        )),
        margin = list(l = 35, r = 15, t = 15, b = 35)
      )
    p
  })

  # Tab 2: Abundance Trajectories
  output$pop_abundance_plot <- renderPlotly({
    h <- sim_history()
    tick()
    if (nrow(h) > 0) {
      plot_ly(data = h, x = ~times) %>%
        add_lines(y = ~log10(pmax(1, foraging + resting + molting + reproducing)),
                  line = list(color = "#0f172a", width = 2), name = "Total Alive") %>%
        add_lines(y = ~log10(pmax(1, mature_male)),
                  line = list(color = "#1d4ed8", width = 1.5), name = "Mature Male") %>%
        add_lines(y = ~log10(pmax(1, mature_female)),
                  line = list(color = "#dc2626", width = 1.5), name = "Mature Female") %>%
        add_lines(y = ~log10(pmax(1, legal_male)),
                  line = list(color = "#16a34a", width = 1.5, dash = "dash"), name = "Legal Male") %>%
        add_lines(y = ~log10(pmax(1, dead)),
                  line = list(color = "#94a3b8", width = 1.2, dash = "dot"), name = "Dead") %>%
        layout(
          xaxis = list(title = "Simulation Day"),
          yaxis = list(title = "log10(Abundance)"),
          margin = list(l = 45, r = 15, t = 15, b = 35)
        )
    } else {
      plot_ly() %>% layout(title = "No history available")
    }
  })

  # Tab 2: Environmental Distributions — sampled from the real bathymetry and
  # seasonal temperature/HSI models inside the AOI box (no synthetic noise).
  output$env_distribution_plot <- renderPlotly({
    aoi_lon <- as.numeric(input$inp_lon)
    aoi_lat <- as.numeric(input$inp_lat)
    span <- 0.5
    in_box <- bathy_raster_df$lon >= aoi_lon - span / 2 & bathy_raster_df$lon <= aoi_lon + span / 2 &
      bathy_raster_df$lat >= aoi_lat - span / 2 & bathy_raster_df$lat <= aoi_lat + span / 2
    box_df <- bathy_raster_df[in_box & !is.na(bathy_raster_df$depth), ]
    if (nrow(box_df) > 400) box_df <- box_df[sample.int(nrow(box_df), 400), ]

    temp_off <- as.numeric(input$ctrl_temp_offset)
    temp_mult <- as.numeric(input$ctrl_temp_mult)
    hsi_off <- as.numeric(input$ctrl_hsi_offset)
    hsi_mult <- as.numeric(input$ctrl_hsi_mult)
    t_vals <- vapply(box_df$depth, function(d) {
      seasonal_bottom_temperature_cycle(0, d, aoi_lat, temp_off, temp_mult)$temps[1]
    }, numeric(1))
    h_vals <- vapply(seq_along(box_df$depth), function(k) {
      aoi_habitat_suitability(t_vals[k], box_df$depth[k], hsi_off, hsi_mult)
    }, numeric(1))

    p1 <- plot_ly(x = box_df$depth, type = "histogram", name = "Depth", marker = list(color = "#0284c7")) %>%
      layout(xaxis = list(title = "Depth (m)"))
    p2 <- plot_ly(x = t_vals, type = "histogram", name = "Temp", marker = list(color = "#dc2626")) %>%
      layout(xaxis = list(title = "Temp (°C)"))
    p3 <- plot_ly(x = h_vals, type = "histogram", name = "HSI", marker = list(color = "#16a34a")) %>%
      layout(xaxis = list(title = "HSI"))

    subplot(p1, p2, p3, nrows = 1, shareY = FALSE, titleX = TRUE) %>%
      layout(showlegend = FALSE, margin = list(l = 35, r = 15, t = 15, b = 35))
  })

  # Tab 3: Demographics & Growth — stage structure from live agent sizes.
  # Clicking a stage (lollipop head) selects that level: it is highlighted
  # and its detail is shown; clicking it again deselects.
  instar_mean_cw <- c(2.7, 5.9, 8.9, 12.0, 15.4, 19.2, 23.4, 28.0, 33.2, 39.5, 47.5, 65.0)
  stage_edges <- (head(instar_mean_cw, -1) + tail(instar_mean_cw, -1)) / 2
  stage_labels <- c(paste("Instar", 1:11), "Instar 12+")
  crab_stage_from_cw <- function(cw) as.integer(cut(cw, breaks = c(-Inf, stage_edges, Inf), labels = FALSE))

  stage_click <- event_data("plotly_click", source = "stage_src")
  observeEvent(stage_click, {
    cd <- stage_click
    if (is.null(cd) || nrow(cd) == 0L) return()
    idx <- suppressWarnings(as.integer(cd$customdata[1]))
    if (is.na(idx)) return()
    selected_stage(if (!identical(selected_stage(), idx)) idx else NULL)
  })

  output$stage_structure_plot <- renderPlotly({
    tick()
    pos <- sim_positions()
    if (nrow(pos) > 0) {
      pos$stage_f <- factor(crab_stage_from_cw(pos$size), levels = seq_along(stage_labels))
      counts <- table(pos$stage_f, pos$sex)
      male_counts <- as.integer(counts[, "M"])
      female_counts <- as.integer(counts[, "F"])
    } else {
      male_counts <- rep(0L, length(stage_labels))
      female_counts <- rep(0L, length(stage_labels))
    }
    x_pos <- seq_along(stage_labels)
    offset <- 0.2
    sel <- selected_stage()

    male_color <- ifelse(x_pos == sel, "#f59e0b", "#0284c7")
    female_color <- ifelse(x_pos == sel, "#f59e0b", "#dc2626")

    p <- plot_ly() %>%
      add_segments(x = x_pos - offset, xend = x_pos - offset, y = 0, yend = male_counts,
                   line = list(color = male_color, width = 2), name = "Males", showlegend = FALSE) %>%
      add_markers(x = x_pos - offset, y = male_counts, customdata = x_pos,
                  marker = list(color = male_color, size = 8, line = list(width = 0.5, color = "#ffffff")),
                  name = "Males") %>%
      add_segments(x = x_pos + offset, xend = x_pos + offset, y = 0, yend = female_counts,
                   line = list(color = female_color, width = 2), name = "Females", showlegend = FALSE) %>%
      add_markers(x = x_pos + offset, y = female_counts, customdata = x_pos,
                  marker = list(color = female_color, size = 8, line = list(width = 0.5, color = "#ffffff")),
                  name = "Females") %>%
      event_register("plotly_click", source = "stage_src") %>%
      layout(
        xaxis = list(tickvals = x_pos, ticktext = stage_labels, tickangle = -45),
        yaxis = list(title = "Abundance"),
        margin = list(l = 45, r = 15, t = 15, b = 65)
      )
    if (!is.null(sel)) {
      p <- p %>%
        layout(title = sprintf("%s — %d males, %d females (click again to clear)",
                               stage_labels[sel], male_counts[sel], female_counts[sel]))
    }
    p
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
    tick()
    pos <- sim_positions()
    if (nrow(pos) == 0) {
      plot_ly() %>% layout(title = "No agents yet — run the simulation")
    } else {
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
    }
  })

  output$cw_distribution_plot <- renderPlotly({
    tick()
    pos <- sim_positions()
    if (nrow(pos) == 0) {
      plot_ly() %>% layout(title = "No agents yet — run the simulation")
    } else {
      plot_ly(data = pos, x = ~size, type = "histogram", marker = list(color = "#6366f1")) %>%
        layout(
          xaxis = list(title = "Carapace Width (mm)"),
          yaxis = list(title = "Frequency"),
          margin = list(l = 45, r = 15, t = 15, b = 35)
        )
    }
  })

  # Tab 4: Macro Models — driven by the simulation history
  output$stock_recruitment_plot <- renderPlotly({
    h <- sim_history()
    tick()
    if (nrow(h) < 2) {
      plot_ly() %>% layout(title = "No history available")
    } else {
      plot_ly(data = h, x = ~times) %>%
        add_lines(y = ~mature_female, line = list(color = "#dc2626", width = 2),
                  name = "Spawning stock (mature females)") %>%
        add_lines(y = ~mature_male, line = list(color = "#1d4ed8", width = 1.5, dash = "dash"),
                  name = "Mature males") %>%
        layout(
          xaxis = list(title = "Simulation Day"),
          yaxis = list(title = "Abundance"),
          annotations = list(list(
            text = "Recruitment not modelled (closed cohort)",
            xref = "paper", yref = "paper", x = 0.02, y = 0.95,
            showarrow = FALSE, font = list(size = 11, color = "#64748b")
          )),
          margin = list(l = 45, r = 15, t = 15, b = 35)
        )
    }
  })

  output$macro_fits_plot <- renderPlotly({
    h <- sim_history()
    tick()
    if (nrow(h) < 2) {
      plot_ly() %>% layout(title = "No history available")
    } else {
      pos <- sim_positions()
      mean_mass <- if (nrow(pos) > 0) {
        mean(mapply(calculate_crab_mass_kg, pos$size, pos$sex, pos$mature == "1"))
      } else 0.24
      h$year <- ceiling(h$times / 365)
      annual <- do.call(rbind, lapply(split(h, h$year), function(g) {
        g[nrow(g), , drop = FALSE]
      }))
      plot_ly(data = annual, x = ~year) %>%
        add_lines(y = ~(foraging + resting + molting + reproducing) * mean_mass,
                  line = list(color = "#1d4ed8", width = 2), name = "Survey biomass (kg)") %>%
        add_bars(y = ~harvested * mean_mass, marker = list(color = "rgba(220, 38, 38, 0.6)"),
                 name = "Catch (kg)") %>%
        layout(
          xaxis = list(title = "Year", tickvals = annual$year, tickformat = "d"),
          yaxis = list(title = "Mass (kg)"),
          margin = list(l = 45, r = 15, t = 15, b = 35)
        )
    }
  })

  output$annual_summary_table <- renderDT({
    h <- sim_history()
    tick()
    if (nrow(h) < 2) {
      df <- data.frame(Year = character(0), End_Day = numeric(0), Alive = integer(0),
                       Harvested = integer(0), Mature_Females = integer(0),
                       Legal_Males = integer(0), Exploit_Rate = character(0))
    } else {
      h$year <- ceiling(h$times / 365)
      rows <- lapply(split(h, h$year), function(g) {
        last <- g[nrow(g), ]
        alive <- last$foraging + last$resting + last$molting + last$reproducing
        exploited <- if (alive + last$harvested > 0) {
          sprintf("%.1f%%", 100 * last$harvested / (alive + last$harvested))
        } else "0.0%"
        data.frame(Year = paste("Year", last$year), End_Day = last$times,
                   Alive = alive, Harvested = last$harvested,
                   Mature_Females = last$mature_female, Legal_Males = last$legal_male,
                   Exploit_Rate = exploited)
      })
      df <- do.call(rbind, rows)
    }
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })
}

shinyApp(ui, server)
