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
    fx <- 1.0 + (lon - etopo_grid$lons[1]) / (etopo_grid$lons[nx] - etopo_grid$lons[1]) * (nx - 1L)
    fy <- 1.0 + (lat - etopo_grid$lats[1]) / (etopo_grid$lats[ny] - etopo_grid$lats[1]) * (ny - 1L)
    fx <- pmax(1.0, pmin(as.numeric(nx), fx))
    fy <- pmax(1.0, pmin(as.numeric(ny), fy))
    ix0 <- pmin(nx - 1L, as.integer(floor(fx)))
    iy0 <- pmin(ny - 1L, as.integer(floor(fy)))
    tx <- fx - ix0
    ty <- fy - iy0
    z00 <- as.numeric(etopo_grid$z[cbind(ix0, iy0)])
    z10 <- as.numeric(etopo_grid$z[cbind(ix0 + 1L, iy0)])
    z01 <- as.numeric(etopo_grid$z[cbind(ix0, iy0 + 1L)])
    z11 <- as.numeric(etopo_grid$z[cbind(ix0 + 1L, iy0 + 1L)])
    z_val <- (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 +
      (1 - tx) * ty * z01 + tx * ty * z11
    return(pmax(10.0, pmin(600.0, z_val)))
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
    lon_seq <- seq(-66.5, -56.5, length.out = 180)
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

    # Tab 1: Controls & AOI Setup (80% width, no scroll)
    tabPanel("1. Controls & AOI Setup",
      br(),
      div(style = "width: 80%; margin: 0 auto;",
        fluidRow(
          column(6,
            div(class = "compact-card",
              div(class = "card-title",
                  "Scotian Shelf Regional Bathymetry (Click Map to Select AOI)"),
              plotOutput("loc_map", height = "300px", click = "map_click"),
              div(style = "margin-top: 6px; background: #f8fafc; padding: 6px 8px; border-radius: 6px; border: 1px solid #cbd5e1;",
                div(style = "display: flex; justify-content: space-between; align-items: center; font-size: 11.5px; margin-bottom: 4px;",
                  span("Selected: ", strong(textOutput("coord_display", inline = TRUE))),
                  span(textOutput("env_summary", inline = TRUE))
                ),
                div(style = "display: flex; gap: 8px; align-items: center; flex-wrap: wrap; font-size: 11px;",
                  span("Coordinates:"),
                  numericInput("inp_lon", "Lon (°W)", -59.0, min = -66.5, max = -56.5, step = 0.1, width = "80px"),
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
              plotlyOutput("temp_climatology_plot", height = "355px")
            )
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
            div(style = "display: flex; justify-content: space-between; align-items: center; gap: 8px; flex-wrap: wrap; margin-top: 6px; padding: 4px 8px; background: #f8fafc; border-radius: 4px; border: 1px solid #e2e8f0; font-size: 11px;",
              div(style = "display: flex; align-items: center; gap: 6px;",
                span(style = "font-weight: 600; color: #334155;", "Background Layer:"),
                selectInput("ctrl_layer", NULL,
                            choices = c("Seabed Depth (m)" = "depth",
                                        "Bottom Temp (°C)" = "temp",
                                        "Habitat Suitability (HSI)" = "hsi"),
                            selected = "depth", width = "170px")
              ),
              div(style = "display: flex; align-items: center; gap: 8px;",
                span(textOutput("exclusion_readout", inline = TRUE)),
                actionButton("btn_clear_mpa", "Clear exclusion zone", class = "preset-btn")
              )
            )
          )
        ),

        column(6,
          div(class = "compact-card",
            div(class = "card-title", "Population Abundance Trajectories"),
            plotlyOutput("pop_abundance_plot", height = "250px")
          ),
          div(class = "compact-card",
            div(class = "card-title", "Environmental Distributions"),
            plotlyOutput("env_distribution_plot", height = "210px")
          )
        )
      ),

      # Full-Width 4-Column Grouped Controls at Bottom of Tab 2 with Hover Tooltips
      fluidRow(
        column(12,
          div(class = "compact-card", style = "background: #f8fafc;",
            fluidRow(
              column(3,
                div(style = "background: #ffffff; border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px;",
                  div(style = "font-size: 12px; font-weight: 700; color: #1e3a8a; border-bottom: 2px solid #cbd5e1; padding-bottom: 4px; margin-bottom: 8px;",
                      "Grid & Spatial Habitat"),
                  div(title = "Number of grid cells along each axis of the spatial simulation domain.", style = "cursor: help;",
                    numericInput("ctrl_grid", "Grid Resolution (cells)", 1000, min = 50, max = 2000, step = 50)
                  ),
                  div(title = "Physical width of each grid cell in meters; scales movement step size and domain area.", style = "cursor: help;",
                    numericInput("ctrl_unit", "Grid Unit / Cell Size (m)", 1.0, min = 0.1, max = 100.0, step = 0.5)
                  ),
                  div(title = "Optimal seabed depth (m) for benthic snow crab habitat suitability (HSI = 1 peak; open-ended deeper).", style = "cursor: help;",
                    numericInput("ctrl_depth_opt", "Depth Opt (m)", 200.0, min = 25.0, max = 1000.0, step = 10.0)
                  ),
                  div(title = "Gaussian width (m) of the shallow depth suitability taper below Depth Opt.", style = "cursor: help;",
                    numericInput("ctrl_depth_sigma", "Depth Sigma (m)", 100.0, min = 10.0, max = 500.0, step = 10.0)
                  ),
                  div(title = "Maximum annual movement budget (m/yr) converted into daily habitat taxis step distance.", style = "cursor: help;",
                    numericInput("ctrl_max_move", "Max Movement (m/yr)", 7300.0, min = 0.0, max = 50000.0, step = 500.0)
                  )
                )
              ),
              column(3,
                div(style = "background: #ffffff; border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px;",
                  div(style = "font-size: 12px; font-weight: 700; color: #1e3a8a; border-bottom: 2px solid #cbd5e1; padding-bottom: 4px; margin-bottom: 8px;",
                      "Thermal Environment"),
                  div(title = "Constant temperature shift (°C) added to the seasonal bottom temperature climatology.", style = "cursor: help;",
                    numericInput("ctrl_temp_offset", "Base Bottom Temp Offset (°C)", 0.0, min = -5.0, max = 5.0, step = 0.5)
                  ),
                  div(title = "Scaling factor applied to the seasonal bottom temperature amplitude.", style = "cursor: help;",
                    numericInput("ctrl_temp_mult", "Seasonal Temp Multiplier", 1.0, min = 0.1, max = 3.0, step = 0.1)
                  ),
                  div(title = "Lower temperature bound (°C) where thermal habitat suitability tapers to zero.", style = "cursor: help;",
                    numericInput("ctrl_temp_pref_min", "Temp Pref Min (°C)", -1.5, min = -5.0, max = 15.0, step = 0.5)
                  ),
                  div(title = "Optimal bottom temperature (°C) for peak thermal habitat suitability.", style = "cursor: help;",
                    numericInput("ctrl_temp_pref_opt", "Temp Pref Opt (°C)", 5.0, min = -2.0, max = 20.0, step = 0.5)
                  ),
                  div(title = "Upper exclusion temperature (°C) above which thermal habitat suitability is zero.", style = "cursor: help;",
                    numericInput("ctrl_temp_pref_max", "Temp Pref Max (°C)", 10.0, min = 0.0, max = 25.0, step = 0.5)
                  )
                )
              ),
              column(3,
                div(style = "background: #ffffff; border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px;",
                  div(style = "font-size: 12px; font-weight: 700; color: #1e3a8a; border-bottom: 2px solid #cbd5e1; padding-bottom: 4px; margin-bottom: 8px;",
                      "Habitat & Carrying Capacity"),
                  div(title = "Additive offset applied to computed cell Habitat Suitability Index (HSI).", style = "cursor: help;",
                    numericInput("ctrl_hsi_offset", "Habitat Suitability Offset", 0.0, min = -1.0, max = 1.0, step = 0.05)
                  ),
                  div(title = "Multiplicative scaling factor applied to cell Habitat Suitability Index (HSI).", style = "cursor: help;",
                    numericInput("ctrl_hsi_mult", "Habitat Multiplier", 1.0, min = 0.1, max = 3.0, step = 0.1)
                  ),
                  div(title = "Scaling multiplier for benthic food availability driving bioenergetic intake and taxis.", style = "cursor: help;",
                    numericInput("ctrl_food", "Food Availability Multiplier", 1.0, min = 0.1, max = 5.0, step = 0.05)
                  ),
                  div(title = "Reference benthic carrying capacity (Total No./km²) for the empirically measured size fraction (instars 8+, stage ≥ 11, CW ≥ 35 mm, and mature crabs); regulates density-dependent baseline mortality when abundance exceeds K.", style = "cursor: help;",
                    numericInput("ctrl_capacity", "Carrying Capacity (Total No./km²)", 6000.0, min = 500, max = 30000, step = 500)
                  ),
                  div(title = "Local crowding box width B in grid cells (0 = auto 1/5 of grid resolution). Triggers elevated crowding mortality and dispersal when local density exceeds carrying capacity.", style = "cursor: help;",
                    numericInput("ctrl_crowding_box", "Crowding Box B (cells, 0=G/5)", 0, min = 0, max = 2000, step = 10)
                  )
                )
              ),
              column(3,
                div(style = "background: #ffffff; border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px;",
                  div(style = "font-size: 12px; font-weight: 700; color: #1e3a8a; border-bottom: 2px solid #cbd5e1; padding-bottom: 4px; margin-bottom: 8px;",
                      "Commercial Fishery"),
                  div(title = "Day of year (1–365) when the commercial snow crab fishing season opens.", style = "cursor: help;",
                    numericInput("ctrl_fish_start", "Fishing Start (Day)", 90, min = 1, max = 365, step = 1)
                  ),
                  div(title = "Day of year (1–365) when the commercial snow crab fishing season closes.", style = "cursor: help;",
                    numericInput("ctrl_fish_end", "Fishing End (Day)", 150, min = 1, max = 365, step = 1)
                  ),
                  div(title = "Target finite annual harvest rate (0–0.99), rescaled internally to a daily fishing hazard over the active fishing season.", style = "cursor: help;",
                    numericInput("ctrl_fish_f", "Annual Harvest Rate", 0.0, min = 0.0, max = 0.99, step = 0.01)
                  ),
                  div(title = "Minimum carapace width (mm) at 50% gear selectivity for commercial male harvest.", style = "cursor: help;",
                    numericInput("ctrl_fish_size", "Legal Male Size (mm)", 95.0, min = 70.0, max = 120.0, step = 1.0)
                  ),
                  div(style = "font-size: 12px; font-weight: 700; color: #1e3a8a; border-bottom: 2px solid #cbd5e1; padding-bottom: 4px; margin: 10px 0 6px 0;",
                      "Biological Processes"),
                  div(title = "Enable Kleiber metabolic maintenance, Holling Type II intake, and mass-conserving ecdysis energetics.", style = "cursor: help;",
                    checkboxInput("tog_energetics", "Bioenergetics", value = TRUE)
                  ),
                  div(title = "Enable elevated mortality and ecdysial inhibition above thermal stress thresholds.", style = "cursor: help;",
                    checkboxInput("tog_thermal", "Thermal Stress", value = TRUE)
                  ),
                  div(title = "Enable size-dependent cannibalism mortality on small benthic juveniles (< 35 mm).", style = "cursor: help;",
                    checkboxInput("tog_cannibalism", "Cannibalism", value = FALSE)
                  ),
                  div(title = "Enable exponential post-terminal-molt mortality after 3 years of shell age.", style = "cursor: help;",
                    checkboxInput("tog_senescence", "Senescence", value = TRUE)
                  ),
                  div(title = "Enable directed softmax habitat taxis and localized crowding dispersal.", style = "cursor: help;",
                    checkboxInput("tog_movement", "Taxis Movement", value = TRUE)
                  ),
                  div(title = "Enable annual female egg extrusion and fecundity senescence.", style = "cursor: help;",
                    checkboxInput("tog_repro", "Reproduction", value = TRUE)
                  )
                )
              )
            )
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

  # Auto-advance simulation while running (timer-driven; isolate prevents self-invalidation)
  auto_timer <- reactiveTimer(250)
  observe({
    auto_timer()
    isolate({
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
  })
  observeEvent(input$btn_pause, { is_running(FALSE) })
  observeEvent(input$btn_stop, {
    is_running(FALSE)
    current_day(0)
    mpa <- exclusion_box()
    sim <- create_crab_simulation(120, list(
      benthic_carrying_capacity_density = as.numeric(input$ctrl_capacity %||% 6000.0),
      local_crowding_box_size = as.integer(input$ctrl_crowding_box %||% 0L),
      grid_size = as.numeric(input$ctrl_grid %||% 1000),
      unit_meters = as.numeric(input$ctrl_unit %||% 1.0),
      depth_opt = as.numeric(input$ctrl_depth_opt %||% 200.0),
      depth_sigma = as.numeric(input$ctrl_depth_sigma %||% 100.0),
      temp_pref_min = as.numeric(input$ctrl_temp_pref_min %||% -1.5),
      temp_pref_opt = as.numeric(input$ctrl_temp_pref_opt %||% 5.0),
      temp_pref_max = as.numeric(input$ctrl_temp_pref_max %||% 10.0),
      max_move_m_per_year = as.numeric(input$ctrl_max_move %||% 7300.0),
      food_availability = as.numeric(input$ctrl_food %||% 1.0),
      fishing_season_start = as.numeric(input$ctrl_fish_start %||% 90.0),
      fishing_season_end = as.numeric(input$ctrl_fish_end %||% 150.0),
      fishing_harvest_rate_annual = as.numeric(input$ctrl_fish_f %||% 0.0),
      legal_size_male = as.numeric(input$ctrl_fish_size %||% 95.0),
      use_energetics = isTRUE(input$tog_energetics),
      use_thermal_stress = isTRUE(input$tog_thermal),
      use_cannibalism = isTRUE(input$tog_cannibalism),
      use_senescence = isTRUE(input$tog_senescence),
      use_movement = isTRUE(input$tog_movement),
      use_reproduction = isTRUE(input$tog_repro),
      mpa = mpa
    ))
    sim_model(sim)
    df <- run_crab_simulation(sim, 0:1)
    sim_history(df)
    sim_positions(crab_agent_positions(sim))
    current_day(1)
    annual_records(data.frame())
  })
  observeEvent(input$btn_rewind, {
    is_running(FALSE)
    h <- sim_history()
    if (nrow(h) > 1) {
      rw <- as.numeric(input$rewind_days %||% 30)
      target_t <- max(1, tail(h$times, 1) - rw)
      keep <- h$times <= target_t
      if (!any(keep)) keep[1] <- TRUE
      h_sub <- h[keep, , drop = FALSE]
      sim_history(h_sub)
      current_day(tail(h_sub$times, 1))
      tick(tick() + 1)
    }
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
    sim_date <- format_sim_date_r(day)
    if (running) {
      paste0("Simulation running — day ", day, " of ", horizon, " (", sim_date, ")")
    } else {
      paste0("Ready — day ", day, " of ", horizon, " (", sim_date, ")")
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
      temp_info$mean_temp, depth,
      offset = as.numeric(input$ctrl_hsi_offset %||% 0.0),
      mult = as.numeric(input$ctrl_hsi_mult %||% 1.0),
      temp_pref_min = as.numeric(input$ctrl_temp_pref_min %||% -1.5),
      temp_pref_opt = as.numeric(input$ctrl_temp_pref_opt %||% 5.0),
      temp_pref_max = as.numeric(input$ctrl_temp_pref_max %||% 10.0),
      depth_opt = as.numeric(input$ctrl_depth_opt %||% 200.0),
      depth_sigma = as.numeric(input$ctrl_depth_sigma %||% 100.0)
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
      scale_x_continuous(limits = c(-66.5, -56.5), expand = c(0, 0)) +
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

    mpa <- exclusion_box()
    shapes_list <- list()
    if (!is.null(mpa) && isTRUE(mpa$active)) {
      mpa_lon <- c(aoi_lon + (mpa$xmin - 0.5) * span, aoi_lon + (mpa$xmax - 0.5) * span)
      mpa_lat <- c(aoi_lat + (mpa$ymin - 0.5) * span, aoi_lat + (mpa$ymax - 0.5) * span)
      shapes_list <- list(list(
        type = "rect", xref = "x", yref = "y",
        x0 = mpa_lon[1], x1 = mpa_lon[2], y0 = mpa_lat[1], y1 = mpa_lat[2],
        line = list(color = "#16a34a", width = 2, dash = "dash"), fillcolor = "rgba(22, 163, 74, 0.1)"
      ))
    }

    layer_choice <- tolower(as.character(input$ctrl_layer %||% "depth"))
    n_bg <- 90L
    lon_bg <- seq(aoi_lon - span / 2, aoi_lon + span / 2, length.out = n_bg)
    lat_bg <- seq(aoi_lat - span / 2, aoi_lat + span / 2, length.out = n_bg)
    bg_df <- expand.grid(lon = lon_bg, lat = lat_bg)
    bg_df$depth <- scotian_shelf_depth(bg_df$lon, bg_df$lat)
    bg_z <- bg_df$depth
    habitat_colorscale <- list(
      list(0.00, "#15385E"),
      list(0.25, "#00A6CA"),
      list(0.50, "#C1EEF1"),
      list(0.75, "#FFE14F"),
      list(1.00, "#FE7958")
    )
    t_off <- as.numeric(input$ctrl_temp_offset %||% 0.0)
    t_mult <- as.numeric(input$ctrl_temp_mult %||% 1.0)
    doy_now <- (current_day() %% 365) + 1
    if (grepl("temp", layer_choice)) {
      bg_z <- vapply(seq_along(bg_df$depth), function(i) {
        seasonal_bottom_temperature_cycle(doy_now, bg_df$depth[i], bg_df$lat[i], t_off, t_mult)$temps[1]
      }, numeric(1))
    } else if (grepl("hsi|habitat|suit", layer_choice)) {
      t_bg <- vapply(seq_along(bg_df$depth), function(i) {
        seasonal_bottom_temperature_cycle(doy_now, bg_df$depth[i], bg_df$lat[i], t_off, t_mult)$temps[1]
      }, numeric(1))
      bg_z <- vapply(seq_along(bg_df$depth), function(i) {
        aoi_habitat_suitability(
          t_bg[i], bg_df$depth[i],
          offset = as.numeric(input$ctrl_hsi_offset %||% 0.0),
          mult = as.numeric(input$ctrl_hsi_mult %||% 1.0),
          temp_pref_min = as.numeric(input$ctrl_temp_pref_min %||% -1.5),
          temp_pref_opt = as.numeric(input$ctrl_temp_pref_opt %||% 5.0),
          temp_pref_max = as.numeric(input$ctrl_temp_pref_max %||% 10.0),
          depth_opt = as.numeric(input$ctrl_depth_opt %||% 200.0),
          depth_sigma = as.numeric(input$ctrl_depth_sigma %||% 100.0)
        )
      }, numeric(1))
    }
    bg_df$layer_val <- bg_z

    p <- plot_ly(source = "spatial_map") %>%
      add_heatmap(
        data = bg_df,
        x = ~lon, y = ~lat, z = ~layer_val,
        zsmooth = "best",
        colorscale = habitat_colorscale,
        reversescale = !grepl("temp", layer_choice),
        showscale = TRUE, opacity = 0.85
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
        dragmode = "select",
        xaxis = list(title = "Longitude (°W)",
                     range = c(aoi_lon - span / 2, aoi_lon + span / 2)),
        yaxis = list(title = "Latitude (°N)",
                     range = c(aoi_lat - span / 2, aoi_lat + span / 2)),
        shapes = shapes_list,
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

    temp_off <- as.numeric(input$ctrl_temp_offset %||% 0.0)
    temp_mult <- as.numeric(input$ctrl_temp_mult %||% 1.0)
    hsi_off <- as.numeric(input$ctrl_hsi_offset %||% 0.0)
    hsi_mult <- as.numeric(input$ctrl_hsi_mult %||% 1.0)
    doy_now <- (current_day() %% 365) + 1
    t_vals <- vapply(box_df$depth, function(d) {
      seasonal_bottom_temperature_cycle(doy_now, d, aoi_lat, temp_off, temp_mult)$temps[1]
    }, numeric(1))
    h_vals <- vapply(seq_along(box_df$depth), function(k) {
      aoi_habitat_suitability(
        t_vals[k], box_df$depth[k],
        offset = hsi_off, mult = hsi_mult,
        temp_pref_min = as.numeric(input$ctrl_temp_pref_min %||% -1.5),
        temp_pref_opt = as.numeric(input$ctrl_temp_pref_opt %||% 5.0),
        temp_pref_max = as.numeric(input$ctrl_temp_pref_max %||% 10.0),
        depth_opt = as.numeric(input$ctrl_depth_opt %||% 200.0),
        depth_sigma = as.numeric(input$ctrl_depth_sigma %||% 100.0)
      )
    }, numeric(1))

    d_min <- min(box_df$depth, na.rm = TRUE)
    d_max <- max(box_df$depth, na.rm = TRUE)
    if (d_max <= d_min) d_max <- d_min + 10.0
    clim_cycle <- seasonal_bottom_temperature_cycle(seq(1, 365, by = 2), mean(box_df$depth, na.rm = TRUE), aoi_lat, temp_off, temp_mult)
    t_cmin <- min(c(clim_cycle$temps, t_vals), na.rm = TRUE)
    t_cmax <- max(c(clim_cycle$temps, t_vals), na.rm = TRUE)
    if (t_cmax <= t_cmin) t_cmax <- t_cmin + 1.0

    p1 <- plot_ly(
      x = box_df$depth, type = "histogram", name = "Depth",
      autobinx = FALSE,
      xbins = list(start = d_min, end = d_max, size = (d_max - d_min) / 24.0),
      marker = list(color = "#1e3a8a")
    ) %>% layout(xaxis = list(title = "Depth (m)"))
    p2 <- plot_ly(
      x = t_vals, type = "histogram", name = "Temp",
      autobinx = FALSE,
      xbins = list(start = t_cmin, end = t_cmax, size = (t_cmax - t_cmin) / 24.0),
      marker = list(color = "#dc2626")
    ) %>% layout(xaxis = list(title = "Temp (°C)", range = c(t_cmin, t_cmax)))
    p3 <- plot_ly(
      x = h_vals, type = "histogram", name = "HSI",
      autobinx = FALSE,
      xbins = list(start = 0.0, end = 1.0, size = 1.0 / 24.0),
      marker = list(color = "#0f766e")
    ) %>% layout(xaxis = list(title = "HSI [0, 1]", range = c(0.0, 1.0)))

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

  observeEvent(event_data("plotly_click", source = "stage_src"), {
    cd <- event_data("plotly_click", source = "stage_src")
    if (is.null(cd) || nrow(cd) == 0L) return()
    idx <- suppressWarnings(as.integer(cd$customdata[1]))
    if (is.na(idx)) return()
    selected_stage(if (!identical(selected_stage(), idx)) idx else NULL)
  }, ignoreNULL = TRUE)

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

    male_color <- if (is.null(sel)) rep("#0284c7", length(x_pos)) else ifelse(x_pos == sel, "#f59e0b", "#0284c7")
    female_color <- if (is.null(sel)) rep("#dc2626", length(x_pos)) else ifelse(x_pos == sel, "#f59e0b", "#dc2626")

    p <- plot_ly(source = "stage_src") %>%
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
      event_register("plotly_click") %>%
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
    tick()
    sizes <- seq(10, 150, by = 2)
    prob_m <- compute_maturity_probability(sizes, "M", 4.0)
    prob_f <- compute_maturity_probability(sizes, "F", 4.0)
    sel_curve <- compute_selectivity(sizes, s50 = as.numeric(input$ctrl_fish_size %||% 95.0), sigma = 5.0)
    p <- plot_ly() %>%
      add_lines(x = sizes, y = prob_m, line = list(color = "#1d4ed8", width = 2), name = "Male Theory (L50=91.8mm)") %>%
      add_lines(x = sizes, y = prob_f, line = list(color = "#dc2626", width = 2), name = "Female Theory (L50=51.6mm)") %>%
      add_lines(x = sizes, y = sel_curve, line = list(color = "#16a34a", width = 1.8, dash = "dash"), name = "Gear Selectivity")
    pos <- sim_positions()
    living <- pos[!pos$state %in% c("dead", "harvested"), , drop = FALSE]
    if (nrow(living) > 3) {
      m_sub <- living[living$sex == "M", , drop = FALSE]
      if (nrow(m_sub) > 2) {
        emp_m <- compute_empirical_ogive_r(m_sub$size, m_sub$mature == "1", bin_width = 10.0)
        if (nrow(emp_m) > 0) {
          p <- p %>% add_markers(x = emp_m$size, y = emp_m$prop,
                                 marker = list(color = "#0284c7", size = 7, symbol = "diamond"),
                                 name = "Empirical Male Maturity")
        }
      }
      f_sub <- living[living$sex == "F", , drop = FALSE]
      if (nrow(f_sub) > 2) {
        emp_f <- compute_empirical_ogive_r(f_sub$size, f_sub$mature == "1", bin_width = 10.0)
        if (nrow(emp_f) > 0) {
          p <- p %>% add_markers(x = emp_f$size, y = emp_f$prop,
                                 marker = list(color = "#e11d48", size = 7, symbol = "diamond"),
                                 name = "Empirical Female Maturity")
        }
      }
    }
    p %>% layout(
      xaxis = list(title = "Carapace Width (mm)"),
      yaxis = list(title = "Probability", range = c(0, 1.05)),
      margin = list(l = 45, r = 15, t = 15, b = 35)
    )
  })

  output$size_frequency_plot <- renderPlotly({
    tick()
    pos <- sim_positions()
    living <- pos[!pos$state %in% c("dead", "harvested"), , drop = FALSE]
    if (nrow(living) == 0) {
      plot_ly() %>% layout(title = "No agents yet — run the simulation")
    } else {
      n_tot <- max(1, nrow(living))
      s_im_m <- living$size[living$sex == "M" & living$mature != "1"]
      s_ma_m <- living$size[living$sex == "M" & living$mature == "1"]
      s_im_f <- living$size[living$sex == "F" & living$mature != "1"]
      s_ma_f <- living$size[living$sex == "F" & living$mature == "1"]
      km_imm <- kde_density_r(s_im_m); km_imm$y <- km_imm$y * (length(s_im_m) / n_tot)
      km_mat <- kde_density_r(s_ma_m); km_mat$y <- km_mat$y * (length(s_ma_m) / n_tot)
      kf_imm <- kde_density_r(s_im_f); kf_imm$y <- kf_imm$y * (length(s_im_f) / n_tot)
      kf_mat <- kde_density_r(s_ma_f); kf_mat$y <- kf_mat$y * (length(s_ma_f) / n_tot)
      plot_ly() %>%
        add_lines(x = km_imm$x, y = km_imm$y, line = list(color = "#38bdf8", width = 2), name = "Immature Male") %>%
        add_lines(x = km_mat$x, y = km_mat$y, line = list(color = "#1d4ed8", width = 2.2), name = "Mature Male") %>%
        add_lines(x = kf_imm$x, y = kf_imm$y, line = list(color = "#fb923c", width = 2), name = "Immature Female") %>%
        add_lines(x = kf_mat$x, y = kf_mat$y, line = list(color = "#dc2626", width = 2.2), name = "Mature Female") %>%
        layout(
          xaxis = list(title = "Carapace Width (mm)", range = c(0, 165)),
          yaxis = list(title = "KDE Density"),
          shapes = list(
            list(type = "line", x0 = as.numeric(input$ctrl_fish_size %||% 95), x1 = as.numeric(input$ctrl_fish_size %||% 95),
                 y0 = 0, y1 = 1, yref = "paper", line = list(color = "#16a34a", width = 2, dash = "dash"))
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
    if (nrow(h) < 4) {
      plot_ly() %>% layout(title = "Run >= 4 steps to populate Stock-Recruitment")
    } else {
      sp <- pmax(1, h$mature_female[-nrow(h)])
      rec <- pmax(1, (h$foraging + h$resting + h$molting + h$reproducing - h$mature_female - h$mature_male)[-1])
      sr <- fit_stock_recruitment_r(sp, rec)
      p <- plot_ly() %>%
        add_markers(x = sp, y = rec, marker = list(color = "#1d4ed8", size = 7), name = "Observed")
      if (isTRUE(sr$isFitted)) {
        sg <- seq(min(sp), max(sp), length.out = 60)
        bh <- (sr$beverton_holt$alpha * sg) / (1 + sr$beverton_holt$beta * sg)
        rk <- sr$ricker$alpha * sg * exp(-sr$ricker$beta * sg)
        p <- p %>%
          add_lines(x = sg, y = bh, line = list(color = "#16a34a", width = 2), name = sprintf("Beverton-Holt (R2=%.2f)", sr$beverton_holt$r2)) %>%
          add_lines(x = sg, y = rk, line = list(color = "#f59e0b", width = 2, dash = "dash"), name = sprintf("Ricker (R2=%.2f)", sr$ricker$r2))
      }
      p %>% layout(
        xaxis = list(title = "Spawning Stock (Mature Females)"),
        yaxis = list(title = "Recruits / Juveniles"),
        margin = list(l = 45, r = 15, t = 15, b = 35)
      )
    }
  })

  output$macro_fits_plot <- renderPlotly({
    input$btn_refit_macro
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
      bio <- (annual$foraging + annual$resting + annual$molting + annual$reproducing) * mean_mass
      cat_kg <- annual$harvested * mean_mass
      p <- plot_ly(x = annual$year) %>%
        add_lines(y = bio, line = list(color = "#1d4ed8", width = 2), name = "Survey biomass (kg)") %>%
        add_bars(y = cat_kg, marker = list(color = "rgba(220, 38, 38, 0.5)"), name = "Catch (kg)")
      if (length(bio) >= 4) {
        sch <- fit_schaefer_r(bio, cat_kg)
        dd <- fit_delay_difference_r(bio, cat_kg)
        p <- p %>%
          add_lines(y = sch$predicted, line = list(color = "#16a34a", width = 2, dash = "dash"),
                    name = sprintf("Schaefer (R2=%.2f)", sch$R2)) %>%
          add_lines(y = dd$predicted, line = list(color = "#9333ea", width = 2, dash = "dot"),
                    name = sprintf("Delay-Diff (R2=%.2f)", dd$R2))
      }
      p %>% layout(
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
                       Legal_Males = integer(0), Harvest_Rate = character(0))
    } else {
      h$year <- ceiling(h$times / 365)
      sp_list <- split(h, h$year)
      prev_cum_harv <- 0L
      rows <- lapply(sp_list, function(g) {
        last <- g[nrow(g), ]
        alive <- last$foraging + last$resting + last$molting + last$reproducing
        ann_harv <- max(0L, as.integer(last$harvested) - prev_cum_harv)
        prev_cum_harv <<- as.integer(last$harvested)
        avail_legal <- last$legal_male + ann_harv
        harvest_rate <- if (avail_legal > 0) {
          sprintf("%.1f%%", 100 * ann_harv / avail_legal)
        } else "0.0%"
        data.frame(Year = paste("Year", last$year), End_Day = last$times,
                   Alive = alive, Harvested = last$harvested,
                   Mature_Females = last$mature_female, Legal_Males = last$legal_male,
                   Harvest_Rate = harvest_rate)
      })
      df <- do.call(rbind, rows)
    }
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })
}

shinyApp(ui, server)
