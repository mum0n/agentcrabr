# Crab Agent-Based Model (TethysBio-faithful)
# All sizes in mm (carapace width), matching TethysBio reference implementation.

library(ABM)

# --- Helper functions (defined outside for scope) ---

calculate_mass <- function(size, sex, is_mature) {
  if (sex == "M") {
    alpha <- if (is_mature) 0.0001830 else 0.0002665
    beta <- if (is_mature) 3.160 else 3.098
  } else {
    alpha <- if (is_mature) 0.0004500 else 0.0003440
    beta <- if (is_mature) 2.920 else 2.980
  }
  alpha * (size ^ beta)
}

clamp <- function(x, lo, hi) max(lo, min(hi, x))

#' Create a crab simulation with TethysBio-faithful biological rules
create_crab_simulation <- function(n_crabs = 50, params = list()) {
  p <- list(
    # --- Grid / Spatial ---
    grid_nx = 100, grid_ny = 100, world_size = 1.0, no_take_proportion = 0.0,

    # --- General / Life History ---
    energy_init = 10.0, energy_max = 20.0,
    size_init_min = 20.0, size_init_max = 40.0,
    legal_size_male = 95.0, maturity_size_gate = 35.0,
    juvenile_max_stage = 6, max_stage_male = 19, max_stage_female = 16,
    adult_size_threshold = 68.0,

    # --- Sex-specific growth (males grow faster/larger) ---
    growth_rate_male = 0.15, growth_rate_female = 0.10,
    max_size_male = 160.0, max_size_female = 95.0,

    # --- Sex-specific maturity (L50, logistic ogive, mm) ---
    maturity_l50_male = 91.8, maturity_l50_female = 51.6,
    maturity_slope_male = 6.9957, maturity_slope_female = 12.9812,

    # --- Fecundity (females only) ---
    fecundity_intercept = 2.0, fecundity_slope = 2.5,

    # --- Fishing (males only) ---
    fishing_probability = 0.05, selectivity_s50 = 95.0, selectivity_sigma = 5.0,

    # --- Energetics (Kleiber + Holling II) ---
    energy_gain_foraging = 0.3, energy_cost_resting = 0.05,
    energy_cost_molting = 3,
    metabolic_rate_base = 0.01, metabolic_mass_exponent = 0.75,
    metabolic_temperature_coefficient = 0.08,
    metabolic_optimal_temperature = 3.0, metabolic_optimal_sigma = 2.5,
    assimilation_efficiency = 0.75, specific_dynamic_action = 0.15,
    half_saturation_food = 2.0, intake_rate_base = 0.05,

    # --- Molting ---
    molt_growth_male = 0.12, molt_growth_female = 0.08,
    molting_rate = 0.05, min_intermolt_days = 30,
    molt_season_start = 15, molt_season_end = 140,
    molt_thermal_onset = 10.0, molt_thermal_bandwidth = 3.0,
    degree_day_base = 2.0,

    # --- Mortality ---
    mortality_rate = 0.02, predation_rate = 0.05,
    starvation_rate = 0.20, thermal_rate = 0.05,
    senescence_rate = 0.05, cannibalism_rate = 0.05,
    repro_rate = 0.1,
    carrying_capacity = 12000, carrying_capacity_exponent = 2.0,

    # --- Handling mortality ---
    handling_mortality_days = 365, handling_mortality_rate = 0.05
  )
  p <- modifyList(p, params)

  # Ensure all parameters are numeric
  for (nm in names(p)) {
    if (is.numeric(p[[nm]])) {
      p[[nm]] <- as.numeric(p[[nm]])
    }
  }

  get_max_size <- function(org) {
    if (org$sex == "M") p$max_size_male else p$max_size_female
  }

  get_maturity_prob <- function(org) {
    if (org$size < p$maturity_size_gate) return(0.0)
    l50 <- if (org$sex == "M") p$maturity_l50_male else p$maturity_l50_female
    slope <- if (org$sex == "M") p$maturity_slope_male else p$maturity_slope_female
    1.0 / (1.0 + exp(-slope * (org$size - l50)))
  }

  is_legal_male <- function(org) {
    org$sex == "M" && org$mature == "1" && org$size >= p$legal_size_male && org$legal == "1"
  }

  # --- Energetics: Kleiber metabolism + Holling Type II feeding ---
  calc_energy_change <- function(org, temp, food) {
    mass <- org$mass
    mass_scaling <- mass ^ p$metabolic_mass_exponent
    metabolic_temp_factor <- exp(p$metabolic_temperature_coefficient * temp)
    t_opt <- p$metabolic_optimal_temperature
    t_sigma <- p$metabolic_optimal_sigma
    thermal_modifier <- exp(-((temp - t_opt)^2) / (2.0 * t_sigma^2))
    maintenance_cost <- p$metabolic_rate_base * mass_scaling *
                        metabolic_temp_factor * thermal_modifier

    intake_demand <- p$intake_rate_base * (food / (food + p$half_saturation_food)) * mass_scaling
    actual_intake <- min(intake_demand, max(0, food))
    assim_eff <- p$assimilation_efficiency
    sda <- p$specific_dynamic_action
    usable_fraction <- 1.0 - sda
    net_assimilation <- (actual_intake * assim_eff * usable_fraction) - maintenance_cost

    return(net_assimilation)
  }

  # --- Molting probability ---
  molt_probability <- function(org, temp, day_of_year) {
    if (org$mature == "1") return(0.0)

    days_since_molt <- org$age_days - org$last_molt_time
    min_intermolt <- p$min_intermolt_days
    if (days_since_molt < min_intermolt) return(0.0)

    min_dd <- p$min_intermolt_days * p$degree_day_base
    if (org$degree_days < min_dd) return(0.0)

    if (org$stage >= p$juvenile_max_stage) {
      in_season <- p$molt_season_start <= day_of_year && day_of_year <= p$molt_season_end
      if (!in_season) return(0.0)
    }

    req_e <- max(p$energy_cost_molting * 0.1, p$energy_cost_molting * 0.05 * org$size)
    if (org$energy_reserve < req_e) return(0.0)

    onset <- p$molt_thermal_onset
    bw <- p$molt_thermal_bandwidth
    therm_inhib <- if (temp > onset && bw > 0) {
      max(0, 1 - ((temp - onset) / bw)^2)
    } else 1.0
    if (therm_inhib <= 0) return(0.0)

    progress <- min(org$degree_days / min_dd, 1.0)
    lambda <- 0.15 * progress * therm_inhib
    return(1 - exp(-lambda))
  }

  # --- Mortality components ---
  mortality_components <- function(org, temp, local_density, day_of_year) {
    dt <- 1.0
    dpy <- 365.0

    q_mort <- org$tq_mortality
    individual_mort_factor <- exp(0.25 * qnorm(q_mort))

    base_rate_by_stage <- c(
      0.35, 0.35, 0.35, 0.35,
      0.15, 0.15, 0.15, 0.15,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
      0.06
    )
    base_rate <- base_rate_by_stage[min(org$stage, length(base_rate_by_stage))]
    base_daily <- (base_rate / dpy) * individual_mort_factor * 1.0 * dt

    size_excess <- max(0, org$size - 5.0)
    predation_risk <- ((0.05 / dpy) * exp(-0.03 * size_excess)) * dt

    days_since_molt <- max(0, org$age_days - org$last_molt_time)
    soft_shell_risk <- if (days_since_molt > 30) 0.0 else
      ((0.12 / dpy) * exp(-days_since_molt / 15)) * dt

    starvation_risk <- if (org$energy_reserve < 0) {
      deficit <- -org$energy_reserve
      crit_pool <- max(0.5, org$mass * 0.05)
      ((0.20 / dpy) * (1 - exp(-deficit / crit_pool))) * dt
    } else 0.0

    thermal_risk <- if (temp > 8.0) {
      raw <- ((0.05 / dpy) * ((temp - 8.0)^1.5)) * dt
      clamp(raw, 0, 1)
    } else 0.0

    terminal <- if (org$sex == "M") p$max_stage_male else p$max_stage_female
    senescence_risk <- if (org$mature == "1" && org$stage >= terminal) {
      age_years <- max(0, (org$age_days - 180 * (org$stage - 6))) / 365
      thresh <- 3.0
      if (age_years > thresh) {
        lambda_per_year <- 0.05 * exp(0.80 * (age_years - thresh))
        raw <- -expm1(-lambda_per_year * dt / dpy)
        clamp(raw, 0, 1)
      } else 0.0
    } else 0.0

    cannibalism_risk <- if (org$size < 35 && local_density > 0) {
      ((0.05 / dpy) *
       min(local_density / 10.0, 2.0)) * dt
    } else 0.0

    raw_total <- base_daily + predation_risk + soft_shell_risk +
      starvation_risk + thermal_risk + senescence_risk + cannibalism_risk
    total_prob <- clamp(raw_total, 0, 1)

    return(list(
      baseline = base_daily, predation = predation_risk,
      soft_shell = soft_shell_risk, starvation = starvation_risk,
      thermal = thermal_risk, senescence = senescence_risk,
      cannibalism = cannibalism_risk, total = total_prob
    ))
  }

  sim <- Simulation$new(as.integer(n_crabs), function(i) {
    sex <- sample(c("M", "F"), 1)[1]
    stage <- as.integer(p$juvenile_max_stage + sample.int(5, 1)[1])
    max_stage <- if (sex == "M") p$max_stage_male else p$max_stage_female
    stage <- as.integer(min(stage, max_stage))
    size_init <- as.numeric(runif(1, p$size_init_min, p$size_init_max)[1])
    if (sex == "M") size_init <- as.numeric(min(size_init * 1.2, p$max_size_male))
    else size_init <- as.numeric(min(size_init, p$max_size_female))

    tq_growth <- as.numeric(runif(1, 0.01, 0.99)[1])
    tq_maturity <- as.numeric(runif(1, 0.01, 0.99)[1])
    tq_fecundity <- as.numeric(runif(1, 0.01, 0.99)[1])
    tq_mortality <- as.numeric(runif(1, 0.01, 0.99)[1])

    l50 <- if (sex == "M") p$maturity_l50_male else p$maturity_l50_female
    slope <- if (sex == "M") p$maturity_slope_male else p$maturity_slope_female
    mat_prob <- if (size_init < p$maturity_size_gate) 0.0 else
      1.0 / (1.0 + exp(-slope * (size_init - l50)))
    is_mature <- if (stage >= max_stage) "1"
      else if (runif(1)[1] < mat_prob) "1" else "0"

    legal <- if (sex == "M" && is_mature == "1" && size_init >= p$legal_size_male) "1" else "0"

    mass <- calculate_mass(size_init, sex, is_mature == "1")

    init_dd <- as.numeric(runif(1, p$degree_day_base * 25, p$degree_day_base * 150)[1])

    list(
      state = "foraging",
      energy = as.numeric(p$energy_init),
      size = as.numeric(size_init),
      stage = as.integer(stage),
      mature = as.character(is_mature),
      sex = as.character(sex),
      x = as.numeric(runif(1)[1]),
      y = as.numeric(runif(1)[1]),
      age_days = as.integer(0),
      degree_days = as.numeric(init_dd),
      energy_reserve = as.numeric(p$energy_init),
      last_molt_time = as.integer(0),
      egg_count = as.integer(0),
      sperm_clutches = as.integer(0),
      hatch_day = as.integer(0),
      is_mating_ready = as.character(if (is_mature == "1") "1" else "0"),
      legal = as.character(legal),
      mass = as.numeric(mass),
      tq_growth = as.numeric(tq_growth),
      tq_maturity = as.numeric(tq_maturity),
      tq_fecundity = as.numeric(tq_fecundity),
      tq_mortality = as.numeric(tq_mortality)
    )
  })

  # --- Transitions ---

  sim$addTransition(
    list(state = "foraging") -> list(state = "resting"),
    newExpWaitingTime(1 / p$resting_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      setState(agent, list(energy = as.numeric(max(0, st$energy - p$energy_cost_resting))))
    }
  )

  sim$addTransition(
    list(state = "resting") -> list(state = "foraging"),
    newExpWaitingTime(1 / p$foraging_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      setState(agent, list(energy = as.numeric(min(p$energy_max, st$energy + p$energy_gain_foraging))))
    }
  )

  sim$addTransition(
    list(state = "foraging") -> list(state = "molting"),
    newExpWaitingTime(1 / p$molting_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      day_of_year <- as.integer(time %% 365)
      mp <- molt_probability(st, 8.0, day_of_year)
      if (runif(1)[1] < mp) {
        growth <- if (st$sex == "M") p$molt_growth_male else p$molt_growth_female
        max_sz <- get_max_size(st)
        new_size <- as.numeric(min(max_sz, st$size + growth))
        new_energy <- as.numeric(max(0, st$energy - p$energy_cost_molting))
        mat_prob <- get_maturity_prob(list(size = new_size, sex = st$sex, mature = st$mature))
        new_mature <- if (st$mature == "1") "1" else if (runif(1)[1] < mat_prob) "1" else "0"
        new_legal <- if (st$sex == "M" && new_mature == "1" && new_size >= p$legal_size_male) "1" else "0"
        new_stage <- as.integer(min(st$stage + 1, if (st$sex == "M") p$max_stage_male else p$max_stage_female))
        new_mass <- as.numeric(calculate_mass(new_size, st$sex, new_mature == "1"))
        setState(agent, list(
          size = new_size, energy = new_energy, mature = new_mature,
          legal = new_legal, stage = new_stage, mass = new_mass,
          last_molt_time = as.integer(st$age_days), degree_days = as.numeric(0)
        ))
      } else {
        setState(agent, list(state = "foraging"))
      }
    }
  )

  sim$addTransition(
    list(state = "molting") -> list(state = "foraging"),
    newExpWaitingTime(2),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      gain <- if (st$sex == "M") p$energy_gain_foraging * 0.5 else p$energy_gain_foraging * 0.3
      setState(agent, list(energy = as.numeric(min(p$energy_max, st$energy + gain))))
    }
  )

  sim$addTransition(
    list(state = "foraging") -> list(state = "dead"),
    newExpWaitingTime(1 / p$mortality_rate),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      st$energy <= 0
    }
  )

  sim$addTransition(
    list(state = "resting") -> list(state = "dead"),
    newExpWaitingTime(1 / (p$mortality_rate * 0.5)),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      st$energy <= 0
    }
  )

  sim$addTransition(
    list(state = "foraging", legal = "1") -> list(state = "harvested"),
    newExpWaitingTime(1 / p$fishing_probability),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      is_legal_male(st)
    }
  )

  sim$addTransition(
    list(state = "foraging", sex = "F", mature = "1") -> list(state = "reproducing"),
    newExpWaitingTime(1 / p$repro_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      fec <- as.numeric(2.0 * (st$size ^ 2.725) * 10^0.7311)
      fec <- as.numeric(max(1000, fec))
      setState(agent, list(egg_count = as.integer(round(fec))))
    }
  )

  sim$addTransition(
    list(state = "reproducing") -> list(state = "foraging"),
    newExpWaitingTime(2),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      setState(agent, list(egg_count = as.integer(0)))
    }
  )

  # --- Loggers (all string values for mature/legal) ---
  sim$addLogger(newCounter("foraging", list(state = "foraging")))
  sim$addLogger(newCounter("resting", list(state = "resting")))
  sim$addLogger(newCounter("molting", list(state = "molting")))
  sim$addLogger(newCounter("dead", list(state = "dead")))
  sim$addLogger(newCounter("harvested", list(state = "harvested")))
  sim$addLogger(newCounter("reproducing", list(state = "reproducing")))
  sim$addLogger(newCounter("mature_male", list(sex = "M", mature = "1")))
  sim$addLogger(newCounter("mature_female", list(sex = "F", mature = "1")))
  sim$addLogger(newCounter("legal_male", list(legal = "1")))

  sim
}

#' Run crab simulation
run_crab_simulation <- function(sim, times = seq(0, 50, by = 1)) {
  res <- sim$run(as.numeric(times))
  df <- as.data.frame(res)
  positions <- data.frame(x = numeric(0), y = numeric(0), state = character(0),
                          sex = character(0), size = numeric(0), mature = character(0),
                          mass = numeric(0), stringsAsFactors = FALSE)
  for (i in seq_len(sim$size)) {
    ag <- sim$agent(i)
    st <- getState(ag)
    positions <- rbind(positions, data.frame(
      x = as.numeric(st$x), y = as.numeric(st$y), state = as.character(st$state),
      sex = as.character(st$sex), size = as.numeric(st$size), mature = as.character(st$mature),
      mass = as.numeric(st$mass), stringsAsFactors = FALSE
    ))
  }
  assign("last_positions", positions, envir = .GlobalEnv)
  df
}

compute_maturity_ogive <- function(sizes, l50, slope) {
  1.0 / (1.0 + exp(-slope * (sizes - l50)))
}

compute_selectivity <- function(sizes, s50, sigma) {
  1.0 / (1.0 + exp(-(sizes - s50) / sigma))
}

aoi_bottom_temperature <- function(depth, lat) {
  z <- depth; phi <- lat
  2.0 + 2.5 * exp(-((z - 150.0) / 60.0)^2) + 0.3 * (45.0 - phi)
}

aoi_habitat_suitability <- function(temp, depth, mult = 1.0, offset = 0.0) {
  t_opt <- exp(-0.5 * ((temp - 2.0) / 2.2)^2)
  z_opt <- exp(-0.5 * ((depth - 150.0) / 75.0)^2)
  raw_hsi <- t_opt * z_opt
  realized <- raw_hsi * mult + offset
  max(0, min(1, realized))
}

aoi_domain_bounds <- function(lon, lat, span_meters = 1000) {
  half_m <- span_meters / 2.0
  lat_deg_len <- 111139.0
  lon_deg_len <- max(1000.0, 111139.0 * cos(lat * pi / 180))
  dlat <- half_m / lat_deg_len
  dlon <- half_m / lon_deg_len
  list(lon = c(lon - dlon, lon + dlon), lat = c(lat - dlat, lat + dlat))
}
