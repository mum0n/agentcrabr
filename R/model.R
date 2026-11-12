# Crab Agent-Based Model (TethysBio Parity Architecture)
# agentcrabr: Agent-based model for snow crab population dynamics
# Biological rules, allometry, and phenology matching TethysBio (Chionoecetes opilio)

library(ABM)

#' Determine if an organism belongs to the empirically measured size fraction
#'
#' Evaluates whether a crab agent belongs to the surveyed benthic size fraction
#' (Instar 8 and older, carapace width >= 35.0 mm, stage >= 11, or mature).
#'
#' @param org A list or object containing crab attributes (stage, size, mature).
#' @return Logical scalar indicating if the crab is empirically measured.
#' @export
is_empirically_measured <- function(org) {
  if (!is.null(org$stage) && as.numeric(org$stage) >= 11) {
    return(TRUE)
  }
  if (!is.null(org$size) && as.numeric(org$size) >= 35.0) {
    return(TRUE)
  }
  if (!is.null(org$mature) && as.character(org$mature) == "1") {
    return(TRUE)
  }
  return(FALSE)
}

#' Calculate benthic carrying capacity K for empirically measured crabs
#'
#' Derives carrying capacity K from user-specified density (default 6000 / km^2)
#' and physical domain area:
#' \deqn{K = \text{density} \times \text{Area}_{\text{km}^2}}
#' where \eqn{\text{Area}_{\text{km}^2} = (N_{\text{grid}} \times u_{\text{meters}} / 1000)^2}.
#'
#' @param density Target density in individuals per km^2 (default 6000.0).
#' @param grid_size Number of cells along domain dimension (default 50).
#' @param unit_meters Physical cell length in meters (default 1.0).
#' @return Integer carrying capacity K (at least 1).
#' @export
benthic_carrying_capacity <- function(density = 6000.0, grid_size = 50,
                                      unit_meters = 1.0) {
  if (density < 0.0) {
    stop("benthic_carrying_capacity: density must be non-negative.")
  }
  if (grid_size <= 0 || unit_meters <= 0.0) {
    stop("benthic_carrying_capacity: grid_size and unit_meters must be positive.")
  }
  area_km2 <- ((as.numeric(grid_size) * as.numeric(unit_meters)) / 1000.0)^2
  max(1L, round(as.numeric(density) * area_km2))
}

#' Compute wet mass in kilograms using sex- and maturity-specific allometry
#'
#' Computes individual wet mass via bipartite allometry:
#' \deqn{W = \frac{\alpha \cdot \text{CW}^\beta}{1000} \quad (\text{kg})}
#' using parameter estimates from Sainte-Marie and DFO assessments:
#' \itemize{
#'   \item Male immature: \eqn{\alpha = 0.0002665, \beta = 3.098}
#'   \item Male mature: \eqn{\alpha = 0.000183, \beta = 3.16}
#'   \item Female immature: \eqn{\alpha = 0.000344, \beta = 2.98}
#'   \item Female mature: \eqn{\alpha = 0.00045, \beta = 2.92}
#' }
#'
#' @param cw Carapace width in millimeters.
#' @param sex Character string, "M" for male or "F" for female.
#' @param is_mature Logical or character ("1"/"0") indicating morphometric maturity.
#' @return Numeric mass in kilograms.
#' @export
calculate_crab_mass_kg <- function(cw, sex, is_mature = FALSE) {
  if (cw <= 0.0) {
    stop("calculate_crab_mass_kg: cw must be positive.")
  }
  mat_flag <- isTRUE(is_mature) || as.character(is_mature) == "1"
  if (sex == "M") {
    if (mat_flag) {
      alpha <- 0.000183
      beta <- 3.16
    } else {
      alpha <- 0.0002665
      beta <- 3.098
    }
  } else if (sex == "F") {
    if (mat_flag) {
      alpha <- 0.00045
      beta <- 2.92
    } else {
      alpha <- 0.000344
      beta <- 2.98
    }
  } else {
    stop(paste("calculate_crab_mass_kg: unrecognized sex", sex))
  }
  (alpha * (cw^beta)) / 1000.0
}

#' Compute post-molt carapace width using Hiatt logarithmic formulation
#'
#' Predicts post-molt size following the sex-specific Hiatt model with maturity
#' penalty and thermal sensitivity:
#' \deqn{\ln(\text{CW}_{t+1}) = \beta_0 + \beta_1 \ln(\text{CW}_t) + \delta_{\text{mat}} I_{\text{mat}} + \gamma_T (T - T_{\text{ref}})}
#' Parameter values match TethysBio:
#' \itemize{
#'   \item Male: \eqn{\beta_0 = 0.30153, \beta_1 = 1.00027, \delta_{\text{mat}} = -0.13564, \text{max} = 160.0} mm
#'   \item Female: \eqn{\beta_0 = 0.1931, \beta_1 = 1.0372, \delta_{\text{mat}} = -0.131, \text{max} = 95.0} mm
#'   \item Thermal coefficient: \eqn{\gamma_T = -0.004, T_{\text{ref}} = 4.0^\circ\text{C}}
#' }
#'
#' @param cw Current carapace width in mm.
#' @param sex Character string, "M" or "F".
#' @param is_mature Logical indicating morphometric maturity.
#' @param temp Ambient bottom temperature in degrees Celsius (default 4.0).
#' @return Predicted post-molt carapace width in mm.
#' @export
hiatt_post_growth_cw <- function(cw, sex, is_mature = FALSE, temp = 4.0) {
  if (cw <= 0.0) {
    stop("hiatt_post_growth_cw: cw must be positive.")
  }
  mat_flag <- if (isTRUE(is_mature) || as.character(is_mature) == "1") 1.0 else 0.0
  if (sex == "M") {
    beta0 <- 0.30153
    beta1 <- 1.00027
    delta_mat <- -0.13564
    max_cw <- 160.0
  } else if (sex == "F") {
    beta0 <- 0.1931
    beta1 <- 1.0372
    delta_mat <- -0.131
    max_cw <- 95.0
  } else {
    stop(paste("hiatt_post_growth_cw: unrecognized sex", sex))
  }
  gamma_t <- -0.004
  t_ref <- 4.0
  ln_next <- beta0 + beta1 * log(cw) + delta_mat * mat_flag + gamma_t * (temp - t_ref)
  min(max_cw, exp(ln_next))
}

#' Compute terminal molt maturity probability with 35mm pre-molt size gate
#'
#' Calculates logistic terminal maturity probability with thermal centroid shift:
#' \deqn{P(\text{CW}) = \frac{1}{1 + \exp\left(-\text{slope} \cdot \ln\left(\frac{\text{CW}}{\text{CW}_{50}}\right)\right)}}
#' where \eqn{\Delta\text{CW}_{50} = -2.0 (T - 4.0^\circ\text{C})}. Crabs below
#' 35.0 mm pre-molt carapace width cannot undergo terminal molt (\eqn{P = 0}).
#'
#' @param cw Current carapace width in mm.
#' @param sex Character string, "M" or "F".
#' @param temp Ambient bottom temperature in degrees Celsius (default 4.0).
#' @return Probability of transitioning to terminal morphometric maturity [0, 1].
#' @export
compute_maturity_probability <- function(cw, sex, temp = 4.0) {
  if (cw < 35.0) {
    return(0.0)
  }
  t_ref <- 4.0
  delta_t <- -2.0 * (temp - t_ref)
  if (sex == "M") {
    l50 <- 91.8 + delta_t
    slope <- 6.9957
  } else if (sex == "F") {
    l50 <- 51.6 + delta_t
    slope <- 12.9812
  } else {
    stop(paste("compute_maturity_probability: unrecognized sex", sex))
  }
  if (l50 <= 0.0) {
    stop("compute_maturity_probability: effective L50 must be positive.")
  }
  1.0 / (1.0 + exp(-slope * log(cw / l50)))
}

#' Compute gear selectivity probability
#'
#' Logistic selectivity curve for commercial snow crab traps:
#' \deqn{S(\text{CW}) = \frac{1}{1 + \exp\left(-\frac{\text{CW} - S_{50}}{\sigma}\right)}}
#'
#' @param sizes Numeric vector of carapace widths in mm.
#' @param s50 Fifty-percent selectivity retention width (default 95.0 mm).
#' @param sigma Slope parameter (default 5.0 mm).
#' @return Numeric vector of gear retention probabilities.
#' @export
compute_selectivity <- function(sizes, s50 = 95.0, sigma = 5.0) {
  if (sigma <= 0.0) {
    stop("compute_selectivity: sigma must be positive.")
  }
  1.0 / (1.0 + exp(-(sizes - s50) / sigma))
}

#' Compute benthic habitat suitability index (HSI)
#'
#' Evaluates composite suitability from bottom temperature and seabed depth:
#' \deqn{\text{HSI} = \exp\left(-0.5 \left(\frac{T - 2.5}{1.5}\right)^2\right) \times \exp\left(-0.5 \left(\frac{Z - 150}{75}\right)^2\right)}
#'
#' @param temp Ambient bottom temperature in degrees Celsius.
#' @param depth Seabed bathymetric depth in meters.
#' @param offset Additive offset applied to HSI (default 0.0).
#' @param mult Multiplicative scaling factor applied to HSI (default 1.0).
#' @return Numeric HSI value bounded in [0.0, 1.0].
#' @export
aoi_habitat_suitability <- function(temp, depth, offset = 0.0, mult = 1.0) {
  t_val <- as.numeric(temp)
  z_val <- as.numeric(depth)
  thermal <- exp(-0.5 * ((t_val - 2.5) / 1.5)^2)
  depth_fit <- exp(-0.5 * ((z_val - 150.0) / 75.0)^2)
  raw <- thermal * depth_fit * mult + offset
  max(0.0, min(1.0, raw))
}

#' Compute seasonal bottom temperature cycle
#'
#' Generates sinusoidal annual climatology at a designated depth and latitude:
#' \deqn{T(d) = T_{\text{mean}} + A \cos\left(2\pi \frac{d - 80.0}{365.25}\right)}
#' where \eqn{T_{\text{mean}} = 2.5 + 3.2 e^{-Z/85} - 0.35 (\text{Lat} - 44.5) + \text{offset}}
#' and \eqn{A = \max(0.2, 3.5 e^{-Z/120}) \times \text{multiplier}}.
#'
#' @param days Numeric vector of simulation days (0 to 365).
#' @param depth Seabed depth in meters (default 150.0).
#' @param lat Latitude in degrees North (default 45.0).
#' @param offset Additive temperature offset in degrees Celsius (default 0.0).
#' @param mult Seasonal amplitude multiplier (default 1.0).
#' @return A list containing `days`, `temps`, `mean_temp`, and `amplitude`.
#' @export
seasonal_bottom_temperature_cycle <- function(days = seq(0, 365, by = 2),
                                              depth = 150.0, lat = 45.0,
                                              offset = 0.0, mult = 1.0) {
  base_mean <- 2.5 + 3.2 * exp(-depth / 85.0) - 0.35 * (lat - 44.5)
  mean_t <- base_mean + offset
  amp <- max(0.2, 3.5 * exp(-depth / 120.0)) * mult
  temps <- mean_t + amp * cos(2 * pi * (days - 80.0) / 365.25)
  list(days = days, temps = temps, mean_temp = mean_t, amplitude = amp)
}

#' Create a crab simulation with sex-specific biological rules
#'
#' Instantiates an agent-based simulation adhering to TethysBio contracts:
#' Hiatt molt increments, size-gated maturity ogives, multi-hazard mortality
#' with quadratic crowding feedback on measured instars 8+, and commercial
#' fishing season restriction to mature legal males (CW >= 95 mm) outside
#' active exclusion areas.
#'
#' @param n_crabs Initial number of agents (default 50).
#' @param params List of parameter overrides.
#' @return An ABM `Simulation` instance.
#' @export
create_crab_simulation <- function(n_crabs = 50, params = list()) {
  p <- list(
    # Spatial domain & Carrying capacity
    unit_meters = 1.0,
    grid_size = 50,
    benthic_carrying_capacity_density = 6000.0,
    benthic_carrying_capacity_exponent = 2.0,
    food_availability = 1.0,
    temp_offset = 0.0,
    temp_multiplier = 1.0,
    # Initial size range
    size_init_min = 20.0,
    size_init_max = 40.0,
    legal_size_male = 95.0,
    # Energy
    energy_init = 10.0,
    energy_max = 20.0,
    energy_gain_foraging = 0.3,
    energy_cost_resting = 0.05,
    energy_cost_molting = 3.0,
    # Rates
    foraging_rate = 0.5,
    resting_rate = 0.3,
    molting_rate = 0.05,
    mortality_rate = 0.002,
    repro_rate = 0.1,
    # Fishing
    fishing_probability = 0.002,
    fishing_season_start = 90,
    fishing_season_end = 150,
    # Spatial MPA refuge and size-structure adjustments
    cw_shift = 0.0,
    mpa = list(xmin = 0.25, xmax = 0.55, ymin = 0.25, ymax = 0.55, active = TRUE)
  )
  p <- modifyList(p, params)

  K <- benthic_carrying_capacity(
    p$benthic_carrying_capacity_density,
    p$grid_size,
    p$unit_meters
  )

  sim <- Simulation$new(n_crabs, function(i) {
    sex <- sample(c("M", "F"), 1)
    size_init <- if (sex == "M") runif(1, p$size_init_min, p$size_init_max * 1.2)
                 else runif(1, p$size_init_min, p$size_init_max)
    size_init <- max(2.5, size_init + as.numeric(p$cw_shift))
    is_measured <- if (size_init >= 35.0) "1" else "0"
    list(
      state = "foraging",
      energy = p$energy_init,
      size = size_init,
      mature = "0",
      sex = sex,
      x = runif(1),
      y = runif(1),
      age = 0,
      egg_count = 0,
      legal = "0",
      measured = is_measured
    )
  })

  # Foraging -> Resting
  sim$addTransition(
    list(state = "foraging") -> list(state = "resting"),
    newExpWaitingTime(1 / p$resting_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      new_energy <- max(0, st$energy - p$energy_cost_resting)
      mat_prob <- compute_maturity_probability(st$size, st$sex, 4.0)
      new_mature <- if (st$mature == "1") "1" else if (runif(1) < mat_prob) "1" else "0"
      new_legal <- if (st$sex == "M" && new_mature == "1" && st$size >= p$legal_size_male) "1" else "0"
      new_measured <- if (st$size >= 35.0 || new_mature == "1") "1" else "0"
      setState(agent, list(
        energy = new_energy,
        mature = new_mature,
        legal = new_legal,
        measured = new_measured
      ))
    }
  )

  # Resting -> Foraging
  sim$addTransition(
    list(state = "resting") -> list(state = "foraging"),
    newExpWaitingTime(1 / p$foraging_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      effective_gain <- p$energy_gain_foraging * as.numeric(p$food_availability)
      new_energy <- min(p$energy_max, st$energy + effective_gain)
      setState(agent, list(energy = new_energy))
    }
  )

  # Foraging -> Molting (Hiatt growth formulation)
  sim$addTransition(
    list(state = "foraging") -> list(state = "molting"),
    newExpWaitingTime(1 / p$molting_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      if (st$mature == "1") {
        # Terminally mature adults do not molt
        return()
      }
      new_size <- hiatt_post_growth_cw(st$size, st$sex, is_mature = FALSE, temp = 4.0)
      new_energy <- max(0, st$energy - p$energy_cost_molting)
      mat_prob <- compute_maturity_probability(new_size, st$sex, 4.0)
      new_mature <- if (runif(1) < mat_prob) "1" else "0"
      new_legal <- if (st$sex == "M" && new_mature == "1" && new_size >= p$legal_size_male) "1" else "0"
      new_measured <- if (new_size >= 35.0 || new_mature == "1") "1" else "0"
      setState(agent, list(
        size = new_size,
        energy = new_energy,
        mature = new_mature,
        legal = new_legal,
        measured = new_measured
      ))
    }
  )

  # Molting -> Foraging
  sim$addTransition(
    list(state = "molting") -> list(state = "foraging"),
    newExpWaitingTime(2),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      gain_base <- if (st$sex == "M") p$energy_gain_foraging * 0.5 else p$energy_gain_foraging * 0.3
      effective_gain <- gain_base * as.numeric(p$food_availability)
      new_energy <- min(p$energy_max, st$energy + effective_gain)
      setState(agent, list(energy = new_energy))
    }
  )

  # Foraging -> Dead (starvation, baseline mortality, carrying capacity feedback)
  sim$addTransition(
    list(state = "foraging") -> list(state = "dead"),
    newExpWaitingTime(1 / p$mortality_rate),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      if (st$energy <= 0) return(TRUE)
      # Density-dependent quadratic crowding feedback on measured fraction
      if (identical(st$measured, "1") && K > 0) {
        if (exists("last_positions", envir = .GlobalEnv)) {
          pos <- get("last_positions", envir = .GlobalEnv)
          n_meas <- sum(pos$measured == "1" | pos$measured == 1, na.rm = TRUE)
          crowd_ratio <- (n_meas / K)^p$benthic_carrying_capacity_exponent
          if (crowd_ratio > 1.0 && runif(1) < min(0.8, (crowd_ratio - 1.0) / crowd_ratio)) {
            return(TRUE)
          }
        }
      }
      FALSE
    }
  )

  # Resting -> Dead
  sim$addTransition(
    list(state = "resting") -> list(state = "dead"),
    newExpWaitingTime(1 / (p$mortality_rate * 0.5)),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      st$energy <= 0
    }
  )

  # Harvest (commercial fishing) restricted to mature legal males outside active MPA
  sim$addTransition(
    list(state = "foraging", legal = "1") -> list(state = "harvested"),
    newExpWaitingTime(1 / p$fishing_probability),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      # Check fishing season window (day of year 90 to 150)
      doy <- as.integer(time) %% 365
      if (doy < p$fishing_season_start || doy > p$fishing_season_end) {
        return(FALSE)
      }
      # Check spatial exclusion refuge
      if (!is.null(p$mpa) && isTRUE(p$mpa$active)) {
        if (st$x >= p$mpa$xmin && st$x <= p$mpa$xmax &&
            st$y >= p$mpa$ymin && st$y <= p$mpa$ymax) {
          return(FALSE)
        }
      }
      st$sex == "M" && st$mature == "1" && st$size >= p$legal_size_male && st$legal == "1"
    }
  )

  # Female reproduction (egg production)
  sim$addTransition(
    list(state = "foraging", sex = "F", mature = "1") -> list(state = "reproducing"),
    newExpWaitingTime(1 / p$repro_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      fec <- 0.7311 * (st$size^2.725)
      setState(agent, list(egg_count = round(fec)))
    }
  )

  # Reproducing -> Foraging
  sim$addTransition(
    list(state = "reproducing") -> list(state = "foraging"),
    newExpWaitingTime(2),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      setState(agent, list(egg_count = 0))
    }
  )

  # Loggers
  sim$addLogger(newCounter("foraging", list(state = "foraging")))
  sim$addLogger(newCounter("resting", list(state = "resting")))
  sim$addLogger(newCounter("molting", list(state = "molting")))
  sim$addLogger(newCounter("dead", list(state = "dead")))
  sim$addLogger(newCounter("harvested", list(state = "harvested")))
  sim$addLogger(newCounter("reproducing", list(state = "reproducing")))
  sim$addLogger(newCounter("mature_male", list(sex = "M", mature = "1")))
  sim$addLogger(newCounter("mature_female", list(sex = "F", mature = "1")))
  sim$addLogger(newCounter("legal_male", list(legal = "1")))
  sim$addLogger(newCounter("measured", list(measured = "1")))

  sim
}

#' Run crab simulation and record agent positions
#'
#' Advances an ABM simulation over designated time steps, returning summary
#' time-series and updating `.GlobalEnv$last_positions`.
#'
#' @param sim An ABM `Simulation` instance.
#' @param times Numeric vector of evaluation time points.
#' @return A `data.frame` of logged state counters across time.
#' @export
run_crab_simulation <- function(sim, times = seq(0, 50, by = 1)) {
  res <- sim$run(times)
  df <- as.data.frame(res)
  positions <- data.frame(
    x = numeric(0), y = numeric(0), state = character(0),
    sex = character(0), size = numeric(0), mature = character(0),
    legal = character(0), measured = character(0),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(sim$size)) {
    ag <- sim$agent(i)
    st <- getState(ag)
    positions <- rbind(positions, data.frame(
      x = st$x, y = st$y, state = st$state,
      sex = st$sex, size = st$size, mature = as.character(st$mature),
      legal = as.character(st$legal), measured = as.character(st$measured),
      stringsAsFactors = FALSE
    ))
  }
  assign("last_positions", positions, envir = .GlobalEnv)
  df
}
