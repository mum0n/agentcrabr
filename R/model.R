# Crab Agent-Based Model (TethysBio Parity Architecture)
# agentcrabr: Agent-based model for snow crab population dynamics
# Biological rules, allometry, and phenology matching TethysBio (Chionoecetes opilio)

library(ABM)

# Side-channel for slowly-varying population quantities (live measured-crab
# count) shared with transition callbacks. It is refreshed at every
# run_crab_simulation() snapshot, avoiding per-event reentrant agent reads.
.crab_cache <- new.env(parent = emptyenv())
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

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
#'   \item Male immature: \eqn{\alpha = 0.000183, \beta = 3.16}
#'   \item Male mature: \eqn{\alpha = 0.0002665, \beta = 3.098}
#'   \item Female immature: \eqn{\alpha = 0.000344, \beta = 2.98}
#'   \item Female mature: \eqn{\alpha = 0.00045, \beta = 2.92}
#' }
#'
#' @param cw Carapace width in millimeters.
#' @param sex Character string, "M" for male or "F" for female.
#' @param is_mature Logical or character ("1"/"0") indicating morphometric maturity.
#' @param egg_count Optional brooded egg count for ovigerous females (default 0).
#' @return Numeric mass in kilograms.
#' @export
calculate_crab_mass_kg <- function(cw, sex, is_mature = FALSE, egg_count = 0.0) {
  if (cw <= 0.0) {
    stop("calculate_crab_mass_kg: cw must be positive.")
  }
  mat_flag <- isTRUE(is_mature) || as.character(is_mature) == "1"
  if (sex == "M") {
    if (mat_flag) {
      alpha <- 0.0002665
      beta <- 3.098
    } else {
      alpha <- 0.000183
      beta <- 3.16
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
  mass_g <- alpha * (cw^beta)
  if (sex == "F" && egg_count > 0.0) {
    mass_g <- mass_g + as.numeric(egg_count) * 7.0e-5
  }
  mass_g / 1000.0
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
hiatt_post_growth_cw <- function(cw, sex, is_mature = FALSE, temp = 4.0, condition = 1.0) {
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
  ln_cw <- log(cw)
  ln_next <- beta0 + beta1 * ln_cw + delta_mat * mat_flag +
    0.05 * (condition - 1.0) + gamma_t * (temp - t_ref)
  min_exp <- if (mat_flag > 0.5) 0.06 else 0.12
  ln_next <- max(ln_next, ln_cw + log(1.0 + min_exp))
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
  t_ref <- 4.0
  delta_t <- -2.0 * (temp - t_ref)
  if (sex == "M") {
    l50 <- max(35.0, 91.8 + delta_t)
    slope <- 6.9957
  } else if (sex == "F") {
    l50 <- max(35.0, 51.6 + delta_t)
    slope <- 12.9812
  } else {
    stop(paste("compute_maturity_probability: unrecognized sex", sex))
  }
  if (l50 <= 0.0) {
    stop("compute_maturity_probability: effective L50 must be positive.")
  }
  ifelse(cw < 35.0, 0.0, 1.0 / (1.0 + exp(-slope * log(pmax(1.0, cw) / l50))))
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
#' Evaluates composite suitability from bottom temperature and seabed depth
#' matching TethysBio (open-ended deep tail at depth >= depth_opt, shallow
#' Gaussian + cosine taper below depth_opt, and bounded thermal window
#' [temp_pref_min, temp_pref_max]):
#'
#' @param temp Ambient bottom temperature in degrees Celsius.
#' @param depth Seabed bathymetric depth in meters.
#' @param offset Additive offset applied to HSI (default 0.0).
#' @param mult Multiplicative scaling factor applied to HSI (default 1.0).
#' @param temp_pref_min Minimum viable temperature (default -1.5).
#' @param temp_pref_opt Optimal temperature inflection (default 5.0).
#' @param temp_pref_max Upper thermal limit (default 10.0).
#' @param depth_opt Optimal depth in meters (default 200.0).
#' @param depth_sigma Shallow Gaussian taper width in meters (default 100.0).
#' @return Numeric HSI value bounded in [0.0, 1.0].
#' @export
aoi_habitat_suitability <- function(temp, depth, offset = 0.0, mult = 1.0,
                                    temp_pref_min = -1.5, temp_pref_opt = 5.0,
                                    temp_pref_max = 10.0, depth_opt = 200.0,
                                    depth_sigma = 100.0) {
  t_val <- as.numeric(temp)
  z_val <- as.numeric(depth)
  if (t_val < temp_pref_min || t_val >= temp_pref_max) {
    s_t <- 0.0
  } else {
    s_t <- 1.0 / (1.0 + exp(1.2 * (t_val - temp_pref_opt)))
    if (t_val < temp_pref_min + 1.5 && temp_pref_opt > temp_pref_min) {
      w_low <- min(1.5, max(1e-6, temp_pref_opt - temp_pref_min))
      if (t_val < temp_pref_min + w_low) {
        s_t <- s_t * 0.5 * (1.0 - cos(pi * (t_val - temp_pref_min) / w_low))
      }
    }
    if (t_val > temp_pref_opt) {
      span_hi <- max(1e-6, temp_pref_max - temp_pref_opt)
      s_t <- s_t * 0.5 * (1.0 + cos(pi * (t_val - temp_pref_opt) / span_hi))
    }
  }
  d_opt <- max(25.0, as.numeric(depth_opt))
  d_sig <- max(10.0, as.numeric(depth_sigma))
  if (z_val <= 20.0) {
    s_z <- 0.0
  } else if (z_val >= d_opt) {
    s_z <- 1.0
  } else {
    s_z <- exp(-0.5 * ((z_val - d_opt) / d_sig)^2)
    if (z_val < 50.0) {
      s_z <- s_z * 0.5 * (1.0 - cos(pi * (z_val - 20.0) / 30.0))
    }
  }
  raw <- s_t * s_z * mult + offset
  max(0.0, min(1.0, raw))
}

#' Compute base-10 Moncton fecundity with parity and senescence
#'
#' @param cw Carapace width in mm.
#' @param shell_age_years Shell age in years (default 1.0).
#' @param minimum_clutch_size Minimum viable clutch size (default 1000).
#' @return Expected egg count.
#' @export
calculate_fecundity_moncton <- function(cw, shell_age_years = 1.0,
                                        minimum_clutch_size = 1000.0) {
  if (cw <= 0.0) {
    stop("calculate_fecundity_moncton: cw must be positive.")
  }
  if (shell_age_years <= 2.0) {
    intercept <- -0.7229
    slope <- 2.94
  } else {
    intercept <- 0.0158
    slope <- 2.5811
  }
  eggs <- 10.0^(intercept + slope * log10(cw))
  if (shell_age_years > 3.0) {
    disc <- max(0.6, 1.0 - 0.08 * (shell_age_years - 3.0))
    eggs <- eggs * disc
  }
  max(minimum_clutch_size, eggs)
}

#' Fit dynamic Schaefer surplus-production model to annual observations
#'
#' @param observed Numeric vector of annual biomass or abundance.
#' @param catch Numeric vector of annual catch (default zeros).
#' @return List with `r`, `K`, `MSY`, `R2`, `rmse`, `aic`, `predicted`.
#' @export
fit_schaefer_r <- function(observed, catch = NULL) {
  obs <- as.numeric(observed)
  n <- length(obs)
  if (n < 3) {
    return(list(r = 0, K = 0, MSY = 0, R2 = 0, rmse = Inf, aic = Inf, predicted = obs))
  }
  c_vec <- if (is.null(catch)) rep(0.0, n) else as.numeric(catch)
  sc <- max(1e-6, mean(obs))
  obs_s <- obs / sc
  c_s <- c_vec / sc
  sim_fn <- function(n0, r_v, k_v) {
    p <- numeric(n)
    p[1] <- max(1e-6, n0)
    for (t in seq_len(n - 1)) {
      p[t + 1] <- max(1e-6, p[t] + r_v * p[t] * (1.0 - p[t] / max(1e-6, k_v)) - c_s[t])
    }
    p
  }
  opt <- optim(
    par = c(0.0, log(2.0), log(max(0.1, obs_s[1]))),
    fn = function(th) {
      r_v <- 0.01 + 2.5 / (1.0 + exp(-th[1]))
      k_v <- exp(th[2])
      n0_v <- exp(th[3])
      sum((obs_s - sim_fn(n0_v, r_v, k_v))^2)
    }
  )
  r_hat <- 0.01 + 2.5 / (1.0 + exp(-opt$par[1]))
  k_hat <- exp(opt$par[2]) * sc
  n0_hat <- exp(opt$par[3]) * sc
  pred <- sim_fn(n0_hat / sc, r_hat, k_hat / sc) * sc
  rss <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2) + 1e-12
  list(
    r = r_hat,
    K = k_hat,
    MSY = r_hat * k_hat / 4.0,
    R2 = 1.0 - rss / ss_tot,
    rmse = sqrt(rss / n),
    aic = n * log(max(1e-12, rss / n)) + 6.0,
    predicted = pred
  )
}

#' Fit dynamic Beverton-Holt delay-difference model
#'
#' @param observed Numeric vector of annual biomass or abundance.
#' @param catch Numeric vector of annual catch (default zeros).
#' @param lag Recruitment lag in years (default 1).
#' @return List with `s`, `M`, `alpha`, `beta`, `R2`, `rmse`, `aic`, `predicted`.
#' @export
fit_delay_difference_r <- function(observed, catch = NULL, lag = 1L) {
  obs <- as.numeric(observed)
  n <- length(obs)
  if (n < 3) {
    return(list(s = 0, M = 0, alpha = 0, beta = 0, R2 = 0, rmse = Inf, aic = Inf, predicted = obs))
  }
  c_vec <- if (is.null(catch)) rep(0.0, n) else as.numeric(catch)
  sc <- max(1e-6, mean(obs))
  obs_s <- obs / sc
  c_s <- c_vec / sc
  sim_fn <- function(s_v, a_v, b_v) {
    p <- numeric(n)
    p[1] <- max(1e-6, obs_s[1])
    for (t in seq_len(n - 1)) {
      esc <- max(0.0, p[t] - c_s[t])
      sp <- p[max(1L, t + 1L - lag)]
      p[t + 1] <- max(1e-6, s_v * esc + (a_v * sp) / (1.0 + b_v * sp))
    }
    p
  }
  opt <- optim(
    par = c(0.8, log(0.6), log(0.5)),
    fn = function(th) {
      s_v <- 0.05 + 0.93 / (1.0 + exp(-th[1]))
      a_v <- exp(th[2])
      b_v <- exp(th[3])
      sum((obs_s - sim_fn(s_v, a_v, b_v))^2)
    }
  )
  s_hat <- 0.05 + 0.93 / (1.0 + exp(-opt$par[1]))
  a_hat <- exp(opt$par[2])
  b_hat <- exp(opt$par[3])
  pred <- sim_fn(s_hat, a_hat, b_hat) * sc
  rss <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2) + 1e-12
  list(
    s = s_hat,
    M = -log(max(1e-6, s_hat)),
    alpha = a_hat,
    beta = b_hat / sc,
    R2 = 1.0 - rss / ss_tot,
    rmse = sqrt(rss / n),
    aic = n * log(max(1e-12, rss / n)) + 6.0,
    predicted = pred
  )
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
                                              offset = 0.0, mult = 1.0,
                                              phase_day = 80.0) {
  d_val <- max(10.0, as.numeric(depth))
  base_mean <- 2.0 + 2.5 * exp(-((d_val - 150.0) / 60.0)^2) + 0.3 * (45.0 - as.numeric(lat))
  mean_t <- base_mean + as.numeric(offset)
  amp <- max(0.5, 5.0 * exp(-d_val / 100.0)) * as.numeric(mult)
  temps <- mean_t + amp * sin(2.0 * pi * (days - phase_day) / 365.0)
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
    local_crowding_box_size = 0L,
    food_availability = 1.0,
    temp_offset = 0.0,
    temp_multiplier = 1.0,
    habitat_offset = 0.0,
    habitat_multiplier = 1.0,
    temp_pref_min = -1.5,
    temp_pref_opt = 5.0,
    temp_pref_max = 10.0,
    depth_opt = 200.0,
    depth_sigma = 100.0,
    max_move_m_per_year = 7300.0,
    # Process toggles
    use_energetics = TRUE,
    use_thermal_stress = TRUE,
    use_cannibalism = FALSE,
    use_senescence = TRUE,
    use_movement = TRUE,
    use_reproduction = TRUE,
    # Initial size range
    size_init_min = 20.0,
    size_init_max = 40.0,
    legal_size_male = 95.0,
    # Energy
    energy_init = 10.0,
    energy_max = 20.0,
    energy_gain_foraging = 0.6,
    energy_cost_resting = 0.05,
    energy_cost_molting = 1.0,
    molt_energy_retention = 0.1,
    # Rates
    foraging_rate = 0.5,
    resting_rate = 0.3,
    molting_rate = 0.05,
    mortality_rate = 0.002,
    repro_rate = 0.1,
    move_step = 0.03,
    # Fishing
    fishing_probability = 0.0,
    fishing_harvest_rate_annual = NULL,
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
  .crab_cache$grid_size <- p$grid_size
  .crab_cache$unit_meters <- p$unit_meters
  .crab_cache$k_density <- p$benthic_carrying_capacity_density
  .crab_cache$box_size <- p$local_crowding_box_size

  f_start <- as.numeric(p$fishing_season_start)
  f_end <- as.numeric(p$fishing_season_end)
  season_days <- if (f_start <= f_end) (f_end - f_start + 1.0) else ((365.0 - f_start + 1.0) + f_end)
  f_daily_rate <- if (!is.null(p$fishing_harvest_rate_annual)) {
    h_ann <- min(0.999999, max(0.0, as.numeric(p$fishing_harvest_rate_annual)))
    -log(1.0 - h_ann) / max(1.0, season_days)
  } else {
    max(0.0, as.numeric(p$fishing_probability))
  }

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

  # ABM events are scheduled on state entry and NOT cancelled on state exit,
  # and rules are not re-checked when a stale event fires. Every transition
  # therefore carries a source-state guard so stale events are ignored.
  still_in <- function(expected) {
    force(expected)
    function(time, agent) identical(getState(agent)$state, expected)
  }

  # Foraging -> Resting (rate semantics: newExpWaitingTime takes a per-day rate)
  sim$addTransition(
    list(state = "foraging") -> list(state = "resting"),
    newExpWaitingTime(p$resting_rate),
    to_change_callback = still_in("foraging"),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      new_energy <- max(0, st$energy - p$energy_cost_resting)
      setState(agent, list(energy = new_energy))
    }
  )

  # Resting -> Foraging (feeding gain, random-walk movement, mortality roll)
  sim$addTransition(
    list(state = "resting") -> list(state = "foraging"),
    newExpWaitingTime(p$foraging_rate),
    to_change_callback = still_in("resting"),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      # Mortality is rolled HERE (once per actual feeding resumption) rather
      # than through dedicated death transitions: ABM events scheduled on
      # state entry are never cancelled, so a death *transition* would stack
      # one Exp event per foraging visit and fire repeatedly. Rolling inside
      # a guarded callback keeps the hazard exactly at the parameterized rate.
      q_nat <- 1.0 - exp(-as.numeric(p$mortality_rate) / as.numeric(p$foraging_rate))
      # Density-dependent crowding on the measured fraction (global + C^0 local B x B)
      q_crowd <- 0.0
      loc_r <- local_crowding_ratio_r(st$x, st$y)
      if (identical(st$measured, "1") && K > 0) {
        n_meas <- if (is.null(.crab_cache$n_meas)) 0L else .crab_cache$n_meas
        glob_ratio <- (n_meas / K)^p$benthic_carrying_capacity_exponent
        loc_ratio <- if (loc_r > 1.0) loc_r^p$benthic_carrying_capacity_exponent else 0.0
        crowd_ratio <- max(glob_ratio, loc_ratio)
        if (crowd_ratio > 1.0) {
          q_crowd <- min(0.8, (crowd_ratio - 1.0) / crowd_ratio)
        }
      }
      if (runif(1) < q_nat || runif(1) < q_crowd) {
        decrement_local_crowding_r(st$x, st$y, identical(st$measured, "1"))
        setState(agent, list(state = "dead"))
        return()
      }
      new_energy <- if (isTRUE(p$use_energetics)) {
        effective_gain <- p$energy_gain_foraging * as.numeric(p$food_availability)
        min(p$energy_max, st$energy + effective_gain)
      } else st$energy

      new_x <- as.numeric(st$x)
      new_y <- as.numeric(st$y)
      move <- as.numeric(p$move_step)
      if (isTRUE(p$use_movement) && move > 0) {
        dxs <- c(0, move, -move, 0, 0, move, move, -move, -move)
        dys <- c(0, 0, 0, move, -move, move, -move, move, -move)
        cxs <- pmin(1.0, pmax(0.0, new_x + dxs))
        cys <- pmin(1.0, pmax(0.0, new_y + dys))
        gamma_c <- as.numeric(p$benthic_carrying_capacity_exponent)
        wts <- vapply(seq_along(cxs), function(k) {
          cr <- local_crowding_ratio_r(cxs[k], cys[k])
          sc <- if (cr > 1.0) 1.0 / (1.0 + (cr - 1.0)^gamma_c) else 1.0
          exp(sc / 0.5)
        }, numeric(1))
        pick <- sample.int(length(cxs), 1L, prob = wts)
        new_x <- cxs[pick]
        new_y <- cys[pick]
      }
      setState(agent, list(energy = new_energy, x = new_x, y = new_y))
    }
  )

  # Foraging -> Molting (Hiatt growth formulation)
  sim$addTransition(
    list(state = "foraging") -> list(state = "molting"),
    newExpWaitingTime(p$molting_rate),
    to_change_callback = still_in("foraging"),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      if (st$mature == "1") {
        # Terminally mature adults do not molt
        return()
      }
      mat_prob <- compute_maturity_probability(st$size, st$sex, 4.0)
      is_maturing <- runif(1) < mat_prob
      prev_cw <- as.numeric(st$size)
      prev_mass_g <- calculate_crab_mass_kg(prev_cw, st$sex, FALSE) * 1000.0
      target_cw <- hiatt_post_growth_cw(prev_cw, st$sex, is_mature = is_maturing, temp = 4.0)
      target_mass_g <- calculate_crab_mass_kg(target_cw, st$sex, is_maturing) * 1000.0
      delta_mass_g <- max(0.0, target_mass_g - prev_mass_g)
      pre_reserve <- as.numeric(st$energy)
      if (isTRUE(p$use_energetics) && delta_mass_g > 0.0 && pre_reserve < delta_mass_g) {
        avail_mass_g <- prev_mass_g + max(0.0, pre_reserve)
        beta_approx <- if (st$sex == "F") 2.95 else 3.12
        scaled_cw <- target_cw * ((avail_mass_g / max(1e-9, target_mass_g))^(1.0 / beta_approx))
        min_exp <- if (is_maturing) 0.06 else 0.12
        target_cw <- max(prev_cw * (1.0 + min_exp), min(target_cw, scaled_cw))
        target_mass_g <- calculate_crab_mass_kg(target_cw, st$sex, is_maturing) * 1000.0
        delta_mass_g <- max(0.0, target_mass_g - prev_mass_g)
      }
      new_size <- target_cw
      new_energy <- if (isTRUE(p$use_energetics)) {
        overhead_cap <- pre_reserve * as.numeric(p$molt_energy_retention %||% 0.1)
        max(0.0, min(overhead_cap, pre_reserve - delta_mass_g))
      } else pre_reserve
      new_mature <- if (is_maturing) "1" else "0"
      new_legal <- if (st$sex == "M" && new_size >= p$legal_size_male && new_size >= 68.0) "1" else "0"
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

  # Molting -> Foraging (recovery ~2 days -> rate 0.5/day)
  sim$addTransition(
    list(state = "molting") -> list(state = "foraging"),
    newExpWaitingTime(0.5),
    to_change_callback = still_in("molting"),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      gain_base <- if (st$sex == "M") p$energy_gain_foraging * 0.5 else p$energy_gain_foraging * 0.3
      effective_gain <- gain_base * as.numeric(p$food_availability)
      new_energy <- min(p$energy_max, st$energy + effective_gain)
      setState(agent, list(energy = new_energy))
    }
  )

  # Foraging -> Dead and Resting -> Dead transitions are intentionally ABSENT:
  # stacked never-cancelled events would multiply the death hazard (see the
  # mortality roll inside the Resting -> Foraging callback above).

  # Harvest (commercial fishing) restricted to mature legal males outside active MPA
  sim$addTransition(
    list(state = "foraging", legal = "1") -> list(state = "harvested"),
    newExpWaitingTime(f_daily_rate),
    to_change_callback = function(time, agent) {
      if (!identical(getState(agent)$state, "foraging")) return(FALSE)
      if (f_daily_rate <= 0.0) return(FALSE)
      st <- getState(agent)
      doy <- as.integer(time) %% 365
      in_season <- if (p$fishing_season_start <= p$fishing_season_end) {
        doy >= p$fishing_season_start && doy <= p$fishing_season_end
      } else {
        doy >= p$fishing_season_start || doy <= p$fishing_season_end
      }
      if (!in_season) return(FALSE)
      if (!is.null(p$mpa) && isTRUE(p$mpa$active)) {
        if (st$x >= p$mpa$xmin && st$x <= p$mpa$xmax &&
            st$y >= p$mpa$ymin && st$y <= p$mpa$ymax) {
          return(FALSE)
        }
      }
      sel <- compute_selectivity(st$size, s50 = p$legal_size_male, sigma = 5.0)
      ok <- st$sex == "M" && st$size >= p$legal_size_male && st$size >= 68.0 && runif(1) < sel
      if (ok) decrement_local_crowding_r(st$x, st$y, identical(st$measured, "1"))
      ok
    }
  )

  # Female reproduction (egg production, base-10 Moncton)
  sim$addTransition(
    list(state = "foraging", sex = "F", mature = "1") -> list(state = "reproducing"),
    newExpWaitingTime(p$repro_rate),
    to_change_callback = still_in("foraging"),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      fec <- calculate_fecundity_moncton(st$size, shell_age_years = 1.0)
      setState(agent, list(egg_count = round(fec)))
    }
  )

  # Reproducing -> Foraging (hatch ~2 days -> rate 0.5/day)
  sim$addTransition(
    list(state = "reproducing") -> list(state = "foraging"),
    newExpWaitingTime(0.5),
    to_change_callback = still_in("reproducing"),
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

  # Seed the live measured-crab count cache used by the crowding feedback.
  .crab_cache$n_meas <- sum(vapply(seq_len(n_crabs), function(i) {
    identical(getState(sim$agent(i))$measured, "1")
  }, logical(1)))

  sim
}

#' Extract the current state of every agent in the population
#'
#' @param sim An ABM `Simulation` instance.
#' @return A `data.frame` with columns x, y, state, sex, size, mature, legal,
#'   measured, energy, age — one row per agent.
#' @export
crab_agent_positions <- function(sim) {
  n <- sim$size
  out <- vector("list", n)
  for (i in seq_len(n)) {
    st <- getState(sim$agent(i))
    out[[i]] <- data.frame(
      x = as.numeric(st$x), y = as.numeric(st$y),
      state = as.character(st$state), sex = as.character(st$sex),
      size = as.numeric(st$size), mature = as.character(st$mature),
      legal = as.character(st$legal), measured = as.character(st$measured),
      energy = as.numeric(st$energy), age = as.numeric(st$age),
      mass_kg = calculate_crab_mass_kg(as.numeric(st$size), as.character(st$sex), identical(st$mature, "1")),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Summarize agent positions into dashboard state counters
#'
#' Rebuilds the historical ABM-logger column set (times, dead, foraging, ...,
#' measured) from a live agent dump. The ABM package's own counters accumulate
#' across repeated `run()` calls, so under the dashboard's chunked stepping
#' they are not usable as a time series; this recomputes clean snapshots.
#'
#' @param positions A `data.frame` from [crab_agent_positions()].
#' @param time The simulation time of the snapshot.
#' @return A one-row `data.frame` of state counters.
#' @export
crab_state_summary <- function(positions, time) {
  living <- !positions$state %in% c("dead", "harvested")
  has_mass <- "mass_kg" %in% names(positions)
  biomass_kg <- if (has_mass) sum(positions$mass_kg[living], na.rm = TRUE) else sum(living) * 0.35
  catch_biomass_kg <- if (has_mass) sum(positions$mass_kg[positions$state == "harvested"], na.rm = TRUE) else sum(positions$state == "harvested") * 0.45
  data.frame(
    times = as.numeric(time),
    dead = sum(positions$state == "dead"),
    foraging = sum(positions$state == "foraging"),
    harvested = sum(positions$state == "harvested"),
    legal_male = sum(positions$legal == "1" & living),
    mature_female = sum(positions$sex == "F" & positions$mature == "1" & living),
    mature_male = sum(positions$sex == "M" & positions$mature == "1" & living),
    measured = sum(positions$measured == "1" & living),
    molting = sum(positions$state == "molting"),
    reproducing = sum(positions$state == "reproducing"),
    resting = sum(positions$state == "resting"),
    biomass_kg = biomass_kg,
    catch_biomass_kg = catch_biomass_kg,
    mean_cw = if (any(living)) mean(positions$size[living], na.rm = TRUE) else 0,
    stringsAsFactors = FALSE
  )
}

#' Run crab simulation to the requested times and snapshot the population
#'
#' Advances an ABM simulation to `max(times)` and returns ONE summary row
#' aggregated from live agent states. The agent dump is also stored in
#' `.GlobalEnv$last_positions` for backward compatibility, but Shiny outputs
#' should prefer [crab_agent_positions()] called directly on the simulation.
#'
#' @param sim An ABM `Simulation` instance.
#' @param times Numeric vector of evaluation time points.
#' @return A one-row `data.frame` of population state counters.
#' @export
run_crab_simulation <- function(sim, times = seq(0, 50, by = 1)) {
  invisible(sim$run(times))
  positions <- crab_agent_positions(sim)
  assign("last_positions", positions, envir = .GlobalEnv)
  # Refresh the cached live measured count and B x B crowding grid.
  .crab_cache$n_meas <- sum(positions$measured == "1" &
                              !positions$state %in% c("dead", "harvested"))
  refresh_local_crowding_r(
    positions,
    grid_size = .crab_cache$grid_size %||% 50,
    unit_meters = .crab_cache$unit_meters %||% 1.0,
    k_density = .crab_cache$k_density %||% 6000.0,
    box_size = .crab_cache$box_size %||% 0L
  )
  crab_state_summary(positions, max(times))
}

# ==============================================================================
# TethysBio Parity Extensions (Spatial, Bioenergetics, Analysis, Export)
# ==============================================================================

#' Stage labels matching TethysBio 20-stage scheme (Egg, z1, z2, m, instars 1..16)
#' @export
cohort_stage_labels_r <- function() {
  c("Egg", "z1", "z2", "m", as.character(1:16))
}

#' Compute geographic bounding box centered at (lon, lat) for a domain span (m)
#' @export
aoi_domain_bounds_r <- function(lon, lat, span_meters = 1000.0) {
  half_m <- as.numeric(span_meters) / 2.0
  lat_deg_len <- 111139.0
  lon_deg_len <- max(1000.0, 111139.0 * cos(as.numeric(lat) * pi / 180.0))
  dlat <- half_m / lat_deg_len
  dlon <- half_m / lon_deg_len
  list(lon = c(lon - dlon, lon + dlon), lat = c(lat - dlat, lat + dlat))
}

#' Check whether (x, y) lies inside an active ExclusionBox list
#' @export
is_in_exclusion_zones_r <- function(x, y, boxes) {
  if (is.null(boxes) || length(boxes) == 0) return(FALSE)
  if (!is.null(boxes$xmin)) boxes <- list(boxes)
  for (b in boxes) {
    if (isTRUE(b$active) && x >= b$xmin && x <= b$xmax && y >= b$ymin && y <= b$ymax) {
      return(TRUE)
    }
  }
  FALSE
}

#' Build 2D spatial depth, temperature, food, and HSI fields
#' @export
build_environmental_fields_r <- function(nx = 30L, ny = 30L, temp_offset = 0.0,
                                         food_max = 10.0, hsi_offset = 0.0, hsi_mult = 1.0) {
  xs <- seq(0, 1, length.out = nx)
  ys <- seq(0, 1, length.out = ny)
  grid <- expand.grid(ix = seq_len(nx), iy = seq_len(ny))
  xv <- xs[grid$ix]
  yv <- ys[grid$iy]
  depth <- pmin(380, pmax(25, 140 + 95 * sin(2.5 * xv) * cos(2.0 * yv) + 45 * (1 - yv)))
  temp <- 2.8 - 0.004 * (depth - 140) + temp_offset
  hsi <- vapply(seq_along(depth), function(i) {
    aoi_habitat_suitability(temp[i], depth[i], mult = hsi_mult, offset = hsi_offset)
  }, numeric(1))
  list(
    nx = nx, ny = ny,
    depth = matrix(depth, nrow = nx, ncol = ny),
    temperature = matrix(temp, nrow = nx, ncol = ny),
    food = matrix(food_max * 0.5, nrow = nx, ncol = ny),
    hsi = matrix(hsi, nrow = nx, ncol = ny)
  )
}

#' Format simulation day count into Year-Month-Day starting at Year 0 (0-01-01)
#' @export
format_sim_date_r <- function(day = 0) {
  d <- max(0L, as.integer(floor(as.numeric(day %||% 0))))
  yr <- d %/% 365L
  rem <- d %% 365L
  mdays <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
  mo <- 1L
  for (m in seq_along(mdays)) {
    if (rem < mdays[m]) {
      mo <- m
      break
    }
    rem <- rem - mdays[m]
  }
  sprintf("%d-%02d-%02d", yr, mo, rem + 1L)
}

#' Build and query C^0 bilinear B x B local crowding ratio grid in .crab_cache
#' @export
refresh_local_crowding_r <- function(positions, grid_size = 1000, unit_meters = 1.0,
                                     k_density = 6000.0, box_size = 0L) {
  gsz <- max(1, as.numeric(grid_size))
  raw_b <- as.integer(box_size %||% 0L)
  b_cells <- max(1L, min(as.integer(gsz), if (raw_b > 0L) raw_b else as.integer(round(gsz / 5.0))))
  n_tiles <- max(1L, as.integer(ceiling(gsz / b_cells)))
  counts <- matrix(0.0, nrow = n_tiles, ncol = n_tiles)
  if (!is.null(positions) && nrow(positions) > 0) {
    living_meas <- !positions$state %in% c("dead", "harvested") & positions$measured == "1"
    if (any(living_meas)) {
      xs <- pmin(1.0, pmax(0.0, as.numeric(positions$x[living_meas])))
      ys <- pmin(1.0, pmax(0.0, as.numeric(positions$y[living_meas])))
      tx <- pmin(n_tiles, pmax(1L, as.integer(floor(xs * n_tiles)) + 1L))
      ty <- pmin(n_tiles, pmax(1L, as.integer(floor(ys * n_tiles)) + 1L))
      for (k in seq_along(tx)) {
        counts[tx[k], ty[k]] <- counts[tx[k], ty[k]] + 1.0
      }
    }
  }
  tile_area_km2 <- max(1e-9, ((b_cells * max(0.01, as.numeric(unit_meters)))^2) / 1e6)
  k_den <- max(1e-6, as.numeric(k_density))
  ratios <- (counts / tile_area_km2) / k_den
  .crab_cache$crowding_ratios <- ratios
  .crab_cache$crowding_ntiles <- n_tiles
  .crab_cache$crowding_delta <- (1.0 / tile_area_km2) / k_den
  invisible(ratios)
}

#' Query C^0 bilinear local crowding ratio at normalized coordinates (x, y) in [0, 1]
#' @export
local_crowding_ratio_r <- function(x, y) {
  ratios <- .crab_cache$crowding_ratios
  if (is.null(ratios)) return(0.0)
  nt <- as.integer(.crab_cache$crowding_ntiles %||% nrow(ratios))
  fx <- min(as.numeric(nt), max(1.0, as.numeric(x) * nt + 0.5))
  fy <- min(as.numeric(nt), max(1.0, as.numeric(y) * nt + 0.5))
  x0 <- as.integer(floor(fx))
  y0 <- as.integer(floor(fy))
  x1 <- min(nt, x0 + 1L)
  y1 <- min(nt, y0 + 1L)
  wx <- fx - x0
  wy <- fy - y0
  (1.0 - wx) * ((1.0 - wy) * ratios[x0, y0] + wy * ratios[x0, y1]) +
    wx * ((1.0 - wy) * ratios[x1, y0] + wy * ratios[x1, y1])
}

decrement_local_crowding_r <- function(x, y, is_measured = TRUE) {
  if (!isTRUE(is_measured)) return(invisible(NULL))
  ratios <- .crab_cache$crowding_ratios
  if (is.null(ratios)) return(invisible(NULL))
  nt <- as.integer(.crab_cache$crowding_ntiles %||% nrow(ratios))
  delta <- as.numeric(.crab_cache$crowding_delta %||% 0.0)
  tx <- min(nt, max(1L, as.integer(floor(pmin(1.0, pmax(0.0, as.numeric(x))) * nt)) + 1L))
  ty <- min(nt, max(1L, as.integer(floor(pmin(1.0, pmax(0.0, as.numeric(y))) * nt)) + 1L))
  ratios[tx, ty] <- max(0.0, ratios[tx, ty] - delta)
  .crab_cache$crowding_ratios <- ratios
  invisible(NULL)
}

#' Kleiber + Holling Type II daily net energy change matching TethysBio
#' @export
calculate_energy_change_r <- function(mass_g, temp, food_avail = 5.0,
                                      half_sat = 2.0, assim_eff = 0.75,
                                      sda = 0.15, metab_base = 0.01,
                                      intake_base = 0.05, t_opt = 3.0,
                                      t_sigma = 2.5) {
  m_use <- max(0.01, as.numeric(mass_g))
  m_scale <- m_use ^ 0.75
  f_pos <- max(0.0, as.numeric(food_avail))
  therm_mod <- exp(-((as.numeric(temp) - t_opt)^2) / (2.0 * (t_sigma^2)))
  intake_demand <- intake_base * (f_pos / (half_sat + f_pos)) * m_scale * therm_mod
  actual_intake <- min(intake_demand, f_pos)
  assim <- actual_intake * assim_eff * (1.0 - sda)
  basal <- metab_base * m_scale * exp(0.08 * as.numeric(temp))
  assim - basal
}

#' Competing-risks mortality components matching TethysBio::mortality_components
#' @export
mortality_components_r <- function(size_mm, temp = 4.0, density = 0.0, k_cap = 6000.0,
                                   is_soft_shell = FALSE, shell_age_days = 0.0,
                                   energy = 10.0, is_mature = FALSE, sex = "M",
                                   dt = 1.0, use_thermal = TRUE,
                                   use_senescence = TRUE, use_cannibalism = FALSE,
                                   local_crowding_ratio = 0.0) {
  dpy <- 365.0
  global_fb <- 1.0 + (max(0.0, density) / max(1.0, k_cap)) ^ 2.0
  local_fb <- if (local_crowding_ratio > 1.0) 1.0 + (local_crowding_ratio ^ 2.0) else 1.0
  feedback <- max(global_fb, local_fb)
  h_base <- (0.08 / dpy) * feedback * dt
  h_pred <- (0.05 / dpy) * exp(-0.03 * max(0.0, size_mm - 5.0)) * dt
  h_soft <- if (is_soft_shell && shell_age_days <= 30.0) (0.12 / dpy) * exp(-shell_age_days / 15.0) * dt else 0.0
  mass_g <- max(0.01, calculate_crab_mass_kg(max(1.0, size_mm), sex, is_mature) * 1000.0)
  crit_pool <- max(0.5, mass_g * 0.05)
  h_starve <- if (energy < 0.0) (0.20 / dpy) * (1.0 - exp(energy / crit_pool)) * dt else 0.0
  h_thermal <- if (use_thermal && temp > 8.0) (0.05 / dpy) * ((temp - 8.0) ^ 1.5) * dt else 0.0
  h_senesc <- if (use_senescence && is_mature && (shell_age_days / dpy) > 3.0) {
    (0.05 / dpy) * exp(0.8 * ((shell_age_days / dpy) - 3.0)) * dt
  } else 0.0
  h_cannib <- if (use_cannibalism && size_mm < 35.0) (0.05 / dpy) * min(2.0, max(1.0, density / 1000.0)) * dt else 0.0
  list(
    base = h_base, predation = h_pred, soft_shell = h_soft,
    starvation = h_starve, thermal = h_thermal,
    senescence = h_senesc, cannibalism = h_cannib
  )
}

#' Advance pelagic/early-juvenile AggregateGroup cohorts
#' @export
step_aggregate_groups_r <- function(groups, temp = 4.0, dt = 1.0, juv_max_stage = 12L) {
  if (length(groups) == 0) return(list(groups = list(), settled = 0))
  out <- list()
  settled <- 0
  for (g in groups) {
    g$age_days <- g$age_days + dt
    g$dd <- g$dd + max(0.0, temp + 1.5) * dt
    surv <- exp(-0.015 * dt)
    g$abundance <- g$abundance * surv
    if (g$abundance < 0.5) next
    p_trans <- 1.0 - exp(-0.03 * dt)
    n_tr <- round(g$abundance * p_trans)
    if (n_tr > 0) {
      g$abundance <- max(0, g$abundance - n_tr)
      nxt <- g$stage_idx + 1L
      if (nxt >= juv_max_stage) {
        settled <- settled + n_tr
      } else {
        out[[length(out) + 1L]] <- list(stage_idx = nxt, abundance = n_tr, age_days = 0, dd = 0)
      }
    }
    if (g$abundance >= 0.5) out[[length(out) + 1L]] <- g
  }
  list(groups = out, settled = settled)
}

#' Fit Ricker and Beverton-Holt Stock-Recruitment models
#' @export
fit_stock_recruitment_r <- function(spawners, recruits) {
  s <- as.numeric(spawners)
  r <- as.numeric(recruits)
  ok <- is.finite(s) & is.finite(r) & s > 0 & r > 0
  s <- s[ok]; r <- r[ok]
  if (length(s) < 3) {
    return(list(isFitted = FALSE, best_model = "beverton_holt",
                beverton_holt = list(alpha = 1, beta = 0.01, r2 = 0, aic = Inf),
                ricker = list(alpha = 1, beta = 0.01, r2 = 0, aic = Inf)))
  }
  # Ricker: log(R/S) = log(a) - b*S
  fit_rk <- lm(log(r / s) ~ s)
  a_rk <- exp(coef(fit_rk)[1]); b_rk <- -coef(fit_rk)[2]
  pred_rk <- a_rk * s * exp(-b_rk * s)
  sst <- sum((r - mean(r)) ^ 2)
  r2_rk <- if (sst > 0) 1 - sum((r - pred_rk) ^ 2) / sst else 0
  aic_rk <- length(r) * log(max(1e-12, mean((r - pred_rk) ^ 2))) + 4

  # Beverton-Holt: S/R = 1/a + (b/a)*S
  fit_bh <- lm((s / r) ~ s)
  a_bh <- 1 / max(1e-6, coef(fit_bh)[1])
  b_bh <- max(0, coef(fit_bh)[2] * a_bh)
  pred_bh <- (a_bh * s) / (1 + b_bh * s)
  r2_bh <- if (sst > 0) 1 - sum((r - pred_bh) ^ 2) / sst else 0
  aic_bh <- length(r) * log(max(1e-12, mean((r - pred_bh) ^ 2))) + 4

  list(
    isFitted = TRUE,
    best_model = if (aic_bh <= aic_rk) "beverton_holt" else "ricker",
    beverton_holt = list(alpha = as.numeric(a_bh), beta = as.numeric(b_bh), r2 = r2_bh, aic = aic_bh, predicted = pred_bh),
    ricker = list(alpha = as.numeric(a_rk), beta = as.numeric(b_rk), r2 = r2_rk, aic = aic_rk, predicted = pred_rk)
  )
}

#' Fit Generalized Leslie stage-projection matrix and dominant eigenvalue
#' @export
fit_generalized_leslie_model_r <- function(stage_matrix) {
  mat <- as.matrix(stage_matrix)
  if (ncol(mat) < 2 || nrow(mat) < 2) {
    return(list(isFitted = FALSE, lambda1 = 1.0, r2 = 0.0, predicted = colSums(mat)))
  }
  k <- nrow(mat)
  X <- mat[, -ncol(mat), drop = FALSE]
  Y <- mat[, -1, drop = FALSE]
  # Ridge-regularized projection L: Y ~ L %*% X
  Xt <- t(X)
  reg <- diag(1e-3, k)
  L <- Y %*% Xt %*% solve(X %*% Xt + reg)
  pred_mat <- cbind(mat[, 1], L %*% X)
  obs_tot <- colSums(mat)
  pred_tot <- colSums(pred_mat)
  sst <- sum((obs_tot - mean(obs_tot)) ^ 2)
  r2 <- if (sst > 0) 1 - sum((obs_tot - pred_tot) ^ 2) / sst else 0
  ev <- eigen(L, only.values = TRUE)$values
  lam <- max(Re(ev))
  list(isFitted = TRUE, matrix = L, lambda1 = as.numeric(lam), r2 = as.numeric(r2), predicted = pred_tot)
}

#' Compare Schaefer, Delay-Difference, and optional Leslie macro models
#' @export
compare_macro_models_r <- function(observed, catches = NULL, stage_matrix = NULL) {
  obs <- as.numeric(observed)
  cat <- if (is.null(catches)) rep(0, length(obs)) else as.numeric(catches)
  sch <- fit_schaefer_r(obs, cat)
  dd <- fit_delay_difference_r(obs, cat)
  les <- if (!is.null(stage_matrix)) fit_generalized_leslie_model_r(stage_matrix) else NULL
  best <- if ((sch$r2 %||% 0) >= (dd$r2 %||% 0)) "schaefer" else "delay_difference"
  list(schaefer = sch, delay_difference = dd, leslie = les, best_model = best)
}

#' Compute Mohn's rho retrospective bias statistic
#' @export
mohn_rho_r <- function(peel_vals, full_vals) {
  if (length(peel_vals) == 0 || length(peel_vals) != length(full_vals)) return(0.0)
  ok <- full_vals != 0
  if (!any(ok)) return(0.0)
  mean((peel_vals[ok] - full_vals[ok]) / full_vals[ok])
}

#' Perform retrospective peels of Schaefer fit and return Mohn's rho
#' @export
retrospective_analysis_r <- function(observed, catches = NULL, n_peels = 4L) {
  obs <- as.numeric(observed)
  cat <- if (is.null(catches)) rep(0, length(obs)) else as.numeric(catches)
  full <- fit_schaefer_r(obs, cat)
  peels <- numeric(0); refs <- numeric(0)
  for (p in seq_len(n_peels)) {
    if (length(obs) - p < 4) break
    sub_f <- fit_schaefer_r(obs[seq_len(length(obs) - p)], cat[seq_len(length(cat) - p)])
    peels <- c(peels, tail(sub_f$predicted, 1))
    refs <- c(refs, full$predicted[length(obs) - p])
  }
  list(mohn_rho = mohn_rho_r(peels, refs), peel_endpoints = peels)
}

#' Compute binned empirical maturity ogive
#' @export
compute_empirical_ogive_r <- function(sizes, is_mature, bin_width = 5.0) {
  s <- as.numeric(sizes)
  m <- as.numeric(is_mature)
  if (length(s) == 0) return(data.frame(size = numeric(0), prop = numeric(0), n = integer(0)))
  lo <- floor(min(s) / bin_width) * bin_width
  hi <- ceiling(max(s) / bin_width) * bin_width + bin_width
  brks <- seq(lo, hi, by = bin_width)
  bins <- cut(s, breaks = brks, include.lowest = TRUE, right = FALSE)
  df <- aggregate(m, by = list(bin = bins), FUN = mean)
  cnt <- aggregate(m, by = list(bin = bins), FUN = length)
  mids <- 0.5 * (brks[-length(brks)] + brks[-1])
  idx <- as.integer(df$bin)
  data.frame(size = mids[idx], prop = df$x, n = cnt$x)
}

#' Compute Gaussian KDE density over a carapace width grid
#' @export
kde_density_r <- function(x, from = 5, to = 160, n = 128, bw = 3.5) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2) {
    return(data.frame(x = seq(from, to, length.out = n), y = rep(0, n)))
  }
  d <- density(vals, bw = bw, from = from, to = to, n = n)
  data.frame(x = d$x, y = d$y)
}

#' Evaluate multi-attribute population calibration criteria
#' @export
evaluate_population_criteria_r <- function(history_df, positions_df = NULL, k_cap = 6000.0) {
  n_tot <- if ("foraging" %in% names(history_df)) {
    history_df$foraging + history_df$resting + history_df$molting + history_df$reproducing
  } else numeric(0)
  pers <- if (length(n_tot) > 0) min(n_tot) >= 5 else FALSE
  no_exp <- if (length(n_tot) > 0) max(n_tot) <= 5.0 * k_cap else TRUE
  drift <- if ("mean_cw" %in% names(history_df) && nrow(history_df) > 2) {
    tail_cw <- tail(history_df$mean_cw, min(10, nrow(history_df)))
    (tail(tail_cw, 1) - head(tail_cw, 1)) / max(1, head(tail_cw, 1))
  } else 0.0
  mat_rmse <- 0.0
  if (!is.null(positions_df) && nrow(positions_df) > 0) {
    m_sub <- positions_df[positions_df$sex == "M", , drop = FALSE]
    if (nrow(m_sub) > 5) {
      emp <- compute_empirical_ogive_r(m_sub$size, m_sub$mature == "1")
      if (nrow(emp) > 0) {
        pred <- vapply(emp$size, function(sz) compute_maturity_probability(sz, "M"), numeric(1))
        mat_rmse <- sqrt(mean((emp$prop - pred) ^ 2))
      }
    }
  }
  list(
    persistence = pers,
    no_explosion = no_exp,
    size_drift = drift,
    maturity_rmse = mat_rmse
  )
}

#' Export scenario parameters to a TOML string and optionally write to disk
#' @export
export_scenario_toml_r <- function(params, path = NULL) {
  lines <- c(
    "[simulation]",
    sprintf("n_crabs = %d", as.integer(params$n_crabs %||% 200)),
    "",
    "[environment]",
    sprintf("temp_offset = %.4f", as.numeric(params$temp_offset %||% 0.0)),
    sprintf("temp_multiplier = %.4f", as.numeric(params$temp_multiplier %||% 1.0)),
    sprintf("benthic_carrying_capacity = %.2f", as.numeric(params$benthic_carrying_capacity %||% 6000.0)),
    "",
    "[fishery]",
    sprintf("fishing_probability = %.4f", as.numeric(params$fishing_probability %||% 0.0)),
    sprintf("legal_size_male = %.2f", as.numeric(params$legal_size_male %||% 95.0)),
    sprintf("fishing_season_start = %.1f", as.numeric(params$fishing_season_start %||% 90.0)),
    sprintf("fishing_season_end = %.1f", as.numeric(params$fishing_season_end %||% 150.0))
  )
  txt <- paste(lines, collapse = "\n")
  if (!is.null(path)) writeLines(txt, con = path)
  txt
}

#' Save simulation summary and agent snapshot to disk (CSV / RDS)
#' @export
save_simulation_results_r <- function(history_df, positions_df = NULL, output_dir = "results") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  h_path <- file.path(output_dir, "annual_production.csv")
  write.csv(history_df, h_path, row.names = FALSE)
  paths <- list(history = h_path)
  if (!is.null(positions_df)) {
    p_path <- file.path(output_dir, "agent_snapshot.csv")
    write.csv(positions_df, p_path, row.names = FALSE)
    paths$positions <- p_path
  }
  paths
}
