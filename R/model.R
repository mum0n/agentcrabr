# Crab Agent-Based Model (TethysBio-faithful)
# agentcrabr: Agent-based model for crab population dynamics
# Biological rules aligned with TethysBio (Chionoecetes opilio)

library(ABM)

#' Create a crab simulation with sex-specific biological rules
create_crab_simulation <- function(n_crabs = 50, params = list()) {
  p <- list(
    # General
    energy_init = 10,
    energy_max = 20,
    size_init_min = 20,
    size_init_max = 40,
    legal_size_male = 95,
    maturity_size_gate = 40,
    # Sex-specific growth (males grow faster/larger)
    growth_rate_male = 0.15,
    growth_rate_female = 0.10,
    max_size_male = 160,
    max_size_female = 95,
    # Sex-specific maturity (L50)
    maturity_l50_male = 91.8,
    maturity_l50_female = 51.6,
    maturity_slope_male = 6.9957,
    maturity_slope_female = 12.9812,
    # Fecundity (females only)
    fecundity_intercept = 2.0,
    fecundity_slope = 2.5,
    # Fishing (males only)
    fishing_probability = 0.05,
    selectivity_s50 = 95,
    selectivity_sigma = 5,
    # Energy
    energy_gain_foraging = 0.3,
    energy_cost_resting = 0.05,
    energy_cost_molting = 3,
    molt_growth_male = 0.12,
    molt_growth_female = 0.08,
    # Rates
    foraging_rate = 0.5,
    resting_rate = 0.3,
    molting_rate = 0.05,
    mortality_rate = 0.02,
    repro_rate = 0.1
  )
  p <- modifyList(p, params)

  sim <- Simulation$new(n_crabs, function(i) {
    sex <- sample(c("M", "F"), 1)
    size_init <- if (sex == "M") runif(1, p$size_init_min, p$size_init_max * 1.2)
                 else runif(1, p$size_init_min, p$size_init_max)
    list(
      state = "foraging",
      energy = p$energy_init,
      size = size_init,
      mature = 0,
      sex = sex,
      x = runif(1),
      y = runif(1),
      age = 0,
      egg_count = 0,
      legal = 0
    )
  })

  # Sex-specific growth increment
  get_growth <- function(org) {
    if (org$sex == "M") p$growth_rate_male else p$growth_rate_female
  }

  # Maturity ogive (sex-specific)
  get_maturity_prob <- function(org) {
    if (org$size < p$maturity_size_gate) return(0.0)
    l50 <- if (org$sex == "M") p$maturity_l50_male else p$maturity_l50_female
    slope <- if (org$sex == "M") p$maturity_slope_male else p$maturity_slope_female
    prob <- 1.0 / (1.0 + exp(-slope * log(max(1.0, org$size) / l50)))
    return(prob)
  }

  # Selectivity (fishing gear)
  get_selectivity <- function(org) {
    if (org$sex != "M") return(0.0)  # Only males
    diff <- (org$size - p$selectivity_s50) / p$selectivity_sigma
    if (diff > 40) return(1.0)
    if (diff < -40) return(0.0)
    return(1.0 / (1.0 + exp(-diff)))
  }

  # Is legal male
  is_legal_male <- function(org) {
    org$sex == "M" && org$mature == 1 && org$size >= p$legal_size_male
  }

  # --- Transitions ---

  # Foraging -> Resting
  sim$addTransition(
    list(state = "foraging") -> list(state = "resting"),
    newExpWaitingTime(1 / p$resting_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      new_energy <- max(0, st$energy - p$energy_cost_resting)
      # Update maturity
      mat_prob <- get_maturity_prob(st)
      new_mature <- if (st$mature == 1) 1 else if (runif(1) < mat_prob) 1 else 0
      # Update legal status
      new_legal <- if (st$sex == "M" && new_mature == 1 && st$size >= p$legal_size_male) 1 else 0
      setState(agent, list(energy = new_energy, mature = new_mature, legal = new_legal))
    }
  )

  # Resting -> Foraging
  sim$addTransition(
    list(state = "resting") -> list(state = "foraging"),
    newExpWaitingTime(1 / p$foraging_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      new_energy <- min(p$energy_max, st$energy + p$energy_gain_foraging)
      setState(agent, list(energy = new_energy))
    }
  )

  # Foraging -> Molting (sex-specific growth)
  sim$addTransition(
    list(state = "foraging") -> list(state = "molting"),
    newExpWaitingTime(1 / p$molting_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      growth <- if (st$sex == "M") p$molt_growth_male else p$molt_growth_female
      max_sz <- if (st$sex == "M") p$max_size_male else p$max_size_female
      new_size <- min(max_sz, st$size * (1 + growth))
      new_energy <- max(0, st$energy - p$energy_cost_molting)
      # Update maturity after growth
      mat_prob <- get_maturity_prob(list(size = new_size, sex = st$sex, mature = st$mature))
      new_mature <- if (st$mature == 1) 1 else if (runif(1) < mat_prob) 1 else 0
      new_legal <- if (st$sex == "M" && new_mature == 1 && new_size >= p$legal_size_male) 1 else 0
      setState(agent, list(size = new_size, energy = new_energy, mature = new_mature, legal = new_legal))
    }
  )

  # Molting -> Foraging
  sim$addTransition(
    list(state = "molting") -> list(state = "foraging"),
    newExpWaitingTime(2),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      gain <- if (st$sex == "M") p$energy_gain_foraging * 0.5 else p$energy_gain_foraging * 0.3
      new_energy <- min(p$energy_max, st$energy + gain)
      setState(agent, list(energy = new_energy))
    }
  )

  # Foraging -> Dead (starvation or predation)
  sim$addTransition(
    list(state = "foraging") -> list(state = "dead"),
    newExpWaitingTime(1 / p$mortality_rate),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      st$energy <= 0
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

  # Harvest (fishing) on legal males
  sim$addTransition(
    list(state = "foraging", legal = TRUE) -> list(state = "harvested"),
    newExpWaitingTime(1 / p$fishing_probability),
    to_change_callback = function(time, agent) {
      st <- getState(agent)
      st$sex == "M" && st$mature == 1 && st$size >= p$legal_size_male && st$legal == 1
    }
  )

  # Female reproduction (egg production)
  sim$addTransition(
    list(state = "foraging", sex = "F", mature = 1) -> list(state = "reproducing"),
    newExpWaitingTime(1 / p$repro_rate),
    changed_callback = function(time, agent) {
      st <- getState(agent)
      # Fecundity: size-dependent
      fec <- p$fecundity_intercept * (st$size ^ p$fecundity_slope)
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
  sim$addLogger(newCounter("mature_male", list(sex = "M", mature = 1)))
  sim$addLogger(newCounter("mature_female", list(sex = "F", mature = 1)))
  sim$addLogger(newCounter("legal_male", list(legal = 1)))

  sim
}

#' Run crab simulation
run_crab_simulation <- function(sim, times = seq(0, 50, by = 1)) {
  res <- sim$run(times)
  as.data.frame(res)
}

#' Compute maturity ogive
compute_maturity_ogive <- function(sizes, l50, slope) {
  1.0 / (1.0 + exp(-slope * log(pmax(1.0, sizes) / l50)))
}

#' Compute selectivity curve
compute_selectivity <- function(sizes, s50, sigma) {
  1.0 / (1.0 + exp(-(sizes - s50) / sigma))
}

#' Compute habitat suitability from bottom temperature and depth
aoi_habitat_suitability <- function(temp, depth) {
  thermal <- exp(-0.5 * ((as.numeric(temp) - 2.0) / 2.2)^2)
  depth_fit <- exp(-0.5 * ((as.numeric(depth) - 150.0) / 75.0)^2)
  max(0.0, min(1.0, thermal * depth_fit))
}
