# Unit tests for agentcrabr core model functions and biological rules

source("C:/Users/Nathalie/projects/agentcrabr/R/model.R")

test_mass_allometry <- function() {
  m_imm <- calculate_crab_mass_kg(50.0, "M", FALSE)
  m_mat <- calculate_crab_mass_kg(95.0, "M", TRUE)
  f_imm <- calculate_crab_mass_kg(40.0, "F", FALSE)
  f_mat <- calculate_crab_mass_kg(70.0, "F", TRUE)

  stopifnot(m_imm > 0.04 && m_imm < 0.06)
  stopifnot(m_mat > 0.30 && m_mat < 0.35)
  stopifnot(f_imm > 0.015 && f_imm < 0.025)
  stopifnot(f_mat > 0.10 && f_mat < 0.13)
  cat("test_mass_allometry: PASS\n")
}

test_hiatt_growth <- function() {
  next_m <- hiatt_post_growth_cw(80.0, "M", FALSE, 4.0)
  next_f <- hiatt_post_growth_cw(60.0, "F", FALSE, 4.0)
  stopifnot(next_m > 105.0 && next_m < 112.0)
  stopifnot(next_f > 75.0 && next_f < 85.0)

  # Check maximum size caps (160mm male, 95mm female)
  cap_m <- hiatt_post_growth_cw(155.0, "M", FALSE, 4.0)
  cap_f <- hiatt_post_growth_cw(90.0, "F", FALSE, 4.0)
  stopifnot(cap_m <= 160.0)
  stopifnot(cap_f <= 95.0)
  cat("test_hiatt_growth: PASS\n")
}

test_maturity_probability <- function() {
  # 35mm size gate
  p_small <- compute_maturity_probability(34.0, "M", 4.0)
  stopifnot(p_small == 0.0)

  # Centroid at L50 = 91.8 mm
  p_l50_m <- compute_maturity_probability(91.8, "M", 4.0)
  stopifnot(abs(p_l50_m - 0.5) < 1e-4)

  # Female centroid at L50 = 51.6 mm
  p_l50_f <- compute_maturity_probability(51.6, "F", 4.0)
  stopifnot(abs(p_l50_f - 0.5) < 1e-4)
  cat("test_maturity_probability: PASS\n")
}

test_benthic_capacity_and_measured <- function() {
  k <- benthic_carrying_capacity(6000.0, 50, 1.0)
  stopifnot(k >= 1)

  org_small <- list(size = 20.0, stage = 5, mature = "0")
  org_meas <- list(size = 40.0, stage = 11, mature = "0")
  org_mat <- list(size = 30.0, stage = 8, mature = "1")

  stopifnot(!is_empirically_measured(org_small))
  stopifnot(is_empirically_measured(org_meas))
  stopifnot(is_empirically_measured(org_mat))
  cat("test_benthic_capacity_and_measured: PASS\n")
}

test_seasonal_temp <- function() {
  cycle <- seasonal_bottom_temperature_cycle(seq(0, 365, by = 10), 150.0, 45.0)
  stopifnot(length(cycle$temps) == length(cycle$days))
  stopifnot(cycle$mean_temp > 2.0 && cycle$mean_temp < 4.0)
  cat("test_seasonal_temp: PASS\n")
}

test_simulation_run <- function() {
  sim <- create_crab_simulation(15)
  df <- run_crab_simulation(sim, 0:3)
  # one clean snapshot row per call, evaluated at the final time
  stopifnot(nrow(df) == 1)
  stopifnot(df$times == 3)
  stopifnot("foraging" %in% names(df))
  alive <- df$foraging + df$resting + df$molting + df$reproducing
  stopifnot(alive >= 0 && alive <= 15)
  pos <- crab_agent_positions(sim)
  stopifnot(nrow(pos) == 15)
  cat("test_simulation_run: PASS\n")
}

test_mass_allometry()
test_hiatt_growth()
test_maturity_probability()
test_benthic_capacity_and_measured()
test_seasonal_temp()
test_simulation_run()
cat("\nALL AGENTCRABR TESTS PASSED SUCCESSFULLY!\n")

