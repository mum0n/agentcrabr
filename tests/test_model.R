# Unit tests for agentcrabr core model functions and biological rules

if (file.exists("R/model.R")) {
  source("R/model.R")
} else {
  source("../R/model.R")
}

test_mass_allometry <- function() {
  m_imm <- calculate_crab_mass_kg(50.0, "M", FALSE)
  m_mat <- calculate_crab_mass_kg(95.0, "M", TRUE)
  f_imm <- calculate_crab_mass_kg(40.0, "F", FALSE)
  f_mat <- calculate_crab_mass_kg(70.0, "F", TRUE)

  stopifnot(abs(m_imm - (0.000183 * 50.0^3.16) / 1000.0) < 1e-8)
  stopifnot(abs(m_mat - (0.0002665 * 95.0^3.098) / 1000.0) < 1e-8)
  stopifnot(abs(f_imm - (0.000344 * 40.0^2.98) / 1000.0) < 1e-8)
  stopifnot(abs(f_mat - (0.00045 * 70.0^2.92) / 1000.0) < 1e-8)
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
  stopifnot(cycle$mean_temp > 2.0 && cycle$mean_temp <= 5.0)
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

test_fecundity_hsi_and_macro <- function() {
  fec <- calculate_fecundity_moncton(70.0, shell_age_years = 1.0)
  stopifnot(abs(fec - 10^(-0.7229 + 2.94 * log10(70.0))) < 1e-6)
  stopifnot(aoi_habitat_suitability(2.0, 200.0) > 0.9)
  stopifnot(aoi_habitat_suitability(2.0, 450.0) > 0.9)
  stopifnot(aoi_habitat_suitability(11.0, 200.0) == 0.0)
  stopifnot(aoi_habitat_suitability(2.0, 15.0) == 0.0)

  obs <- c(100, 120, 138, 150, 158, 162)
  catches <- c(5, 6, 7, 8, 8, 8)
  sch <- fit_schaefer_r(obs, catches)
  dd <- fit_delay_difference_r(obs, catches)
  stopifnot(sch$R2 > 0.5 && sch$MSY > 0.0)
  stopifnot(dd$R2 > 0.5 && dd$M > 0.0)
  cat("test_fecundity_hsi_and_macro: PASS\n")
}

test_fecundity_hsi_and_macro()

test_extended_r_parity <- function() {
  stopifnot(length(cohort_stage_labels_r()) == 20)
  b <- aoi_domain_bounds_r(-60.0, 45.0, 2000.0)
  stopifnot(b$lon[1] < -60.0 && b$lon[2] > -60.0)
  mpa <- list(xmin = 0.2, xmax = 0.5, ymin = 0.2, ymax = 0.5, active = TRUE)
  stopifnot(is_in_exclusion_zones_r(0.3, 0.3, mpa))
  stopifnot(!is_in_exclusion_zones_r(0.8, 0.8, mpa))

  flds <- build_environmental_fields_r(10, 10)
  stopifnot(all(dim(flds$hsi) == c(10, 10)))
  e_delta <- calculate_energy_change_r(200.0, 2.0, food_avail = 5.0)
  stopifnot(is.finite(e_delta))
  comps <- mortality_components_r(50.0, temp = 10.0, is_soft_shell = TRUE, shell_age_days = 5.0)
  stopifnot(comps$thermal > 0 && comps$soft_shell > 0)

  ag <- step_aggregate_groups_r(list(list(stage_idx = 11L, abundance = 100, age_days = 0, dd = 0)), temp = 4.0, dt = 1.0)
  stopifnot(length(ag$groups) >= 1)

  sp <- c(20, 30, 40, 50, 60, 70)
  rec <- c(25, 38, 48, 55, 58, 56)
  sr <- fit_stock_recruitment_r(sp, rec)
  stopifnot(isTRUE(sr$isFitted))

  smat <- rbind(c(100, 95, 92, 90, 88), c(40, 45, 48, 50, 52))
  les <- fit_generalized_leslie_model_r(smat)
  stopifnot(isTRUE(les$isFitted) && les$lambda1 > 0.5)

  ret <- retrospective_analysis_r(c(100, 120, 138, 150, 158, 162, 165), c(5, 6, 7, 8, 8, 8, 8), n_peels = 2)
  stopifnot(is.finite(ret$mohn_rho))

  og <- compute_empirical_ogive_r(c(30, 50, 95, 110), c(0, 0, 1, 1))
  stopifnot(nrow(og) >= 2)
  kd <- kde_density_r(c(40, 60, 80, 100))
  stopifnot(max(kd$y) > 0)

  toml_txt <- export_scenario_toml_r(list(n_crabs = 200, fishing_probability = 0.0))
  stopifnot(grepl("n_crabs = 200", toml_txt))
  stopifnot(format_sim_date_r(0) == "0-01-01")
  stopifnot(format_sim_date_r(600) == "1-08-24")
  stopifnot(is.finite(local_crowding_ratio_r(0.5, 0.5)))
  cat("test_extended_r_parity: PASS\n")
}

test_extended_r_parity()
cat("\nALL AGENTCRABR TESTS PASSED SUCCESSFULLY!\n")

