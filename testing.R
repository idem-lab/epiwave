## header
library(lowerGPreff)
library(dplyr)
library(greta)

study_seq <- seq(from = as.Date('2021-06-01'),
                 to = as.Date('2021-12-31'), 'days')

## user specific folder with not synced data
not_synced_folder <- '../../BDSS_governance/BDSS/data'

## notification counts
notif_file <- paste0(not_synced_folder, '/COVID_live_cases.rds')

notif_dat <- notif_file |>
  readRDS() |>
  dplyr::rename(value = new_cases) |>
  dplyr::filter(date %in% study_seq)
class(notif_dat) <- c('lowerGPreff_fixed_timeseries',
                      'lowerGPreff_timeseries',
                      class(notif_dat))

## hospitalisation counts
hosp_file <- paste0(not_synced_folder, '/COVID_live_cases_in_hospital.rds')

hosp_occupancy <- hosp_file |>
  readRDS() |>
  dplyr::filter(date %in% study_seq)
hosp_dat <- hosp_occupancy |>
  dplyr::mutate(value = cases_in_hospital/3)
class(hosp_dat) <- c('lowerGPreff_fixed_timeseries',
                     'lowerGPreff_timeseries',
                     class(hosp_dat))


# jurisdictions
jurisdictions <- unique(notif_dat$jurisdiction)
n_jurisdictions <- length(jurisdictions)

## delay data
oldnotif_delay <- readRDS(paste0(not_synced_folder, '/ECDF_delay_constant_PCR.rds'))
min_delay <- 0
max_delay <- 41
notif_delay_ecdf <- create_lowerGPreff_massfun(
  min_delay, max_delay, oldnotif_delay, normalise = TRUE)

## days of infection timeseries
infection_days <- seq(from = as.Date('2021-05-03'),
                      to = as.Date('2022-01-29'), 'days')

n_days_infection <- length(infection_days)

## proportions
car_constant <- 0.75
car <- create_lowerGPreff_fixed_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = car_constant) # 1 for sero prop

# create IHR
ihr <- create_lowerGPreff_fixed_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = list(car_constant * greta::uniform(0, 1))
)

## incubation period
incub_dist <- distributional::dist_weibull(
  shape = 1.83, scale = 4.93)
incubation_period_distribution <- parametric_dist_to_distribution(incub_dist)

## formatted delay distribution
notif_delay_dist <-
  add(notif_delay_ecdf, incubation_period_distribution)
notif_full_delay_dist <- create_lowerGPreff_massfun_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = notif_delay_dist)


hosp_dist <- distributional::dist_weibull(shape = 2.51, scale = 10.17)
hosp_delay_ecdf <- parametric_dist_to_distribution(hosp_dist)
hosp_full_delay_dist <- create_lowerGPreff_massfun_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = hosp_delay_ecdf)


## generation interval distribution
gi_distribution_data <- readr::read_csv(
  file = paste0(not_synced_folder,'/nishiura_samples.csv'),
  col_types = readr::cols(param1 = readr::col_double(),
                          param2 = readr::col_double()))

gi_dist <- distributional::dist_lognormal(
  mu = mean(gi_distribution_data$param1),
  sigma = mean(gi_distribution_data$param2))
generation_interval_distribution <- parametric_dist_to_distribution(gi_dist)

## optional day-of-week effect model
dow_model <- lowerGPreff::create_dow_priors(
  n_jurisdictions)

## infection timeseries model
infection_model_objects <- lowerGPreff::create_infection_timeseries(
  n_days_infection,
  n_jurisdictions,
  effect_type = 'growth_rate')

## observation models
notif_observation_model_objects <- create_observation_model(
  infection_timeseries = infection_model_objects$infection_timeseries,
  delay_distribution = notif_full_delay_dist,
  proportion_observed = car,
  count_data = notif_dat,
  dow_model = dow_model,
  data_id = 'notif')

## PROBLEM: here with the IHR being lists.
# hosp_observation_model_objects <- create_observation_model(
#   infection_timeseries = infection_model_objects$infection_timeseries,
#   delay_distribution = hosp_full_delay_dist,
#   proportion_observed = ihr,
#   count_data = hosp_dat,
#   data_id = 'hosp')

## Reff
reff_model_objects <- lowerGPreff::estimate_reff(
  infection_timeseries = infection_model_objects$infection_timeseries,
  generation_interval_mass_fxns = generation_interval_distribution)

## greta model fit
m <- greta::model(infection_model_objects$infection_timeseries,
                  reff_model_objects$reff)

plot(m)

fit <- GPreff::fit_model(
  model = m,
  n_chains = 2,
  max_convergence_tries = 1,
  warmup = 100,
  n_samples = 100,
  n_extra_samples = 100)

# long is classic output, then map to infection traj.
infections_out <- GPreff::generate_long_estimates(
  infection_model_objects$infection_timeseries,
  fit,
  infection_days,
  jurisdictions)

reff_out <- GPreff::generate_long_estimates(
  reff_model_objects$reff,
  fit,
  infection_days,
  jurisdictions)

infection_traj <- GPreff::build_trajectories(
  param = infection_model_objects$infection_timeseries,
  infection_days,
  fit,
  nsim = 1000,
  jurisdictions)

GPreff::plot_reff_interval_curves(
  'reff2.png',
  reff_out,
  dates = infection_days,
  start_date = min(study_seq),
  end_date = max(study_seq),
  jurisdictions = jurisdictions)

infection_sims <- greta::calculate(infection_model_objects$infection_timeseries,
                                   values = fit,
                                   nsim = 1000)

GPreff::plot_timeseries_sims(
  'infection_timeseries2.png',
  infection_sims[[1]],
  type = "infection",
  dates = infection_days,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = FALSE)

