## header
# library(lowerGPreff)
# library(GPreff)
source('R/pmf_dev.R')
library(dplyr)
library(greta)

study_seq <- seq(from = as.Date('2021-06-01'),
                 to = as.Date('2021-12-31'), 'days')

## notification counts
notif_file <- '../../BDSS_governance/BDSS/data/COVID_live_cases.rds'

notif_dat <- notif_file |>
  readRDS() |>
  dplyr::rename(count = new_cases) |>
  dplyr::filter(date %in% study_seq)

## hospitalisation counts
hosp_file <- '../../BDSS_governance/BDSS/data/COVID_live_cases_in_hospital.rds'

hosp_occupancy <- hosp_file |>
  readRDS() |>
  dplyr::filter(date %in% study_seq)
hosp_dat <- hosp_occupancy |>
  dplyr::mutate(count = cases_in_hospital/3)

# jurisdictions
jurisdictions <- unique(notif_dat$jurisdiction)
n_jurisdictions <- length(jurisdictions)

## delay data
oldnotif_delay <- readRDS('../../BDSS_governance/BDSS/data/ECDF_delay_constant_PCR.rds')
## created massfun version using the following. to update with pmf from data
min_delay <- 0
max_delay <- 41
notif_delay_ecdf <- make_massfun(min_delay, max_delay, oldnotif_delay, normalise = TRUE)

## days of infection timeseries
infection_days <- seq(from = as.Date('2021-05-03'),
                      to = as.Date('2022-01-29'), 'days')

n_days_infection <- length(infection_days)

## proportions
car <- create_lowerGPreff_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  constant_val = 0.75) #will error because we removed col_name from function

# we may pass in timevarying or greta arrays in the "second class"
ihr <- lowerGPreff::create_ihr_prior(car)

## incubation period
incub_dist <- distributional::dist_weibull(
  shape = 1.83, scale = 4.93)
# incubation_period_distribution <- lowerGPreff::make_cdf(
#   option = "Delta")
incubation_period_distribution <- params_as_delay_dist(incub_dist)

## formatted delay distribution
# notif_delay_dist <- GPreff::expand_constant_value(
#   dates = infection_days,
#   jurisdictions = jurisdictions,
#   constant_val = list(notif_delay_ecdf),
#   col_name = 'delay_fxn')
# notif_full_delay_dist <- lowerGPreff::extend_delay_data(
#   notif_delay_dist,
#   incubation_period_distribution)
notif_delay_dist <-
  combine_massfuns(notif_delay_ecdf, incubation_period_distribution)
notif_full_delay_dist <- create_lowerGPreff_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  constant_val = list(notif_delay_dist))


# hosp_delay_ecdf <- make_cdf(option = 'None',
#                             weibull_shape = 2.51,
#                             weibull_scale = 10.17)
# hosp_delay_dist <- GPreff::expand_constant_value(
#   dates = infection_days,
#   jurisdictions = jurisdictions,
#   constant_val = list(hosp_delay_ecdf),
#   col_name = 'delay_fxn')
# hosp_full_delay_dist <- lowerGPreff::extend_delay_data(
#   hosp_delay_dist)
hosp_dist <- distributional::dist_weibull(shape = 2.51, scale = 10.17)
hosp_delay_ecdf <- params_as_delay_dist(hosp_dist)
hosp_full_delay_dist <- create_lowerGPreff_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  constant_val = list(hosp_delay_ecdf))


## generation interval distribution
gi_distribution_data <- readr::read_csv(
  file = '../../BDSS_governance/BDSS/data/nishiura_samples.csv',
  col_types = readr::cols(param1 = readr::col_double(),
                          param2 = readr::col_double()))

# generation_interval_distribution <- lowerGPreff::make_generation_interval_density(
#   gi_distribution_data)
gi_dist <- distributional::dist_lognormal(
  mu = mean(gi_distribution_data$param1),
  sigma = mean(gi_distribution_data$param2))
generation_interval_distribution <- params_as_delay_dist(gi_dist)

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

hosp_observation_model_objects <- lowerGPreff::create_observation_model(
  infection_timeseries = infection_model_objects$infection_timeseries,
  delay_distribution = hosp_full_delay_dist,
  proportion_observed = ihr,
  count_data = hosp_dat,
  data_id = 'hosp')

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
  n_chains = 4,
  max_convergence_tries = 1,
  warmup = 1000,
  n_samples = 1000,
  n_extra_samples = 1000)

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

infection_traj <- build_trajectories(
  param = infection_model_objects$infection_timeseries,
  infection_days,
  fit,
  nsim = 1000,
  jurisdictions)

GPreff::plot_reff_interval_curves(
  'outputs/reff2.png',
  reff_out,
  dates = infection_days,
  start_date = min(study_seq),
  end_date = max(study_seq),
  jurisdictions = jurisdictions)

infection_sims <- greta::calculate(infection_model_objects$infection_timeseries,
                                   values = fit,
                                   nsim = 1000)

GPreff::plot_timeseries_sims(
  'outputs/infection_timeseries2.png',
  infection_sims[[1]],
  type = "infection",
  dates = infection_days,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = FALSE)

