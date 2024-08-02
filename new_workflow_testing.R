## header
library(epiwave)
library(epiwave.params)
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

# notif_dat <- notif_dat[notif_dat$jurisdiction == 'NSW',]
class(notif_dat) <- c('epiwave_fixed_timeseries',
                      'epiwave_timeseries',
                      class(notif_dat))

## hospitalisation counts
hosp_file <- paste0(not_synced_folder, '/COVID_live_cases_in_hospital.rds')

hosp_occupancy <- hosp_file |>
  readRDS() |>
  dplyr::filter(date %in% study_seq)
hosp_dat <- hosp_occupancy |>
  dplyr::mutate(value = cases_in_hospital/3)
class(hosp_dat) <- c('epiwave_fixed_timeseries',
                     'epiwave_timeseries',
                     class(hosp_dat))

# jurisdictions
jurisdictions <- unique(notif_dat$jurisdiction)
n_jurisdictions <- length(jurisdictions)

## days of infection timeseries
infection_days <- seq(from = as.Date('2021-05-03'),
                      to = as.Date('2022-01-29'), 'days')

n_days_infection <- length(infection_days)

## proportions
car_constant <- 0.75
car <- create_epiwave_fixed_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = car_constant) # 1 for sero prop

# create IHR
ihr <- greta::uniform(0, 1)
ihr_timeseries <- create_epiwave_fixed_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = list(car_constant * ihr)
)

### new

incubation <- readRDS('tests/test_distributions/incubation.rds')
gi <- readRDS('tests/test_distributions/gi.rds')
onset_to_notification <- readRDS('tests/test_distributions/onset_to_notification.rds')
# notification_to_hospitalisation <- lowerGPReff::data_to_distribution(delay_hospitalisation_timeseries)
hosp_dist <- distributional::dist_weibull(shape = 2.51, scale = 10.17)
hosp_delay_ecdf <- parametric_dist_to_distribution(hosp_dist)
hosp_full_delay_dist <- create_epiwave_massfun_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = hosp_delay_ecdf)


wavylistything <- fit_waves(
  observations = define_observation_model(
    cases = prepare_observation_data(
      timeseries_data = notif_dat,
      delay = add_distributions(incubation,
                                onset_to_notification),
      proportion_observed = car,
      type = "count",
      dow_model = create_dow_priors(n_jurisdictions)) # make an on/off?
    ,
    hospitalisations = prepare_observation_data(
      timeseries_data = hosp_dat,
      delay = hosp_full_delay_dist,
      proportion_observed = car,
      type = "count",
      ihr_correction = ihr)
  ),

  # proportion = ,
  target_infection_dates = infection_days,
  n_chains = 3,
  max_convergence_tries = 1,
  warmup = 100,
  n_samples = 100,
  n_extra_samples = 100
)

fitted_reff <- compute_reff(wavylistything, gi)


# outptus
infections_out <- GPreff::generate_long_estimates(
  wavylistything$infection_model$infection_timeseries,
  wavylistything$fit,
  infection_days,
  jurisdictions)

infection_traj <- GPreff::build_trajectories(
  param = wavylistything$infection_model$infection_timeseries,
  infection_days,
  wavylistything$fit,
  nsim = 1000,
  jurisdictions)

GPreff::plot_reff_interval_curves(
  'reff2.png',
  fitted_reff,
  dates = infection_days,
  start_date = min(study_seq),
  end_date = max(study_seq),
  jurisdictions = jurisdictions)

infection_sims <- greta::calculate(
  wavylistything$infection_model$infection_timeseries,
  values = wavylistything$fit,
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

