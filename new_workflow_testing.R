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

### new

incubation <- readRDS('tests/test_distributions/incubation.rds')
gi <- readRDS('tests/test_distributions/gi.rds')
onset_to_notification <- readRDS('tests/test_distributions/onset_to_notification.rds')
# notification_to_hospitalisation <- lowerGPReff::data_to_distribution(delay_hospitalisation_timeseries)


wavylistything <- fit_waves(
  observations = define_observation_model(
    cases = prepare_observation_data(
      timeseries_data = notif_dat,
      delay = add_distributions(incubation,
                                onset_to_notification),
      proportion_observed = car,
      type = "count", # check with August
      dow = create_dow_priors(n_jurisdictions)) # make an on/off?
    # ,
    # hospitalisations =
  ),

  target_infection_dates = infection_days
)

fitted_reff <- compute_reff(wavylistything, gi)





### old

# long is classic output, then map to infection traj.
infections_out <- GPreff::generate_long_estimates(
  infection_model_objects$infection_timeseries,
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

