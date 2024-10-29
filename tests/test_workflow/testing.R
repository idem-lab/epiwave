## header
library(epiwave)
library(epiwave.params)
library(dplyr)
library(greta)


## user modified
study_seq <- seq(from = as.Date('2021-08-01'),
                 to = as.Date('2021-09-15'), 'days')
infection_days <- seq(from = as.Date('2021-07-03'),
                      to = as.Date('2021-10-12'), 'days')
jurisdictions <- 'NSW'

# specific folder with not synced data
not_synced_folder <- '../../BDSS_governance/BDSS/data'


## data prep
# notification counts
notif_file <- paste0(not_synced_folder, '/COVID_live_cases.rds')
notif_dat <- notif_file |>
  readRDS() |>
  dplyr::rename(value = new_cases) |>
  dplyr::filter(date %in% study_seq)
if (!is.na(jurisdictions)) {
  notif_dat <- notif_dat[notif_dat$jurisdiction == jurisdictions,]
}
class(notif_dat) <- c('epiwave_fixed_timeseries',
                      'epiwave_timeseries',
                      class(notif_dat))

# hospitalisation counts
hosp_file <- paste0(not_synced_folder, '/COVID_live_cases_in_hospital.rds')
hosp_dat <- hosp_file |>
  readRDS() |>
  dplyr::filter(date %in% study_seq) |>
  dplyr::mutate(value = cases_in_hospital/3)
if (!is.na(jurisdictions)) {
  hosp_dat <- hosp_dat[hosp_dat$jurisdiction == jurisdictions,]
}
class(hosp_dat) <- c('epiwave_fixed_timeseries',
                     'epiwave_timeseries',
                     class(hosp_dat))

# jurisdictions
if (is.na(jurisdictions)) {
  jurisdictions <- unique(notif_dat$jurisdiction)
}
n_jurisdictions <- length(jurisdictions)

# days of infection timeseries
n_days_infection <- length(infection_days)

# proportions
car_constant <- 0.75
car <- create_epiwave_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = car_constant) # 1 for sero prop

# create IHR
chr <- greta::uniform(0, 1)
# ihr <- chr * car
# wrapper for ihr specific flow
ihr <- create_epiwave_greta_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  car = car,
  chr_prior = chr)


## delays
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


fit_object <- fit_waves(
  observations = define_observation_model(
    cases = prepare_observation_data(
      timeseries_data = notif_dat,
      delay = epiwave.params::add_distributions(incubation,
                                                onset_to_notification),
      proportion_infections = car,
      type = "count",
      dow_model = create_dow_priors(n_jurisdictions)) # make an on/off?
    ,
    hospitalisations = prepare_observation_data(
      timeseries_data = hosp_dat,
      delay = hosp_full_delay_dist,
      proportion_infections = ihr,
      type = "count")
    # ,
    # other
  ),

  infection_model = 'flat_prior', #'gp_growth_rate',
  target_infection_dates = infection_days#,
  # n_chains = 2,
  # max_convergence_tries = 2,
  # warmup = 100,
  # n_samples = 100,
  # n_extra_samples = 100
)

fitted_reff <- compute_reff(fit_object, gi)


# outputs
infections_out <- epiwave.pipelines::generate_long_estimates(
  fit_object$infection_model$infection_timeseries,
  fit_object$fit,
  infection_days,
  jurisdictions)

infection_traj <- epiwave.pipelines::build_trajectories(
  param = fit_object$infection_model$infection_timeseries,
  infection_days,
  fit_object$fit,
  nsim = 1000,
  jurisdictions)

epiwave.pipelines::plot_reff_interval_curves(
  'gp_both_reff_9Oct.pdf',
  fitted_reff,
  dates = infection_days,
  start_date = min(study_seq),
  end_date = max(study_seq),
  jurisdictions = jurisdictions)

infection_sims <- greta::calculate(
  fit_object$infection_model$infection_timeseries,
  values = fit_object$fit,
  nsim = 1000)

epiwave.pipelines::plot_timeseries_sims(
  'output/gp_both_infections_Kate-9Oct.png',
  infection_sims[[1]],
  type = "infection",
  dates = infection_days,
  states = jurisdictions,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = as.Date('2021-12-08'))

