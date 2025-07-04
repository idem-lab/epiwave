## header
library(epiwave)
library(epiwave.params)
library(dplyr)
library(greta)
# devtools::load_all()
## user modified
infection_days <- seq(from = as.Date('2021-04-01'),
                 to = as.Date('2022-01-01'), 'days')
study_seq <- seq(from = as.Date('2021-06-01'),
                      to = as.Date('2021-12-01'), 'days')
jurisdictions <- 'VIC' #c('NSW', 'VIC')

# specific folder with not synced data
not_synced_folder <- '../data'


## data prep
# notification counts
notif_file <- paste0(not_synced_folder, '/COVID_live_cases.rds')
notif_dat <- notif_file |>
  readRDS() |>
  dplyr::rename(value = new_cases) |>
  dplyr::filter(date %in% study_seq)
if (exists('jurisdictions')) {
  notif_dat <- notif_dat[notif_dat$jurisdiction %in% jurisdictions,]
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
if (exists('jurisdictions')) {
  hosp_dat <- hosp_dat[hosp_dat$jurisdiction %in% jurisdictions,]
}
class(hosp_dat) <- c('epiwave_fixed_timeseries',
                     'epiwave_timeseries',
                     class(hosp_dat))

# jurisdictions
if (!exists('jurisdictions')) {
  jurisdictions <- unique(notif_dat$jurisdiction)
}
n_jurisdictions <- length(jurisdictions)

# days of infection timeseries
n_days_infection <- length(infection_days)

# proportions
car <- 0.75

# create IHR
chr <- greta::uniform(0, 1)
ihr <- create_epiwave_greta_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  car = car,
  chr_prior = chr)

## delays
incubation <- readRDS('tests/test_distributions/incubation.rds')
gi <- readRDS('tests/test_distributions/gi.rds')
onset_to_notification <- readRDS('tests/test_distributions/onset_to_notification.rds')
hosp_dist <- distributional::dist_weibull(shape = 2.51, scale = 10.17)
hosp_delay_ecdf <- parametric_dist_to_distribution(hosp_dist)

observation_models <- define_observation_model(

  target_infection_dates = infection_days,
  target_jurisdictions = jurisdictions,

  cases = define_observation_data(
    timeseries_data = notif_dat,
    delay_from_infection =
      epiwave.params::add_distributions(
      incubation,
      onset_to_notification),
    proportion_infections = car,
    dow_model = TRUE),

  hospitalisations = define_observation_data(
    timeseries_data = hosp_dat,
    delay_from_infection = hosp_delay_ecdf,
    proportion_infections = ihr)
  # ,
  # other
)

fit_object <- fit_waves(
  observations = observation_models,
  infection_model = 'gp_growth_rate'# 'flat_prior'#,define_infection_model() #
)

rhats <- coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)
max(rhats$psrf[, 1], na.rm = TRUE)
bayesplot::mcmc_trace(fit_object$fit, pars = 'incidence[170,1]')

coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)

posterior_area_plots <- function(draws, select_pars, label = "Posterior") {
  plot_title <- ggtitle(
    paste0(label, " distributions"),
    "with medians and 90% intervals"
  )
  mcmc_areas(draws, pars = select_pars, prob = 0.9, point_est = "median") +
    plot_title +
    labs(y = "parameter")
}


#
# fitted_reff <- compute_reff(fit_object, gi)
#
#
# # outputs
# infections_out <- epiwave.pipelines::generate_long_estimates(
#   fit_object$infection_model$infection_timeseries,
#   fit_object$fit,
#   infection_days,
#   jurisdictions)

infection_traj <- epiwave.pipelines::build_trajectories(
  param = fit_object$infection_model,
  infection_days,
  fit_object$fit,
  nsim = 1000,
  jurisdictions)

infection_traj_diagnostic_vis(infection_traj)

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

#XXX
#
# infection_traj <- epiwave.pipelines::build_trajectories(
#   param = fit_object$infection_model$infection_timeseries,
#   infection_days,
#   fit_object$fit,
#   nsim = 1000,
#   jurisdictions)
#
# epiwave.pipelines::plot_reff_interval_curves(
#   'gp_both_reff_9Oct.pdf',
#   fitted_reff,
#   dates = infection_days,
#   start_date = min(study_seq),
#   end_date = max(study_seq),
#   jurisdictions = jurisdictions)
# infection_sims <- greta::calculate(
#   fit_object$infection_model$infection_timeseries,
#   values = fit_object$fit,
#   nsim = 1000)
#
# epiwave.pipelines::plot_timeseries_sims(
#   'outputs/gp_reff.png',
#   infection_sims[[1]],
#   type = "reff",
#   dates = infection_days,
#   states = jurisdictions,
#   start_date = study_seq[1],
#   end_date = study_seq[length(study_seq)],
#   dim_sim = "2",
#   infection_nowcast = TRUE,
#   nowcast_start = as.Date('2021-12-08'))
#
#
# epiwave.pipelines::plot_timeseries_sims(
#   'outputs/gp_both_infections_Kate-9Oct.png',
#   infection_sims[[1]],
#   type = "infection",
#   dates = infection_days,
#   states = jurisdictions,
#   start_date = study_seq[1],
#   end_date = study_seq[length(study_seq)],
#   dim_sim = "2",
#   infection_nowcast = TRUE,
#   nowcast_start = as.Date('2021-12-08'))
