## header
# library(epiwave)
library(epiwave.params)
library(dplyr)
library(greta)
devtools::load_all()
## user modified

study_seq <- seq(from = as.Date('2026-07-19'),
                      to = as.Date('2027-07-17'), 'days')
infection_days <- seq(from = as.Date('2026-06-01'),
                       to = as.Date('2027-07-17'), 'days')
jurisdictions <- 'study_area'


## data prep
# notification counts
notif_dat <- 'simdata/sim_study_cases.csv' |>
  read.csv() |>
  t() |>
  as.data.frame()
notif_dat$date <- study_seq
colnames(notif_dat) <- c('study_area', 'date')
notif_dat <- notif_dat |>
  tidyr::pivot_longer(!date, names_to = "jurisdiction", values_to = "value")
class(notif_dat) <- c('epiwave_fixed_timeseries',
                      'epiwave_timeseries',
                      class(notif_dat))

# hospitalisation counts
hosp_dat <- 'simdata/sim_study_hosp.csv' |>
  read.csv() |>
  t() |>
  as.data.frame()
hosp_dat$date <- study_seq
colnames(hosp_dat) <- c('study_area', 'date')
hosp_dat <- hosp_dat |>
  tidyr::pivot_longer(!date, names_to = "jurisdiction", values_to = "value")
class(hosp_dat) <- c('epiwave_fixed_timeseries',
                      'epiwave_timeseries',
                      class(hosp_dat))

# jurisdictions
jurisdictions <- unique(notif_dat$jurisdiction)
n_jurisdictions <- length(jurisdictions)

# proportions
car <- 0.42 # minimum viable product scenario
# ihr <- 0.005
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
#Exp(mean=2)
incubation <- distributional::dist_exponential(rate = 1/2) |>
  parametric_dist_to_distribution()
onset_to_notification <- distributional::dist_gamma(shape = 2, rate = 1) |>
  parametric_dist_to_distribution()
# Gamma(shape=2, scale =1)
# mean(rgamma(1000000, shape = 2,rate = 1))
onset_to_hospitalisation <- distributional::dist_gamma(shape = 5, rate = 1) |>
  parametric_dist_to_distribution()
# Gamma(shape = 5, scale = 1)

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
    dow_model = FALSE),

  hospitalisations = define_observation_data(
    timeseries_data = hosp_dat,
    delay_from_infection =
      epiwave.params::add_distributions(
        incubation, onset_to_hospitalisation),
    proportion_infections = ihr)
)


fit_object <- fit_waves(
  observations = observation_models,
  infection_model_type = 'gp_growth_rate')#,# define_infection_model()'gp_growth_rate' #
  # n_chains = 4,
  # max_convergence_tries = 2,
  # warmup = 1000,
  # n_samples = 1000,
  # # n_extra_samples = 100
#)

rhats <- coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)
max(rhats$psrf[, 1], na.rm = TRUE)
bayesplot::mcmc_trace(fit_object$fit, pars = 'incidence[157,1]')

coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)



ihr_out <- greta::calculate(ihr$ihr, values = fit_object$fit, nsim = 1000)

mean(ihr_out$`ihr$ihr`)

mean_ihr <- data.frame(apply(ihr_out$`ihr$ihr`[,,1], 2, mean))

# greta::calculate()

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
  fit_object$infection_model,
  values = fit_object$fit,
  nsim = 1000)

epiwave.pipelines::plot_timeseries_sims(
  'infections_gp_growth_rate_cases.png',
  infection_sims[[1]],
  type = "infection",
  dates = infection_days,
  states = jurisdictions,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)] - 14,
  dim_sim = "2",
  infection_nowcast = FALSE)

