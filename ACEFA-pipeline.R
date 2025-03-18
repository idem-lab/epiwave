## header
library(epiwave)
library(epiwave.params)
library(epiwave.pipelines)
library(dplyr)
library(greta)
library(readr)
library(lubridate)

origin_date <- as.Date("2024-08-27")#set for each week
forecast_date <- as.Date("2024-08-29")

study_seq <- seq(from = as.Date('2024-01-01'),
                 to = origin_date, 'days')

## user specific folder with not synced data
not_synced_folder <- '~/not_synced/ACEFA'

## notification counts
notif_file <- paste0(not_synced_folder, '/SARSCOV2-case-count-', forecast_date, ".csv")

notif_dat <- notif_file |>
  read_csv() |>
  dplyr::rename(value = cases, jurisdiction = location, date = notification_date) |>
  dplyr::select(!(test_type)) |>
  dplyr::filter(date %in% study_seq)

class(notif_dat) <- c('epiwave_fixed_timeseries',
                      'epiwave_timeseries',
                      class(notif_dat))
# jurisdictions
jurisdictions <- unique(notif_dat$jurisdiction)
n_jurisdictions <- length(jurisdictions)

#delay distribution
oldnotif_delay <- readRDS(paste0(not_synced_folder, '/ECDF_delay_constant_PCR.rds'))
min_delay <- 0
max_delay <- 41
notif_delay_ecdf <- create_epiwave_massfun(
  min_delay, max_delay, oldnotif_delay, normalise = TRUE)

## days of infection timeseries
infection_days <- seq(from = as.Date(min(study_seq) - days(28)),
                      to = as.Date(origin_date + days(28)), 'days')

n_days_infection <- length(infection_days)

## proportions
car_constant <- 1
car <- create_epiwave_fixed_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = car_constant) # 1 for sero prop

## incubation period
incub_dist <- distributional::dist_weibull(
  shape = 1.83, scale = 4.93)
incubation_period_distribution <- parametric_dist_to_distribution(incub_dist)

## formatted delay distribution
notif_delay_dist <-
  add(notif_delay_ecdf, incubation_period_distribution) #add function currently not working, load in manually
notif_full_delay_dist <- create_epiwave_massfun_timeseries(
  dates = infection_days,
  jurisdictions = jurisdictions,
  value = notif_delay_dist)


## generation interval distribution
#gi_distribution_data <- readr::read_csv(
 # file = paste0(not_synced_folder,'/nishiura_samples.csv'),
  #col_types = readr::cols(param1 = readr::col_double(),
   #                       param2 = readr::col_double()))

#gi_dist <- distributional::dist_lognormal(
 # mu = mean(gi_distribution_data$param1),
  #sigma = mean(gi_distribution_data$param2))
#generation_interval_distribution <- parametric_dist_to_distribution(gi_dist)

gi_dist <-  distributional::dist_gamma(shape = 0.89, rate = 0.27)
generation_interval_distribution <- parametric_dist_to_distribution(gi_dist)

## optional day-of-week effect model
dow_model <- epiwave::create_dow_priors(
  n_jurisdictions)

## infection timeseries model
infection_model_objects <- epiwave::create_infection_timeseries(
  n_days_infection,
  n_jurisdictions,
  effect_type = 'growth_rate_deriv')

## observation models
notif_observation_model_objects <- create_observation_model(
  infection_timeseries = infection_model_objects$infection_timeseries,
  delay_distribution = notif_full_delay_dist,
  proportion_observed = car,
  count_data = notif_dat,
  dow_model = dow_model,
  data_id = 'notif')

## Reff
reff_model_objects <- epiwave::estimate_reff(
  infection_timeseries = infection_model_objects$infection_timeseries,
  generation_interval_mass_fxns = generation_interval_distribution)

## greta model fit
m <- greta::model(infection_model_objects$infection_timeseries,
                  reff_model_objects$reff)

plot(m)

z_raw <- attr(infection_model_objects$gp, "gp_info")$v
n_chains <- 10
inits <- replicate(n_chains,
                   initials(
                     z_raw = rep(0, length(z_raw))
                   ),
                   simplify = FALSE)

fit <- greta::mcmc(
  model = m,
  sampler = hmc(Lmin = 35, Lmax = 40),
  chains = n_chains,
  warmup = 2e3,
  n_samples = 2e3,
  initial_values = inits,
  one_by_one = TRUE)

r_hats <- coda::gelman.diag(fit, autoburnin = FALSE,
                            multivariate = FALSE)$psrf[, 1]
# infection R hat witout nowcast phase
summary(r_hats[1:(length(study_seq)+28)])
# reff R hat witout nowcast phase
summary(r_hats[(1:(length(study_seq)+28))+ n_days_infection])
# # long is classic output, then map to infection traj.
# infections_out <- epiwave.pipelines::generate_long_estimates(
#   infection_model_objects$infection_timeseries,
#   fit,
#   infection_days,
#   jurisdictions)
#
# infectiousness_out <- epiwave.pipelines::generate_long_estimates(
#   reff_model_objects$infectiousness,
#   fit,
#   infection_days,
#   jurisdictions)
#
# infection_traj <- epiwave.pipelines::build_trajectories(
#   param = infection_model_objects$infection_timeseries,
#   infection_days,
#   fit,
#   nsim = 1000,
#   jurisdictions)
#
#
# #plot_reff <- plot_timeseries_sims(
#  # reff_out,
#   #type = "reff",
#   #days = infection_days)
#
#
# epiwave.pipelines::plot_reff_interval_curves(
#   'reff.png',
#   reff_out,
#   dates = infection_days,
#   start_date = min(study_seq),
#   end_date = max(study_seq),
#   jurisdictions = jurisdictions)
# growth_rate <- (greta::apply(infection_model_objects$gp, 2, 'cumsum'))
# #growth_rate <- infection_model_objects$gp
#
# sim_check(growth_rate,
#           title = "prior growth rate")
#
# sim_check(growth_rate,
#           fitted_values = fit,
#           title = "posterior growth rate")
#
# sim_check(infection_model_objects$infection_timeseries,
#           fitted_values = fit,
#           title = "posterior infections")

infection_sims <- greta::calculate(infection_model_objects$infection_timeseries,
                                   values = fit,
                                   nsim = 1000)

nowcast_shift <-notif_delay_dist$delay[which((cumsum(notif_delay_dist$mass) > 0.95))[1]]

nowcast_date <- origin_date - nowcast_shift

plot_timeseries_sims(
  'infection_timeseries2.png',
  infection_sims[[1]],
  type = "infection",
  dates = infection_days,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = nowcast_date)

reff_sims <- greta::calculate(reff_model_objects$reff,
                                   values = fit,
                                   nsim = 1000)

plot_timeseries_sims(
  'reff_timeseries2.png',
  reff_sims[[1]],
  type = "reff",
  dates = infection_days,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = nowcast_date)

plot_timeseries_sims(
  'reff_timeseries2_1month.png',
  reff_sims[[1]],
  type = "reff",
  dates = infection_days,
  start_date = origin_date - days(30),
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = nowcast_date)

calculate(infection_model_objects$gp_lengthscale,
          values = fit,
          nsim = 1000)[[1]] %>% summary
calculate(infection_model_objects$gp_variance,
          values = fit,
          nsim = 1000)[[1]] %>% summary

calculate(dow_model$dow_dist,
          values = fit,
          nsim = 1000)[[1]] %>% apply(3,mean)

(1/(calculate(notif_observation_model_objects$notif_size[1],
          values = fit,
          nsim = 1000)[[1]] %>% apply(3,mean)))^2

case_sims <- greta::calculate(notif_observation_model_objects$notif_case_mat,
                              values = fit,
                              nsim = 1000)

residual_diag_cases(case_sims,
                    cases_true = notif_dat$value)

plot_timeseries_sims(
  'cases_timeseries2.png',
  case_sims[[1]],
  type = "notification",
  dates = study_seq,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = FALSE,
  #nowcast_start = nowcast_date,
  case_validation_data = notif_dat %>% rename("count" = "value"))


case_w_forecast_sims <- greta::negative_binomial(notif_observation_model_objects$notif_size,
                                                 notif_observation_model_objects$notif_prob)
case_w_forecast_sims <- greta::calculate(case_w_forecast_sims,
                              values = fit,
                              nsim = 1000)

plot_timeseries_sims(
  'cases_timeseries2_w_forecast.png',
  case_w_forecast_sims[[1]],
  type = "notification",
  dates = study_seq,
  start_date = origin_date - days(30),
  states = jurisdictions,
  dim_sim = "2",
  case_forecast =  TRUE,
  infection_nowcast = TRUE,
  nowcast_start = origin_date,
  case_validation_data = notif_dat %>% rename("count" = "value"))
