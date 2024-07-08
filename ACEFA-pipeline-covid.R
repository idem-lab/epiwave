## header
library(epiwave)
library(epiwave.params)
library(epiwave.pipelines)
library(dplyr)
library(greta)
library(readr)
library(lubridate)
library(tidyr)
library(ggplot2)

#set forecast date for the week and then create output folder and date sequence from that
forecast_date <- as.Date("2024-06-27") #set for each week
origin_date <- as.Date(forecast_date - days(2))

outpath <- paste0("outputs/", forecast_date)
dir.create(outpath)

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
oldnotif_delay <- readRDS(paste0(not_synced_folder, '/delays/ECDF_delay_constant_PCR.rds'))
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
  effect_type = 'growth_rate')

## observation models
notif_observation_model_objects <- create_observation_model(
  infection_timeseries = infection_model_objects$infection_timeseries,
  delay_distribution = notif_full_delay_dist,
  proportion_observed = car,
  count_data = notif_dat,
  dow_model = dow_model,
  data_id = 'notif')

## Reff
reff_model_objects <- estimate_reff(
  infection_timeseries = infection_model_objects$infection_timeseries,
  generation_interval_mass_fxns = generation_interval_distribution)

## greta model fit
m <- greta::model(infection_model_objects$infection_timeseries,
                  reff_model_objects$reff)

plot(m)

fit <- epiwave.pipelines::fit_model(
  model = m,
  n_chains = 10,
  max_convergence_tries = 1,
  warmup = 1000,
  n_samples = 1000,
  n_extra_samples = 100)


#calculate nowcast period - currently 50% like earlier in COVID pandemic
#(according to Freya) they only shifted to 95% because there were right truncation or data completeness issues issues

nowcast_shift <- notif_delay_dist$delay[which((cumsum(notif_delay_dist$mass) > 0.5))[1]]
nowcast_start <- origin_date - nowcast_shift

#plot reff

reff_sims <- greta::calculate(reff_model_objects$reff,
                              values = fit,
                              nsim = 1000)
#reff one month
plot_timeseries_sims(
  paste0(outpath,'/covid_reff_one_month.pdf'),
  reff_sims[[1]],
  type = "reff",
  dates = infection_days,
  start_date = origin_date - months(1),
  end_date = origin_date,
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = nowcast_start)

#reff six month
plot_timeseries_sims(
  paste0(outpath,'/covid_reff_six_months.pdf'),
  reff_sims[[1]],
  type = "reff",
  dates = infection_days,
  start_date = origin_date - months(6),
  end_date = origin_date,
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = nowcast_start)

#plot infection timeseries
infection_sims <- greta::calculate(infection_model_objects$infection_timeseries,
                                   values = fit,
                                   nsim = 1000)

plot_timeseries_sims(
  paste0(outpath,'/covid_infection_timeseries.pdf'),
  infection_sims[[1]],
  type = "infection",
  dates = infection_days,
  start_date = study_seq[1],
  end_date = study_seq[length(study_seq)],
  states = jurisdictions,
  dim_sim = "2",
  infection_nowcast = TRUE,
  nowcast_start = nowcast_start)

#estimate of reff on last day before nowcast and last of nowcast
mean <- apply(reff_sims[[1]], 2:3, FUN = "mean", na.rm = TRUE)
ci_95_lo <- apply(reff_sims[[1]], 2:3, quantile, c(0.025), na.rm = TRUE)
ci_95_hi <- apply(reff_sims[[1]], 2:3, quantile, c(0.975), na.rm = TRUE)
reff_est <- data.frame(infection_days, mean, ci_95_lo, ci_95_hi)

nowcast_50 <- nowcast_start
nowcast_95 <- as_date(origin_date - days(14))

reff_dates <- data.frame(infection_days = c(nowcast_50, nowcast_95, origin_date), date_label = c("nowcast_50", "nowcast_95", "origin"))
reff_estimates <- reff_est |>
  filter(infection_days %in% reff_dates$infection_days)

reff_all <- full_join(reff_estimates, reff_dates)
reff_all <- reff_all[,c(5,1:4)]

write_csv(reff_all, file = paste0(outpath,"/reff_estimates_2024-07-04.csv"))


##END - below is just alternative plotting code

#code for alternative reff plots based on curvewise intervals
# long is classic output, then map to infection traj.
infections_out <- epiwave.pipelines::generate_long_estimates(
  infection_model_objects$infection_timeseries,
  fit,
  infection_days,
  jurisdictions)

infectiousness_out <- epiwave.pipelines::generate_long_estimates(
  reff_model_objects$infectiousness,
  fit,
  infection_days,
  jurisdictions)

reff_out <- epiwave.pipelines::generate_long_estimates(
  reff_model_objects$reff,
  fit,
  infection_days,
  jurisdictions)

#saveRDS(reff_out, file = paste0(outpath,"/reff_out_",forecast_date,".RDS"))

infection_traj <- epiwave.pipelines::build_trajectories(
  param = infection_model_objects$infection_timeseries,
  infection_days,
  fit,
  nsim = 1000,
  jurisdictions)

#plot reff one month - change axis scale breaks manually to date_breaks = "1 week"
#change line alpha to 0.02
plot_reff_interval_curves(
  'covid_reff_1month.pdf',
  reff_out,
  dates = infection_days,
  start_date = max(study_seq) - months(1),
  end_date = max(study_seq))

#plot reff six month - change axis scale breaks manually to date_breaks = "4 weeks"
#change line alpha to 0.002
plot_reff_interval_curves(
  'covid_reff_6months.pdf',
  reff_out,
  dates = infection_days,
  start_date = min(study_seq),
  #start_date = max(study_seq) - months(1),
  end_date = max(study_seq))
